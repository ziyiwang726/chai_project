#!/usr/bin/env python3
import os
import time
import re
import glob
import json
import pandas as pd
import requests
from xml.etree import ElementTree as ET

from llm_provider_client import (
    assert_provider_ready,
    call_llm_json,
    ensure_parent_dir,
    get_provider_file_tag,
    get_provider_model,
    normalize_provider,
    pricing_for_model,
    runtime_path,
    summarize_usage,
    with_provider_tag,
)
from ncbi_request_utils import BASE_EUTILS as NCBI_DB, NCBIClient, PMC_IDCONV
from taxon_result_utils import (
    build_target_payloads,
    is_result_specific_to_taxon,
    load_taxon_metadata,
    lookup_result_for_target,
)

# --- Configuration ---
INPUT_FILE = "taxon_article_fulltext_matrix.csv"
OUTPUT_FILE = os.getenv("LLM_OUTPUT_FILE", with_provider_tag("taxon_article_matrix_filled.csv"))
FALLBACK_PREV_FILE = os.getenv("LLM_FALLBACK_PREV_FILE", with_provider_tag("taxon_article_matrix_filled_fixed.csv"))
TAXON_INFO_FILE = "taxa_ids_filtered.csv"
LLM_PROVIDER = normalize_provider()
MODEL_NAME = get_provider_model(LLM_PROVIDER)
REASONING_EFFORT = os.getenv("OPENAI_REASONING_EFFORT", "high")
USAGE_TRACKER_FILE = os.getenv(
    "OPENAI_COST_TRACKER_FILE",
    runtime_path("cost", f"{LLM_PROVIDER}_cost_tracker.json"),
)
RESET_USAGE_TRACKER = os.getenv("OPENAI_COST_TRACKER_RESET", "0") == "1"

try:
    assert_provider_ready(LLM_PROVIDER)
except Exception as exc:
    print(f"Error: {exc}")
    raise SystemExit(1)

# Define contexts for marker column
CONTEXT_DEFINITIONS = {
    "A": "Melanoma RESPONSE to immunotherapy (RECIST/ORR etc.)",
    "B": "Melanoma SURVIVAL (OS/PFS)",
    "C": "Other cancers + any checkpoint or immunotherapy evidence (exclude melanoma)",
}
CONTEXT_PROMPT = "\n".join([f"- {k}: {v}" for k, v in CONTEXT_DEFINITIONS.items()])

REVIEW_PUB_TYPES = {"Review", "Systematic Review", "Scoping Review"}
META_PUB_TYPES = {"Meta-Analysis", "Systematic Review"}
ANIMAL_MESH_MARKERS = {
    "animals", "mice", "mouse", "rats", "rat", "dogs", "dog", "swine",
    "xenograft", "murinae", "murine",
}
ANIMAL_TEXT_MARKERS = {
    " mouse ", " mice ", " murine ", " rat ", " rats ", " xenograft ",
    " canine ", " porcine ", " macaque ",
}
HUMAN_TEXT_MARKERS = {
    " patient ", " patients ", " human ", " humans ",
    " clinical trial ", " cohort ", " prospective ",
}
LOWER_LEVEL_MARKERS = (
    " otu", "otus", " asv", "asvs", " sgb", "sgbs", " strain",
    " amplicon", " sequence variant", " metagenome-assembled genome",
)


class TokenTracker:
    def __init__(self):
        if RESET_USAGE_TRACKER and os.path.exists(USAGE_TRACKER_FILE):
            try:
                os.remove(USAGE_TRACKER_FILE)
            except OSError:
                pass
        state = self._load_state()
        self.prompt_tokens = int(state.get("prompt_tokens", 0) or 0)
        self.cached_prompt_tokens = int(state.get("cached_prompt_tokens", 0) or 0)
        self.completion_tokens = int(state.get("completion_tokens", 0) or 0)
        self.estimated_cost_usd = float(state.get("estimated_cost_usd", 0.0) or 0.0)

    def update(self, usage):
        usage_summary = summarize_usage(usage, provider=LLM_PROVIDER, model_name=MODEL_NAME)
        self.prompt_tokens += usage_summary["prompt_tokens"]
        self.cached_prompt_tokens += usage_summary["cached_prompt_tokens"]
        self.completion_tokens += usage_summary["completion_tokens"]
        self.estimated_cost_usd += usage_summary["estimated_cost_usd"]
        self._persist_state()
        return usage_summary

    def _load_state(self):
        if not os.path.exists(USAGE_TRACKER_FILE):
            return {}
        try:
            with open(USAGE_TRACKER_FILE, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            return {}

    def _persist_state(self):
        payload = {
            "model": MODEL_NAME,
            "prompt_tokens": self.prompt_tokens,
            "cached_prompt_tokens": self.cached_prompt_tokens,
            "completion_tokens": self.completion_tokens,
            "estimated_cost_usd": self.estimated_cost_usd,
        }
        with open(ensure_parent_dir(USAGE_TRACKER_FILE), "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)


token_tracker = TokenTracker()
ncbi = NCBIClient("find_pvalues")

taxon_metadata = load_taxon_metadata(TAXON_INFO_FILE)
CHECKPOINT_TAG = get_provider_file_tag(LLM_PROVIDER)
CHECKPOINT_DIR = os.getenv("LLM_CHECKPOINT_DIR") or runtime_path("checkpoints", is_dir=True)
os.makedirs(CHECKPOINT_DIR, exist_ok=True)
MAX_TAXA_PER_REQUEST = max(int(os.getenv("LLM_MAX_TAXA_PER_REQUEST", "8")), 1)


def checkpoint_glob():
    suffix = f"_{CHECKPOINT_TAG}" if CHECKPOINT_TAG else ""
    return os.path.join(CHECKPOINT_DIR, f"checkpoint_col_*{suffix}.csv")


def checkpoint_path(article_id):
    suffix = f"_{CHECKPOINT_TAG}" if CHECKPOINT_TAG else ""
    return os.path.join(CHECKPOINT_DIR, f"checkpoint_col_{article_id}{suffix}.csv")


def checkpoint_article_id(path):
    base = os.path.basename(path)
    pattern = r"^checkpoint_col_(.+?)(?:_" + re.escape(CHECKPOINT_TAG) + r")?\.csv$" if CHECKPOINT_TAG else r"^checkpoint_col_(.+)\.csv$"
    match = re.match(pattern, base)
    return match.group(1) if match else None


def chunked(values, chunk_size):
    for start in range(0, len(values), chunk_size):
        yield values[start:start + chunk_size]


def combine_usage_summaries(usage_summaries):
    summaries = [summary for summary in usage_summaries if summary]
    if not summaries:
        return None
    return {
        "prompt_tokens": sum(summary["prompt_tokens"] for summary in summaries),
        "cached_prompt_tokens": sum(summary["cached_prompt_tokens"] for summary in summaries),
        "completion_tokens": sum(summary["completion_tokens"] for summary in summaries),
        "estimated_cost_usd": sum(summary["estimated_cost_usd"] for summary in summaries),
        "long_context_applied": any(summary["long_context_applied"] for summary in summaries),
        "pricing_known": all(summary["pricing_known"] for summary in summaries),
    }


def retry_sleep_seconds(exc, attempt, base=4):
    msg = str(exc).lower()
    if (
        "429" in msg
        or "529" in msg
        or "too many requests" in msg
        or "rate limit" in msg
        or "server error" in msg
        or "internal server error" in msg
    ):
        return min(60, 15 * attempt)
    return base * attempt


def _text_flags(title, abstract=""):
    txt = f" {title or ''} {abstract or ''} ".lower()
    has_animal = any(k in txt for k in ANIMAL_TEXT_MARKERS)
    has_human = any(k in txt for k in HUMAN_TEXT_MARKERS)
    return has_animal, has_human


def find_first_text_by_localname(root, localname):
    """
    Namespace-agnostic XML text lookup by local tag name.
    """
    if root is None:
        return ""
    for elem in root.iter():
        if elem.tag.split("}")[-1] == localname:
            txt = "".join(elem.itertext()).strip()
            if txt:
                return txt
    return ""


def pmcid_to_pmid(pmcid):
    try:
        resp = ncbi.get(
            PMC_IDCONV,
            params=ncbi.default_params({"ids": pmcid, "format": "json"}),
            timeout=20,
        )
        resp.raise_for_status()
        recs = resp.json().get("records", [])
        if recs:
            pmid = str(recs[0].get("pmid") or "").strip()
            return pmid if pmid else None
    except Exception:
        return None
    return None


def fetch_pubmed_metadata_by_pmid(pmid):
    """Metadata screening: review/meta/pure-animal (text + pubmed indexing)."""
    out = {
        "valid": False,
        "title": "",
        "is_review": False,
        "is_meta": False,
        "is_animal_only": False,
        "reason": "",
        "title_abstract": "",
    }
    try:
        resp = ncbi.post(
            f"{NCBI_DB}/efetch.fcgi",
            data=ncbi.default_params({"db": "pubmed", "id": str(pmid), "retmode": "xml"}),
            timeout=30,
        )
        resp.raise_for_status()
        root = ET.fromstring(resp.content)
        art = root.find(".//PubmedArticle")
        if art is None:
            out["reason"] = "invalid"
            return out

        title = (art.findtext(".//ArticleTitle") or "").strip()
        abs_texts = art.findall(".//Abstract/AbstractText")
        abstract = " ".join([(t.text or "") for t in abs_texts]).strip()
        out["title"] = title
        out["title_abstract"] = f"TITLE: {title}\nABSTRACT: {abstract}".strip()
        out["valid"] = bool(title)

        pub_types = {
            (pt.text or "").strip()
            for pt in art.findall(".//PublicationType")
            if (pt.text or "").strip()
        }
        mesh_terms = {
            (mh.findtext("DescriptorName") or "").strip().lower()
            for mh in art.findall(".//MeshHeading")
        }

        is_review = bool(REVIEW_PUB_TYPES & pub_types) or (
            re.search(r"\breview\b", title, flags=re.IGNORECASE) is not None
        )
        is_meta = bool(META_PUB_TYPES & pub_types) or (
            re.search(r"\bmeta[- ]analysis\b|\bsystematic review\b", title, flags=re.IGNORECASE) is not None
        )

        has_humans_mesh = "humans" in mesh_terms
        has_animals_mesh = any(
            (m in ANIMAL_MESH_MARKERS) or any(k in m for k in ANIMAL_MESH_MARKERS)
            for m in mesh_terms
        )
        has_animal_text, has_human_text = _text_flags(title, abstract)
        is_animal_only = (
            (has_animals_mesh and not has_humans_mesh) or
            (has_animal_text and not has_human_text and not has_humans_mesh)
        )

        out["is_review"] = is_review
        out["is_meta"] = is_meta
        out["is_animal_only"] = is_animal_only

        if is_review:
            out["reason"] = "review"
        elif is_meta:
            out["reason"] = "meta_analysis"
        elif is_animal_only:
            out["reason"] = "animal_only"
        else:
            out["reason"] = "ok"
    except Exception:
        out["reason"] = "metadata_fetch_error"
    return out


def screen_article(doc_id):
    """
    Pre-LLM NCBI screen. Uses PubMed metadata when possible.
    Returns dict with skip flag and reason.
    """
    if doc_id.upper().startswith("PMID"):
        pmid = doc_id.replace("PMID", "")
        meta = fetch_pubmed_metadata_by_pmid(pmid)
        return {
            "skip": (not meta["valid"]) or meta["is_review"] or meta["is_meta"] or meta["is_animal_only"],
            "reason": meta.get("reason", ""),
            "meta": meta,
        }

    # PMCID path: map to PMID then screen via PubMed metadata
    pmid = pmcid_to_pmid(doc_id)
    if pmid:
        meta = fetch_pubmed_metadata_by_pmid(pmid)
        return {
            "skip": (not meta["valid"]) or meta["is_review"] or meta["is_meta"] or meta["is_animal_only"],
            "reason": meta.get("reason", ""),
            "meta": meta,
        }

    # No PMID mapping: do conservative title-only checks from PMC
    try:
        clean_id = doc_id.replace("PMC", "")
        resp = ncbi.post(
            f"{NCBI_DB}/efetch.fcgi",
            data=ncbi.default_params({"db": "pmc", "id": clean_id, "retmode": "xml"}),
            timeout=30,
        )
        resp.raise_for_status()
        root = ET.fromstring(resp.content)
        title = find_first_text_by_localname(root, "article-title")
        title_l = title.lower()
        is_review = "review" in title_l
        is_meta = ("meta-analysis" in title_l) or ("meta analysis" in title_l)
        return {
            "skip": (not bool(title)) or is_review or is_meta,
            "reason": "review" if is_review else ("meta_analysis" if is_meta else ("ok" if title else "invalid")),
            "meta": {"title": title},
        }
    except Exception:
        return {"skip": True, "reason": "metadata_fetch_error", "meta": {}}


def get_paper_text(doc_id):
    """Fetch full text (PMC) or title+abstract (PMID)."""
    is_pmc = doc_id.upper().startswith("PMC")
    clean_id = doc_id.replace("PMID", "").replace("PMC", "")
    text = ""
    try:
        db = "pmc" if is_pmc else "pubmed"
        resp = ncbi.post(
            f"{NCBI_DB}/efetch.fcgi",
            data=ncbi.default_params({"db": db, "id": clean_id, "retmode": "xml"}),
            timeout=40,
        )
        resp.raise_for_status()
        root = ET.fromstring(resp.content)

        if is_pmc:
            text = "".join(root.itertext())
        else:
            art = root.find(".//PubmedArticle")
            if art is not None:
                title = art.findtext(".//ArticleTitle") or ""
                abs_texts = art.findall(".//Abstract/AbstractText")
                abstract = " ".join([t.text or "" for t in abs_texts])
                text = f"TITLE: {title}\nABSTRACT: {abstract}"
    except Exception:
        return None

    text = (text or "").strip()
    if not text:
        return None
    return text


def analyze_paper_batch(paper_id, paper_text, taxa_list):
    """Analyze one paper for multiple taxa in one provider request."""
    if not paper_text:
        return {}, None

    system_prompt = (
        "You extract statistically significant taxon-associated outcomes from biomedical paper text. "
        "Return strict JSON only. Never include markdown."
    )
    target_payloads = build_target_payloads(taxa_list, taxon_metadata)

    user_prompt = f"""
Paper ID: {paper_id}

Target taxa JSON list:
{json.dumps(target_payloads, ensure_ascii=True)}

Contexts (marker):
{CONTEXT_PROMPT}

Rules:
1) Return ONLY taxa with a statistically significant result.
2) For each taxon, return ONE best context marker in {{A,B,C}}.
3) Direction must encode patient benefit or harm, not the raw coefficient sign.
4) Use '+' when the reported association benefits patients and '-' when it harms patients. If benefit vs harm cannot be determined from text, set '?'.
5) For response-type outcomes, '+' means associated with response, responders, clinical benefit, disease control, CR/PR, or other favorable response outcomes; '-' means associated with non-response, resistance, progression, or other unfavorable response outcomes.
6) For survival or time-to-event outcomes, '+' means longer survival or lower hazard/risk of progression, recurrence, or death; '-' means shorter survival or higher hazard/risk of progression, recurrence, or death.
7) Apply rule 6 to Cox models, Fine-Gray models, cause-specific hazard models, and related survival models. HR/subdistribution HR/cause-specific HR < 1 implies '+', HR/subdistribution HR/cause-specific HR > 1 implies '-'. A negative coefficient in a hazard model implies '+', and a positive coefficient implies '-'.
8) For marker C, still use '+' for patient-beneficial associations and '-' for patient-harmful associations. If the reported outcome is response-like, apply rule 5; if it is survival/time-to-event, apply rules 6-7.
9) p-value: if significant but not numerically reported, set "0.05".
10) Extract sample size as integer n_obs when available (total analyzed observations); else null.
11) Do not infer from figures/plots. Use text only.
12) For a species target, OTU/ASV/SGB/strain evidence will count when it ties to that exact species.
13) For any target above species, do not treat species, OTUs, ASVs, SGBs, strains, or any finer descendants as evidence for the target.
14) Use the exact provided target taxon string as the JSON key.

Output schema (JSON object):
{{
  "TaxonName": {{
    "p": "0.003",
    "dir": "+",
    "marker": "A",
    "eff": "HR 0.72",
    "n_obs": 128
  }}
}}

Paper text:
{paper_text}
""".strip()

    for attempt in range(1, 4):
        try:
            parsed, usage, _ = call_llm_json(
                system_prompt,
                user_prompt,
                provider=LLM_PROVIDER,
                reasoning_effort=REASONING_EFFORT,
            )
            usage_summary = token_tracker.update(usage)
            if isinstance(parsed, list):
                parsed = parsed[0] if parsed else {}
            if not isinstance(parsed, dict):
                return {}, usage_summary
            return parsed, usage_summary
        except Exception as e:
            print(f"  [Attempt {attempt}/3] {LLM_PROVIDER} error on {paper_id}: {e}")
            if attempt < 3:
                time.sleep(retry_sleep_seconds(e, attempt))
    return {}, None


def parse_p_numeric(raw):
    s = str(raw or "").strip()
    if not s:
        return "0.05"
    s = s.lower().replace("p", "")
    s = re.sub(r"[=\s]", "", s)

    # handle threshold styles conservatively
    if s.startswith("<"):
        s = s[1:]
    if s.startswith("<="):
        s = s[2:]

    m = re.search(r"\d*\.?\d+(?:e-?\d+)?", s)
    if not m:
        return "0.05"

    try:
        p = float(m.group(0))
    except Exception:
        return "0.05"

    if not (0 < p <= 1):
        return "0.05"
    return f"{p:.6g}"


def parse_n_obs(raw):
    if raw is None:
        return None
    if isinstance(raw, (int, float)):
        n = int(raw)
        return n if n > 0 else None
    s = str(raw).strip()
    m = re.search(r"\d+", s)
    if not m:
        return None
    n = int(m.group(0))
    return n if n > 0 else None


def build_cell_value(res_obj):
    if not isinstance(res_obj, dict):
        return ""

    direction = str(res_obj.get("dir", "")).strip()
    if direction not in {"+", "-"}:
        # Explicitly exclude significant results with missing/unknown direction
        return ""

    marker = str(res_obj.get("marker", "")).strip().upper()
    if marker not in {"A", "B", "C"}:
        return ""

    p_val = parse_p_numeric(res_obj.get("p"))
    eff = res_obj.get("eff")
    n_obs = parse_n_obs(res_obj.get("n_obs"))

    n_part = f" {{n={n_obs}}}" if n_obs is not None else ""
    eff_part = f" [{eff}]" if isinstance(eff, str) and eff.strip() else ""
    return f"{direction}{p_val} {marker}{n_part}{eff_part}".strip()


def merge_existing_results(df, taxon_col, article_cols):
    """
    Reuse existing extracted results so only new/unfilled cells are sent to the active provider.
    Priority: OUTPUT_FILE, then FALLBACK_PREV_FILE.
    """
    prev_path = None
    for candidate in [OUTPUT_FILE, FALLBACK_PREV_FILE]:
        if os.path.exists(candidate):
            prev_path = candidate
            break
    if prev_path is None:
        return df, 0, None, set()

    prev = pd.read_csv(prev_path).fillna("")
    if taxon_col not in prev.columns:
        return df, 0, prev_path, set()

    prev_idx = prev.set_index(taxon_col, drop=False)
    reused = 0
    scanned_cols = set()

    for col in article_cols:
        if col not in prev_idx.columns:
            continue
        col_nonempty = prev_idx[col].astype(str).str.strip().ne("").any()
        if col_nonempty:
            scanned_cols.add(col)
        old_vals = prev_idx.reindex(df[taxon_col])[col].fillna("").astype(str)
        keep_mask = (old_vals.str.strip() != "").to_numpy()
        if keep_mask.any():
            df.loc[keep_mask, col] = old_vals.to_numpy()[keep_mask]
            reused += int(keep_mask.sum())

    return df, reused, prev_path, scanned_cols


def main():
    print("[CHECKPOINT] Loading input matrix...")
    checkpoint_files = glob.glob(checkpoint_glob())
    df = pd.read_csv(INPUT_FILE)

    # Keep non-target cells empty strings
    df = df.replace([0, "0"], "")
    df = df.fillna("")

    taxon_col = df.columns[0]
    article_cols = [c for c in df.columns if ("PMC" in c or "PMID" in c)]

    # Reuse previously extracted results to avoid re-scanning old articles.
    df, reused_cells, reused_from, scanned_cols = merge_existing_results(df, taxon_col, article_cols)
    if reused_from:
        print(f"[CHECKPOINT] Reused {reused_cells} existing extracted cells from {reused_from}")
        print(f"[CHECKPOINT] Will skip {LLM_PROVIDER} re-scan for {len(scanned_cols)} already-scanned article columns")

    start_index = 0
    if checkpoint_files:
        latest_file = max(checkpoint_files, key=os.path.getctime)
        print(f"[CHECKPOINT] Found checkpoint: {latest_file}. Resuming...")
        df = pd.read_csv(latest_file).fillna("")
        last_id = checkpoint_article_id(latest_file)
        if last_id:
            try:
                start_index = article_cols.index(last_id) + 1
                print(f"[CHECKPOINT] Resume index: {start_index}")
            except ValueError:
                print("[CHECKPOINT] Checkpoint column not found in current input. Restarting from 0.")

    excluded_counts = {"review": 0, "meta_analysis": 0, "animal_only": 0, "invalid": 0, "metadata_fetch_error": 0}
    pricing = pricing_for_model(LLM_PROVIDER, MODEL_NAME)
    cached_rate = pricing["cached_input"]
    cached_rate_text = f"${cached_rate:.2f}/1M cached-input" if cached_rate is not None else "no cached-input rate"
    print(
        f"[CHECKPOINT] provider={LLM_PROVIDER}, model={MODEL_NAME}, reasoning_effort={REASONING_EFFORT}, "
        f"pricing=${pricing['input']:.2f}/1M input, {cached_rate_text}, ${pricing['output']:.2f}/1M output"
    )

    for idx in range(start_index, len(article_cols)):
        col = article_cols[idx]

        # Never re-scan article columns that were already processed in previous outputs.
        if col in scanned_cols:
            flagged = df[col].astype(str).str.strip().isin(["1", "1.0"])
            if flagged.any():
                # Clear remaining flags so future runs also skip these cells.
                df.loc[flagged, col] = ""
            continue

        mask = df[col].astype(str).str.strip().isin(["1", "1.0"])
        target_indices = df.index[mask].tolist()
        if not target_indices:
            continue

        target_taxa = df.loc[target_indices, taxon_col].tolist()
        print(f"[{idx+1}/{len(article_cols)}] {col}: {len(target_taxa)} taxa flagged")

        # 1) Pre-LLM NCBI screening
        screening = screen_article(col)
        if screening["skip"]:
            reason = screening.get("reason", "invalid") or "invalid"
            excluded_counts[reason] = excluded_counts.get(reason, 0) + 1
            df.loc[target_indices, col] = ""
            print(f"  -> Skipped before LLM ({reason})")
        else:
            # 2) Fetch text
            text = get_paper_text(col)
            if not text:
                df.loc[target_indices, col] = ""
                print("  -> Skipped (no retrievable text)")
            else:
                text_source = "PMC full text" if col.upper().startswith("PMC") else "PubMed title+abstract"
                print(f"  -> Retrieved {len(text):,} characters from {text_source}")
                # 3) LLM extraction
                results = {}
                usage_parts = []
                total_batches = (len(target_taxa) + MAX_TAXA_PER_REQUEST - 1) // MAX_TAXA_PER_REQUEST
                for batch_num, taxa_batch in enumerate(chunked(target_taxa, MAX_TAXA_PER_REQUEST), start=1):
                    batch_results, batch_usage = analyze_paper_batch(col, text, taxa_batch)
                    if batch_results:
                        results.update(batch_results)
                    if batch_usage:
                        usage_parts.append(batch_usage)
                    if total_batches > 1:
                        print(f"  -> Batch {batch_num}/{total_batches} processed {len(taxa_batch)} taxa")
                usage_summary = combine_usage_summaries(usage_parts)
                if usage_summary:
                    cost_bits = [
                        f"prompt={usage_summary['prompt_tokens']:,}",
                        f"completion={usage_summary['completion_tokens']:,}",
                        f"article=${usage_summary['estimated_cost_usd']:.4f}",
                        f"cumulative=${token_tracker.estimated_cost_usd:.4f}",
                    ]
                    if usage_summary["cached_prompt_tokens"]:
                        cost_bits.insert(1, f"cached={usage_summary['cached_prompt_tokens']:,}")
                    print(f"  -> Cost: {', '.join(cost_bits)}")
                    if usage_summary["long_context_applied"]:
                        print("     Long-context pricing applied (>272,000 input tokens).")

                for row_idx, taxon_name in zip(target_indices, target_taxa):
                    res_obj = lookup_result_for_target(results, taxon_name, taxon_metadata)
                    if not is_result_specific_to_taxon(taxon_name, res_obj, taxon_metadata):
                        res_obj = None
                    cell = build_cell_value(res_obj)
                    df.at[row_idx, col] = cell

        # 4) Column-level checkpoint
        new_ckpt = checkpoint_path(col)
        df.to_csv(new_ckpt, index=False)
        for f in glob.glob(checkpoint_glob()):
            if f != new_ckpt:
                try:
                    os.remove(f)
                except OSError:
                    pass

        print("  [CHECKPOINT] column saved")
        time.sleep(0.8)

    df.to_csv(OUTPUT_FILE, index=False)
    print(f"[CHECKPOINT] Processing complete. Saved to {OUTPUT_FILE}")
    print(
        f"[CHECKPOINT] Token usage: prompt={token_tracker.prompt_tokens}, "
        f"cached_prompt={token_tracker.cached_prompt_tokens}, "
        f"completion={token_tracker.completion_tokens}"
    )
    print(f"[CHECKPOINT] Estimated cumulative {LLM_PROVIDER} cost: ${token_tracker.estimated_cost_usd:.4f}")
    print(f"[CHECKPOINT] Exclusion summary: {excluded_counts}")

    for f in glob.glob(checkpoint_glob()):
        try:
            os.remove(f)
        except OSError:
            pass


if __name__ == "__main__":
    main()
