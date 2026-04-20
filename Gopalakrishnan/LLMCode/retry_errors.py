#!/usr/bin/env python3
import os
import time
import re
import json
import pandas as pd
import requests
from xml.etree import ElementTree as ET

from llm_provider_client import (
    assert_provider_ready,
    call_llm_json,
    ensure_parent_dir,
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
LLM_PROVIDER = normalize_provider()
MODEL_NAME = get_provider_model(LLM_PROVIDER)
REASONING_EFFORT = os.getenv("OPENAI_REASONING_EFFORT", "high")
USAGE_TRACKER_FILE = os.getenv(
    "OPENAI_COST_TRACKER_FILE",
    runtime_path("cost", f"{LLM_PROVIDER}_cost_tracker.json"),
)
ORIGINAL_FILE = os.getenv("LLM_ORIGINAL_FILE", "taxon_article_fulltext_matrix.csv")
FILLED_FILE = os.getenv("LLM_FILLED_FILE", with_provider_tag("taxon_article_matrix_filled.csv"))
OUTPUT_FILE = os.getenv("LLM_OUTPUT_FILE", FILLED_FILE)
TAXON_INFO_FILE = "taxa_ids_filtered.csv"

# Optional explicit retry list (kept from prior run); if empty, script auto-detects stale columns.
RETRY_IDS = [
    "PMC9724183",
    "PMID1341252",
    "PMC10619027",
    "PMC11073541",
    "PMID34088571",
    "PMC11868406",
    "PMC11959079",
    "PMC10245822",
    "PMC7890436",
]


def _split_retry_tokens(text):
    parts = re.split(r"[\s,]+", text.strip())
    return [part for part in parts if part]


def load_custom_retry_ids():
    env_value = os.getenv("LLM_RETRY_IDS", "").strip()
    file_value = os.getenv("LLM_RETRY_IDS_FILE", "").strip()

    ids = []
    if env_value:
        ids.extend(_split_retry_tokens(env_value))

    if file_value:
        try:
            with open(file_value, "r", encoding="utf-8") as fh:
                ids.extend(_split_retry_tokens(fh.read()))
        except OSError as exc:
            raise RuntimeError(f"Could not read LLM_RETRY_IDS_FILE={file_value}: {exc}") from exc

    if not ids:
        return None

    deduped = []
    seen = set()
    for doc_id in ids:
        if doc_id not in seen:
            seen.add(doc_id)
            deduped.append(doc_id)
    return deduped

try:
    assert_provider_ready(LLM_PROVIDER)
except Exception as exc:
    print(f"Error: {exc}")
    raise SystemExit(1)
ncbi = NCBIClient("retry_errors")
MAX_TAXA_PER_REQUEST = max(int(os.getenv("LLM_MAX_TAXA_PER_REQUEST", "8")), 1)

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


def save_filled_outputs(df):
    df.to_csv(OUTPUT_FILE, index=False)
    if os.path.abspath(FILLED_FILE) != os.path.abspath(OUTPUT_FILE):
        df.to_csv(FILLED_FILE, index=False)


class TokenTracker:
    def __init__(self):
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
taxon_metadata = load_taxon_metadata(TAXON_INFO_FILE)


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
    out = {
        "valid": False,
        "is_review": False,
        "is_meta": False,
        "is_animal_only": False,
        "reason": "",
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

        out["is_review"] = bool(REVIEW_PUB_TYPES & pub_types) or (
            re.search(r"\breview\b", title, flags=re.IGNORECASE) is not None
        )
        out["is_meta"] = bool(META_PUB_TYPES & pub_types) or (
            re.search(r"\bmeta[- ]analysis\b|\bsystematic review\b", title, flags=re.IGNORECASE) is not None
        )

        has_humans_mesh = "humans" in mesh_terms
        has_animals_mesh = any(
            (m in ANIMAL_MESH_MARKERS) or any(k in m for k in ANIMAL_MESH_MARKERS)
            for m in mesh_terms
        )
        has_animal_text, has_human_text = _text_flags(title, abstract)
        out["is_animal_only"] = (
            (has_animals_mesh and not has_humans_mesh) or
            (has_animal_text and not has_human_text and not has_humans_mesh)
        )

        if out["is_review"]:
            out["reason"] = "review"
        elif out["is_meta"]:
            out["reason"] = "meta_analysis"
        elif out["is_animal_only"]:
            out["reason"] = "animal_only"
        else:
            out["reason"] = "ok"
    except Exception:
        out["reason"] = "metadata_fetch_error"
    return out


def screen_article(doc_id):
    if doc_id.upper().startswith("PMID"):
        pmid = doc_id.replace("PMID", "")
        meta = fetch_pubmed_metadata_by_pmid(pmid)
        return (not meta["valid"]) or meta["is_review"] or meta["is_meta"] or meta["is_animal_only"], meta["reason"]

    pmid = pmcid_to_pmid(doc_id)
    if pmid:
        meta = fetch_pubmed_metadata_by_pmid(pmid)
        return (not meta["valid"]) or meta["is_review"] or meta["is_meta"] or meta["is_animal_only"], meta["reason"]

    # fallback title-only for PMCID without PMID mapping
    try:
        clean_id = doc_id.replace("PMC", "")
        resp = ncbi.post(
            f"{NCBI_DB}/efetch.fcgi",
            data=ncbi.default_params({"db": "pmc", "id": clean_id, "retmode": "xml"}),
            timeout=30,
        )
        resp.raise_for_status()
        root = ET.fromstring(resp.content)
        title = find_first_text_by_localname(root, "article-title").lower()
        if not title:
            return True, "invalid"
        if "review" in title:
            return True, "review"
        if "meta-analysis" in title or "meta analysis" in title:
            return True, "meta_analysis"
        return False, "ok"
    except Exception:
        return True, "metadata_fetch_error"


def get_paper_text_robust(doc_id):
    is_pmc = doc_id.upper().startswith("PMC")
    clean_id = doc_id.replace("PMID", "").replace("PMC", "")
    db = "pmc" if is_pmc else "pubmed"

    for attempt in range(1, 4):
        try:
            resp = ncbi.post(
                f"{NCBI_DB}/efetch.fcgi",
                data=ncbi.default_params({"db": db, "id": clean_id, "retmode": "xml"}),
                timeout=45,
            )
            resp.raise_for_status()
            root = ET.fromstring(resp.content)

            if is_pmc:
                text = "".join(root.itertext())
            else:
                art = root.find(".//PubmedArticle")
                if art is None:
                    return None
                title = art.findtext(".//ArticleTitle") or ""
                abs_texts = art.findall(".//Abstract/AbstractText")
                abstract = " ".join([t.text or "" for t in abs_texts])
                text = f"TITLE: {title}\nABSTRACT: {abstract}"

            text = (text or "").strip()
            if len(text) > 50:
                return text
        except Exception as e:
            print(f"  [Attempt {attempt}/3] Fetch failed for {doc_id}: {e}")
            time.sleep(retry_sleep_seconds(e, attempt, base = 2))

    return None


def parse_p_numeric(raw):
    s = str(raw or "").strip().lower().replace("p", "")
    s = re.sub(r"[=\s]", "", s)
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
    m = re.search(r"\d+", str(raw).strip())
    if not m:
        return None
    n = int(m.group(0))
    return n if n > 0 else None


def build_cell_value(res_obj):
    if not isinstance(res_obj, dict):
        return ""
    d = str(res_obj.get("dir", "")).strip()
    if d not in {"+", "-"}:
        return ""

    marker = str(res_obj.get("marker", "")).strip().upper()
    if marker not in {"A", "B", "C"}:
        return ""

    p_val = parse_p_numeric(res_obj.get("p"))
    n_obs = parse_n_obs(res_obj.get("n_obs"))
    eff = res_obj.get("eff")

    n_part = f" {{n={n_obs}}}" if n_obs is not None else ""
    eff_part = f" [{eff}]" if isinstance(eff, str) and eff.strip() else ""
    return f"{d}{p_val} {marker}{n_part}{eff_part}".strip()


def analyze_paper_batch_robust(paper_id, paper_text, taxa_list):
    if not paper_text:
        return {}, None

    system_prompt = (
        "You extract statistically significant taxon-associated outcomes from biomedical paper text. "
        "Return strict JSON only."
    )
    target_payloads = build_target_payloads(taxa_list, taxon_metadata)
    user_prompt = f"""
Paper ID: {paper_id}
Target taxa: {json.dumps(target_payloads, ensure_ascii=True)}
Contexts:
{CONTEXT_PROMPT}

Return only taxa with significant findings. Output JSON object:
{{"Taxon": {{"p": "0.05", "dir": "+/-/?", "marker": "A/B/C", "eff": "string or null", "n_obs": 100}}}}
Rules:
- Direction must encode patient benefit or harm, not the raw coefficient sign.
- Use '+' when the reported association benefits patients and '-' when it harms patients. If benefit vs harm cannot be determined from text, return dir='?'.
- For response-type outcomes, '+' means associated with response, responders, clinical benefit, disease control, CR/PR, or other favorable response outcomes; '-' means associated with non-response, resistance, progression, or other unfavorable response outcomes.
- For survival or time-to-event outcomes, '+' means longer survival or lower hazard/risk of progression, recurrence, or death; '-' means shorter survival or higher hazard/risk of progression, recurrence, or death.
- Apply the survival rule to Cox models, Fine-Gray models, cause-specific hazard models, and related survival models. HR/subdistribution HR/cause-specific HR < 1 implies '+', HR/subdistribution HR/cause-specific HR > 1 implies '-'. A negative coefficient in a hazard model implies '+', and a positive coefficient implies '-'.
- For marker C, still use '+' for patient-beneficial associations and '-' for patient-harmful associations. If the reported outcome is response-like, apply the response rule; if it is survival/time-to-event, apply the survival rule.
- If significant but p-value missing, return p="0.05".
- Do not infer from figures/plots; use text only.
- For a species target, OTU/ASV/SGB/strain evidence will count when it ties to that exact species.
- For any target above species, do not treat species, OTUs, ASVs, SGBs, strains, or any finer descendants as evidence for the target.
- Use the exact provided target taxon string as the JSON key.

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
            if isinstance(parsed, dict):
                return parsed, usage_summary
            return {}, usage_summary
        except Exception as e:
            print(f"  [Attempt {attempt}/3] {LLM_PROVIDER} error on {paper_id}: {e}")
            time.sleep(retry_sleep_seconds(e, attempt))

    return {}, None


def main():
    print(f"[CHECKPOINT] Loading {FILLED_FILE}...")
    df_filled = pd.read_csv(FILLED_FILE).fillna("")

    print(f"[CHECKPOINT] Loading {ORIGINAL_FILE}...")
    df_orig = pd.read_csv(ORIGINAL_FILE).fillna("")

    taxon_col = df_orig.columns[0]

    # Auto-detect unresolved columns from FILLED file only.
    # This avoids re-scanning all previously processed articles.
    auto_retry = []
    for c in df_filled.columns[1:]:
        if not (c.startswith("PMC") or c.startswith("PMID")):
            continue
        m = df_filled[c].astype(str).str.strip().isin(["1", "1.0"])
        if m.any():
            auto_retry.append(c)

    custom_retry_ids = load_custom_retry_ids()
    if custom_retry_ids is not None:
        retry_ids = custom_retry_ids
        print(f"[CHECKPOINT] Starting retry over {len(retry_ids)} custom IDs")
    else:
        retry_ids = list(dict.fromkeys(RETRY_IDS + auto_retry))
        print(f"[CHECKPOINT] Starting retry over {len(retry_ids)} IDs")
    pricing = pricing_for_model(LLM_PROVIDER, MODEL_NAME)
    cached_rate = pricing["cached_input"]
    cached_rate_text = f"${cached_rate:.2f}/1M cached-input" if cached_rate is not None else "no cached-input rate"
    print(
        f"[CHECKPOINT] provider={LLM_PROVIDER}, model={MODEL_NAME}, reasoning_effort={REASONING_EFFORT}, "
        f"pricing=${pricing['input']:.2f}/1M input, {cached_rate_text}, ${pricing['output']:.2f}/1M output"
    )

    for i, doc_id in enumerate(retry_ids, start=1):
        print(f"\n[{i}/{len(retry_ids)}] Retrying {doc_id}")

        if doc_id not in df_orig.columns:
            print("  -> Not present in original matrix.")
            continue

        mask = df_orig[doc_id].astype(str).str.strip().isin(["1", "1.0"])
        target_indices = df_orig.index[mask].tolist()
        if not target_indices:
            print("  -> No target taxa in original matrix.")
            continue

        target_taxa = df_orig.loc[target_indices, taxon_col].tolist()
        print(f"  -> {len(target_taxa)} taxa to check")

        skip, reason = screen_article(doc_id)
        if skip:
            print(f"  -> Skipped before LLM ({reason})")
            continue

        text = get_paper_text_robust(doc_id)
        if not text:
            print("  -> Unable to retrieve article text")
            continue
        text_source = "PMC full text" if doc_id.upper().startswith("PMC") else "PubMed title+abstract"
        print(f"  -> Retrieved {len(text):,} characters from {text_source}")

        results = {}
        usage_parts = []
        total_batches = (len(target_taxa) + MAX_TAXA_PER_REQUEST - 1) // MAX_TAXA_PER_REQUEST
        for batch_num, taxa_batch in enumerate(chunked(target_taxa, MAX_TAXA_PER_REQUEST), start=1):
            batch_results, batch_usage = analyze_paper_batch_robust(doc_id, text, taxa_batch)
            if batch_results:
                results.update(batch_results)
            if batch_usage:
                usage_parts.append(batch_usage)
            if total_batches > 1:
                print(f"  -> Batch {batch_num}/{total_batches} processed {len(taxa_batch)} taxa")
        usage_summary = combine_usage_summaries(usage_parts)
        if not results:
            print(f"  -> {LLM_PROVIDER} returned empty result")
            continue
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

        updates = 0
        for row_idx, taxon_name in zip(target_indices, target_taxa):
            res_obj = lookup_result_for_target(results, taxon_name, taxon_metadata)
            if not is_result_specific_to_taxon(taxon_name, res_obj, taxon_metadata):
                res_obj = None
            val = build_cell_value(res_obj)
            if doc_id in df_filled.columns:
                existing = str(df_filled.at[row_idx, doc_id] or "").strip()
                if val or not existing:
                    df_filled.at[row_idx, doc_id] = val
                if val and val != existing:
                    updates += 1

        print(f"  -> Updated {updates} entries")
        save_filled_outputs(df_filled)
        print("  [CHECKPOINT] Saved partial output")
        time.sleep(0.8)

    print(f"\n[CHECKPOINT] Retry complete. Saved to {OUTPUT_FILE}")
    print(
        f"[CHECKPOINT] Token usage: prompt={token_tracker.prompt_tokens}, "
        f"cached_prompt={token_tracker.cached_prompt_tokens}, "
        f"completion={token_tracker.completion_tokens}"
    )
    print(f"[CHECKPOINT] Estimated cumulative {LLM_PROVIDER} cost: ${token_tracker.estimated_cost_usd:.4f}")


if __name__ == "__main__":
    main()
