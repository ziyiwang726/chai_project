#!/usr/bin/env python3
import json
import os
import re
import time

import pandas as pd

from llm_provider_client import (
    assert_provider_ready,
    call_llm_json,
    ensure_parent_dir,
    get_provider_model,
    normalize_provider,
    runtime_path,
    summarize_usage,
    with_provider_tag,
)
from taxon_result_utils import get_taxon_metadata, load_taxon_metadata, normalize_taxon_text


LLM_PROVIDER = normalize_provider()
MODEL_NAME = get_provider_model(LLM_PROVIDER)
REASONING_EFFORT = os.getenv("OPENAI_REASONING_EFFORT", "high")
FILLED_FILE = os.getenv("LLM_FILLED_FILE", with_provider_tag("taxon_article_matrix_filled.csv"))
OUTPUT_FILE = os.getenv("LLM_OUTPUT_FILE", FILLED_FILE)
TAXON_INFO_FILE = "taxa_ids_filtered.csv"
BATCH_SIZE = max(int(os.getenv("BRACKET_CLEANUP_BATCH_SIZE", "20")), 1)
USAGE_TRACKER_FILE = os.getenv(
    "BRACKET_CLEANUP_COST_TRACKER_FILE",
    runtime_path("cost", f"{LLM_PROVIDER}_bracket_cleanup_cost.json"),
)
LOG_FILE = os.getenv(
    "BRACKET_CLEANUP_LOG_FILE",
    runtime_path("logs", with_provider_tag("bracket_cleanup_log.csv")),
)

BRACKET_RE = re.compile(r"\[(.*)\]\s*$")
MOUSE_RE = re.compile(r"\b(mouse|mice|murine)\b", flags=re.IGNORECASE)
BINOMIAL_RE = re.compile(r"\b([A-Z][a-z]+)\s+([a-z][A-Za-z0-9._-]+)\b")


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


def extract_bracket_text(cell_value):
    match = BRACKET_RE.search(str(cell_value or ""))
    return (match.group(1) or "").strip() if match else ""


def fallback_decision(item):
    bracket_text = item["bracket_text"]
    if MOUSE_RE.search(bracket_text):
        return {"decision": "DELETE", "genus": "", "notes": "mouse-model mention"}

    genera = []
    for genus, species in BINOMIAL_RE.findall(bracket_text):
        if species.lower() in {"sp", "spp", "cf", "aff"}:
            continue
        norm = normalize_taxon_text(genus)
        if norm and norm not in genera:
            genera.append(norm)

    if len(genera) == 1:
        return {"decision": "DELETE", "genus": genera[0], "notes": "explicit species-level mention"}
    return {"decision": "KEEP", "genus": "", "notes": "fallback keep"}


def batched(items, size):
    for start in range(0, len(items), size):
        yield items[start:start + size]


def build_genus_row_map(df, taxon_col, metadata_store):
    genus_row_map = {}
    row_meta = {}
    for row_idx, taxon_name in df[taxon_col].items():
        meta = get_taxon_metadata(taxon_name, metadata_store)
        row_meta[row_idx] = meta
        if meta["taxon_rank"] != "genus":
            continue
        for candidate in (taxon_name, meta.get("taxon", ""), meta.get("current_scientific_name", "")):
            norm = normalize_taxon_text(candidate)
            if norm and norm not in genus_row_map:
                genus_row_map[norm] = row_idx
    return genus_row_map, row_meta


def build_candidates(df, taxon_col, row_meta):
    candidates = []
    for row_idx, row in df.iterrows():
        taxon_name = row[taxon_col]
        rank = row_meta[row_idx]["taxon_rank"]
        for col in df.columns:
            if col == taxon_col:
                continue
            value = str(row[col] or "").strip()
            if not value:
                continue
            bracket_text = extract_bracket_text(value)
            if not bracket_text:
                continue
            candidates.append(
                {
                    "id": f"{row_idx}::{col}",
                    "row_idx": int(row_idx),
                    "article_id": col,
                    "taxon_name": taxon_name,
                    "taxon_rank": rank,
                    "cell_value": value,
                    "bracket_text": bracket_text,
                }
            )
    return candidates


def classify_batch(batch):
    system_prompt = (
        "You are curating a microbiome evidence matrix. For each item, decide whether the bracketed note "
        "should stay on the current row, be deleted, or be moved to a genus row.\n"
        "Rules:\n"
        "1. If the bracket explicitly mentions mouse, mice, or murine evidence, return DELETE.\n"
        "2. If the bracket explicitly names a bacterial species or a lower taxon tied to a species "
        "(for example 'Prevotella copri ASVs'), return DELETE. This is species-level evidence and must "
        "not be attributed to a genus row.\n"
        "3. If no explicit bacterial species is named, use KEEP.\n"
        "Return strict JSON: {\"items\":[{\"id\":\"...\",\"decision\":\"KEEP|DELETE\","
        "\"genus\":\"\", \"notes\":\"\"}, ...]}"
    )
    user_prompt = json.dumps(
        {
            "task": "Classify bracketed evidence notes.",
            "items": [
                {
                    "id": item["id"],
                    "current_taxon": item["taxon_name"],
                    "current_rank": item["taxon_rank"],
                    "bracket_text": item["bracket_text"],
                }
                for item in batch
            ],
        },
        ensure_ascii=True,
    )

    parsed, usage, _ = call_llm_json(
        system_prompt,
        user_prompt,
        provider=LLM_PROVIDER,
        reasoning_effort=REASONING_EFFORT,
    )

    items = parsed.get("items") if isinstance(parsed, dict) else None
    mapped = {}
    if isinstance(items, list):
        for item in items:
            if not isinstance(item, dict):
                continue
            item_id = str(item.get("id", "")).strip()
            if not item_id:
                continue
            mapped[item_id] = {
                "decision": str(item.get("decision", "KEEP")).strip().upper() or "KEEP",
                "genus": normalize_taxon_text(item.get("genus", "")),
                "notes": str(item.get("notes", "") or "").strip(),
            }
    return mapped, usage


def main():
    try:
        assert_provider_ready(LLM_PROVIDER)
    except Exception as exc:
        print(f"Error: {exc}")
        raise SystemExit(1)

    print(f"[CHECKPOINT] Loading {FILLED_FILE}")
    df = pd.read_csv(FILLED_FILE, dtype=str).fillna("")
    taxon_col = df.columns[0]
    metadata_store = load_taxon_metadata(TAXON_INFO_FILE)
    genus_row_map, row_meta = build_genus_row_map(df, taxon_col, metadata_store)
    candidates = build_candidates(df, taxon_col, row_meta)
    token_tracker = TokenTracker()

    print(
        f"[CHECKPOINT] provider={LLM_PROVIDER}, model={MODEL_NAME}, "
        f"candidate_bracket_cells={len(candidates)}, batch_size={BATCH_SIZE}"
    )

    log_rows = []
    deleted = 0
    moved = 0
    kept = 0
    conflicts = 0

    for batch_idx, batch in enumerate(batched(candidates, BATCH_SIZE), start=1):
        print(f"[{batch_idx}] Cleaning {len(batch)} bracketed cells")
        try:
            decisions, usage = classify_batch(batch)
        except Exception as exc:
            print(f"  -> LLM batch failed: {exc}. Falling back to regex-based decisions.")
            decisions = {}
            usage = None

        if usage is not None:
            usage_summary = token_tracker.update(usage)
            cost_bits = [
                f"prompt={usage_summary['prompt_tokens']:,}",
                f"completion={usage_summary['completion_tokens']:,}",
                f"batch=${usage_summary['estimated_cost_usd']:.4f}",
                f"cumulative=${token_tracker.estimated_cost_usd:.4f}",
            ]
            if usage_summary["cached_prompt_tokens"]:
                cost_bits.insert(1, f"cached={usage_summary['cached_prompt_tokens']:,}")
            print(f"  -> Cost: {', '.join(cost_bits)}")

        for item in batch:
            row_idx = item["row_idx"]
            col = item["article_id"]
            current_value = str(df.at[row_idx, col]).strip()
            if not current_value:
                continue

            decision = decisions.get(item["id"]) or fallback_decision(item)
            action = decision["decision"]
            genus_norm = normalize_taxon_text(decision.get("genus", ""))
            outcome = "kept"

            if action == "DELETE":
                df.at[row_idx, col] = ""
                deleted += 1
                outcome = "deleted"
            elif action == "MOVE_TO_GENUS":
                df.at[row_idx, col] = ""
                deleted += 1
                outcome = "deleted_species_level"
            else:
                kept += 1

            log_rows.append(
                {
                    "provider": LLM_PROVIDER,
                    "article_id": col,
                    "source_taxon": item["taxon_name"],
                    "source_rank": item["taxon_rank"],
                    "decision": action,
                    "genus": genus_norm,
                    "outcome": outcome,
                    "notes": decision.get("notes", ""),
                    "cell_value": current_value,
                    "bracket_text": item["bracket_text"],
                }
            )

        df.to_csv(OUTPUT_FILE, index=False)
        time.sleep(0.3)

    pd.DataFrame(log_rows).to_csv(ensure_parent_dir(LOG_FILE), index=False)
    print(f"[CHECKPOINT] Saved cleaned matrix to {OUTPUT_FILE}")
    print(f"[CHECKPOINT] Saved audit log to {LOG_FILE}")
    print(
        f"[CHECKPOINT] Summary: kept={kept}, moved={moved}, deleted={deleted}, "
        f"target_conflicts={conflicts}, cumulative_cost=${token_tracker.estimated_cost_usd:.4f}"
    )


if __name__ == "__main__":
    main()
