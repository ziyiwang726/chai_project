#!/usr/bin/env python3
import csv
import re
from pathlib import Path
from statistics import NormalDist


SCRIPT_DIR = Path(__file__).resolve().parent
PROVIDER_INPUTS = {
    "claude": {
        "summary": SCRIPT_DIR.parent / "sideCov/claudeGenerated/weighted/taxon_ABCD_probit_stoufferZ_weighted_claude.csv",
        "matrix": SCRIPT_DIR / "taxon_article_matrix_filled_claude.csv",
    },
    "gemini": {
        "summary": SCRIPT_DIR.parent / "sideCov/geminiGenerated/weighted/taxon_ABCD_probit_stoufferZ_weighted_gemini.csv",
        "matrix": SCRIPT_DIR / "taxon_article_matrix_filled_gemini.csv",
    },
    "openai": {
        "summary": SCRIPT_DIR.parent / "sideCov/openAIGenerated/weighted/taxon_ABCD_probit_stoufferZ_weighted_openai.csv",
        "matrix": SCRIPT_DIR / "taxon_article_matrix_filled_openai.csv",
    },
}
OUTPUT_FILE = SCRIPT_DIR / "conflicting_weighted_taxa_article_contributions.csv"
CONTEXT_LABELS = {
    "A": "Melanoma response to immunotherapy",
    "B": "Melanoma survival",
    "C": "Other cancers + immunotherapy",
}

CELL_RE = re.compile(r"^\s*([+-])\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*([ABCD])\b")
CELL_RE_FALLBACK = re.compile(r"^\s*([+-])\s*([ABCD])\b")
N_OBS_RE_1 = re.compile(r"(?i)\{\s*n\s*=\s*([0-9]+)\s*\}")
N_OBS_RE_2 = re.compile(r"(?i)\bn\s*[=:]\s*([0-9]+)")
BRACKET_RE = re.compile(r"\[(.*)\]\s*$")
NORM = NormalDist()
EPS = 1e-15


def parse_float(value):
    text = str(value or "").strip()
    if not text or text.upper() == "NA":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def sign_of(value):
    if value is None:
        return 0
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def sign_label(sign_num):
    return {1: "positive", -1: "negative", 0: "zero_or_missing"}[sign_num]


def load_summary(path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for line_number, row in enumerate(reader, start=2):
            rows.append(
                {
                    "taxon": (row.get("taxon") or "").strip(),
                    "line_number": line_number,
                    "A": parse_float(row.get("A")),
                    "B": parse_float(row.get("B")),
                    "C": parse_float(row.get("C")),
                }
            )
    return rows


def load_matrix(path):
    matrix = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return matrix
        article_cols = [col for col in reader.fieldnames if col != "taxon"]
        for line_number, row in enumerate(reader, start=2):
            taxon = (row.get("taxon") or "").strip()
            matrix[taxon] = {
                "line_number": line_number,
                "articles": {col: (row.get(col) or "").strip() for col in article_cols},
            }
    return matrix


def parse_cell(raw_value):
    text = str(raw_value or "").strip()
    if not text:
        return None

    effect_match = BRACKET_RE.search(text)
    effect_text = effect_match.group(1) if effect_match else ""
    core = re.sub(r"\s*\[.*$", "", text)

    match = CELL_RE.match(core)
    sign_chr = None
    p_value = None
    group = None
    if match:
        sign_chr, p_text, group = match.groups()
        p_value = float(p_text)
    else:
        fallback = CELL_RE_FALLBACK.match(core)
        if fallback:
            sign_chr, group = fallback.groups()
            p_value = 0.05

    if sign_chr is None or group is None:
        return None

    if group == "D":
        group = "C"

    n_match = N_OBS_RE_1.search(text) or N_OBS_RE_2.search(text)
    n_obs = float(n_match.group(1)) if n_match else None
    weight = n_obs if n_obs and n_obs > 0 else 1.0
    p_clamped = min(max(p_value, EPS), 1 - EPS)
    sign_num = 1 if sign_chr == "+" else -1
    article_z = sign_num * NORM.inv_cdf(1 - p_clamped / 2)

    return {
        "raw_value": text,
        "effect_text": effect_text,
        "group": group,
        "article_direction": sign_chr,
        "article_p": p_value,
        "article_n": n_obs,
        "article_weight": weight,
        "article_z": article_z,
        "article_weighted_contribution": weight * article_z,
    }


def mixed_sign_taxa(summary_rows):
    selected = []
    for row in summary_rows:
        signs = [sign_of(row[group]) for group in ("A", "B", "C")]
        if any(sign > 0 for sign in signs) and any(sign < 0 for sign in signs):
            selected.append(row)
    return selected


def build_rows():
    output_rows = []

    for provider, paths in PROVIDER_INPUTS.items():
        summary_rows = mixed_sign_taxa(load_summary(paths["summary"]))
        matrix_rows = load_matrix(paths["matrix"])

        for summary in summary_rows:
            taxon = summary["taxon"]
            matrix_entry = matrix_rows.get(taxon)
            if matrix_entry is None:
                continue

            positive_groups = [group for group in ("A", "B", "C") if sign_of(summary[group]) > 0]
            negative_groups = [group for group in ("A", "B", "C") if sign_of(summary[group]) < 0]

            parsed_entries = []
            for article, raw_value in matrix_entry["articles"].items():
                parsed = parse_cell(raw_value)
                if parsed is None:
                    continue
                parsed["article"] = article
                parsed_entries.append(parsed)

            group_stats = {}
            for group in ("A", "B", "C"):
                group_entries = [entry for entry in parsed_entries if entry["group"] == group]
                group_entries.sort(
                    key=lambda entry: (-abs(entry["article_weighted_contribution"]), entry["article"])
                )
                positive_articles = sum(entry["article_direction"] == "+" for entry in group_entries)
                negative_articles = sum(entry["article_direction"] == "-" for entry in group_entries)
                group_stats[group] = {
                    "entries": group_entries,
                    "count": len(group_entries),
                    "positive_articles": positive_articles,
                    "negative_articles": negative_articles,
                    "has_mixed_signs": positive_articles > 0 and negative_articles > 0,
                }

            conflict_kind = (
                "cross_context_and_within_group"
                if any(stats["has_mixed_signs"] for stats in group_stats.values())
                else "cross_context_only"
            )

            for group in ("A", "B", "C"):
                stats = group_stats[group]
                if not stats["entries"]:
                    continue

                group_summary_z = summary[group]
                group_summary_sign = sign_of(group_summary_z)
                opposite_groups = [
                    other_group
                    for other_group in ("A", "B", "C")
                    if other_group != group and sign_of(summary[other_group]) == -group_summary_sign
                ]

                for entry in stats["entries"]:
                    article_sign = 1 if entry["article_direction"] == "+" else -1
                    output_rows.append(
                        {
                            "provider": provider,
                            "taxon": taxon,
                            "conflict_kind": conflict_kind,
                            "summary_row": summary["line_number"],
                            "matrix_row": matrix_entry["line_number"],
                            "summary_A": summary["A"],
                            "summary_B": summary["B"],
                            "summary_C": summary["C"],
                            "summary_positive_groups": ",".join(positive_groups),
                            "summary_negative_groups": ",".join(negative_groups),
                            "group": group,
                            "context_label": CONTEXT_LABELS[group],
                            "group_summary_z": group_summary_z,
                            "group_summary_sign": sign_label(group_summary_sign),
                            "opposite_direction_summary_groups": ",".join(opposite_groups),
                            "group_article_count": stats["count"],
                            "group_positive_articles": stats["positive_articles"],
                            "group_negative_articles": stats["negative_articles"],
                            "group_has_mixed_article_signs": int(stats["has_mixed_signs"]),
                            "article": entry["article"],
                            "raw_value": entry["raw_value"],
                            "effect_text": entry["effect_text"],
                            "article_direction": entry["article_direction"],
                            "article_p": entry["article_p"],
                            "article_n": entry["article_n"],
                            "article_weight": entry["article_weight"],
                            "article_z": entry["article_z"],
                            "article_weighted_contribution": entry["article_weighted_contribution"],
                            "article_opposes_group_summary": int(
                                group_summary_sign != 0 and article_sign != group_summary_sign
                            ),
                        }
                    )

    output_rows.sort(
        key=lambda row: (
            row["provider"],
            row["taxon"],
            row["group"],
            -abs(row["article_weighted_contribution"]),
            row["article"],
        )
    )
    return output_rows


def main():
    rows = build_rows()
    fieldnames = [
        "provider",
        "taxon",
        "conflict_kind",
        "summary_row",
        "matrix_row",
        "summary_A",
        "summary_B",
        "summary_C",
        "summary_positive_groups",
        "summary_negative_groups",
        "group",
        "context_label",
        "group_summary_z",
        "group_summary_sign",
        "opposite_direction_summary_groups",
        "group_article_count",
        "group_positive_articles",
        "group_negative_articles",
        "group_has_mixed_article_signs",
        "article",
        "raw_value",
        "effect_text",
        "article_direction",
        "article_p",
        "article_n",
        "article_weight",
        "article_z",
        "article_weighted_contribution",
        "article_opposes_group_summary",
    ]

    with OUTPUT_FILE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
