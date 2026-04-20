#!/usr/bin/env python3
import csv
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_FILES = {
    "openai": SCRIPT_DIR / "taxon_article_matrix_filled_openai.csv",
    "claude": SCRIPT_DIR / "taxon_article_matrix_filled_claude.csv",
    "gemini": SCRIPT_DIR / "taxon_article_matrix_filled_gemini.csv",
}
OUTPUT_FILE = SCRIPT_DIR / "article_disagreement_ranked.csv"

BRACKET_RE = re.compile(r"\s*\[[^\]]*\]\s*$")
SIGNED_P_RE = re.compile(r"^\s*([+-])\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)")


def load_matrix(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {}, []

    taxon_col = rows[0].keys().__iter__().__next__()
    article_cols = [col for col in rows[0].keys() if col != taxon_col]
    matrix = {}
    for row in rows:
        taxon = row[taxon_col]
        matrix[taxon] = {col: (row.get(col, "") or "").strip() for col in article_cols}
    return matrix, article_cols


def strip_brackets(text):
    return BRACKET_RE.sub("", str(text or "")).strip()


def extract_signed_p(text):
    core = strip_brackets(text)
    if not core:
        return ""
    match = SIGNED_P_RE.match(core)
    if not match:
        return core
    return f"{match.group(1)}{match.group(2)}"


def compare_cell(values):
    present = {provider: bool(raw) for provider, raw in values.items()}
    normalized = {provider: extract_signed_p(raw) for provider, raw in values.items()}

    nonempty_norm = [normalized[p] for p in normalized if normalized[p]]
    nonempty_signs = [val[0] for val in nonempty_norm if val and val[0] in "+-"]
    nonempty_pvals = [val[1:] for val in nonempty_norm if len(val) > 1]

    presence_diff = len(set(present.values())) > 1
    sign_diff = len(set(nonempty_signs)) > 1
    pvalue_diff = len(set(nonempty_pvals)) > 1

    if not any((presence_diff, sign_diff, pvalue_diff)):
        return None

    diff_parts = []
    if presence_diff:
        diff_parts.append("presence_mismatch")
    if sign_diff:
        diff_parts.append("sign_mismatch")
    if pvalue_diff:
        diff_parts.append("pvalue_mismatch")

    return {
        "presence_diff": int(presence_diff),
        "sign_diff": int(sign_diff),
        "pvalue_diff": int(pvalue_diff),
        "difference_type": ";".join(diff_parts),
        "normalized": normalized,
        "raw": values,
    }


def build_rows():
    provider_data = {}
    article_union = set()
    taxon_union = set()

    for provider, path in DEFAULT_FILES.items():
        matrix, article_cols = load_matrix(path)
        provider_data[provider] = matrix
        article_union.update(article_cols)
        taxon_union.update(matrix.keys())

    article_summaries = {}
    article_rows = {}

    for article in sorted(article_union):
        diff_rows = []
        for taxon in sorted(taxon_union):
            values = {
                provider: provider_data.get(provider, {}).get(taxon, {}).get(article, "")
                for provider in DEFAULT_FILES
            }
            if not any(values.values()):
                continue

            cell_diff = compare_cell(values)
            if cell_diff is None:
                continue

            diff_rows.append(
                {
                    "article": article,
                    "taxon": taxon,
                    "sign_diff": cell_diff["sign_diff"],
                    "pvalue_diff": cell_diff["pvalue_diff"],
                    "presence_diff": cell_diff["presence_diff"],
                    "difference_type": cell_diff["difference_type"],
                    "openai_value": cell_diff["normalized"]["openai"],
                    "claude_value": cell_diff["normalized"]["claude"],
                    "gemini_value": cell_diff["normalized"]["gemini"],
                    "openai_raw": cell_diff["raw"]["openai"],
                    "claude_raw": cell_diff["raw"]["claude"],
                    "gemini_raw": cell_diff["raw"]["gemini"],
                }
            )

        if not diff_rows:
            continue

        article_summary = {
            "informative_taxa": len(diff_rows),
            "sign_diff": sum(row["sign_diff"] for row in diff_rows),
            "pvalue_diff": sum(row["pvalue_diff"] for row in diff_rows),
            "presence_diff": sum(row["presence_diff"] for row in diff_rows),
        }
        article_summary["total_diff"] = (
            article_summary["sign_diff"] +
            article_summary["pvalue_diff"] +
            article_summary["presence_diff"]
        )
        article_summaries[article] = article_summary
        article_rows[article] = sorted(
            diff_rows,
            key=lambda row: (
                -row["sign_diff"],
                -row["pvalue_diff"],
                -row["presence_diff"],
                row["taxon"],
            ),
        )

    ranked_articles = sorted(
        article_summaries,
        key=lambda article: (
            -article_summaries[article]["sign_diff"],
            -article_summaries[article]["pvalue_diff"],
            -article_summaries[article]["presence_diff"],
            -article_summaries[article]["informative_taxa"],
            article,
        ),
    )

    output_rows = []
    for rank, article in enumerate(ranked_articles, start=1):
        summary = article_summaries[article]
        for row in article_rows[article]:
            output_rows.append(
                {
                    "article_rank": rank,
                    "article": article,
                    "informative_taxa": summary["informative_taxa"],
                    "sign_diff": summary["sign_diff"],
                    "pvalue_diff": summary["pvalue_diff"],
                    "presence_diff": summary["presence_diff"],
                    "total_diff": summary["total_diff"],
                    **row,
                }
            )
    return output_rows


def main():
    rows = build_rows()
    fieldnames = [
        "article_rank",
        "article",
        "informative_taxa",
        "sign_diff",
        "pvalue_diff",
        "presence_diff",
        "total_diff",
        "taxon",
        "difference_type",
        "openai_value",
        "claude_value",
        "gemini_value",
        "openai_raw",
        "claude_raw",
        "gemini_raw",
    ]

    with OUTPUT_FILE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
