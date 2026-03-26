#!/usr/bin/env python3
import csv
import re
from collections import Counter
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILES = {
    "openai": SCRIPT_DIR / "taxon_article_matrix_filled_openai.csv",
    "claude": SCRIPT_DIR / "taxon_article_matrix_filled_claude.csv",
    "gemini": SCRIPT_DIR / "taxon_article_matrix_filled_gemini.csv",
}
OUTPUT_FILE = SCRIPT_DIR / "taxon_article_matrix_filled_combined.csv"
BRACKET_RE = re.compile(r"\s*\[[^\]]*\]\s*$")


def load_matrix(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {}, []

    taxon_col = next(iter(rows[0].keys()))
    article_cols = [col for col in rows[0].keys() if col != taxon_col]
    matrix = {}
    for row in rows:
        taxon = row[taxon_col]
        matrix[taxon] = {col: (row.get(col, "") or "").strip() for col in article_cols}
    return matrix, article_cols


def strip_brackets(value):
    return BRACKET_RE.sub("", str(value or "")).strip()


def representative_value(raw_values, core_value):
    agreeing = [raw for raw in raw_values if raw and strip_brackets(raw) == core_value]
    if not agreeing:
        return ""
    agreeing.sort(key=lambda raw: (-len(raw), raw))
    return agreeing[0]


def main():
    provider_data = {}
    article_union = set()
    taxon_union = set()

    for provider, path in INPUT_FILES.items():
        matrix, article_cols = load_matrix(path)
        provider_data[provider] = matrix
        article_union.update(article_cols)
        taxon_union.update(matrix.keys())

    fieldnames = ["taxon"] + sorted(article_union)
    rows_out = []
    kept_cells = 0

    for taxon in sorted(taxon_union):
        row = {"taxon": taxon}
        for article in fieldnames[1:]:
            raw_values = [
                provider_data.get(provider, {}).get(taxon, {}).get(article, "")
                for provider in ("openai", "claude", "gemini")
            ]
            core_values = [strip_brackets(raw) for raw in raw_values if strip_brackets(raw)]
            if not core_values:
                row[article] = ""
                continue

            counts = Counter(core_values)
            consensus_core = next((core for core, count in counts.items() if count >= 2), None)
            if consensus_core is None:
                row[article] = ""
                continue

            row[article] = representative_value(raw_values, consensus_core)
            if row[article]:
                kept_cells += 1

        rows_out.append(row)

    with OUTPUT_FILE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote {OUTPUT_FILE} with {len(rows_out)} taxa rows and {kept_cells} consensus cells")


if __name__ == "__main__":
    main()
