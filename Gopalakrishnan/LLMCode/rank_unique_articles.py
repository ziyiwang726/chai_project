#!/usr/bin/env python3
"""
Build a ranked list of unique article IDs from taxa_ids_filtered.csv.

Usage:
  python rank_unique_articles.py --in taxa_ids_filtered.csv --out unique_articles_ranked.csv
"""

import argparse
from collections import Counter
import pandas as pd

ID_COLS_DEFAULT = [
    "melanoma_response_ids",
    "melanoma_survival_ids",
    "other_cancer_io_ids",
]


def parse_cell(cell):
    if pd.isna(cell):
        return []
    s = str(cell).strip()
    if not s:
        return []
    return [x.strip() for x in s.split(";") if x.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_csv", default="taxa_ids_filtered.csv")
    ap.add_argument("--out", dest="out_csv", default="unique_articles_ranked.csv")
    args = ap.parse_args()

    print(f"[CHECKPOINT] Reading {args.in_csv}")
    df = pd.read_csv(args.in_csv)

    id_cols = [c for c in ID_COLS_DEFAULT if c in df.columns]
    if not id_cols:
        id_cols = [c for c in df.columns if c.endswith("_ids")]
    if not id_cols:
        id_cols = list(df.columns[-3:])

    print(f"[CHECKPOINT] ID columns: {id_cols}")

    counter = Counter()
    for col in id_cols:
        for cell in df[col].tolist():
            ids = parse_cell(cell)
            counter.update(ids)

    ranked = sorted(counter.items(), key=lambda x: x[1], reverse=True)
    out = pd.DataFrame(ranked, columns=["id", "count"])
    out.to_csv(args.out_csv, index=False)

    print(f"[CHECKPOINT] Wrote {len(out)} unique IDs to {args.out_csv}")


if __name__ == "__main__":
    main()
