#!/usr/bin/env python3
import os
import re
import subprocess
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_THRESHOLD = float(os.getenv("LOW_ABUNDANCE_PHYLA_THRESHOLD", "0.10"))


def normalize_taxon_name(text):
    text = str(text or "")
    text = re.sub(r"^[a-z]__", "", text, flags=re.IGNORECASE)
    text = text.replace("[", " ").replace("]", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def taxonomy_csv_path():
    return Path(os.getenv("LOW_ABUNDANCE_PHYLA_TAXONOMY_CSV", SCRIPT_DIR / "d1Taxonomy.csv"))


def otu_table_csv_path():
    env_path = os.getenv("LOW_ABUNDANCE_PHYLA_OTU_CSV", "").strip()
    if env_path:
        return Path(env_path)

    cmd = [
        "Rscript",
        "-e",
        "cat(system.file('extdata','d1OTUtable.csv', package='CATMicrobiome'))",
    ]
    resolved = subprocess.check_output(cmd, text=True).strip()
    if not resolved:
        raise RuntimeError("Could not resolve CATMicrobiome d1OTUtable.csv via Rscript.")
    return Path(resolved)


def compute_phylum_abundance(threshold=DEFAULT_THRESHOLD):
    tax_path = taxonomy_csv_path()
    otu_path = otu_table_csv_path()

    taxonomy = pd.read_csv(tax_path, dtype=str)
    otu = pd.read_csv(otu_path, index_col=0)

    if "OTUname" not in taxonomy.columns or "phylum" not in taxonomy.columns:
        raise ValueError("d1Taxonomy.csv must contain OTUname and phylum columns.")

    merged = taxonomy[["OTUname", "phylum"]].copy()
    merged["phylum_clean"] = merged["phylum"].map(normalize_taxon_name)
    merged = merged.merge(otu, left_on="OTUname", right_index=True, how="inner")

    sample_cols = [c for c in otu.columns if c in merged.columns]
    if not sample_cols:
        raise ValueError("No OTU abundance columns were found after joining taxonomy to OTU table.")

    phylum_totals = merged.groupby("phylum_clean", dropna=False)[sample_cols].sum(numeric_only=True).sum(axis=1)
    phylum_totals = phylum_totals[phylum_totals.index.notna() & (phylum_totals.index != "")]

    total_abundance = float(phylum_totals.sum())
    if total_abundance <= 0:
        raise ValueError("Total phylum abundance is zero; cannot compute relative abundance.")

    out = (
        pd.DataFrame(
            {
                "phylum": phylum_totals.index,
                "total_abundance": phylum_totals.values,
                "relative_abundance": phylum_totals.values / total_abundance,
            }
        )
        .sort_values(["relative_abundance", "phylum"], ascending=[False, True])
        .reset_index(drop=True)
    )
    out["is_low_abundance"] = out["relative_abundance"] < float(threshold)
    return out


def compute_low_abundance_phyla(threshold=DEFAULT_THRESHOLD):
    abundance = compute_phylum_abundance(threshold=threshold)
    return {
        normalize_taxon_name(name)
        for name in abundance.loc[abundance["is_low_abundance"], "phylum"].tolist()
        if normalize_taxon_name(name)
    }


if __name__ == "__main__":
    table = compute_phylum_abundance()
    print(table.to_string(index=False))
