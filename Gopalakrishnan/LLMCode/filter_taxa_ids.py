#!/usr/bin/env python3
"""
filter_taxa_ids.py

Post-process a taxa IDs CSV:

- Read an input CSV (default: .runtime/results_melanoma_taxa_updated.csv).
- Detect `*_ids` columns containing ID lists (semicolon-separated),
  with IDs like:
      PMC1234567
      PMID12345678
- For each ID:

  1. If it's a PMCID (starts with 'PMC'):
       - Check validity via NCBI E-utilities (db=pmc, efetch).
       - If invalid (no article or fetch error), DROP it.

  2. For ALL IDs (PMC and PMID):
       - Fetch PubMed/PMC metadata.
       - Drop retracted papers.
       - Drop papers whose PubMed citation is not indexed for MEDLINE.
       - Drop review/meta-analysis/pure-animal studies.

- Keep all rows, just clean up the ID lists.
- Write a new CSV (default: taxa_ids_filtered.csv).

Usage:
    python filter_taxa_ids.py --in .runtime/results_melanoma_taxa_updated.csv --out taxa_ids_filtered.csv

Environment variables:
    NCBI_EMAIL   (optional but recommended)
    NCBI_API_KEY (optional; increases rate limits)
"""

import os
import re
import argparse
import pandas as pd
from xml.etree import ElementTree as ET

from ncbi_request_utils import BASE_EUTILS, PMC_IDCONV, NCBIClient

client = NCBIClient("taxa_ids_filter")
DEFAULT_INPUT_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    ".runtime",
    "results_melanoma_taxa_updated.csv",
)


def eutil_params(extra=None):
    return client.default_params(extra)

REVIEW_PUB_TYPES = {
    "Review",
    "Systematic Review",
    "Scoping Review",
}
META_PUB_TYPES = {
    "Meta-Analysis",
    "Systematic Review",
}
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

# Deepak original paper (exclude as non-auxiliary input for this workflow)
EXCLUDED_IDS = {"PMC5827966"}

# -------------------- ID parsing helpers --------------------

def parse_id(raw_id):
    """
    Parse an ID string into (kind, canonical_id).

    Returns:
      kind: 'pmcid', 'pmid', or None
      canonical_id:
        - for pmcid: 'PMC' + digits
        - for pmid: digits only
        - for None: original string
    """
    if not raw_id:
        return None, raw_id

    s = raw_id.strip()
    if not s:
        return None, raw_id

    u = s.upper()

    # PMCID: PMC1234567 or variations
    if u.startswith("PMC"):
        digits = re.sub(r"\D", "", u[3:])
        if digits:
            return "pmcid", f"PMC{digits}"
        else:
            return "pmcid", u  # weird, but keep as-is

    # PMID: PMID123456 or bare digits
    if u.startswith("PMID"):
        digits = re.sub(r"\D", "", s[4:])
        if digits:
            return "pmid", digits
        else:
            return "pmid", s

    if s.isdigit():
        return "pmid", s

    # Unknown type; treat as PMID-esque, but mark type None
    return None, s

# -------------------- Fetch title / validity --------------------

pmcid_info = {}   # canonical_pmcid -> dict metadata
pmid_info  = {}   # canonical_pmid  -> dict metadata
idconv_map = {}   # canonical_pmcid -> pmid or None


def _base_info():
    return {
        "valid": False,
        "title": "",
        "journal": "",
        "medline_status": "",
        "is_medline": False,
        "is_retracted": False,
        "is_review": False,
        "is_meta": False,
        "is_animal_only": False,
        "exclude_reason": "",
    }


def _text_flags(title, abstract=""):
    txt = f" {title or ''} {abstract or ''} ".lower()
    has_animal = any(k in txt for k in ANIMAL_TEXT_MARKERS)
    has_human = any(k in txt for k in HUMAN_TEXT_MARKERS)
    return has_animal, has_human


def classify_pubmed_article(article):
    """
    Infer retraction / MEDLINE / review / meta / pure-animal status from
    PubMed XML article node.
    """
    medline_citation = article.find("MedlineCitation")
    medline_status = (medline_citation.attrib.get("Status") or "").strip() if medline_citation is not None else ""
    title = (article.findtext(".//ArticleTitle") or "").strip()
    journal = (
        (article.findtext(".//MedlineJournalInfo/MedlineTA") or "").strip()
        or (article.findtext(".//Journal/Title") or "").strip()
    )
    abs_texts = article.findall(".//Abstract/AbstractText")
    abstract = " ".join([(t.text or "") for t in abs_texts]).strip()

    pub_types = {
        (pt.text or "").strip()
        for pt in article.findall(".//PublicationType")
        if (pt.text or "").strip()
    }
    retraction_ref_types = {
        (cc.attrib.get("RefType") or "").strip().lower()
        for cc in article.findall(".//CommentsCorrections")
        if (cc.attrib.get("RefType") or "").strip()
    }
    mesh_terms = {
        (mh.findtext("DescriptorName") or "").strip().lower()
        for mh in article.findall(".//MeshHeading")
    }

    is_review = (
        bool(REVIEW_PUB_TYPES & pub_types) or
        re.search(r"\breview\b", title, flags=re.IGNORECASE) is not None
    )
    is_meta = (
        bool(META_PUB_TYPES & pub_types) or
        re.search(r"\bmeta[- ]analysis\b|\bsystematic review\b", title, flags=re.IGNORECASE) is not None
    )
    is_retracted = (
        "Retracted Publication" in pub_types
        or "Retraction of Publication" in pub_types
        or "retractionin" in retraction_ref_types
    )

    has_humans_mesh = "humans" in mesh_terms
    has_animals_mesh = any(
        (m in ANIMAL_MESH_MARKERS) or any(k in m for k in ANIMAL_MESH_MARKERS)
        for m in mesh_terms
    )
    has_animal_text, has_human_text = _text_flags(title, abstract)

    # Pure animal = clear animal signal with no human signal
    is_animal_only = (
        (has_animals_mesh and not has_humans_mesh) or
        (has_animal_text and not has_human_text and not has_humans_mesh)
    )

    return {
        "title": title,
        "journal": journal,
        "medline_status": medline_status,
        "is_medline": medline_status.upper() == "MEDLINE",
        "is_retracted": is_retracted,
        "is_review": is_review,
        "is_meta": is_meta,
        "is_animal_only": is_animal_only,
    }


def pmcid_to_pmid(pmcid):
    if pmcid in idconv_map:
        return idconv_map[pmcid]
    try:
        r = client.get(
            PMC_IDCONV,
            params={"ids": pmcid, "format": "json"},
            timeout=20,
        )
        js = r.json()
        recs = js.get("records", [])
        if recs:
            pmid = str(recs[0].get("pmid") or "").strip()
            idconv_map[pmcid] = pmid if pmid else None
            return idconv_map[pmcid]
    except Exception:
        pass
    idconv_map[pmcid] = None
    return None

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

def fetch_pmcid_info(pmcid):
    """
    Fetch article title for a PMCID via NCBI (db=pmc, efetch).
    Returns dict with title + exclusion flags.
    """
    if pmcid in pmcid_info:
        return pmcid_info[pmcid]

    info = _base_info()

    try:
        r = client.get(
            f"{BASE_EUTILS}/efetch.fcgi",
            params=eutil_params({"db": "pmc", "id": pmcid, "retmode": "xml"}),
            timeout=30,
        )
        root = ET.fromstring(r.content)

        # Try to find <article-title> with namespace-agnostic lookup
        title = find_first_text_by_localname(root, "article-title")

        if title:
            info["valid"] = True
            info["title"] = title
            info["is_review"] = is_review_title(title)
            info["is_meta"] = re.search(
                r"\bmeta[- ]analysis\b|\bsystematic review\b",
                title,
                flags=re.IGNORECASE,
            ) is not None
        else:
            info["valid"] = False
            info["exclude_reason"] = "invalid_id"
    except Exception:
        info["valid"] = False
        info["exclude_reason"] = "invalid_id"

    # Try to classify via linked PMID metadata for better review/meta/animal detection
    if info["valid"]:
        linked_pmid = pmcid_to_pmid(pmcid)
        if linked_pmid:
            pmid_meta = fetch_pmid_info(linked_pmid)
            if pmid_meta.get("valid"):
                info["journal"] = pmid_meta.get("journal", "")
                info["medline_status"] = pmid_meta.get("medline_status", "")
                info["is_medline"] = pmid_meta.get("is_medline", False)
                info["is_retracted"] = pmid_meta.get("is_retracted", False)
                info["is_review"] = info["is_review"] or pmid_meta.get("is_review", False)
                info["is_meta"] = info["is_meta"] or pmid_meta.get("is_meta", False)
                info["is_animal_only"] = pmid_meta.get("is_animal_only", False)
                if not info["title"]:
                    info["title"] = pmid_meta.get("title", "")
        else:
            # PMC-only content is outside MEDLINE indexing and should not seed article columns.
            info["is_medline"] = False

        if info["is_retracted"]:
            info["exclude_reason"] = "retracted"
        elif not info["is_medline"]:
            info["exclude_reason"] = "non_medline_journal"
        elif info["is_review"]:
            info["exclude_reason"] = "review"
        elif info["is_meta"]:
            info["exclude_reason"] = "meta_analysis"
        elif info["is_animal_only"]:
            info["exclude_reason"] = "animal_only"

    pmcid_info[pmcid] = info
    return info

def fetch_pmid_info(pmid):
    """
    Fetch article title for a PMID via NCBI (db=pubmed, efetch).
    Returns dict with title + exclusion flags.
    """
    if pmid in pmid_info:
        return pmid_info[pmid]

    info = _base_info()

    try:
        r = client.get(
            f"{BASE_EUTILS}/efetch.fcgi",
            params=eutil_params({"db": "pubmed", "id": pmid, "retmode": "xml"}),
            timeout=30,
        )
        root = ET.fromstring(r.content)

        art = root.find(".//PubmedArticle")
        if art is not None:
            cls = classify_pubmed_article(art)
            title = cls["title"]
            if title:
                info["valid"] = True
                info["title"] = title
                info["journal"] = cls["journal"]
                info["medline_status"] = cls["medline_status"]
                info["is_medline"] = cls["is_medline"]
                info["is_retracted"] = cls["is_retracted"]
                info["is_review"] = cls["is_review"]
                info["is_meta"] = cls["is_meta"]
                info["is_animal_only"] = cls["is_animal_only"]
                if info["is_retracted"]:
                    info["exclude_reason"] = "retracted"
                elif not info["is_medline"]:
                    info["exclude_reason"] = "non_medline_journal"
                elif info["is_review"]:
                    info["exclude_reason"] = "review"
                elif info["is_meta"]:
                    info["exclude_reason"] = "meta_analysis"
                elif info["is_animal_only"]:
                    info["exclude_reason"] = "animal_only"
    except Exception:
        info["valid"] = False
        info["exclude_reason"] = "invalid_id"

    pmid_info[pmid] = info
    return info

def is_review_title(title):
    """
    Return True if 'review' appears as a whole word in the title (case-insensitive).
    """
    if not title:
        return False
    return re.search(r"\breview\b", title, flags=re.IGNORECASE) is not None

# -------------------- Main processing --------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_csv", default=DEFAULT_INPUT_CSV,
                    help="Input CSV file (default: .runtime/results_melanoma_taxa_updated.csv)")
    ap.add_argument("--out", dest="out_csv", default="taxa_ids_filtered.csv",
                    help="Output CSV file (filtered)")
    ap.add_argument("--dryrun", action="store_true",
                    help="Do not write output, just print summary")
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)

    id_cols = [c for c in df.columns if c.endswith("_ids")]
    if not id_cols:
        id_cols = list(df.columns[-3:])
    print(f"ID columns: {id_cols}")

    # Collect unique IDs across the detected ID columns
    unique_pmcids = set()
    unique_pmids  = set()

    for col in id_cols:
        for val in df[col].astype(str).tolist():
            if pd.isna(val) or val is None:
                continue
            s = str(val).strip()
            if not s:
                continue
            parts = [p.strip() for p in s.split(";") if p.strip()]
            for p in parts:
                kind, cid = parse_id(p)
                if kind == "pmcid":
                    unique_pmcids.add(cid)
                elif kind == "pmid":
                    unique_pmids.add(cid)

    print(f"Unique PMCID count: {len(unique_pmcids)}")
    print(f"Unique PMID  count: {len(unique_pmids)}")

    # Pre-fetch info for all unique PMCs and PMIDs
    print("Fetching PMCID info...")
    for pmcid in sorted(unique_pmcids):
        info = fetch_pmcid_info(pmcid)
        if not info["valid"]:
            print(f"  PMCID invalid or no title: {pmcid}")
        elif info["exclude_reason"]:
            print(f"  PMCID excluded ({info['exclude_reason']}): {pmcid}")
    print("Fetching PMID info...")
    for pmid in sorted(unique_pmids):
        info = fetch_pmid_info(pmid)
        if not info["valid"]:
            print(f"  PMID invalid or no title: {pmid}")
        elif info["exclude_reason"]:
            print(f"  PMID excluded ({info['exclude_reason']}): PMID{pmid}")

    # Now clean each ID column
    dropped_invalid = 0
    dropped_excluded = 0
    dropped_retracted = 0
    dropped_non_medline = 0
    dropped_review = 0
    dropped_meta = 0
    dropped_animal = 0

    for col in id_cols:
        print(f"Processing column: {col}")
        new_vals = []
        for val in df[col].astype(str).tolist():
            if pd.isna(val) or val is None:
                new_vals.append(val)
                continue
            s = str(val).strip()
            if not s:
                new_vals.append("")
                continue

            parts = [p.strip() for p in s.split(";") if p.strip()]
            kept_ids = []
            for p in parts:
                kind, cid = parse_id(p)

                # Unknown ID type: keep as-is (we don't know how to validate / fetch)
                if kind is None:
                    kept_ids.append(p)
                    continue

                if kind == "pmcid":
                    if cid in EXCLUDED_IDS:
                        dropped_excluded += 1
                        continue
                    info = fetch_pmcid_info(cid)
                    # Step 1: drop invalid PMC codes
                    if not info["valid"]:
                        dropped_invalid += 1
                        # skip
                        continue
                else:  # "pmid"
                    info = fetch_pmid_info(cid)
                    if not info["valid"]:
                        dropped_invalid += 1
                        # treat as invalid -> drop
                        continue

                # Step 2: drop review/meta-analysis/pure-animal articles
                if info.get("is_retracted"):
                    dropped_retracted += 1
                    continue
                if not info.get("is_medline", False):
                    dropped_non_medline += 1
                    continue
                if info.get("is_review"):
                    dropped_review += 1
                    continue
                if info.get("is_meta"):
                    dropped_meta += 1
                    continue
                if info.get("is_animal_only"):
                    dropped_animal += 1
                    continue

                # Keep this ID; preserve canonical formatting
                if kind == "pmcid":
                    kept_ids.append(cid)
                else:  # pmid
                    kept_ids.append(f"PMID{cid}")

            new_vals.append("; ".join(kept_ids))

        df[col] = new_vals

    if args.dryrun:
        print("Dry run enabled; not writing output CSV.")
    else:
        df.to_csv(args.out_csv, index=False)
        print(f"Filtered table written to: {args.out_csv}")
        print("Drop summary:")
        print(f"  invalid: {dropped_invalid}")
        print(f"  explicitly excluded IDs: {dropped_excluded}")
        print(f"  retracted: {dropped_retracted}")
        print(f"  non-MEDLINE journal/citation: {dropped_non_medline}")
        print(f"  review: {dropped_review}")
        print(f"  meta-analysis/systematic review: {dropped_meta}")
        print(f"  pure animal: {dropped_animal}")

if __name__ == "__main__":
    main()
