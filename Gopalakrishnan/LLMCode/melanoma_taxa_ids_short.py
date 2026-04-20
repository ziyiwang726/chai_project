#!/usr/bin/env python3
"""
melanoma_taxa_updated.py

For each bacterial taxon in an input CSV:

  1) Resolve it to NCBI Taxonomy (rank-aware).
     - Prefer ranks: phylum > class > order > family > genus > species.
     - Avoid kingdom/superkingdom/domain if any lower rank is available.
     - Apply hard overrides for renamed phyla:
         Firmicutes    -> Bacillota (phylum)
         Bacteroidetes -> Bacteroidota (phylum)

  2) Extract:
       - current_scientific_name (NCBI ScientificName, cleaned)
       - taxon_rank

  3) Build PubMed queries using ONLY:
       - original taxon name (from CSV)
       - current_scientific_name
     No synonyms, no BlastName in queries.

     Q1: melanoma + immunotherapy + RESPONSE
       (bacteria original + current)
       AND melanoma
       AND immunotherapy terms
       → then filter hits by RESPONSE keywords in title/abstract

     Q2: melanoma + immunotherapy + SURVIVAL
       (bacteria original + current)
       AND melanoma
       AND immunotherapy terms
       → then filter hits by SURVIVAL keywords in title/abstract

     Q3/Q4 combined: other cancers + checkpoint or other immunotherapy
       (bacteria original + current)
       AND any cancer (Neoplasms/cancer/tumor*/carcinoma*)
       AND either PD‑1/PD‑L1/checkpoint drug terms or non‑PD1 IO terms
       AND NOT melanoma

     For each query we return a '; '-separated list of IDs, using PMCID if
     available (via pubmed_pmc elink), otherwise PMIDxxxx.

  4) Output CSV columns:

       taxon                       # original name from input CSV
       current_scientific_name     # NCBI current name, cleaned
       taxon_rank                  # NCBI rank (phylum, class, etc.)
       melanoma_response_ids       # Q1 IDs
       melanoma_survival_ids       # Q2 IDs
       other_cancer_io_ids         # Q3/Q4 combined IDs

Default output now points to `.runtime/results_melanoma_taxa_updated.csv` so the
raw pre-filtered table stays out of the working directory. The canonical
visible taxon/ID table is `taxa_ids_filtered.csv`.
"""

import os
import time
import re
import argparse
import pandas as pd
import requests
from xml.etree import ElementTree as ET

from low_abundance_phyla import compute_low_abundance_phyla, normalize_taxon_name

DEFAULT_OUTPUT_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    ".runtime",
    "results_melanoma_taxa_updated.csv",
)

# -------------------- NCBI config --------------------
BASE_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "your.name@institution.edu")
NCBI_KEY   = os.getenv("NCCI_API_KEY") or os.getenv("NCBI_API_KEY")  # support typo or correct

SLEEP = 0.35 if not NCBI_KEY else 0.12  # polite throttle
NCBI_MAX_RETRIES = int(os.getenv("NCBI_MAX_RETRIES", "5"))
NCBI_RETRY_STATUS_CODES = {429, 500, 502, 503, 504}

def eutil_params(extra=None):
    p = {"tool": "melanoma_taxa_ids_rank_aware", "email": NCBI_EMAIL}
    if NCBI_KEY:
        p["api_key"] = NCBI_KEY
    if extra:
        p.update(extra)
    return p

session = requests.Session()
session.headers.update({"User-Agent": f"melanoma_taxa_ids_rank_aware ({NCBI_EMAIL})"})


def ncbi_get(url, params=None, timeout=30, retries=NCBI_MAX_RETRIES):
    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            resp = session.get(url, params=params, timeout=timeout)
            time.sleep(SLEEP)
            if resp.status_code in NCBI_RETRY_STATUS_CODES:
                raise requests.HTTPError(
                    f"{resp.status_code} Server Error for url: {resp.url}",
                    response=resp,
                )
            resp.raise_for_status()
            return resp
        except (requests.Timeout, requests.ConnectionError, requests.HTTPError) as exc:
            last_exc = exc
            status_code = getattr(getattr(exc, "response", None), "status_code", None)
            is_retryable = status_code in NCBI_RETRY_STATUS_CODES or isinstance(
                exc,
                (requests.Timeout, requests.ConnectionError),
            )
            if not is_retryable or attempt == retries:
                raise
            wait_seconds = min(2 ** attempt, 20)
            print(
                f"  [NCBI retry {attempt}/{retries}] GET failed ({exc}). "
                f"Sleeping {wait_seconds}s before retry."
            )
            time.sleep(wait_seconds)
    raise last_exc

# -------------------- Outcome keyword sets (for filtering Q1/Q2) --------------------
KW_RESPONSE = {
    "response", "responder", "responders", "non-responder",
    "nonresponder", "non-response", "nonresponse",
    "objective response", "orr", "recist", "best overall response",
    "bor", "clinical benefit", "no clinical benefit"
}

KW_SURVIVAL = {
    "overall survival", "progression-free survival",
    "pfs", " os ", " os)", "(os ", " time to progression", " ttp",
    "hazard ratio", " hr ", " hr("
}

def has_any_keyword(text, keywords):
    if not text:
        return False
    t = text.lower()
    return any(k in t for k in keywords)

# -------------------- Helpers --------------------
def sanitize_taxon(t):
    t = str(t).strip()
    t = re.sub(r"[\[\]]", "", t)
    t = t.replace("_", " ")
    t = re.sub(r"['\u2019]", "", t)
    t = re.sub(r"(?i)\bunclassified\b", "", t)
    return " ".join(t.split())

def normalize_scientific_name(sci):
    """
    Strip authorship / year and suffixes like 'corrig.' from NCBI ScientificName.

    Heuristics:
      - Remove anything in parentheses: 'Escherichia coli (Migula 1895)...'
      - Remove 'corrig.' and everything after it.
      - Remove trailing 'Author 1990' or 'Author, 1990'.
    """
    if not sci:
        return ""
    name = sci.strip()

    # Drop parenthetical authorship
    name = re.sub(r"\s*\(.*?\)", "", name).strip()

    # Cut at 'corrig.' and following
    name = re.sub(r"\s+corrig\.\b.*$", "", name, flags=re.IGNORECASE).strip()

    # Cut at 'Author 1990'
    m = re.search(r"\s+[A-Z][a-zA-Z\-]+\s+\d{4}", name)
    if m:
        name = name[:m.start()].strip()

    # Cut at 'Author, 1990'
    m2 = re.search(r"\s+[A-Z][a-zA-Z\-]+,\s*\d{4}", name)
    if m2:
        name = name[:m2.start()].strip()

    return name

def normalize_for_match(s):
    """Lowercase, strip, collapse spaces for comparing names."""
    if not s:
        return ""
    s = re.sub(r"\s+", " ", s.strip())
    return s.lower()

# -------------------- Taxonomy (rank-aware selection) --------------------

RANK_SCORE = {
    "species": 4,
    "subspecies": 4,
    "genus": 4,
    "family": 5,
    "order": 6,
    "class": 7,
    "phylum": 8,
    "division": 8,
    "kingdom": 1,
    "superkingdom": 1,
    "domain": 1,
    "no rank": 2,
}

def rank_score(rank):
    return RANK_SCORE.get(rank.lower(), 2)

# Hard overrides for renamed phyla
PHYLUM_OVERRIDES = {
    "firmicutes": "1239",      # Bacillota (phylum)
    "bacteroidetes": "976",    # Bacteroidota (phylum)
}

def fetch_taxon_by_id(taxid, original_name=None):
    """
    Fetch a single taxonomy record by TaxID and basic info.
    Used for overrides and generic fetch.
    """
    try:
        ef = session.get(
            f"{BASE_EUTILS}/efetch.fcgi",
            params=eutil_params({"db": "taxonomy", "id": taxid, "retmode": "xml"}),
            timeout=20,
        )
        time.sleep(SLEEP)
        ef.raise_for_status()
        root = ET.fromstring(ef.content)
        taxon_node = root.find(".//Taxon")
        if taxon_node is None:
            return {"valid": False}

        sci = taxon_node.findtext("ScientificName") or ""
        lineage = taxon_node.findtext("Lineage") or ""
        if "Bacteria" not in lineage:
            return {"valid": False}

        rank = taxon_node.findtext("Rank") or ""
        blast_name = taxon_node.findtext(".//BlastName") or ""

        sci_clean = normalize_scientific_name(sci)

        # aliases are used only internally for scoring (not for PubMed queries)
        aliases = set()
        aliases.add(sci)
        for tag in ["Synonym", "EquivalentName", "GenbankSynonym", "GenbankCommonName"]:
            for node in taxon_node.findall(f".//{tag}"):
                txt = (node.text or "").strip()
                if txt:
                    aliases.add(txt)
        if blast_name:
            aliases.add(blast_name)
        if original_name:
            aliases.add(original_name)

        return {
            "valid": True,
            "ncbi_name": sci,
            "current_name": sci_clean,
            "taxid": str(taxid),
            "rank": rank,
            "blast_name": blast_name,
            "aliases": sorted(aliases),
        }
    except Exception:
        return {"valid": False}

def get_taxonomy_info(name):
    """
    Resolve a bacterial taxon in NCBI Taxonomy and choose best node
    (rank-aware, avoids kingdom-level where possible).

    Returns:
      dict(valid: bool,
           ncbi_name,
           current_name,
           taxid,
           rank,
           blast_name,
           aliases: [str])
    """
    clean_original = str(name).strip()
    if not clean_original:
        return {"valid": False}

    target_norm = normalize_for_match(normalize_scientific_name(clean_original))

    # Hard overrides for Firmicutes / Bacteroidetes
    override_taxid = PHYLUM_OVERRIDES.get(target_norm)
    if override_taxid:
        info = fetch_taxon_by_id(override_taxid, original_name=clean_original)
        if info.get("valid"):
            return info
        # if override fails, fall through to generic search

    # ESearch
    try:
        r = session.get(
            f"{BASE_EUTILS}/esearch.fcgi",
            params=eutil_params({
                "db": "taxonomy",
                "term": f"\"{clean_original}\"[All Names]",
                "retmode": "json",
                "retmax": 50,
            }),
            timeout=20,
        )
        time.sleep(SLEEP)
        r.raise_for_status()
        idlist = r.json().get("esearchresult", {}).get("idlist", [])
    except Exception:
        idlist = []

    if not idlist:
        # fallback: ScientificName only
        try:
            r = session.get(
                f"{BASE_EUTILS}/esearch.fcgi",
                params=eutil_params({
                    "db": "taxonomy",
                    "term": f"\"{clean_original}\"[SCIN]",
                    "retmode": "json",
                    "retmax": 10,
                }),
                timeout=20,
            )
            time.sleep(SLEEP)
            r.raise_for_status()
            idlist = r.json().get("esearchresult", {}).get("idlist", [])
        except Exception:
            idlist = []

    if not idlist:
        return {"valid": False}

    candidates = []

    for taxid in idlist:
        try:
            ef = session.get(
                f"{BASE_EUTILS}/efetch.fcgi",
                params=eutil_params({"db": "taxonomy", "id": taxid, "retmode": "xml"}),
                timeout=20,
            )
            time.sleep(SLEEP)
            ef.raise_for_status()
            root = ET.fromstring(ef.content)
            taxon_node = root.find(".//Taxon")
            if taxon_node is None:
                continue

            sci = taxon_node.findtext("ScientificName") or ""
            lineage = taxon_node.findtext("Lineage") or ""
            if "Bacteria" not in lineage:
                continue

            rank = taxon_node.findtext("Rank") or ""
            blast_name = taxon_node.findtext(".//BlastName") or ""

            sci_clean = normalize_scientific_name(sci)
            sci_norm = normalize_for_match(sci_clean)

            aliases = set()
            aliases.add(sci)
            for tag in ["Synonym", "EquivalentName", "GenbankSynonym", "GenbankCommonName"]:
                for node in taxon_node.findall(f".//{tag}"):
                    txt = (node.text or "").strip()
                    if txt:
                        aliases.add(txt)
            if blast_name:
                aliases.add(blast_name)
            aliases.add(clean_original)

            alias_norms = {normalize_for_match(normalize_scientific_name(a)) for a in aliases}

            match_score = 0
            if target_norm and (target_norm == sci_norm or target_norm in alias_norms):
                match_score += 10

            match_score += rank_score(rank)

            if rank.lower() in {"kingdom", "superkingdom", "domain"}:
                match_score -= 2

            candidates.append({
                "taxid": taxid,
                "sci": sci,
                "sci_clean": sci_clean,
                "rank": rank,
                "blast_name": blast_name,
                "aliases": sorted(aliases),
                "match_score": match_score,
            })
        except Exception:
            continue

    if not candidates:
        return {"valid": False}

    # Prefer non-kingdom/superkingdom/domain if present
    non_superhigh = [
        c for c in candidates
        if c["rank"].lower() not in {"kingdom", "superkingdom", "domain"}
    ]
    if non_superhigh:
        candidates = non_superhigh

    candidates.sort(key=lambda x: x["match_score"], reverse=True)
    best = candidates[0]

    if best["match_score"] < 5:
        return {"valid": False}

    return {
        "valid": True,
        "ncbi_name": best["sci"],
        "current_name": best["sci_clean"],
        "taxid": best["taxid"],
        "rank": best["rank"],
        "blast_name": best["blast_name"],
        "aliases": best["aliases"],  # NOT used in PubMed queries
    }

# -------------------- Bacteria clause for PubMed (ONLY original + current) --------------------

def bacteria_clause(original_name, current_name):
    """
    Build (original_name OR current_scientific_name) across Title/Abstract and MeSH.

    Used for ALL four queries (Q1–Q4). No synonyms, no BlastName, no other aliases.
    """
    names = set()
    if original_name:
        names.add(original_name.strip())
    if current_name:
        names.add(current_name.strip())

    ta_parts = []
    mesh_parts = []
    for nm in names:
        if not nm:
            continue
        nm = nm.replace("_", " ")
        ta_parts.append(f"\"{nm}\"[Title/Abstract]")
        mesh_parts.append(f"\"{nm}\"[MeSH Terms]")

    if not ta_parts:
        return ""
    ta = " OR ".join(ta_parts)
    mesh = " OR ".join(mesh_parts)
    return f"(({ta}) OR ({mesh}))"

# -------------------- Shared disease/therapy clauses --------------------

def melanoma_clause():
    return "(\"melanoma\"[Title/Abstract] OR \"melanoma\"[MeSH Terms])"

def cancer_not_melanoma_clause():
    cancer = "(\"Neoplasms\"[MeSH Terms] OR cancer[Title/Abstract] OR neoplasm*[Title/Abstract] OR tumor*[Title/Abstract] OR tumour*[Title/Abstract] OR carcinoma*[Title/Abstract])"
    not_mel = "NOT (\"melanoma\"[MeSH Terms] OR \"melanoma\"[Title/Abstract])"
    return f"({cancer} AND {not_mel})"

def immunotherapy_clause_general():
    """
    General immunotherapy / checkpoint / IO terms (for melanoma queries).
    """
    terms = [
        "immunotherapy",
        "immune checkpoint",
        "checkpoint inhibitor",
        "PD-1", "PD1",
        "PD-L1", "PDL1",
        "nivolumab", "pembrolizumab",
        "ipilimumab",
        "atezolizumab", "durvalumab", "avelumab",
        "cemiplimab", "relatlimab",
        "IL-2", "interleukin 2",
        "interferon alpha",
        "talimogene laherparepvec", "T-VEC",
    ]
    return "(" + " OR ".join([f"\"{t}\"[Title/Abstract]" for t in terms]) + ")"

def pd1_clause():
    """
    PD1/PDL1 / anti-PD1 antibody terms.
    """
    terms = [
        "PD-1", "PD1", "PDCD1",
        "PD-L1", "PDL1", "CD274",
        "nivolumab", "pembrolizumab",
        "atezolizumab", "durvalumab", "avelumab",
        "cemiplimab", "relatlimab",
    ]
    return "(" + " OR ".join([f"\"{t}\"[Title/Abstract]" for t in terms]) + ")"

def not_pd1_clause():
    """
    Exclusion for PD1-related terms (Title/Abstract).
    """
    return f"NOT {pd1_clause()}"

def other_io_clause():
    """
    Non-PD1 immunotherapy terms: CTLA4, CAR-T, IL-2, IFN, T-VEC, etc.
    (Deliberately no PD1/PD-L1 terms here.)
    """
    terms = [
        "ipilimumab",
        "CTLA-4", "CTLA4",
        "CAR-T", "chimeric antigen receptor",
        "TIL therapy",
        "interleukin 2", "IL-2",
        "interferon alpha",
        "talimogene laherparepvec", "T-VEC",
        "immunotherapy", "immune checkpoint", "checkpoint inhibitor",
    ]
    return "(" + " OR ".join([f"\"{t}\"[Title/Abstract]" for t in terms]) + ")"

# -------------------- Query builders --------------------

def build_query_melanoma_immuno(bact_clause):
    """
    Q1/Q2 base: bacteria (original + current) AND melanoma AND immunotherapy.
    RESPONSE/SURVIVAL is enforced by downstream filtering.
    """
    mel = melanoma_clause()
    immuno = immunotherapy_clause_general()
    return f"({bact_clause}) AND {mel} AND {immuno}"

def build_query_other_cancer_pd1(bact_clause):
    """
    Q3: bacteria (original + current)
        AND any cancer (NOT melanoma)
        AND PD1/PDL1/checkpoint drugs.

    Inclusion:
      - bacteria clause (original + current names)
      - any cancer terms (Neoplasms, cancer, tumor*, carcinoma*)
      - PD1/PDL1/drug terms in Title/Abstract

    Exclusion:
      - melanoma in MeSH or Title/Abstract
    """
    cancer_nm = cancer_not_melanoma_clause()
    pd1 = pd1_clause()
    return f"({bact_clause}) AND {cancer_nm} AND {pd1}"

def build_query_other_cancer_other_io(bact_clause):
    """
    Q4: bacteria (original + current)
        AND any cancer (NOT melanoma)
        AND non-PD1 IO terms
        AND NOT PD1 terms.

    Inclusion:
      - bacteria clause (original + current names)
      - any cancer terms (Neoplasms, cancer, tumor*, carcinoma*)
      - non-PD1 IO terms (CTLA4, CAR-T, IL-2, IFN, T-VEC, ipilimumab...)

    Exclusion:
      - melanoma (MeSH or Title/Abstract)
      - PD1/PDL1/drug terms (Title/Abstract)
    """
    cancer_nm = cancer_not_melanoma_clause()
    other_io = other_io_clause()
    not_pd1 = not_pd1_clause()
    return f"({bact_clause}) AND {cancer_nm} AND {other_io} AND {not_pd1}"

# -------------------- PubMed utilities --------------------
def pubmed_search(term, retmax=100):
    params = eutil_params({
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "retmax": retmax,
        "sort": "relevance",
    })
    try:
        r = ncbi_get(f"{BASE_EUTILS}/esearch.fcgi", params=params, timeout=30)
    except (requests.Timeout, requests.ConnectionError, requests.HTTPError) as exc:
        query_preview = re.sub(r"\s+", " ", term).strip()[:140]
        print(f"  [NCBI warning] PubMed search failed after retries: {exc}")
        print(f"  [NCBI warning] Query preview: {query_preview}")
        return []
    js = r.json()
    return js.get("esearchresult", {}).get("idlist", [])

def pubmed_fetch(pmids):
    if not pmids:
        return []
    out = []
    for chunk_start in range(0, len(pmids), 50):
        chunk = pmids[chunk_start:chunk_start+50]
        params = eutil_params({
            "db": "pubmed",
            "retmode": "xml",
            "id": ",".join(chunk),
        })
        try:
            r = ncbi_get(f"{BASE_EUTILS}/efetch.fcgi", params=params, timeout=60)
        except (requests.Timeout, requests.ConnectionError, requests.HTTPError) as exc:
            print(
                f"  [NCBI warning] PubMed fetch failed after retries for PMID chunk "
                f"{chunk_start + 1}-{chunk_start + len(chunk)}: {exc}"
            )
            continue
        root = ET.fromstring(r.content)
        for art in root.findall(".//PubmedArticle"):
            pmid = (art.findtext(".//PMID") or "").strip()
            title = (art.findtext(".//ArticleTitle") or "").strip()
            ab_parts = [t.text or "" for t in art.findall(".//Abstract/AbstractText")]
            abstract = " ".join(ab_parts).strip()
            out.append({
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
            })
    return out

pmcid_cache = {}

def pmid_to_pmcid(pmid):
    """
    Return PMCID string like PMC1234567 or None (cached).

    Uses elink with linkname=pubmed_pmc so we only get the PMCID
    corresponding to THAT PubMed article (not referenced articles).
    """
    if pmid in pmcid_cache:
        return pmcid_cache[pmid]

    try:
        r = session.get(
            f"{BASE_EUTILS}/elink.fcgi",
            params=eutil_params({
                "dbfrom": "pubmed",
                "db": "pmc",
                "linkname": "pubmed_pmc",
                "retmode": "json",
                "id": pmid,
            }),
            timeout=20,
        )
        time.sleep(SLEEP)
        r.raise_for_status()
        js = r.json()

        pmcid = None
        for linkset in js.get("linksets", []):
            for ldb in (linkset.get("linksetdbs") or []):
                if ldb.get("dbto") != "pmc":
                    continue
                if ldb.get("linkname") != "pubmed_pmc":
                    continue
                for li in ldb.get("links", []):
                    cid = li.get("id") if isinstance(li, dict) else li
                    if cid:
                        cid = str(cid).strip()
                        if not cid.upper().startswith("PMC"):
                            cid = "PMC" + cid
                        pmcid = cid
                        break
            if pmcid:
                break

        pmcid_cache[pmid] = pmcid
        return pmcid
    except Exception:
        pmcid_cache[pmid] = None
        return None

def encode_ids(pmids):
    """
    Given a list of PMIDs, return a '; '-separated string of IDs,
    preferring PMCID if available (PMCxxxx) else PMIDxxxx.
    """
    ids = []
    for pmid in pmids:
        pmcid = pmid_to_pmcid(pmid)
        if pmcid:
            ids.append(pmcid)
        else:
            ids.append(f"PMID{pmid}")
    return "; ".join(ids) if ids else ""


def unique_preserve_order(items):
    seen = set()
    out = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out

# -------------------- Context filter helpers (for Q1/Q2) --------------------
def filter_pmids_by_context(pmids, context):
    """
    context: 'response' | 'survival'
    Uses title + abstract keyword search to enforce:
      - Q1: RESPONSE-related language
      - Q2: SURVIVAL-related language
    """
    if not pmids:
        return []

    papers = pubmed_fetch(pmids)
    kept = []
    for p in papers:
        txt = (p["title"] + " " + p["abstract"]).lower()
        if context == "response":
            if has_any_keyword(txt, KW_RESPONSE):
                kept.append(p["pmid"])
        elif context == "survival":
            if has_any_keyword(txt, KW_SURVIVAL):
                kept.append(p["pmid"])

    pmid_set = set(kept)
    return [p for p in pmids if p in pmid_set]

# -------------------- Main --------------------
def should_skip_rank(rank, skip_ranks, original_name="", current_name="", allowed_low_abundance_phyla=None):
    """
    Skip selected taxonomy ranks before PubMed search/LLM pipeline.
    If kingdom is requested, superkingdom/domain are also skipped.
    """
    r = (rank or "").strip().lower()
    allowed_low_abundance_phyla = allowed_low_abundance_phyla or set()

    if r == "phylum":
        candidates = {
            normalize_taxon_name(original_name).lower(),
            normalize_taxon_name(current_name).lower(),
        }
        candidates.discard("")
        return not bool(candidates & allowed_low_abundance_phyla)

    if r in skip_ranks:
        return True
    if "kingdom" in skip_ranks and r in {"superkingdom", "domain"}:
        return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_csv", default="deepakTaxaUnique.csv",
                    help="Input CSV with a 'taxon' column or taxa in first column")
    ap.add_argument("--out", dest="out_csv", default=DEFAULT_OUTPUT_CSV,
                    help="Output CSV (default: .runtime/results_melanoma_taxa_updated.csv)")
    ap.add_argument("--retmax", type=int, default=150,
                    help="Max PubMed hits per query per taxon (before filtering)")
    ap.add_argument(
        "--skip-ranks",
        default="kingdom",
        help="Comma-separated ranks to skip for literature search (default: kingdom; low-abundance phyla are handled separately)",
    )
    ap.add_argument(
        "--phylum-threshold",
        type=float,
        default=float(os.getenv("LOW_ABUNDANCE_PHYLA_THRESHOLD", "0.10")),
        help="Include only phyla below this overall relative-abundance threshold (default: 0.10)",
    )
    args = ap.parse_args()

    skip_ranks = {x.strip().lower() for x in str(args.skip_ranks).split(",") if x.strip()}
    print(f"Skipping PubMed search for ranks: {sorted(skip_ranks)}")
    allowed_low_abundance_phyla = {
        normalize_taxon_name(name).lower()
        for name in compute_low_abundance_phyla(threshold=args.phylum_threshold)
    }
    print(
        "Including low-abundance phyla only (< "
        f"{args.phylum_threshold:.2f} of total phylum abundance): "
        f"{sorted(allowed_low_abundance_phyla)}"
    )

    df = pd.read_csv(args.in_csv)
    if df.shape[1] > 1:
        lower = [c.lower() for c in df.columns]
        col = df.columns[lower.index("taxon")] if "taxon" in lower else df.columns[0]
        taxa = df[col].astype(str).tolist()
    else:
        taxa = df.iloc[:, 0].astype(str).tolist()

    rows = []
    skipped_by_rank = 0
    for i, raw in enumerate(taxa, 1):
        original_name = str(raw)
        taxon = sanitize_taxon(raw)
        print(f"[{i}/{len(taxa)}] Taxon: {taxon} (original: {original_name})")

        info = get_taxonomy_info(taxon)
        if not info.get("valid"):
            print("  -> not a valid bacterial taxon under current rules, skipping")
            continue

        ncbi_name    = info["ncbi_name"]
        current_name = info["current_name"]
        rank         = info["rank"]
        blast_name   = info.get("blast_name", "")
        # aliases exist but are NOT used for PubMed queries anymore
        aliases      = info["aliases"]

        print(f"  Chosen taxon: {ncbi_name} (current='{current_name}', "
              f"rank={rank}, blast_name={blast_name}); "
              f"{len(aliases)} aliases (for taxonomy scoring only)")

        if should_skip_rank(
            rank,
            skip_ranks,
            original_name=taxon,
            current_name=current_name,
            allowed_low_abundance_phyla=allowed_low_abundance_phyla,
        ):
            skipped_by_rank += 1
            print(f"  -> skipped due to excluded rank: {rank}")
            continue

        # Build bacteria clause (ONLY original + current names)
        bact_clause = bacteria_clause(original_name, current_name)
        if not bact_clause:
            print("  -> no usable names for bacteria clause, skipping")
            continue

        # --- Build queries (using ONLY original name + current_scientific_name) ---

        # Q1/Q2 base: melanoma + immunotherapy
        q_mel_immuno = build_query_melanoma_immuno(bact_clause)

        # Q3: other cancers + PD1 (NOT melanoma)
        q_pd1        = build_query_other_cancer_pd1(bact_clause)

        # Q4: other cancers + non-PD1 IO (NOT melanoma, NOT PD1)
        q_other_io   = build_query_other_cancer_other_io(bact_clause)

        # --- PubMed searches ---
        pmids_mel_immuno = pubmed_search(q_mel_immuno, retmax=args.retmax)
        pmids_pd1        = pubmed_search(q_pd1,        retmax=args.retmax)
        pmids_other_io   = pubmed_search(q_other_io,   retmax=args.retmax)
        pmids_other_any  = unique_preserve_order(pmids_pd1 + pmids_other_io)

        print(f"  Melanoma immuno hits (before filter): {len(pmids_mel_immuno)}")
        print(f"  Other cancer PD1 hits: {len(pmids_pd1)}")
        print(f"  Other cancer other IO hits: {len(pmids_other_io)}")
        print(f"  Other cancer combined IO hits: {len(pmids_other_any)}")

        # --- Split melanoma immuno into RESPONSE vs SURVIVAL via abstract text ---
        pmids_mel_resp = filter_pmids_by_context(pmids_mel_immuno, "response")
        pmids_mel_surv = filter_pmids_by_context(pmids_mel_immuno, "survival")

        print(f"    -> Melanoma RESPONSE (after filter): {len(pmids_mel_resp)}")
        print(f"    -> Melanoma SURVIVAL (after filter): {len(pmids_mel_surv)}")

        row = {
            "taxon": original_name,                         # original string from CSV
            "current_scientific_name": current_name,        # cleaned NCBI name
            "taxon_rank": rank,                             # NCBI rank
            "melanoma_response_ids": encode_ids(pmids_mel_resp),
            "melanoma_survival_ids": encode_ids(pmids_mel_surv),
            "other_cancer_io_ids": encode_ids(pmids_other_any),
        }
        rows.append(row)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out_csv, index=False)
    print(f"\nSaved {len(rows)} taxa to {args.out_csv}")
    print(f"Skipped {skipped_by_rank} taxa due to excluded ranks.")

if __name__ == "__main__":
    main()
