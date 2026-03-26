#!/usr/bin/env python3
import csv
import re


RANK_LEVELS = [
    "superkingdom",
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
]
RANK_INDEX = {rank: idx for idx, rank in enumerate(RANK_LEVELS)}
LOWER_LEVEL_MARKERS = (
    " otu",
    " otus",
    " asv",
    " asvs",
    " sgb",
    " sgbs",
    " strain",
    " strains",
    " sequence variant",
    " sequence variants",
    " amplicon",
    " amplicons",
)
MOUSE_MODEL_MARKERS = (
    " mouse ",
    " mice ",
    " murine ",
)


def normalize_taxon_text(text):
    text = str(text or "")
    text = re.sub(r"^[a-z]__", "", text, flags=re.IGNORECASE)
    text = text.replace("[", " ").replace("]", " ")
    text = re.sub(r"\s+", " ", text).strip().lower()
    return text


def _whole_name_in_text(name, text):
    if not name or not text:
        return False
    pattern = r"(?<![A-Za-z0-9])" + re.escape(name) + r"(?![A-Za-z0-9])"
    return re.search(pattern, text) is not None


def _contains_lower_level_marker(text):
    text = f" {normalize_taxon_text(text)} "
    return any(marker in text for marker in LOWER_LEVEL_MARKERS)


def _mentions_mouse_model(text):
    text = f" {normalize_taxon_text(text)} "
    return any(marker in text for marker in MOUSE_MODEL_MARKERS)


def _extract_species_genera(raw_text):
    genera = set()
    for genus, species in re.findall(r"\b([A-Z][a-z]+)\s+([a-z][A-Za-z0-9._-]+)\b", str(raw_text or "")):
        if species.lower() in {"sp", "spp", "cf", "aff"}:
            continue
        genera.add(normalize_taxon_text(genus))
    return genera


def load_taxon_metadata(path):
    by_name = {}
    known_names = []

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            taxon = row.get("taxon", "")
            current_name = row.get("current_scientific_name", "") or taxon
            rank = normalize_taxon_text(row.get("taxon_rank", ""))
            meta = {
                "taxon": taxon,
                "current_scientific_name": current_name,
                "taxon_rank": rank,
            }

            for candidate in (taxon, current_name):
                norm = normalize_taxon_text(candidate)
                if not norm:
                    continue
                by_name.setdefault(norm, meta)
                known_names.append((norm, rank))

    deduped_names = []
    seen = set()
    for item in known_names:
        if item in seen:
            continue
        seen.add(item)
        deduped_names.append(item)

    return {
        "by_name": by_name,
        "known_names": deduped_names,
    }


def get_taxon_metadata(taxon_name, metadata_store):
    norm = normalize_taxon_text(taxon_name)
    meta = metadata_store.get("by_name", {}).get(norm, {})
    return {
        "taxon": meta.get("taxon", taxon_name),
        "current_scientific_name": meta.get("current_scientific_name", taxon_name),
        "taxon_rank": normalize_taxon_text(meta.get("taxon_rank", "")),
    }


def build_target_payloads(target_taxa, metadata_store):
    payloads = []
    for taxon_name in target_taxa:
        meta = get_taxon_metadata(taxon_name, metadata_store)
        rank = meta["taxon_rank"]
        if rank == "species":
            rule = (
                "Accept exact-species findings. OTU/ASV/SGB/strain findings count only when the text "
                "explicitly ties them to this exact species."
            )
        elif rank == "genus":
            rule = (
                "Accept findings about this genus itself only. If the text explicitly names a species "
                "within the genus, do not count it for the genus. Discard mouse-model findings."
            )
        else:
            rule = (
                "Accept findings only when they are explicitly about this taxon itself. Do not treat "
                "species, OTUs, ASVs, SGBs, strains, or any finer descendants as evidence for this taxon. "
                "Discard mouse-model findings."
            )
        payloads.append(
            {
                "taxon": taxon_name,
                "current_scientific_name": meta["current_scientific_name"],
                "taxon_rank": rank,
                "evidence_rule": rule,
            }
        )
    return payloads


def _mentions_other_species(effect_text_norm, allowed_names, metadata_store):
    for name, rank in metadata_store.get("known_names", []):
        if rank != "species" or name in allowed_names:
            continue
        if _whole_name_in_text(name, effect_text_norm):
            return True
    return False


def _mentions_finer_taxon(effect_text_norm, target_rank, allowed_names, metadata_store):
    target_idx = RANK_INDEX.get(target_rank)
    if target_idx is None:
        return False

    for name, rank in metadata_store.get("known_names", []):
        rank_idx = RANK_INDEX.get(rank)
        if rank_idx is None or rank_idx <= target_idx or name in allowed_names:
            continue
        if _whole_name_in_text(name, effect_text_norm):
            return True
    return False


def is_result_specific_to_taxon(taxon_name, res_obj, metadata_store):
    if not isinstance(res_obj, dict):
        return False

    meta = get_taxon_metadata(taxon_name, metadata_store)
    rank = meta["taxon_rank"]
    effect_text = str(res_obj.get("eff", "") or "")
    effect_text_norm = normalize_taxon_text(effect_text)
    if not effect_text_norm:
        return True
    if _mentions_mouse_model(effect_text):
        return False

    allowed_names = {
        normalize_taxon_text(taxon_name),
        normalize_taxon_text(meta.get("taxon", "")),
        normalize_taxon_text(meta.get("current_scientific_name", "")),
    }
    allowed_names.discard("")
    species_genera = _extract_species_genera(effect_text)

    if rank == "species":
        if _mentions_other_species(effect_text_norm, allowed_names, metadata_store):
            return False
        if _contains_lower_level_marker(effect_text_norm):
            return any(_whole_name_in_text(name, effect_text_norm) for name in allowed_names)
        return True

    if rank == "genus":
        if species_genera:
            return False
        if _contains_lower_level_marker(effect_text_norm):
            return any(_whole_name_in_text(name, effect_text_norm) for name in allowed_names)
        return not _mentions_finer_taxon(effect_text_norm, rank, allowed_names, metadata_store)

    if species_genera:
        return False
    if _contains_lower_level_marker(effect_text_norm):
        return False
    if _mentions_finer_taxon(effect_text_norm, rank, allowed_names, metadata_store):
        return False
    return True


def lookup_result_for_target(results, taxon_name, metadata_store):
    if not isinstance(results, dict):
        return None

    key_map = {
        normalize_taxon_text(key): value
        for key, value in results.items()
        if isinstance(key, str)
    }
    meta = get_taxon_metadata(taxon_name, metadata_store)
    candidate_keys = [
        normalize_taxon_text(taxon_name),
        normalize_taxon_text(meta.get("taxon", "")),
        normalize_taxon_text(meta.get("current_scientific_name", "")),
    ]

    for key in candidate_keys:
        if key and key in key_map:
            return key_map[key]
    return results.get(taxon_name)
