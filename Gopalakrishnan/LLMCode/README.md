# LLMCode Pipeline README

## Purpose
This folder runs a staged evidence pipeline for microbiome taxa and immunotherapy literature.

The streamlined execution paths are:
- `./run_pipeline.sh` (full pipeline)
- `./run_from_filtered_ids.sh` (resume/update from `taxa_ids_filtered.csv`)

## Streamlined Changes Added
- Added `--from-filtered-ids` mode in `run_pipeline.sh`.
- Added `run_from_filtered_ids.sh` as a one-command downstream runner.
- Added incremental OpenAI reuse in `findPValues.py` to avoid rescanning already-processed article columns.
- Limited retry auto-detection in `retry_errors.py` to unresolved entries only.
- Added robust PMC title parsing (namespace-agnostic XML) to avoid false "PMCID invalid or no title" exclusions.

## Recommended Run Modes

### A) Full run
```bash
./run_pipeline.sh --python-only
```

### B) Resume/update after editing `taxa_ids_filtered.csv`
```bash
./run_from_filtered_ids.sh --python-only
```
This skips upstream taxonomy + PubMed ID discovery and starts from downstream files.
The resume path now re-applies `filter_taxa_ids.py` to the edited `taxa_ids_filtered.csv`
before rebuilding `unique_articles_ranked.csv`.

### Cornell proxy switching
If you are using the Cornell LiteLLM proxy, set `OPENAI_BASE_URL=https://api.ai.it.cornell.edu` and switch models with `LLM_PROVIDER` + `LLM_MODEL`.

```bash
export LLM_PROVIDER=openai
export LLM_MODEL=openai.gpt-5.4
```

```bash
export LLM_PROVIDER=gemini
export LLM_MODEL=google.gemini-3.1-pro-preview
```

```bash
export LLM_PROVIDER=claude
export LLM_MODEL=anthropic.claude-4.6-opus
```

With the proxy configured, the same Cornell key can be used for all three providers.

### C) Skip LLM calls
```bash
./run_pipeline.sh --python-only --no-llm
```

## File-by-File Functions, Inputs, Outputs

### Environment + Orchestration

1) `setup_conda_env.sh`
- Function: create/update conda environment with required Python/R dependencies.
- Input: optional env name argument (default `chai-pipeline`).
- Output: conda env with packages installed.

2) `run_pipeline.sh`
- Function: orchestrate the pipeline stages with checkpoints, archiving, and key preflight checks.
- Input:
  - files depending on stage (`deepakTaxonomy.csv`, `deepakTaxaUnique.csv`, `taxa_ids_filtered.csv`, etc.)
  - env vars from `.env` (`OPENAI_API_KEY`, optional `OPENAI_BASE_URL`, optional `LLM_MODEL`, `NCBI_EMAIL`, optional `NCBI_API_KEY`, `OPENAI_MODEL`).
- Output: all stage outputs listed below.
- Key flags:
  - `--python-only`
  - `--from-filtered-ids`
  - `--no-llm`
  - `--skip-retry`
  - `--no-archive`

3) `run_from_filtered_ids.sh`
- Function: convenience wrapper for interrupted workflows.
- Input: existing `taxa_ids_filtered.csv`.
- Output: downstream outputs regenerated from that point.

### Stage 1: Taxa Preprocessing

4) `uniqueTaxa.R`
- Function: deduplicate and clean taxa names from taxonomy table.
- Input: `deepakTaxonomy.csv` (or custom via args).
- Output: `deepakTaxaUnique.csv`.

### Stage 2: Taxonomy Resolution + PubMed ID Discovery

5) `melanoma_taxa_ids_short.py`
- Function: resolve taxa via NCBI taxonomy and generate query-specific PMCID/PMID sets.
- Input: `deepakTaxaUnique.csv`.
- Output: `.runtime/results_melanoma_taxa_updated.csv` (transient raw output).
- Notes:
  - uses only original + current scientific name in query clauses.
  - default rank skip for search: `kingdom,phylum,species` (`--skip-ranks`).
  - the canonical visible taxon/ID table is still `taxa_ids_filtered.csv`.

### Stage 3: ID Filtering

6) `filter_taxa_ids.py`
- Function: clean ID lists, remove invalid/excluded study IDs.
- Input: `.runtime/results_melanoma_taxa_updated.csv` by default, or `taxa_ids_filtered.csv` on resume.
- Output: `taxa_ids_filtered.csv`.
- Note: the pipeline now keeps only `taxa_ids_filtered.csv` in the working directory; the raw pre-filter file lives under `.runtime/`.
- Filters include:
  - invalid PMCID/PMID,
  - retracted papers,
  - papers whose PubMed citation is not indexed for MEDLINE,
  - review/meta-analysis/systematic review,
  - pure animal studies,
  - explicit Deepak paper exclusion (`PMC5827966`).

### Stage 4: Unique Article Ranking

7) `rank_unique_articles.py`
- Function: compile unique PMCID/PMID list and rank by frequency across taxa/query columns.
- Input: `taxa_ids_filtered.csv`.
- Output: `unique_articles_ranked.csv`.

### Stage 5: Taxon-Article Presence Matrix

8) `taxaPaperTable.py`
- Function: fetch article text (PMC full-text or PubMed abstract) and build 0/1 taxon-article presence matrix.
- Input:
  - `unique_articles_ranked.csv`
  - `taxa_ids_filtered.csv`
- Output: `taxon_article_fulltext_matrix.csv`.
- Restriction: skips taxa with ranks in `{kingdom, superkingdom, domain, phylum, species}`.

### Stage 6: LLM Extraction (OpenAI)

9) `findPValues.py`
- Function: extract directional significance, p-values, context marker (A/B/C/D), optional sample size.
- Input: `taxon_article_fulltext_matrix.csv`.
- Output: `taxon_article_matrix_filled.csv`.
- Runtime state: checkpoints, cost trackers, and bracket-cleanup logs now live under `.runtime/`.
- Checkpoints: `.runtime/checkpoints/checkpoint_col_*.csv` during run.
- Token-saving behavior:
  - reuses existing extracted cells from prior `taxon_article_matrix_filled.csv` or `taxon_article_matrix_filled_fixed.csv`;
  - skips OpenAI rescans for already-scanned article columns.

10) `retry_errors.py`
- Function: targeted retries for unresolved/failed extraction columns.
- Input:
  - `taxon_article_fulltext_matrix.csv`
  - `taxon_article_matrix_filled.csv`
- Output: `taxon_article_matrix_filled_fixed.csv`.
- Retry scope:
  - explicit retry ID list + unresolved `1/1.0` columns from `taxon_article_matrix_filled.csv`.

### Stage 7: Aggregation and Final Mapping

11) `cleanTable.R`
- Function:
  - parse extracted cell strings,
  - enforce directional entries,
  - convert p-values to signed z,
  - compute unweighted and weighted Stouffer z.
- Input: `taxon_article_matrix_filled_fixed.csv`.
- Outputs:
  - `taxon_ABCD.csv`
  - `taxon_ABCD_probit_stoufferZ_unweighted.csv`
  - `taxon_ABCD_probit_stoufferZ_weighted.csv`

12) `produceFinalTable.R`
- Function: map taxon-level ABCD scores onto taxonomy tables using nearest-available fallback from the selected target rank to the allowed coarser auxiliary ranks.
- Inputs:
  - `taxon_ABCD_probit_stoufferZ_unweighted.csv`
  - `taxon_ABCD_probit_stoufferZ_weighted.csv`
  - `taxa_ids_filtered.csv`
  - `d1Taxonomy.csv`
- Outputs:
  - `d1Taxonomy_OTU_ABCD_{class|order|family|genus}_and_below_{unweighted|weighted}.csv`
  - `sideCov/openAIGenerated/{unweighted|weighted}/selectionTarget*/d1Taxonomy_*_ABCD_*_and_below_{unweighted|weighted}.csv`

### Utilities

13) `countUniqueStudies.R`
- Utility script for counting PMID mentions in a results table.

## Restrictions Added (Current Behavior)

1. Taxonomic rank exclusions for search/summarization:
- `kingdom`, `phylum`, `species` are skipped by default in PubMed search stage.
- `superkingdom`/`domain` are also treated as excluded in downstream matrix stage.

2. NCBI-level article filtering before LLM:
- remove review/systematic-review/meta-analysis.
- remove pure animal studies.
- remove retracted papers.
- remove papers whose PubMed citation is not indexed for MEDLINE.
- remove explicit non-auxiliary ID `PMC5827966`.

3. Text-only evidence policy:
- extraction prompts enforce text-only reading (no figure/plot inference).

4. Directional requirement:
- significant result entries without clear direction are excluded.

5. Missing p-value fallback:
- if significant but numeric p-value missing, use `p=0.05`.

6. Incremental token-saving:
- previously scanned article columns are reused and skipped in `findPValues.py`.

7. Robust PMCID parsing:
- namespace-agnostic `article-title` extraction prevents false invalid-PMCID drops.

8. Weighted and unweighted outputs:
- both Stouffer variants are generated and carried through OTU-level and `sideCov/openAIGenerated` outputs.

## Notes on Interruptions
- `findPValues.py` stores column-level checkpoints in `.runtime/checkpoints/` and can resume.
- For interrupted downstream runs, prefer:
```bash
./run_from_filtered_ids.sh --python-only
```
- To avoid archiving behavior during debugging:
```bash
./run_pipeline.sh --from-filtered-ids --python-only --no-archive
```
