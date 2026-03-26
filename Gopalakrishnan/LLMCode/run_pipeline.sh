#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Load environment variables from local .env if present.
# This parser tolerates simple quoted values and legacy keys with spaces.
if [[ -f ".env" ]]; then
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    [[ "$line" =~ ^[[:space:]]*$ ]] && continue
    [[ "$line" == *=* ]] || continue

    key="${line%%=*}"
    value="${line#*=}"

    key="${key#"${key%%[![:space:]]*}"}"
    key="${key%"${key##*[![:space:]]}"}"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"

    if [[ "$value" == \"*\" && "$value" == *\" ]]; then
      value="${value:1:${#value}-2}"
    elif [[ "$value" == \'*\' && "$value" == *\' ]]; then
      value="${value:1:${#value}-2}"
    fi

    case "$key" in
      "API Key")
        key="API_KEY"
        ;;
      "Cornell_API_Key")
        key="CORNELL_API_KEY"
        ;;
      "API Base URI")
        key="API_BASE_URI"
        ;;
      "API_Base_URI")
        key="API_BASE_URI"
        ;;
      "API Base URL")
        key="API_BASE_URL"
        ;;
      "API_Base_URL")
        key="API_BASE_URL"
        ;;
    esac

    [[ "$key" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
    export "$key=$value"
  done < ".env"
fi

export OPENAI_API_KEY="${OPENAI_API_KEY:-${API_KEY:-}}"
export OPENAI_API_KEY="${OPENAI_API_KEY:-${CORNELL_API_KEY:-}}"
export OPENAI_BASE_URL="${OPENAI_BASE_URL:-${API_BASE_URI:-}}"
export LLM_MODEL="${LLM_MODEL:-}"

export OPENAI_MODEL="${OPENAI_MODEL:-gpt-5.4}"
export OPENAI_REASONING_EFFORT="${OPENAI_REASONING_EFFORT:-high}"
export LLM_PROVIDER="${LLM_PROVIDER:-openai}"
export LLM_FILE_TAG="${LLM_FILE_TAG:-}"
export LLM_RUNTIME_DIR="${LLM_RUNTIME_DIR:-$SCRIPT_DIR/.runtime}"
export LLM_CHECKPOINT_DIR="${LLM_CHECKPOINT_DIR:-$LLM_RUNTIME_DIR/checkpoints}"
RAW_TAXA_IDS_CSV="${RAW_TAXA_IDS_CSV:-$LLM_RUNTIME_DIR/results_melanoma_taxa_updated.csv}"

PYTHON_BIN="${PYTHON_BIN:-python}"
RETMAX="${RETMAX:-150}"
RUN_R=1
RUN_LLM=1
SKIP_RETRY=0
ARCHIVE_OLD=1
START_FROM="full"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
  log "ERROR: $*"
  exit 1
}

ensure_runtime_layout() {
  mkdir -p "$LLM_RUNTIME_DIR" "$LLM_CHECKPOINT_DIR" "$LLM_RUNTIME_DIR/cost" "$LLM_RUNTIME_DIR/logs"
}

cleanup_root_transient_files() {
  local removed=0
  local path

  local direct_paths=(
    ".DS_Store"
    ".Rhistory"
    "__pycache__"
    "results_melanoma_taxa_updated.csv"
    "taxa_ids.csv"
    "taxa_ids_filtered.refiltered.csv"
    "article_exclusion_report.csv"
  )

  for path in "${direct_paths[@]}"; do
    if [[ -e "$path" ]]; then
      rm -rf "$path"
      removed=$((removed + 1))
    fi
  done

  shopt -s nullglob
  local legacy_matches=(
    checkpoint_col_*.csv
    *_cost_tracker.json
    *_bracket_cleanup_cost.json
    bracket_cleanup_log*.csv
  )
  if (( ${#legacy_matches[@]} > 0 )); then
    rm -f "${legacy_matches[@]}"
    removed=$((removed + ${#legacy_matches[@]}))
  fi
  shopt -u nullglob

  if [[ "$removed" -gt 0 ]]; then
    log "Cleanup step: removed ${removed} transient root artifact(s)."
  else
    log "Cleanup step: no transient root artifacts found."
  fi
}

usage() {
  cat <<'EOF'
Usage:
  ./run_pipeline.sh [options]

Options:
  --python-only   Run only Python steps (skip R scripts)
  --from-filtered-ids
                 Start from taxa_ids_filtered.csv (skip unique taxa + taxonomy search + ID filtering)
  --no-llm        Skip LLM extraction/retry steps
  --skip-retry    Skip retry_errors.py
  --no-archive    Do not move previous outputs to archive/
  --provider P    LLM provider: openai, claude, or gemini
  --retmax N      PubMed retmax for melanoma_taxa_ids_short.py (default: 150)
  -h, --help      Show this help

Environment variables:
  NCBI_EMAIL      Recommended for NCBI requests
  NCBI_API_KEY    Optional NCBI API key
  LLM_PROVIDER    Optional (default: openai)
  OPENAI_API_KEY  Required for provider=openai unless --no-llm
  OPENAI_BASE_URL Optional OpenAI-compatible proxy base URL
                 Example: https://api.ai.it.cornell.edu
  LLM_MODEL       Optional exact model ID override
                 Example: openai.gpt-5.4
  ANTHROPIC_API_KEY
                 Required for provider=claude unless --no-llm
  GEMINI_API_KEY  Required for provider=gemini unless --no-llm
  OPENAI_MODEL    Optional (default in scripts: gpt-5.4)
  OPENAI_REASONING_EFFORT
                 Optional (default: high)
  LLM_RUNTIME_DIR Optional runtime state directory (default: ./.runtime)
  PYTHON_BIN      Python executable (default: python)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python-only)
      RUN_R=0
      shift
      ;;
    --from-filtered-ids)
      START_FROM="filtered_ids"
      shift
      ;;
    --no-llm)
      RUN_LLM=0
      shift
      ;;
    --skip-retry)
      SKIP_RETRY=1
      shift
      ;;
    --no-archive)
      ARCHIVE_OLD=0
      shift
      ;;
    --provider)
      [[ $# -ge 2 ]] || die "--provider requires a value"
      LLM_PROVIDER="$2"
      shift 2
      ;;
    --retmax)
      [[ $# -ge 2 ]] || die "--retmax requires a value"
      RETMAX="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown option: $1"
      ;;
  esac
done

case "${LLM_PROVIDER}" in
  openai)
    export OPENAI_MODEL="${OPENAI_MODEL:-gpt-5.4}"
    export OPENAI_COST_TRACKER_FILE="${OPENAI_COST_TRACKER_FILE:-$LLM_RUNTIME_DIR/cost/openai_cost_tracker.json}"
    ;;
  claude)
    export ANTHROPIC_MODEL="${ANTHROPIC_MODEL:-claude-opus-4-6}"
    export OPENAI_COST_TRACKER_FILE="${OPENAI_COST_TRACKER_FILE:-$LLM_RUNTIME_DIR/cost/claude_cost_tracker.json}"
    ;;
  gemini)
    export GEMINI_MODEL="${GEMINI_MODEL:-gemini-3.1-pro-preview}"
    export OPENAI_COST_TRACKER_FILE="${OPENAI_COST_TRACKER_FILE:-$LLM_RUNTIME_DIR/cost/gemini_cost_tracker.json}"
    ;;
  *)
    die "Unsupported provider: ${LLM_PROVIDER}"
    ;;
esac
export LLM_FILE_TAG="${LLM_FILE_TAG:-${LLM_PROVIDER}}"
export BRACKET_CLEANUP_COST_TRACKER_FILE="${BRACKET_CLEANUP_COST_TRACKER_FILE:-$LLM_RUNTIME_DIR/cost/${LLM_PROVIDER}_bracket_cleanup_cost.json}"
export BRACKET_CLEANUP_LOG_FILE="${BRACKET_CLEANUP_LOG_FILE:-$LLM_RUNTIME_DIR/logs/bracket_cleanup_log_${LLM_FILE_TAG}.csv}"

command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python not found: $PYTHON_BIN"
if [[ "$RUN_R" == "1" ]]; then
  command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH"
fi
if [[ "$RUN_LLM" == "1" ]]; then
  case "$LLM_PROVIDER" in
    openai)
      [[ -n "${OPENAI_API_KEY:-}" ]] || die "OPENAI_API_KEY (or API_KEY) is not set for provider=openai. Use --no-llm to skip."
      ;;
    claude)
      if [[ -n "${OPENAI_BASE_URL:-}" ]]; then
        [[ -n "${OPENAI_API_KEY:-}" ]] || die "OPENAI_API_KEY/CORNELL_API_KEY is not set for provider=claude with OPENAI_BASE_URL."
      else
        [[ -n "${ANTHROPIC_API_KEY:-}" ]] || die "ANTHROPIC_API_KEY is not set (required for provider=claude). Use --no-llm to skip."
      fi
      ;;
    gemini)
      if [[ -n "${OPENAI_BASE_URL:-}" ]]; then
        [[ -n "${OPENAI_API_KEY:-}" ]] || die "OPENAI_API_KEY/CORNELL_API_KEY is not set for provider=gemini with OPENAI_BASE_URL."
      else
        [[ -n "${GEMINI_API_KEY:-}" ]] || die "GEMINI_API_KEY is not set (required for provider=gemini). Use --no-llm to skip."
      fi
      ;;
  esac
fi
if [[ -z "${NCBI_EMAIL:-}" ]]; then
  log "Warning: NCBI_EMAIL is not set. NCBI may throttle requests."
fi

run_cmd() {
  log "RUN: $*"
  "$@"
}

provider_tagged_name() {
  local name="$1"
  local tag="${LLM_FILE_TAG:-}"
  if [[ -z "$tag" ]]; then
    printf '%s\n' "$name"
    return 0
  fi

  local stem="${name%.*}"
  local ext=""
  if [[ "$name" == *.* ]]; then
    ext=".${name##*.}"
  fi
  printf '%s_%s%s\n' "$stem" "$tag" "$ext"
}

archive_previous_outputs() {
  if [[ "$ARCHIVE_OLD" != "1" ]]; then
    log "Archive step skipped (--no-archive)."
    return 0
  fi

  local out_files=(
    "deepakTaxaUnique.csv"
    "taxa_ids_filtered.csv"
    "unique_articles_ranked.csv"
    "taxon_article_fulltext_matrix.csv"
    "taxon_ABCD.csv"
    "taxon_ABCD_probit_stoufferZ_unweighted.csv"
    "taxon_ABCD_probit_stoufferZ_weighted.csv"
    "d1Taxonomy_OTU_ABCD_class_and_below_unweighted.csv"
    "d1Taxonomy_OTU_ABCD_order_and_below_unweighted.csv"
    "d1Taxonomy_OTU_ABCD_family_and_below_unweighted.csv"
    "d1Taxonomy_OTU_ABCD_genus_and_below_unweighted.csv"
    "d1Taxonomy_OTU_ABCD_species_and_below_unweighted.csv"
    "d1Taxonomy_OTU_ABCD_class_and_below_weighted.csv"
    "d1Taxonomy_OTU_ABCD_order_and_below_weighted.csv"
    "d1Taxonomy_OTU_ABCD_family_and_below_weighted.csv"
    "d1Taxonomy_OTU_ABCD_genus_and_below_weighted.csv"
    "d1Taxonomy_OTU_ABCD_species_and_below_weighted.csv"
  )

  local found=0
  local f
  for f in "${out_files[@]}"; do
    if [[ "$RUN_R" != "1" && "$f" == "deepakTaxaUnique.csv" ]]; then
      # In --python-only mode this is a required input file, not an output.
      continue
    fi
    if [[ "$START_FROM" == "filtered_ids" && "$f" == "taxa_ids_filtered.csv" ]]; then
      # In --from-filtered-ids mode this is a required input file, not an output.
      continue
    fi
    if [[ -f "$f" ]]; then
      found=1
      break
    fi
  done

  if [[ "$found" -eq 0 ]]; then
    log "Archive step: no previous output files found."
    return 0
  fi

  local ts
  local arch_dir
  ts="$(date '+%Y%m%d_%H%M%S')"
  arch_dir="archive/${ts}"
  mkdir -p "$arch_dir"

  local moved=0
  for f in "${out_files[@]}"; do
    if [[ "$RUN_R" != "1" && "$f" == "deepakTaxaUnique.csv" ]]; then
      continue
    fi
    if [[ "$START_FROM" == "filtered_ids" && "$f" == "taxa_ids_filtered.csv" ]]; then
      continue
    fi
    if [[ -f "$f" ]]; then
      mv "$f" "$arch_dir/"
      moved=$((moved + 1))
    fi
  done

  log "Archive step: moved ${moved} file(s) to ${arch_dir}"
}

check_ncbi_key() {
  if [[ -z "${NCBI_API_KEY:-}" ]]; then
    log "NCBI key check: NCBI_API_KEY is not set (continuing without API key)."
    return 0
  fi

  local url
  local resp
  url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=melanoma&retmode=json&retmax=1&api_key=${NCBI_API_KEY}"
  resp="$(curl -sS --max-time 20 "$url" || true)"

  if [[ -z "$resp" ]]; then
    log "Warning: NCBI key check failed (empty response). Continuing, but key may be unusable."
    return 0
  fi

  if echo "$resp" | grep -Eiq "api key.*(invalid|not valid|incorrect)|\"error\"|<ERROR>"; then
    log "Warning: NCBI_API_KEY appears invalid/unusable based on E-utilities response."
    log "NCBI response snippet: $(echo "$resp" | tr '\n' ' ' | cut -c1-180)"
    return 0
  fi

  if echo "$resp" | grep -q "\"esearchresult\""; then
    log "NCBI key check: OK (E-utilities accepted the API key)."
    return 0
  fi

  log "Warning: NCBI key check inconclusive. Response did not contain expected fields."
  log "NCBI response snippet: $(echo "$resp" | tr '\n' ' ' | cut -c1-180)"
  return 0
}

check_llm_provider() {
  if [[ "$RUN_LLM" != "1" ]]; then
    log "LLM preflight skipped (--no-llm)."
    return 0
  fi

  log "Checking LLM provider='${LLM_PROVIDER}'..."
  local out
  set +e
  out="$("$PYTHON_BIN" - <<'PY'
import os
import sys

try:
    from llm_provider_client import assert_provider_ready, get_provider_model, normalize_provider
except Exception as e:
    print(f"LLM_PREFLIGHT_FAIL: helper import failed: {e}")
    sys.exit(2)

provider = normalize_provider(os.getenv("LLM_PROVIDER", "openai"))
model = get_provider_model(provider)

try:
    assert_provider_ready(provider)
    print(f"LLM_PREFLIGHT_OK provider={provider} model={model}")
    sys.exit(0)
except Exception as e:
    print("LLM_PREFLIGHT_FAIL: " + str(e))
sys.exit(2)
PY
)"
  local status=$?
  set -e

  if [[ $status -ne 0 ]]; then
    log "$out"
    die "LLM provider preflight failed. Verify provider API key/model and install deps via ./setup_conda_env.sh."
  fi

  log "$out"
  log "LLM provider preflight: OK."
}

log "Pipeline start in: $SCRIPT_DIR"
log "Start mode: ${START_FROM}"
log "LLM configuration: provider=${LLM_PROVIDER}, model=${LLM_MODEL:-${OPENAI_MODEL:-${ANTHROPIC_MODEL:-${GEMINI_MODEL:-}}}}, file_tag=${LLM_FILE_TAG}, reasoning=${OPENAI_REASONING_EFFORT}, cost_tracker=${OPENAI_COST_TRACKER_FILE}"
ensure_runtime_layout
cleanup_root_transient_files
archive_previous_outputs
check_ncbi_key
check_llm_provider

if [[ "$START_FROM" == "full" ]]; then
  if [[ "$RUN_R" == "1" ]]; then
    run_cmd Rscript uniqueTaxa.R deepakTaxonomy.csv deepakTaxaUnique.csv
  else
    log "Skipping uniqueTaxa.R because --python-only was set."
    if [[ ! -f deepakTaxaUnique.csv ]]; then
      die "deepakTaxaUnique.csv not found. Generate it first or run without --python-only."
    fi
  fi

  run_cmd "$PYTHON_BIN" melanoma_taxa_ids_short.py \
    --in deepakTaxaUnique.csv \
    --out "$RAW_TAXA_IDS_CSV" \
    --retmax "$RETMAX"

  run_cmd "$PYTHON_BIN" filter_taxa_ids.py \
    --in "$RAW_TAXA_IDS_CSV" \
    --out taxa_ids_filtered.csv
elif [[ "$START_FROM" == "filtered_ids" ]]; then
  log "Skipping upstream steps; expecting existing taxa_ids_filtered.csv"
  [[ -f taxa_ids_filtered.csv ]] || die "taxa_ids_filtered.csv not found."
  filtered_tmp="taxa_ids_filtered.refiltered.csv"
  run_cmd "$PYTHON_BIN" filter_taxa_ids.py \
    --in taxa_ids_filtered.csv \
    --out "$filtered_tmp"
  mv "$filtered_tmp" taxa_ids_filtered.csv
  log "Refreshed taxa_ids_filtered.csv with current article exclusion rules"
else
  die "Unsupported start mode: ${START_FROM}"
fi

run_cmd "$PYTHON_BIN" rank_unique_articles.py \
  --in taxa_ids_filtered.csv \
  --out unique_articles_ranked.csv

run_cmd "$PYTHON_BIN" taxaPaperTable.py

if [[ "$RUN_LLM" == "1" ]]; then
  run_cmd env OPENAI_COST_TRACKER_RESET=1 "$PYTHON_BIN" findPValues.py
  filled_file="$(provider_tagged_name taxon_article_matrix_filled.csv)"
  fixed_file="$(provider_tagged_name taxon_article_matrix_filled_fixed.csv)"
  if [[ "$SKIP_RETRY" == "0" ]]; then
    run_cmd env OPENAI_COST_TRACKER_RESET=0 "$PYTHON_BIN" retry_errors.py
  else
    log "Skipping retry_errors.py (--skip-retry)."
    if [[ -f "$filled_file" && ! -f "$fixed_file" ]]; then
      cp "$filled_file" "$fixed_file"
      log "Created ${fixed_file} from ${filled_file}"
    fi
  fi
  run_cmd "$PYTHON_BIN" normalize_patient_direction.py
  run_cmd "$PYTHON_BIN" cleanup_bracket_entries.py
else
  log "Skipping LLM steps (--no-llm)."
fi

if [[ "$RUN_R" == "1" ]]; then
  filled_file="$(provider_tagged_name taxon_article_matrix_filled.csv)"
  [[ -f "$filled_file" ]] || die "Missing ${filled_file} for R steps."
  run_cmd Rscript cleanTable.R
  run_cmd Rscript produceFinalTable.R
fi

cleanup_root_transient_files
log "Pipeline completed successfully."
