#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-chai-pipeline}"
PY_VER="${PY_VER:-3.11}"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
  log "ERROR: $*"
  exit 1
}

if ! command -v conda >/dev/null 2>&1; then
  die "conda not found. Install Miniconda/Anaconda first."
fi

# Ensure conda shell functions are available in this script.
eval "$(conda shell.bash hook)"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  log "Conda env already exists: $ENV_NAME"
else
  log "Creating conda env: $ENV_NAME (python=$PY_VER)"
  conda create -y -n "$ENV_NAME" "python=${PY_VER}" pip
fi

log "Installing Python packages..."
conda run -n "$ENV_NAME" python -m pip install --upgrade pip
conda run -n "$ENV_NAME" python -m pip install \
  pandas \
  requests \
  biopython \
  openai \
  python-dotenv \
  google-genai

log "Installing R + R packages..."
conda install -y -n "$ENV_NAME" -c conda-forge \
  r-base \
  r-dplyr \
  r-tidyr \
  r-stringr \
  r-readr

log "Running import checks..."
conda run -n "$ENV_NAME" python - <<'PY'
import importlib
mods = ["pandas", "requests", "Bio", "openai", "dotenv", "google.genai"]
missing = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception:
        missing.append(m)
if missing:
    raise SystemExit(f"Missing modules: {missing}")
print("Python module check: OK")
PY

conda run -n "$ENV_NAME" Rscript -e "library(dplyr); library(tidyr); library(stringr); library(readr); cat('R package check: OK\\n')"

log "Setup complete."
log "Activate with: conda activate ${ENV_NAME}"
log "Then run: ./run_pipeline.sh --python-only"
