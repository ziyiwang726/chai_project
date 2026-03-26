#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Streamlined entrypoint for interrupted workflows:
# starts from taxa_ids_filtered.csv and runs downstream steps.
exec ./run_pipeline.sh --from-filtered-ids "$@"
