#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workers="${1:-4}"

python "${repo_root}/scripts/validate_inputs.py"
python "${repo_root}/scripts/run_enrichment.py" \
  --resource all \
  --workers "${workers}" \
  --quick-input-check

