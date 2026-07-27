#!/bin/bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
ENV_NAME="${BREAKDANCE_ENV_NAME:-breakdance}"
CREATE_CONDA="${CREATE_CONDA:-0}"

if [[ "$CREATE_CONDA" == "1" ]]; then
  conda create -n "$ENV_NAME" python=3.11 -y
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
fi

"$PYTHON_BIN" -m pip install --upgrade pip
"$PYTHON_BIN" -m pip install -r requirements.txt

echo "Breakdance setup complete."
echo "Next steps:"
echo "  python run_harmonize.py --bam-list /path/to/bams_to_run.txt --summary-dir harmonized_summary"
