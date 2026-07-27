#!/bin/bash
set -euo pipefail

python run_harmonize.py \
  --bam-list example_bams_to_run.txt \
  --output-root /path/to/all_sv_calls \
  --summary-dir harmonized_summary \
  --figure-dir harmonized_summary/figures
