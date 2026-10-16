#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
DB="${1:-${PHREEQC_DATABASE:-}}"
if ! command -v phreeqc >/dev/null 2>&1; then
  echo 'phreeqc not found in PATH; run export PATH="$HOME/.local/bin:$PATH"' >&2
  exit 2
fi
if [[ -z "$DB" || ! -f "$DB" ]]; then
  echo 'Usage: bash run_stage6.sh /actual/path/to/phreeqc.dat' >&2
  exit 3
fi
python stage6_spatial_prepare.py
phreeqc phreeqc/spatial_fixed_ph.pqi results_stage6/spatial_phreeqc.out "$DB"
[[ -s results_stage6/phreeqc_selected.tsv ]] || { echo 'Missing selected output' >&2; exit 4; }
python stage6_spatial_compare.py
