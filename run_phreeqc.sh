#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if ! command -v phreeqc >/dev/null 2>&1; then
  echo 'PHREEQC executable not found. Install USGS PHREEQC v3 first.' >&2
  echo 'https://water.usgs.gov/water-resources/software/PHREEQC/' >&2
  exit 2
fi
DB="${1:-${PHREEQC_DATABASE:-}}"
if [[ -z "$DB" ]]; then
  for p in /usr/share/phreeqc/database/phreeqc.dat /usr/share/phreeqc/phreeqc.dat /usr/local/share/phreeqc/database/phreeqc.dat /usr/local/share/doc/phreeqc/database/phreeqc.dat /usr/local/share/doc/iphreeqc/database/phreeqc.dat; do
    if [[ -f "$p" ]]; then DB="$p"; break; fi
  done
fi
if [[ -z "$DB" || ! -f "$DB" ]]; then
  echo 'PHREEQC database not found. Supply path: bash run_phreeqc.sh /path/to/phreeqc.dat' >&2
  exit 3
fi
mkdir -p results_stage5
for mode in charge_balance fixed_ph; do
  echo "PHREEQC: $mode with database $DB"
  phreeqc "phreeqc/benchmark_${mode}.pqi" "results_stage5/phreeqc_${mode}_run.out" "$DB"
  [[ -s "results_stage5/phreeqc_${mode}_selected.tsv" ]] || { echo "No selected output for $mode" >&2; exit 4; }
done
python compare_phreeqc.py
