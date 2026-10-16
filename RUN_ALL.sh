#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python -m unittest discover -s tests -v
python run.py
python run_stage2.py
python run_stage3.py
python run_stage5.py
if [[ $# -ge 1 ]]; then
  export PATH="$HOME/.local/bin:$PATH"
  bash run_phreeqc.sh "$1"
  bash run_stage6.sh "$1"
else
  echo "Python-only runs finished. For PHREEQC, use: bash RUN_ALL.sh /path/to/phreeqc.dat"
fi
