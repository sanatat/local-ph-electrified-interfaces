# Reproducing the calculations on Ubuntu or Linux

The aim of this guide is to help someone else run the same code and inspect the assumptions. **Passing tests does not mean the chemistry is experimentally validated.**

## Python setup

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

The current version contains **33** Python consistency tests. They do not require PHREEQC.

## Run the Python models

```bash
python run.py
python verification_run.py
python run_stage2.py
python run_stage3.py
python run_stage5.py
```

Plots and tables appear under local `results/` and `results_stage*/` folders. Git ignores them so that stale generated outputs are not accidentally published. Selected illustration figures are included in `figures/`.

## Install and locate PHREEQC

The equilibrium comparison needs the standalone [USGS PHREEQC version 3](https://www.usgs.gov/software/phreeqc-version-3/) application and its database. The original Ubuntu installation used PHREEQC 3.8.6; it was installed under `$HOME/.local/bin`, with `phreeqc.dat` under `$HOME/.local/share/doc/phreeqc/database/`.

Your installation location may be different. Check:

```bash
export PATH="$HOME/.local/bin:$PATH"
command -v phreeqc
find "$HOME/.local" -type f -name phreeqc.dat
```

If you run `phreeqc` with no arguments, it starts prompting for input files. Exit that prompt with `Ctrl+C` before entering Linux shell commands.

## Re-run the five-point speciation benchmark

```bash
DB="$HOME/.local/share/doc/phreeqc/database/phreeqc.dat"  # change if needed
python run_stage5.py
bash run_phreeqc.sh "$DB"
cat results_stage5/actual_phreeqc_comparison.csv
```

The fixed-pH and charge-balanced calculations should be interpreted separately. The latter deliberately changes pH to satisfy electroneutrality.

## Generate 160 spatial PHREEQC solutions

```bash
bash run_stage6.sh "$DB"
cat results_stage6/stage6_summary.json
```

Outputs include `results_stage6/spatial_comparison.csv` and figures named `spatial_calcite_omega.png`, `spatial_calcium_speciation.png`, and `spatial_charge_alkalinity_audit.png`. Inspect `results_stage6/spatial_phreeqc.out` for PHREEQC warnings before interpreting the chemistry.

If `phreeqc not found in PATH` appears, repeat the `export PATH=...` command. If `ModuleNotFoundError: pandas` appears, run `python -m pip install -r requirements.txt` inside the active virtual environment.

## Before publishing or citing the results

Check the PHREEQC database and version, Python version, package versions, and input totals. Review the generated charge/alkalinity differences rather than forcing them to zero. Never describe these calculations as measurements or include other researchers' copyrighted figures. Public repositories should exclude private applications, contact details, credentials and unlicensed data.
