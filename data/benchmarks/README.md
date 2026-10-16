# Published benchmark transcription (not experimental data)

The CSV in this folder records a compact set of **five PHREEQC 3.8.6 results reported from an Ubuntu terminal session** in October 2026. It is included so the comparison of ideal carbonate chemistry, approximate Davies activities, and PHREEQC saturation ratios can be traced to numeric values instead of a screenshot.

It is **not an external measurement dataset**, and the PHREEQC executable was not rerun in the environment used to package this repository. Generate your own outputs using `python run_stage5.py` followed by `bash run_phreeqc.sh /path/to/phreeqc.dat` and inspect the database, totals, pH assumptions and charge-balance residuals.
