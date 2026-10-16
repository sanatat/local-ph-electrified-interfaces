"""Generate numerical verification plots and a machine-readable QA summary."""
import json, csv
from dataclasses import replace
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from model import Case, solve_case, ph_from_alkalinity
from validation import (analytic_alkalinity,grid_study,numeric_balance,
                        numerical_surface_ta)

OUT=Path(__file__).parent/'results';OUT.mkdir(exist_ok=True)
base=Case()
rows=grid_study()
with (OUT/'grid_verification.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

fig,axes=plt.subplots(1,2,figsize=(10,4.1))
axes[0].loglog([r['dx_um'] for r in rows],[r['wall_pH_abs_error'] for r in rows],marker='o')
axes[0].set(xlabel='Grid spacing (micrometers)', ylabel='Absolute pH error at wall',
            title='Wall pH vs analytic benchmark')
axes[0].grid(alpha=0.25,which='both')
axes[1].plot([r['n_cells'] for r in rows],[r['wall_pH'] for r in rows],marker='o',label='Numerical wall')
axes[1].axhline(rows[0]['analytic_wall_pH'],linestyle='--',label='Analytical reference')
axes[1].set(xlabel='Number of finite volumes',ylabel='pH at x=0 (180 s)',
            title='Grid convergence')
axes[1].legend(fontsize=8);axes[1].grid(alpha=.25)
fig.suptitle('Numerical verification only — not experimental validation',fontsize=12)
fig.tight_layout();fig.savefig(OUT/'numerical_verification.png',dpi=180);plt.close(fig)

c=replace(base,n_cells=80)
times=tuple(np.linspace(0.,180.,361));s=solve_case(c,times=times)
balance=numeric_balance(s,c)
refwall = analytic_alkalinity(np.array([0.]),180,c)[0]
wallta = numerical_surface_ta(s,c)[-1]

# Heatmap is a *calculated* ideal screening index, not precipitation prediction.
s2=solve_case(base,times=(0.,180.))
summary={
  'status':'Verified against analytic diffusion solution; not experimentally validated',
  'model':'1D effective alkalinity diffusion + fixed-DIC instantaneous carbonate equilibrium',
  'settings':vars(base),
  'analytical_wall_pH_at_180s':float(ph_from_alkalinity(refwall,c.dic_molm3)),
  'numeric_wall_pH_at_180s_N80':float(ph_from_alkalinity(wallta,c.dic_molm3)),
  'discrete_mass_balance_relative_error_max':balance['relative_error'],
  'grid_verification':rows,
  'known_missing_physics':['electromigration','membrane-selectivity','different ionic diffusivities',
  'spatial DIC transport','convective flow','activity corrections','CaCO3 precipitation kinetics',
  'electrode Faradaic efficiency measurements']}
(OUT/'verification_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
