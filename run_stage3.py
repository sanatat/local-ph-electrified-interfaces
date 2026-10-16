"""Generate Stage 3 carbonate/pH/calcite-screening research outputs."""
import csv,json
from dataclasses import replace
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from stage3_carbonate import CarbonateCase,solve_carbonate
HERE=Path(__file__).resolve().parent
OUT=HERE/'results_stage3';OUT.mkdir(exist_ok=True)
base=CarbonateCase()
s=solve_carbonate(base)
fig,axs=plt.subplots(1,2,figsize=(12,4.6))
for k,t in enumerate(s['t_s']):
    axs[0].plot(np.r_[s['x_um'],400],np.r_[s['pH'][k],base.bulk_ph],label=f'{t:g} s')
    axs[1].plot(np.r_[s['x_um'],400],np.r_[s['omega_ideal'][k],s['omega_ideal'][0,-1]],label=f'{t:g} s')
axs[0].set(xlabel='Distance from cathode (µm)',ylabel='pH',title='Carbonate-buffered near-cathode pH')
axs[1].axhline(1.,ls='--',color='black',label='Ω=1 threshold (ideal)')
axs[1].set(xlabel='Distance from cathode (µm)',ylabel='Ideal calcite saturation Ω',title='CaCO₃ thermodynamic screening — NOT precipitation')
for ax in axs:ax.legend(fontsize=8);ax.grid(alpha=.25)
fig.suptitle('Stage 3: carbonate NP, ideal chemistry, NO experimental validation')
fig.tight_layout();fig.savefig(OUT/'carbonate_ph_and_omega.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4.5))
summary=[]
for jj in [0.,0.15,0.5,1.0]:
    c=replace(base,current_Am2=jj)
    ss=solve_carbonate(c,times=(0.,180.))
    ax.plot(np.r_[ss['x_um'],400],np.r_[ss['pH'][-1],base.bulk_ph],label=f'{jj:g} A/m²')
    summary.append({'current_A_m2':jj,'pH_near_cathode_180s':float(ss['pH'][-1,0]),
                    'max_omega_ideal_180s':float(ss['omega_ideal'][-1].max())})
ax.set(xlabel='Distance from cathode (µm)',ylabel='Calculated pH',title='Carbonate-buffered pH sensitivity to cathodic current (180 s)')
ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(OUT/'carbonate_current_sensitivity.png',dpi=170);plt.close(fig)
with (OUT/'carbonate_profiles.csv').open('w',newline='') as file:
    wr=csv.writer(file)
    wr.writerow(['time_s','distance_um','pH','carbonate_molm3','bicarbonate_molm3','DIC_molm3','TA_molm3','Ca_molm3','calcite_Omega_ideal'])
    for k,t in enumerate(s['t_s']):
        for i,x in enumerate(s['x_um']):
            wr.writerow([float(t),float(x),float(s['pH'][k,i]),float(s['co3_molm3'][k,i]),
                         float(s['hco3_molm3'][k,i]),float(s['dic_molm3'][k,i]),
                         float(s['ta_molm3'][k,i]),float(s['ca_molm3'][k,i]),float(s['omega_ideal'][k,i])])
qerr=float(np.max(np.abs(s['face']['q'][1:]+base.current_Am2/96485.33212)))
result={'interpretation':'Ideal theoretical screening; not deposition prediction or experimental validation',
        'case':vars(base),'current_sensitivity':summary,
        'charge_flux_closure_max_abs_mol_m2_s':qerr,
        'bulk_omega_ideal':float(s['omega_ideal'][0,0]),
        'known_missing':['activities and ion-pair complexes','calcite nucleation/growth or precipitation sink',
                         'membrane-selective transport','fluid flow','electrode kinetics','experimental validation']}
(OUT/'stage3_summary.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
