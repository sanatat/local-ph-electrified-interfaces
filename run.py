"""Run demonstration and produce plots + CSV tables. No experimental datasets used."""
from dataclasses import replace
from pathlib import Path
import csv
import matplotlib.pyplot as plt
import numpy as np
from model import Case, solve_case, ph_from_alkalinity, alkalinity_from_ph

OUT=Path(__file__).parent / 'results'
OUT.mkdir(exist_ok=True)
base=Case()
times=(0., 30., 90., 180.)
cases={'unbuffered':replace(base,dic_molm3=0.0),
       'carbonate_buffer':replace(base,dic_molm3=2.0)}
sols={name:solve_case(case, times=times) for name,case in cases.items()}

fig,ax=plt.subplots(figsize=(8,5))
for name,s in sols.items():
    for k,t in enumerate(s['time_s']):
        if t==0: continue
        ax.plot(s['x_um'],s['pH'][k],label=f'{name}: {t:.0f} s')
ax.axhline(base.bulk_ph,linestyle=':',label='bulk pH')
ax.set(xlabel='Distance from cathode (micrometers)',ylabel='Calculated pH',
       title='Illustrative near-cathode pH profiles (not validated)')
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'pH_profiles.png',dpi=170)
plt.close(fig)

# Simple screening across dissolved inorganic carbon and current density.
current=[0.025,0.05,0.1,0.2,0.4]
dic=[0.0,0.5,1.,2.,3.]
rows=[]
for c in dic:
    for j in current:
        result=solve_case(replace(base,dic_molm3=c,current_Am2=j),times=(0.,180.))
        local_ph=float(result['pH'][-1,0])
        rows.append({'DIC_mM':c,'current_A_per_m2':j,'near_wall_pH':local_ph,
                     'delta_pH':local_ph-base.bulk_ph,
                     'ideal_calcite_omega':float(result['omega_calcite_ideal'][-1,0])})
with (OUT/'parameter_sweep.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

fig,ax=plt.subplots(figsize=(8,4.8))
for c in dic:
    data=[r for r in rows if r['DIC_mM']==c]
    ax.plot([r['current_A_per_m2'] for r in data],
            [r['near_wall_pH'] for r in data],marker='o',label=f'DIC={c:g} mM')
ax.set(xlabel='Current density (A/m²)',ylabel='Predicted pH near cathode',
       title='Exploratory buffering–current sensitivity; 180 s')
ax.legend()
fig.tight_layout();fig.savefig(OUT/'buffering_sensitivity.png',dpi=170);plt.close(fig)

for name,s in sols.items():
    with (OUT/f'{name}_profiles.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['time_s','distance_um','pH','alkalinity_mol_per_m3','ideal_calcite_omega'])
        for ti,t in enumerate(s['time_s']):
            for i,x in enumerate(s['x_um']):
                wr.writerow([t,x,s['pH'][ti,i],s['alkalinity_molm3'][ti,i],s['omega_calcite_ideal'][ti,i]])
print('Files created under:',OUT)
for name,s in sols.items():
    print(name, 'pH near wall at 180 s:', round(float(s['pH'][-1,0]),3),
          'ideal CaCO3 index:',round(float(s['omega_calcite_ideal'][-1,0]),2))
