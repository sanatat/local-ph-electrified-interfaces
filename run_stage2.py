"""Produce Stage-2 electrodiffusion diagnostics, no experimental validation."""
from pathlib import Path
from dataclasses import replace
import csv,json
import numpy as np
import matplotlib.pyplot as plt
from stage2_np import NPCase, solve_np, F

HERE=Path(__file__).resolve().parent
OUT=HERE/'results_stage2'; OUT.mkdir(exist_ok=True)
base=NPCase()
sol=solve_np(base, times=(0.,30.,90.,180.))
fig,axes=plt.subplots(1,2,figsize=(12,4.3))
for k,t in enumerate(sol['t_s']):
    axes[0].plot(np.r_[sol['x_um'],base.layer_m*1e6],
                 np.r_[sol['pH'][k],base.bulk_ph],label=f'{t:g} s')
axes[0].axhline(base.bulk_ph, color='gray',ls=':',label='reservoir pH')
axes[0].set(xlabel='Distance from cathode (µm)',ylabel='pH',title='Near-cathode pH — unbuffered NaCl')
axes[0].legend(fontsize=8);axes[0].grid(alpha=.3)
axes[0].set_xlim(0,base.layer_m*1e6)
flux=sol['flux']; xf=np.arange(1,base.n_cells+1)*base.layer_m/base.n_cells*1e6
# exact coordinates of each flux face (last face lies at reservoir at L)
axes[1].plot(xf,flux['a_diffusion'][1:],label='OH/H alkalinity diffusion')
axes[1].plot(xf,flux['a_migration'][1:],label='OH/H alkalinity migration')
axes[1].plot(xf,flux['a'][1:],label='Total alkalinity flux',ls='--')
axes[1].set(xlabel='Distance from cathode (µm)',ylabel='mol / (m² s)',title='Nernst–Planck fluxes at 180 s')
axes[1].legend(fontsize=7);axes[1].grid(alpha=.3)
fig.suptitle('Stage 2: electroneutral Nernst–Planck — idealized, NOT experimentally validated')
fig.tight_layout(); fig.savefig(OUT/'np_pH_and_flux.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(7.8,4.4))
for density in [0.0,0.15,0.5,1.0]:
    c=replace(base,cathodic_current_Am2=density)
    s=solve_np(c,times=(0.,180.))
    ax.plot(np.r_[s['x_um'],c.layer_m*1e6],
            np.r_[s['pH'][-1],c.bulk_ph],label=f'j={density:g} A/m²')
ax.set(xlabel='Distance from cathode (µm)',ylabel='Calculated pH',
       title='Sensitivity to applied cathodic current: unbuffered NaCl (180 s)')
ax.legend();ax.grid(alpha=.3);fig.tight_layout();fig.savefig(OUT/'np_current_sensitivity.png',dpi=170);plt.close(fig)
with (OUT/'np_profiles.csv').open('w',newline='') as fp:
    w=csv.writer(fp);w.writerow(['time_s','distance_um','pH','alkalinity_molm3','cl_molm3','na_molm3'])
    for k,t in enumerate(sol['t_s']):
        for i,x in enumerate(sol['x_um']):
            w.writerow([float(t),float(x),float(sol['pH'][k,i]),float(sol['alkalinity_molm3'][k,i]),
                        float(sol['cl_molm3'][k,i]),float(sol['na_molm3'][k,i])])
# check integrated charge flux across all internal faces is imposed current
max_q_err=float(np.nanmax(np.abs(sol['flux']['q'][1:]+base.cathodic_current_Am2/F)))
summary={'status':'Demonstration of electroneutral Nernst–Planck, not experimental validation',
         'case':vars(base),'wall_cell_center_pH_180s':float(sol['pH'][-1,0]),
         'last_cell_center_pH_180s':float(sol['pH'][-1,-1]),
         'reservoir_pH_at_boundary':base.bulk_ph,
         'max_internal_charge_flux_absolute_error_molm2s':max_q_err,
         'known_missing_physics':['carbonate and other buffers','membranes','faradaic efficiency',
                                   'double layer','fluid flow','nonideal activity','scaling','validation with measured pH']}
(OUT/'stage2_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
