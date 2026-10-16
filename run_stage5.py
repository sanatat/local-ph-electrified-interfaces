"""Stage 5 benchmark material generator. No external PHREEQC execution or data fitting."""
from pathlib import Path
from dataclasses import asdict
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from stage3_carbonate import CarbonateCase,solve_carbonate
from stage5_thermo import Snapshot,davies_reference
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results_stage5';OUT.mkdir(exist_ok=True)
PQC=ROOT/'phreeqc';PQC.mkdir(exist_ok=True)

def make_snapshots():
    snaps=[]
    for j in [0,0.15,0.5,1.0]:
        d=solve_carbonate(CarbonateCase(current_Am2=j),times=(0,180))
        # baseline and cathode center; 0.5 uses initial bulk in a separate snapshot below
        snaps.append(Snapshot(f'near_cathode_j{str(j).replace(".","p")}',float(d['na_molm3'][-1,0]),float(d['cl_molm3'][-1,0]),
                              float(d['ca_molm3'][-1,0]),float(d['dic_molm3'][-1,0]),float(d['pH'][-1,0]),
                              float(d['omega_ideal'][-1,0]),180.,j))
    d=solve_carbonate(CarbonateCase(current_Am2=0),times=(0,180))
    initial=Snapshot('bulk_reference',float(d['na_molm3'][0,0]),float(d['cl_molm3'][0,0]),
                     float(d['ca_molm3'][0,0]),float(d['dic_molm3'][0,0]),
                     float(d['pH'][0,0]),float(d['omega_ideal'][0,0]),0.,0.)
    return [initial]+snaps

def pqc_text(snaps,charge_balanced=True):
    mode='charge_balance' if charge_balanced else 'fixed_ph'
    head=f'''TITLE Stage 5: ideal Stage 3 inventories vs PHREEQC thermodynamics -- {mode}
# Approximation: mol/m3 numerically equals mmol/L, approximately mmol/kg water in dilute solutions.
# Na/Cl/Ca/C(4) totals from Stage 3 near-wall cell; no precipitated solids.
# Fixed-pH deck compares a common numerical pH but need not be charge balanced.
# Charge-balance deck uses pH charge (pH may change from Stage 3).
# DO NOT add equilibrium calcite phases, as this would force saturation and modify inventories.
SELECTED_OUTPUT 1
    -file results_stage5/phreeqc_{mode}_selected.tsv
    -reset false
    -solution true
    -high_precision true
USER_PUNCH 1
    -headings ph_activity si_calcite ionic_strength ca_free_molkg co3_free_molkg
    -start
    10 PUNCH -LA("H+"), SI("Calcite"), MU, MOL("Ca+2"), MOL("CO3-2")
    -end
END
'''
    for idx,s in enumerate(snaps,1):
        ph=f'{s.ph_ideal:.11f}' + (' charge' if charge_balanced else '')
        head+=f'''# {idx} {s.label} current={s.j_Am2} A/m2 time={s.time_s} s
SOLUTION {idx} {s.label}
    temp 25
    units mmol/kgw
    pH {ph}
    pe 4
    Na {s.na_mM:.12g}
    Cl {s.cl_mM:.12g}
    Ca {s.ca_mM:.12g}
    C(4) {s.dic_mM:.12g}
    -water 1
END
'''
    return head

def main():
    ss=make_snapshots()
    (OUT/'stage3_snapshots.json').write_text(json.dumps([asdict(s) for s in ss],indent=2))
    rows=[]
    for s in ss:
        cb=davies_reference(s,True)
        fp=davies_reference(s,False)
        rows.append(dict(snapshot=s.label,current_A_m2=s.j_Am2,model_pH_ideal=s.ph_ideal,
                        davies_pH_charge_balanced=cb['pH_activity'],delta_pH=cb['pH_activity']-s.ph_ideal,
                        ideal_omega=s.omega_ideal,omega_davies_charge_balanced=cb['omega_davies'],
                        omega_ratio_davies_over_ideal=cb['omega_davies']/s.omega_ideal,
                        omega_davies_fixed_ph=fp['omega_davies'],
                        gamma1=cb['gamma1'],gamma2=cb['gamma2'],ionic_strength_M=cb['ionic_strength_M'],
                        charge_residual_M=cb['charge_mol_L']))
    df=pd.DataFrame(rows)
    df.to_csv(OUT/'davies_vs_ideal.csv',index=False)
    for mode in (True,False):
        fn='charge_balance' if mode else 'fixed_ph'
        (PQC/f'benchmark_{fn}.pqi').write_text(pqc_text(ss,mode))
    fig,ax=plt.subplots(figsize=(8.2,4.7))
    currents=np.array([r['current_A_m2'] for r in rows if r['snapshot']!='bulk_reference'])
    ideal=np.array([r['ideal_omega'] for r in rows if r['snapshot']!='bulk_reference'])
    dav=np.array([r['omega_davies_charge_balanced'] for r in rows if r['snapshot']!='bulk_reference'])
    ax.plot(currents,ideal,'o-',label='Stage 3: ideal-concentration Omega')
    ax.plot(currents,dav,'s-',label='Stage 5: approximate Davies-activity Omega')
    ax.axhline(1,ls='--',c='gray',label='Omega = 1')
    ax.set(xlabel='Cathodic current density (A/m²)',ylabel='Calcite supersaturation indicator, Ω',title='Activity correction changes CaCO₃ screening predictions')
    ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'activity_effect.png',dpi=160);plt.close(fig)
    print(df[['snapshot','model_pH_ideal','davies_pH_charge_balanced','ideal_omega','omega_davies_charge_balanced']].to_string(index=False))
    print('IMPORTANT: Davies reference is approximate; PHREEQC decks GENERATED, not executed.')
if __name__=='__main__':main()
