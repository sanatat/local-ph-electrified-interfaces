"""Compare fixed-pH PHREEQC aqueous speciation with Stage-3 ideal spatial profiles."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results_stage6'
COLS=('soln','ph_activity','si_calcite','ionic_strength','ca_free_molkg','co3_free_molkg','caco3_aq_molkg','cahco3_aq_molkg','caoh_aq_molkg')


def compare():
    inp=OUT/'input_stage6_profiles.csv'
    sel=OUT/'phreeqc_selected.tsv'
    if not inp.exists() or not sel.exists():
        raise FileNotFoundError('Generate input with stage6_spatial_prepare.py and run PHREEQC first.')
    model=pd.read_csv(inp)
    ph=pd.read_csv(sel,sep=r'\s+',engine='python')
    ph.columns=[str(c).strip().lower() for c in ph.columns]
    missing=set(COLS)-set(ph.columns)
    if missing:
        raise ValueError(f'Missing PHREEQC selected-output columns: {sorted(missing)}; present: {list(ph.columns)}')
    ph['soln']=pd.to_numeric(ph['soln'],errors='raise').astype(int)
    if ph['soln'].duplicated().any():raise ValueError('Duplicate PHREEQC solution ids')
    df=model.merge(ph,on='soln',validate='one_to_one',how='left',indicator=True)
    if not (df['_merge']=='both').all() or len(df)!=len(model):
        raise RuntimeError('PHREEQC output has missing or extra solution records')
    df.drop(columns=['_merge'],inplace=True)
    numeric_cols=['ph_activity','si_calcite','ca_free_molkg','co3_free_molkg','caco3_aq_molkg']
    if df[numeric_cols].isna().any().any():raise ValueError('NaN values in PHREEQC output')
    df['omega_phreeqc']=10.**df['si_calcite']
    df['ca_free_fraction']=df['ca_free_molkg']*1000/df['ca_total_mmolkg']
    df['caco3_aq_fraction']=df['caco3_aq_molkg']*1000/df['ca_total_mmolkg']
    df['omega_ratio_phreeqc_to_ideal']=df['omega_phreeqc']/df['omega_ideal']
    df['delta_ph']=df['ph_activity']-df['ph_stage3']
    # 'alkalinity' in selected output is eq/kgw; convert to meq/kgw (approximately model mmol/L as eq/L).
    alkcols=[c for c in ph.columns if c.startswith('alk')]
    if alkcols:
        df['phreeqc_alkalinity_meqkg']=df[alkcols[0]]*1000
        df['ta_difference_meqkg']=df['phreeqc_alkalinity_meqkg']-df['ta_stage3_mmolkg']
    df.to_csv(OUT/'spatial_comparison.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.7))
    for cur,g in df.groupby('current_A_m2',sort=True):
        if cur==0: continue
        axs[0].plot(g['distance_um'],g['omega_phreeqc'],label=f'{cur:g} A/m²')
        axs[1].plot(g['distance_um'],g['omega_ideal'],label=f'{cur:g} A/m²')
    for ax in axs:
        ax.axhline(1,ls='--',lw=1,color='black')
        ax.set(xlabel='Distance from cathode (µm)',ylabel='Calcite saturation ratio Ω')
        ax.grid(alpha=.3);ax.legend(fontsize=8)
    axs[0].set_title('PHREEQC: fixed-pH, nonideal speciation')
    axs[1].set_title('Stage 3: ideal speciation')
    fig.suptitle('One-way coupling, no precipitation or transport feedback')
    fig.tight_layout();fig.savefig(OUT/'spatial_calcite_omega.png',dpi=170);plt.close(fig)
    g=df.loc[df['current_A_m2']==0.5].sort_values('distance_um')
    fig,ax=plt.subplots(figsize=(8,4.5))
    ax.plot(g['distance_um'],g['ca_total_mmolkg'],label='Input total Ca')
    ax.plot(g['distance_um'],g['ca_free_molkg']*1000,label='Free Ca²⁺, PHREEQC')
    ax.plot(g['distance_um'],g['caco3_aq_molkg']*1000,label='CaCO₃(aq), PHREEQC')
    ax.set(xlabel='Distance from cathode (µm)',ylabel='Calcium species (mmol/kgw)',
           title='Dissolved calcium speciation at 0.5 A/m², 180 s')
    ax.legend();ax.grid(alpha=.3);fig.tight_layout();fig.savefig(OUT/'spatial_calcium_speciation.png',dpi=170);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5))
    error_cols=[c for c in df.columns if 'pct_err' in c or 'percent_error' in c]
    for cur,g in df.groupby('current_A_m2',sort=True):
        if error_cols:
            axs[0].plot(g['distance_um'],g[error_cols[0]],label=f'{cur:g} A/m²')
        else:
            axs[0].plot(g['distance_um'],g['delta_ph'],label=f'{cur:g} A/m²')
        if 'ta_difference_meqkg' in df.columns:
            axs[1].plot(g['distance_um'],g['ta_difference_meqkg'],label=f'{cur:g} A/m²')
    axs[0].set(xlabel='Distance (µm)',ylabel='Charge balance error (%)' if error_cols else 'ΔpH',title='Fixed-pH charge imbalance' if error_cols else 'pH fixed-input consistency')
    axs[1].set(xlabel='Distance (µm)',ylabel='Δ total alkalinity (meq/kgw)',title='Speciation vs Stage-3 alkalinity')
    for ax in axs:ax.grid(alpha=.3);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'spatial_charge_alkalinity_audit.png',dpi=170);plt.close(fig)
    summary={'n_spatial_solutions':len(df),'method':'One-way fixed-pH speciation, no precipitated solids',
             'per_current_near_wall':[],'caveats':['Model concentrations are approximate mmol/kgw',
             'PHREEQC pH fixed; chemical equilibrium does not feed back to transport',
             'Stage-3 alkalinity and PHREEQC alkalinity may differ',
             'Aqueous CaCO3 complexes are not solid calcite',
             'No experimental validation']}
    for cur,g in df.groupby('current_A_m2',sort=True):
        row=g.sort_values('distance_um').iloc[0]
        summary['per_current_near_wall'].append({
         'current_A_m2':float(cur),'distance_um':float(row['distance_um']),
         'model_ph':float(row['ph_stage3']), 'phreeqc_ph':float(row['ph_activity']),
         'omega_ideal':float(row['omega_ideal']), 'omega_phreeqc':float(row['omega_phreeqc']),
         'ca_free_fraction':float(row['ca_free_fraction']), 'caco3_aq_fraction':float(row['caco3_aq_fraction'])})
    (OUT/'stage6_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
    print('Saved plots and CSV in',OUT)
    return df

if __name__=='__main__':compare()
