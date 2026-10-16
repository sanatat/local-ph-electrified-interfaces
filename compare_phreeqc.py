"""Read *real* PHREEQC USER_PUNCH selected outputs after external execution."""
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parent
DIR=ROOT/'results_stage5'
SNAPS=pd.read_json(DIR/'stage3_snapshots.json')

def read_user_punch(filename):
    if not filename.is_file(): raise FileNotFoundError(f'PHREEQC has NOT been run: {filename}. First run: bash run_phreeqc.sh /path/to/phreeqc.dat')
    df=pd.read_csv(filename,sep=r'\s+',engine='python')
    df.columns=[str(c).strip().lower() for c in df.columns]
    needed=['ph_activity','si_calcite','ionic_strength','ca_free_molkg','co3_free_molkg']
    if any(c not in df.columns for c in needed):
        raise ValueError(f'Expected PHREEQC USER_PUNCH headings {needed}. Found {list(df.columns)}; check selected output and versions.')
    df=df.loc[:,needed].apply(pd.to_numeric,errors='coerce').dropna()
    if len(df)!=len(SNAPS):
        raise ValueError(f'{len(df)} PHREEQC rows vs {len(SNAPS)} input cases; inspect selected output for duplicates or missing rows.')
    return df.reset_index(drop=True)

def main():
    comp=pd.read_csv(DIR/'davies_vs_ideal.csv')
    for mode in ('charge_balance','fixed_ph'):
        ref=read_user_punch(DIR/f'phreeqc_{mode}_selected.tsv')
        comp[f'phreeqc_{mode}_pH']=ref['ph_activity']
        comp[f'phreeqc_{mode}_SI_calcite']=ref['si_calcite']
        comp[f'phreeqc_{mode}_Omega']=10**ref['si_calcite']
        comp[f'phreeqc_{mode}_ionic_strength']=ref['ionic_strength']
    comp['phreeqc_charge_balance_delta_pH']=comp.phreeqc_charge_balance_pH-comp.model_pH_ideal
    comp['phreeqc_charge_balance_omega_ratio']=comp.phreeqc_charge_balance_Omega/comp.ideal_omega
    dest=DIR/'actual_phreeqc_comparison.csv';comp.to_csv(dest,index=False)
    print(comp[['snapshot','model_pH_ideal','davies_pH_charge_balanced','phreeqc_charge_balance_pH','ideal_omega','omega_davies_charge_balanced','phreeqc_charge_balance_Omega']].to_string(index=False))
    print('Wrote',dest)
if __name__=='__main__': main()
