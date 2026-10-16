"""Independent, limited Davies-activity check for Stage 3 local carbonate snapshots.

NOT PHREEQC. Davies approximation; no complexes; fixed soluble inventories.
Uses equilibrium electroneutrality to infer activity-based pH. Temperature 25 C.
"""
from dataclasses import dataclass
from scipy.optimize import brentq
import math

PK1, PK2, PKW = 6.35, 10.33, 14.00
K1, K2, KW = 10**(-PK1), 10**(-PK2), 10**(-PKW)
KSP=3.3e-9  # Same thermodynamic proxy as Stage 3, for transparent comparison

@dataclass
class Snapshot:
    label: str
    na_mM: float
    cl_mM: float
    ca_mM: float
    dic_mM: float
    ph_ideal: float
    omega_ideal: float
    time_s: float
    j_Am2: float


def davies_gamma(z,I):
    if I<0 or I>0.5:raise ValueError('Davies must not be used above 0.5 mol/kg (and is approximate already below).')
    sqrt=math.sqrt(I)
    return 10**(-0.509*z*z*(sqrt/(1+sqrt)-0.3*I))


def speciation_at_ph(s:Snapshot, ph_activity:float, ionic_strength:float):
    # Concentrations approximated as molarity (mol/L), not strict molality.
    g1=davies_gamma(1, ionic_strength);g2=davies_gamma(2, ionic_strength)
    ah=10**(-ph_activity)
    co2=s.dic_mM*1e-3/(1+K1/(ah*g1)+K1*K2/(ah**2*g2))
    hco3=K1*co2/(ah*g1)
    co3=K1*K2*co2/(ah*ah*g2)
    h=ah/g1;oh=KW/(ah*g1)
    na=s.na_mM*1e-3;cl=s.cl_mM*1e-3;ca=s.ca_mM*1e-3
    charge=na+2*ca+h-cl-hco3-2*co3-oh
    I=0.5*(na+4*ca+cl+h+hco3+4*co3+oh)
    omega=(g2*ca)*(g2*co3)/KSP
    return {'pH_activity':ph_activity,'CO2_M':co2,'HCO3_M':hco3,'CO3_M':co3,'H_M':h,'OH_M':oh,
            'charge_mol_L':charge,'ionic_strength_M':I,'gamma1':g1,'gamma2':g2,
            'omega_davies':omega,'si_davies':math.log10(omega) if omega>0 else float('-inf')}


def davies_reference(s, charge_balanced=True):
    if min(s.na_mM,s.cl_mM,s.ca_mM,s.dic_mM)<0: raise ValueError('negative dissolved totals')
    I=(s.na_mM+4*s.ca_mM+s.cl_mM)/2e3
    for k in range(50):
        if charge_balanced:
            def f(ph):return speciation_at_ph(s,ph,I)['charge_mol_L']
            if f(2)*f(13.9)>=0:
                raise ValueError('no charge-balance root for pH between 2 and 13.9')
            ph=brentq(f,2,13.9,xtol=1e-13)
        else: ph=s.ph_ideal
        ans=speciation_at_ph(s,ph,I)
        updated=ans['ionic_strength_M']
        if abs(updated-I)<1e-14:
            ans['iterations']=k+1;return ans
        I=0.5*(I+updated)
    raise RuntimeError('Davies ionic strength did not converge')
