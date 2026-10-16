"""Stage 3 research prototype: coupled 1-D electroneutral Nernst-Planck + carbonate equilibrium.

Conservative FV balances for [Ca], [Cl], dissolved inorganic carbon (DIC) and
alkalinity (TA). [Na] is reconstructed via charge neutrality. Seven ionic species
(Na+, Ca2+, Cl-, H+, OH-, HCO3-, CO3(2-)) and neutral CO2(aq) participate in fluxes.
Water and carbonate reactions are assumed to equilibrate instantly. The model uses
ideal concentrations for equilibria, no Ca complexation, no convection, and no
mineral precipitation: calcite saturation Omega is a *screening indicator* only.
Results are illustrative and NOT experimentally validated.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

F=96485.33212
R=8.314462618
T=298.15
K=F/(R*T)
P_K1=6.35
P_K2=10.33
K1=10**(-P_K1)
K2=10**(-P_K2)
KW=1.0e-14
KSP_CALCITE=3.3e-9       # ideal, mol^2 L^-2; not activity corrected
# Diffusion coefficients (m2/s), illustrative 25 C infinite dilution approximations.
D={'na':1.33e-9,'ca':0.79e-9,'cl':2.03e-9,'h':9.31e-9,
   'oh':5.27e-9,'hco3':1.18e-9,'co3':0.92e-9,'co2':1.91e-9}
Z={'na':1,'ca':2,'cl':-1,'h':1,'oh':-1,'hco3':-1,'co3':-2,'co2':0}

@dataclass(frozen=True)
class CarbonateCase:
    cl_molm3:float=10.0
    ca_molm3:float=0.5
    dic_molm3:float=2.0
    bulk_ph:float=7.5
    current_Am2:float=0.5
    layer_m:float=400e-6
    n_cells:int=40
    time_s:float=180.0


def equilibrium(ta_molm3,dic_molm3):
    """Vectorized pH root from total alkalinity and DIC. Outputs in mol/m3.

    TA=HCO3 + 2 CO3 + OH - H, DIC=CO2 + HCO3 + CO3. pH bracket [2, 13.99].
    Negative DIC represents a physically inadmissible state and is rejected.
    """
    a,c=np.broadcast_arrays(np.asarray(ta_molm3,float), np.asarray(dic_molm3,float))
    if np.any(c < -1e-8):
        raise ValueError('Negative DIC; model failed positivity.')
    c=np.maximum(c,0.)
    lo=np.full(a.shape,2.,dtype=float)
    hi=np.full(a.shape,13.99,dtype=float)
    def spec(ph):
        h_l=10.**(-ph)
        denom=h_l*h_l+K1*h_l+K1*K2
        co2=c*h_l*h_l/denom
        hco3=c*K1*h_l/denom
        co3=c*K1*K2/denom
        h=1000*h_l
        oh=1000*KW/h_l
        tcalc=hco3+2*co3+oh-h
        return tcalc,co2,hco3,co3,h,oh
    tl=spec(lo)[0];th=spec(hi)[0]
    if np.any((a<tl)|(a>th)):
        raise ValueError('Total alkalinity is outside permitted pH bracket.')
    for _ in range(42):
        mid=(lo+hi)/2
        f=spec(mid)[0]
        hi=np.where(f>a,mid,hi)
        lo=np.where(f<=a,mid,lo)
    ph=(lo+hi)/2
    _,co2,hco3,co3,h,oh=spec(ph)
    return {'pH':ph,'co2':co2,'hco3':hco3,'co3':co3,'h':h,'oh':oh}


def total_alkalinity(ph,dic_molm3):
    h_l=10**(-np.asarray(ph,float))
    c=np.asarray(dic_molm3,float)
    den=h_l*h_l+K1*h_l+K1*K2
    return c*(K1*h_l+2*K1*K2)/den+1000*(KW/h_l-h_l)


def ideal_calcite_omega(ca_molm3,co3_molm3):
    return (np.asarray(ca_molm3)/1000)*(np.asarray(co3_molm3)/1000)/KSP_CALCITE


def face_flux(ca_L,cl_L,dic_L,ta_L,ca_R,cl_R,dic_R,ta_R,dist,q):
    """Species fluxes and conserved element/charge flux across FV faces."""
    eq_L=equilibrium(ta_L,dic_L);eq_R=equilibrium(ta_R,dic_R)
    aL=np.asarray(ta_L);aR=np.asarray(ta_R)
    speciesL={'ca':np.asarray(ca_L),'cl':np.asarray(cl_L),
              'na':np.asarray(cl_L)+aL-2*np.asarray(ca_L), **{x:eq_L[x] for x in ('h','oh','hco3','co3','co2')}}
    speciesR={'ca':np.asarray(ca_R),'cl':np.asarray(cl_R),
              'na':np.asarray(cl_R)+aR-2*np.asarray(ca_R), **{x:eq_R[x] for x in ('h','oh','hco3','co3','co2')}}
    den=0.;charge_diff=0.
    for k,z in Z.items():
        if z:
            den=den+z*z*D[k]*(speciesL[k]+speciesR[k])/2
            charge_diff=charge_diff+z*D[k]*(speciesR[k]-speciesL[k])/dist
    if np.any(den<=0):raise ValueError('Nonpositive ionic conductivity')
    gp=-(q+charge_diff)/(K*den)
    jf={}
    for k,z in Z.items():
        cc=(speciesL[k]+speciesR[k])/2
        jf[k]=-D[k]*(speciesR[k]-speciesL[k])/dist-z*D[k]*K*cc*gp
    jt=jf['hco3']+2*jf['co3']+jf['oh']-jf['h']
    jd=jf['co2']+jf['hco3']+jf['co3']
    qq=jf['na']+2*jf['ca']-jf['cl']-jt
    return {'ca':jf['ca'],'cl':jf['cl'],'dic':jd,'ta':jt,
            'na':jf['na'],'q':qq,'grad_phi':gp,'species_flux':jf}


def jac_sparsity(n):
    jac=lil_matrix((4*n,4*n),dtype=int)
    for b in range(4):
        for i in range(n):
            for bb in range(4):
                for k in range(max(0,i-1),min(n,i+2)):
                    jac[b*n+i,bb*n+k]=1
    return jac.tocsr()


def solve_carbonate(case=CarbonateCase(), times=(0.,30.,90.,180.)):
    if case.n_cells<5 or case.layer_m<=0 or case.time_s<=0 or case.current_Am2<0:
        raise ValueError('Invalid grid, domain, time or current')
    if min(case.cl_molm3,case.ca_molm3,case.dic_molm3)<0 or not 3<=case.bulk_ph<=12:
        raise ValueError('Invalid electrolyte or bulk pH')
    ts=np.asarray(times,float)
    if ts.ndim!=1 or len(ts)==0 or ts[0]<0 or ts[-1]>case.time_s or np.any(np.diff(ts)<=0):
        raise ValueError('Invalid output times')
    n=case.n_cells;dx=case.layer_m/n
    ta0=float(total_alkalinity(case.bulk_ph,case.dic_molm3))
    na0=case.cl_molm3+ta0-2*case.ca_molm3
    if na0<=0:raise ValueError('Bulk electrolyte is not electro-neutral with positive Na.')
    values={'ca':case.ca_molm3,'cl':case.cl_molm3,
            'dic':case.dic_molm3,'ta':ta0}
    j=case.current_Am2/F;q=-j
    distances=np.full(n,dx);distances[-1]=dx/2
    def get_faces(y):
        ca,cl,dic,ta=np.split(y,4)
        f=face_flux(ca,cl,dic,ta,
                    np.r_[ca[1:],values['ca']],np.r_[cl[1:],values['cl']],
                    np.r_[dic[1:],values['dic']],np.r_[ta[1:],values['ta']],distances,q)
        result={k:np.r_[j if k=='ta' else 0.,f[k]] for k in ('ca','cl','dic','ta')}
        result['q']=np.r_[q,f['q']]
        result['na']=np.r_[0.,f['na']]
        result['grad_phi']=f['grad_phi']
        return result
    def rhs(t,y):
        face=get_faces(y)
        return np.concatenate([(face[k][:-1]-face[k][1:])/dx for k in ('ca','cl','dic','ta')])
    y0=np.concatenate([np.full(n,values[k]) for k in ('ca','cl','dic','ta')])
    sol=solve_ivp(rhs,(0,case.time_s),y0,t_eval=ts,method='BDF',rtol=3e-7,atol=2e-9,
                  jac_sparsity=jac_sparsity(n),max_step=15.)
    if not sol.success:raise RuntimeError(sol.message)
    ca,cl,dic,ta=np.split(sol.y.T,4,axis=1)
    if np.any(ca<-1e-8) or np.any(cl<-1e-8) or np.any(dic<-1e-8):
        raise RuntimeError('Negative species inventory in computed outputs.')
    eq=equilibrium(ta,dic)
    na=cl+ta-2*ca
    if np.any(na<=0):raise RuntimeError('Nonpositive Na from electroneutral constraint')
    ff=get_faces(sol.y[:,-1])
    return {'x_um':(np.arange(n)+0.5)*dx*1e6,'t_s':ts,
            'ca_molm3':ca,'cl_molm3':cl,'dic_molm3':dic,'ta_molm3':ta,
            'na_molm3':na,'pH':eq['pH'],'co2_molm3':eq['co2'],
            'hco3_molm3':eq['hco3'],'co3_molm3':eq['co3'],
            'omega_ideal':ideal_calcite_omega(ca,eq['co3']),
            'face':ff,'case':case,'bulk_ta_molm3':ta0,
            'bulk_na_molm3':na0}
