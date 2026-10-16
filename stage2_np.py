"""Stage 2: dilute, unbuffered NaCl / NaOH / H+ 1-D electroneutral Nernst-Planck.

All concentrations: mol/m^3; electric current density: A/m^2; potentials: V.
Tracks chloride and total alkalinity A=[OH-]-[H+] in water; imposes instantaneous
water dissociation Kw and electroneutrality [Na+]=[Cl-]+A.

Imposes a constant ionic charge flux i/F=-j/F, consistent with a cathodic
OH-equivalent source J_A=j/F and zero Na/Cl electrode flux. The cathode water
reduction mechanism and faradaic yield are assumptions; no experimental data.
This model does *not* include carbonate buffering, membrane selectivity, EDL,
surface kinetics, flow, activity corrections, precipitation or gas evolution.
The dilute-solution assumption becomes questionable at high ionic strengths.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

F = 96485.33212
R = 8.314462618
T = 298.15
KAPPA = F/(R*T)
KW_M6 = 1.0e-8  # Kw=1e-14 (mol/L)^2 =1e-8 (mol/m^3)^2
D_NA = 1.33e-9
D_CL = 2.03e-9
D_H  = 9.31e-9
D_OH = 5.27e-9


@dataclass(frozen=True)
class NPCase:
    salt_molm3: float = 5.0             # NaCl=5 mM
    bulk_ph: float = 7.5
    cathodic_current_Am2: float = 0.5   # positive magnitude
    layer_m: float = 400e-6
    n_cells: int = 50
    time_s: float = 180.0


def water_species(alkalinity):
    """Ideal water equilibrium [OH]-[H]=A, [OH][H]=Kw; mol/m^3."""
    a = np.asarray(alkalinity, dtype=float)
    root = np.sqrt(a*a+4*KW_M6)
    # Stable for either strong basic or strong acid environments
    oh = np.where(a>=0, (a+root)/2, 2*KW_M6/(root-a))
    h = KW_M6/oh
    return oh, h


def initial_alkalinity(ph):
    oh=1000*10**(float(ph)-14)
    h=KW_M6/oh
    return oh-h


def ph_from_a(a):
    oh, _ = water_species(a)
    return 14+np.log10(oh/1000)


def _face_fluxes(cl_L,a_L,cl_R,a_R,face_dx,q):
    """Compute chloride, alkalinity and Na NP fluxes; q is mol charge/m2/s.

    q := J_Na+J_H-J_Cl-J_OH = J_Na-J_Cl-J_A.
    """
    cl_l=np.asarray(cl_L, dtype=float); a_l=np.asarray(a_L, dtype=float)
    cl_r=np.asarray(cl_R, dtype=float); a_r=np.asarray(a_R, dtype=float)
    dx=np.asarray(face_dx, dtype=float)
    oh_l,h_l=water_species(a_l); oh_r,h_r=water_species(a_r)
    na_l=cl_l+a_l; na_r=cl_r+a_r
    cl=0.5*(cl_l+cl_r); oh=0.5*(oh_l+oh_r)
    h=0.5*(h_l+h_r); na=0.5*(na_l+na_r)
    den = D_NA*na+D_CL*cl+D_H*h+D_OH*oh
    if np.any(den <= 0):
        raise ValueError('Nonpositive ionic conductivity (bad concentration).')
    # Sum_i z_i D_i * dc_i/dx, z=[+1,-1,+1,-1]
    charge_diff=(D_NA*(na_r-na_l)-D_CL*(cl_r-cl_l)
               +D_H*(h_r-h_l)-D_OH*(oh_r-oh_l))/dx
    grad_phi=-(q+charge_diff)/(KAPPA*den)
    # Separate diffusion / migration components for audit
    cl_dif=-D_CL*(cl_r-cl_l)/dx
    cl_mig=+D_CL*KAPPA*cl*grad_phi
    a_dif=(-D_OH*(oh_r-oh_l)+D_H*(h_r-h_l))/dx
    a_mig=KAPPA*(D_OH*oh+D_H*h)*grad_phi
    na_dif=-D_NA*(na_r-na_l)/dx
    na_mig=-D_NA*KAPPA*na*grad_phi
    cl_flux=cl_dif+cl_mig
    a_flux=a_dif+a_mig
    na_flux=na_dif+na_mig
    return {'cl':cl_flux, 'a':a_flux, 'na':na_flux, 'q':na_flux-cl_flux-a_flux,
            'grad_phi':grad_phi, 'a_diffusion':a_dif, 'a_migration':a_mig,
            'cl_diffusion':cl_dif, 'cl_migration':cl_mig}


def _jac_sparsity(n):
    """Each of two FV residuals depends on neighbors in both state variables."""
    jac=lil_matrix((2*n,2*n), dtype=int)
    for i in range(n):
        for k in range(max(0,i-1),min(n,i+2)):
            jac[i,k]=jac[i,n+k]=jac[n+i,k]=jac[n+i,n+k]=1
    return jac.tocsr()


def solve_np(case=NPCase(), times=(0.,30.,90.,180.)):
    if case.salt_molm3 <= 0 or not (5 <= case.bulk_ph <= 12):
        raise ValueError('Require positive background salt and basic/neutral bulk pH.')
    if case.cathodic_current_Am2 < 0 or case.layer_m <= 0 or case.n_cells < 5:
        raise ValueError('Invalid cathodic current, layer length or grid size.')
    ts=np.asarray(times, dtype=float)
    if ts.ndim != 1 or np.any(np.diff(ts)<=0) or ts[0]<0 or ts[-1]>case.time_s:
        raise ValueError('Times must increase and be within 0..time_s.')
    n=case.n_cells
    dx=case.layer_m/n
    a0=initial_alkalinity(case.bulk_ph)
    cl0=case.salt_molm3
    J=case.cathodic_current_Am2/F
    q=-J  # OH-equivalent cathode current, negative charge travels toward +x
    right_dx=np.full(n,dx)
    right_dx[-1]=dx/2

    def faces(y):
        cl=y[:n]; a=y[n:]
        # faces 1..n-1 between cells; face n: last cell -> reservoir
        f=_face_fluxes(cl, a,
                       np.r_[cl[1:], cl0], np.r_[a[1:], a0],
                       right_dx, q)
        # boundary face 0: no Cl crossing the metal; imposed alkalinity flux
        for k in ['cl','a','na','q','grad_phi','a_diffusion','a_migration','cl_diffusion','cl_migration']:
            f[k]=np.r_[np.nan if k not in ('cl','a','na','q') else
                         {'cl':0.0,'a':J,'na':0.0,'q':q}[k], f[k]]
        return f

    def rhs(t,y):
        f=faces(y)
        return np.r_[ (f['cl'][:-1]-f['cl'][1:])/dx,
                      (f['a'][:-1]-f['a'][1:])/dx ]

    y0=np.r_[np.full(n,cl0),np.full(n,a0)]
    sol=solve_ivp(rhs,(0.,case.time_s),y0,t_eval=ts,method='BDF',
                  jac_sparsity=_jac_sparsity(n),rtol=4e-8,atol=1e-9)
    if not sol.success:
        raise RuntimeError(sol.message)
    cl=sol.y[:n].T
    a=sol.y[n:].T
    na=cl+a
    if np.any((cl<=0)|(na<=0)):
        raise RuntimeError('Unphysical negative ions in computed output.')
    last_faces=faces(sol.y[:,-1])
    return {'x_um':(np.arange(n)+0.5)*dx*1e6,
            't_s':ts,'cl_molm3':cl,'alkalinity_molm3':a,
            'na_molm3':na,'pH':ph_from_a(a), 'flux':last_faces,
            'injected_alkalinity_flux_molm2s':J, 'ionic_charge_flux_molm2s':q,
            'initial_alkalinity_molm3':a0,
            'case':case}
