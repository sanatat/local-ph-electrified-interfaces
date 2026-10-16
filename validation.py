"""Independent numerical checks of the baseline 1D alkalinity model.

Analytical benchmark solves constant-coefficient diffusion on [0,L]
with constant inward alkalinity flux at x=0, fixed reservoir TA at x=L,
and initially uniform TA. This verifies the PDE solver numerically; it
cannot establish physical accuracy of the chosen electrochemical assumptions.
"""
from dataclasses import dataclass, replace
import numpy as np
from model import Case, solve_case, ph_from_alkalinity, alkalinity_from_ph, FARADAY


def analytic_alkalinity(x_m, t_s, case=Case(), modes=250):
    """Exact cosine-series solution for TA, for t>0, mol/m^3.

    TA(x,t) = TA_b + J(L-x)/D
        - 2J/(LD) * sum_n cos(k_n*x) * exp(-D*k_n^2*t) / k_n^2
    k_n=(n+1/2)*pi/L, n=0,1,...
    At t=0, TA is exactly bulk TA.
    """
    x = np.asarray(x_m, dtype=float)
    if np.any((x < 0) | (x > case.layer_m)) or t_s < 0:
        raise ValueError('x must be within 0..L, time >= 0')
    bulk = alkalinity_from_ph(case.bulk_ph, case.dic_molm3)
    if t_s == 0 or case.current_Am2 == 0 or case.oh_efficiency == 0:
        return np.full_like(x, bulk, dtype=float)
    J = case.oh_efficiency * case.current_Am2 / FARADAY
    n = np.arange(modes, dtype=float)
    k = (n + 0.5) * np.pi / case.layer_m
    transient = np.sum(np.cos(np.outer(x.ravel(), k)) *
                       (np.exp(-case.d_eff_m2s * k*k * t_s) / (k*k)), axis=1)
    return (bulk + J*(case.layer_m-x)/case.d_eff_m2s
            - 2*J/(case.layer_m*case.d_eff_m2s)*transient.reshape(x.shape))


def numerical_surface_ta(solution, case):
    """Second-order wall extrapolation from FV cell-centre using known gradient."""
    dx = case.layer_m / case.n_cells
    return solution['alkalinity_molm3'][:, 0] + solution['OH_flux_molm2s'] * dx / (2*case.d_eff_m2s)


def numeric_balance(solution, case):
    """Verify discrete conservation: d M/dt=J-(2 D/dx)*(TA_last-TA_bulk).

    Uses quadrature over a fine output time series; errors include finite dt.
    """
    from scipy.integrate import cumulative_trapezoid
    dx = case.layer_m / case.n_cells
    ts = solution['time_s']
    T = solution['alkalinity_molm3']
    M = np.sum(T-solution['bulk_ta_molm3'], axis=1) * dx
    flux_out = 2*case.d_eff_m2s*(T[:,-1]-solution['bulk_ta_molm3'])/dx
    net = solution['OH_flux_molm2s']-flux_out
    accounted = cumulative_trapezoid(net, ts, initial=0)
    residual = M - M[0] - accounted
    scale = max(solution['OH_flux_molm2s']*max(ts), 1e-14)
    return {'residual':residual,'relative_error':float(np.max(np.abs(residual))/scale),
            'M':M,'flux_out':flux_out}


def grid_study(base=Case(), sizes=(20,40,80,160), final_time=180.):
    rows=[]
    for n in sizes:
        c=replace(base,n_cells=n,time_s=final_time)
        s=solve_case(c, times=(0.,final_time))
        x=s['x_um']/1e6
        analytic=analytic_alkalinity(x,final_time,c)
        cell_linf=float(np.max(np.abs(s['alkalinity_molm3'][-1]-analytic)))
        wall_ta=float(numerical_surface_ta(s,c)[-1])
        wall_ref=float(analytic_alkalinity(np.array([0.]),final_time,c)[0])
        rows.append({'n_cells':n, 'dx_um':c.layer_m/n*1e6,
                     'ta_linf_molm3':cell_linf,
                     'wall_pH':ph_from_alkalinity(wall_ta,c.dic_molm3),
                     'analytic_wall_pH':ph_from_alkalinity(wall_ref,c.dic_molm3),
                     'wall_pH_abs_error':abs(ph_from_alkalinity(wall_ta,c.dic_molm3) -
                                              ph_from_alkalinity(wall_ref,c.dic_molm3))})
    return rows
