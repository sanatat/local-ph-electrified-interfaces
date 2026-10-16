"""Exploratory 1D total-alkalinity diffusion model for a cathode boundary layer.

This is deliberately a reduced-order proof-of-concept, NOT Nernst-Planck transport,
membrane electrodialysis, fitted experimental data, or validated electrochemistry.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

FARADAY = 96485.33212  # C/mol
PK1, PK2, PKW = 6.35, 10.33, 14.0  # illustrative ~25 C values
KA1, KA2, KW = 10**(-PK1), 10**(-PK2), 10**(-PKW)
KSP_CALCITE = 3.3e-9  # (mol/L)^2 ideal-activity screening ONLY


def carbonate_species(pH, dic_molm3):
    """Equilibrium carbonate species for fixed dissolved inorganic carbon (mol/m^3)."""
    h = 10.0**(-float(pH))
    den = h*h + KA1*h + KA1*KA2
    c_l = float(dic_molm3) / 1000.0
    return c_l * h*h/den, c_l*KA1*h/den, c_l*KA1*KA2/den


def alkalinity_from_ph(pH, dic_molm3):
    """Total alkalinity proxy mol/m^3, using only carbonate and water equilibria."""
    h = 10.0**(-float(pH))
    _, hco3_l, co3_l = carbonate_species(pH, dic_molm3)
    return 1000.0 * (hco3_l + 2*co3_l + KW/h - h)


def ph_from_alkalinity(ta_molm3, dic_molm3):
    return brentq(lambda ph: alkalinity_from_ph(ph, dic_molm3) - float(ta_molm3),
                  2.0, 13.99, xtol=1e-11)


def omega_calcite(pH, dic_molm3, ca_molm3=0.5):
    """Concentration-based saturation INDICATOR; not activity-corrected precipitation."""
    _, _, co3 = carbonate_species(pH, dic_molm3)
    return (ca_molm3 / 1000.0) * co3 / KSP_CALCITE


@dataclass(frozen=True)
class Case:
    dic_molm3: float = 2.0        # = 2 mM DIC
    bulk_ph: float = 7.5
    current_Am2: float = 0.15      # imposed current density, not measured flux
    oh_efficiency: float = 1.0    # OH equivalents/current: illustrative maximum
    d_eff_m2s: float = 1.2e-9     # effective TA diffusion only
    layer_m: float = 400e-6       # near-electrode stagnant layer thickness
    n_cells: int = 60
    time_s: float = 180.0
    ca_molm3: float = 0.5         # Ca concentration held fixed = 0.5 mM


def solve_case(case=Case(), times=(0., 30., 90., 180.)):
    """FV: dTA/dt = D d2TA/dx2, cathode molar OH-flux = eta*j/F.
    Boundary x=L held at bulk TA; DIC is held spatially constant.
    Return pH, TA and heuristic calcite supersaturation indicators.
    """
    if not (0.0 <= case.oh_efficiency <= 1.0):
        raise ValueError('OH-current efficiency must be in [0,1].')
    if case.current_Am2 < 0 or case.dic_molm3 < 0 or case.d_eff_m2s <= 0:
        raise ValueError('Invalid non-negative physical parameter.')
    n = case.n_cells
    dx = case.layer_m / n
    bulk_ta = alkalinity_from_ph(case.bulk_ph, case.dic_molm3)
    flux = case.oh_efficiency * case.current_Am2 / FARADAY
    # Finite-volume diffusion: left boundary specified flux; right boundary reservoir
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        mat[i, i] = -2.0
        if i > 0: mat[i, i-1] = 1.0
        if i < n-1: mat[i, i+1] = 1.0
    mat[0, 0] = -1.0  # zero diffusive outflow at wall; supplied flux in source
    mat[-1, -1] = -3.0  # Dirichlet reservoir at half a grid cell
    mat *= case.d_eff_m2s / dx**2
    source = np.zeros(n)
    source[0] = flux / dx
    source[-1] = 2 * case.d_eff_m2s * bulk_ta / dx**2

    ts = np.asarray(times, dtype=float)
    if ts.min() < 0 or ts.max() > case.time_s or not np.all(np.diff(ts) > 0):
        raise ValueError('times must be increasing, positive and within run.')
    sol = solve_ivp(lambda t, y: mat @ y + source, (0, case.time_s),
                    np.full(n, bulk_ta), method='BDF', t_eval=ts,
                    jac=lambda t, y: mat, rtol=1e-8, atol=1e-10)
    if not sol.success: raise RuntimeError(sol.message)
    ta = sol.y.T
    ph = np.array([[ph_from_alkalinity(v, case.dic_molm3) for v in row] for row in ta])
    omega = np.array([[omega_calcite(v, case.dic_molm3, case.ca_molm3) for v in row] for row in ph])
    return {'x_um': (np.arange(n) + 0.5) * dx * 1e6,
            'time_s': ts, 'alkalinity_molm3': ta, 'pH': ph,
            'omega_calcite_ideal': omega, 'OH_flux_molm2s': flux,
            'bulk_ta_molm3': bulk_ta}
