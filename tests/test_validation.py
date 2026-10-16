"""Stronger validation tests than the original illustrative four behavior checks."""
import unittest
from dataclasses import replace
import numpy as np
from model import Case, solve_case, ph_from_alkalinity
from validation import analytic_alkalinity, numerical_surface_ta, grid_study, numeric_balance

class ScientificChecks(unittest.TestCase):
    def test_analytic_benchmark(self):
        c=Case(n_cells=80)
        s=solve_case(c,times=(0.,30.,90.,180.))
        for row,t in enumerate(s['time_s'][1:],start=1):
            ref=analytic_alkalinity(s['x_um']*1e-6,t,c)
            self.assertLess(float(np.max(abs(s['alkalinity_molm3'][row]-ref))), 0.00005)

    def test_extrapolated_wall_converges(self):
        r=grid_study(sizes=(20,40,80))
        self.assertGreater(r[0]['wall_pH_abs_error'],r[-1]['wall_pH_abs_error'])
        self.assertLess(r[-1]['wall_pH_abs_error'],0.00001)

    def test_diffusion_mass_balance(self):
        c=Case(n_cells=60)
        s=solve_case(c,times=tuple(np.linspace(0.,180.,361)))
        bal=numeric_balance(s,c)
        self.assertLess(bal['relative_error'], 1e-4)

    def test_linear_alkalinity_response_to_flux(self):
        c=Case()
        a=solve_case(c,times=(0.,180.))
        b=solve_case(replace(c,current_Am2=2*c.current_Am2),times=(0.,180.))
        da=a['alkalinity_molm3'][-1]-a['bulk_ta_molm3']
        db=b['alkalinity_molm3'][-1]-b['bulk_ta_molm3']
        np.testing.assert_allclose(db,2*da,rtol=1e-6,atol=1e-9)

if __name__=='__main__': unittest.main()
