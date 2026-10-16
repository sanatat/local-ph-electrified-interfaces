"""Independent chemical consistency and FV numerical checks for carbonate stage."""
import unittest
from dataclasses import replace
import numpy as np
from stage3_carbonate import (
    CarbonateCase, solve_carbonate, equilibrium, total_alkalinity,
    face_flux, ideal_calcite_omega, F)

class CarbonateChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case=CarbonateCase(n_cells=25)
        cls.sol=solve_carbonate(cls.case,times=(0.,60.,120.,180.))

    def test_equilibrium_mass_and_alkalinity(self):
        s=self.sol
        chem=equilibrium(s['ta_molm3'],s['dic_molm3'])
        np.testing.assert_allclose(chem['co2']+chem['hco3']+chem['co3'],s['dic_molm3'],rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(chem['hco3']+2*chem['co3']+chem['oh']-chem['h'],s['ta_molm3'],rtol=0,atol=3e-9)
        np.testing.assert_allclose(chem['h']*chem['oh'],1e-8,rtol=1e-11)

    def test_electroneutrality_all_species(self):
        s=self.sol
        e=equilibrium(s['ta_molm3'],s['dic_molm3'])
        charge=s['na_molm3']+2*s['ca_molm3']+e['h']-s['cl_molm3']-e['hco3']-2*e['co3']-e['oh']
        self.assertLess(float(np.max(np.abs(charge))),3e-9)

    def test_charged_flux_and_wall_selectivity(self):
        s=self.sol;f=s['face'];j=self.case.current_Am2/F
        np.testing.assert_allclose(f['q'], -j,rtol=0,atol=2e-17)
        for key in ('ca','cl','dic','na'):
            self.assertEqual(f[key][0],0.0)
        self.assertAlmostEqual(f['ta'][0],j)

    def test_zero_current_is_uniform(self):
        c=replace(self.case,current_Am2=0)
        s=solve_carbonate(c,times=(0.,180.))
        np.testing.assert_allclose(s['pH'],c.bulk_ph,rtol=0,atol=1e-8)
        np.testing.assert_allclose(s['omega_ideal'][0],s['omega_ideal'][1],rtol=0,atol=1e-9)

    def test_cathode_ph_and_calcite_saturation(self):
        s=self.sol
        self.assertGreater(float(s['pH'][-1,0]),self.case.bulk_ph)
        self.assertGreater(float(s['omega_ideal'][-1,0]),float(s['omega_ideal'][0,0]))
        self.assertGreater(float(s['omega_ideal'][-1,0]),1.)

    def test_convergence(self):
        vals=[]
        for nc in (15,25,40):
            c=replace(self.case,n_cells=nc)
            vals.append(float(solve_carbonate(c,times=(0.,180.))['pH'][-1,0]))
        self.assertLess(abs(vals[2]-vals[1]),abs(vals[1]-vals[0]))
        self.assertLess(abs(vals[2]-vals[1]),0.035)

    def test_conserved_inventory_balances(self):
        c=replace(self.case,n_cells=25)
        ts=np.linspace(0,c.time_s,61)
        s=solve_carbonate(c,times=tuple(ts))
        dx=c.layer_m/c.n_cells; q=-c.current_Am2/F
        outward={name:[] for name in ('ca','cl','dic','ta')}
        for idx in range(len(ts)):
            f=face_flux(s['ca_molm3'][idx,-1],s['cl_molm3'][idx,-1],
              s['dic_molm3'][idx,-1],s['ta_molm3'][idx,-1],
              c.ca_molm3,c.cl_molm3,c.dic_molm3,s['bulk_ta_molm3'],dx/2,q)
            for name in outward:outward[name].append(float(f[name]))
        vals={'ca':s['ca_molm3'],'cl':s['cl_molm3'],
              'dic':s['dic_molm3'],'ta':s['ta_molm3']}
        init={'ca':c.ca_molm3,'cl':c.cl_molm3,'dic':c.dic_molm3,'ta':s['bulk_ta_molm3']}
        for name in outward:
            inventory=float(np.sum((vals[name][-1]-init[name])*dx))
            injected=c.current_Am2/F if name=='ta' else 0.
            expected=float(np.trapezoid(injected-np.array(outward[name]),ts))
            self.assertLess(abs(inventory-expected),1e-8, f'{name} conservation mismatch')

    def test_omega_indicator_not_negative(self):
        s=self.sol
        self.assertTrue(np.all(s['omega_ideal']>=0))
        self.assertAlmostEqual(float(ideal_calcite_omega(0.5,0.)),0.)

if __name__=='__main__':unittest.main()
