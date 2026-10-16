"""Structural and physical consistency checks for Stage 2."""
import unittest
from dataclasses import replace
import numpy as np
from stage2_np import NPCase, solve_np, initial_alkalinity, water_species, F

class NPChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case=NPCase(n_cells=35)
        cls.sol=solve_np(cls.case,times=(0.,30.,90.,180.))

    def test_water_equilibrium_and_electroneutrality(self):
        s=self.sol
        oh,h=water_species(s['alkalinity_molm3'])
        np.testing.assert_allclose(oh*h,1e-8,rtol=1e-11)
        np.testing.assert_allclose(s['na_molm3']+h,s['cl_molm3']+oh,atol=1e-12)

    def test_zero_current_uniform(self):
        c=replace(self.case,cathodic_current_Am2=0)
        s=solve_np(c,times=(0.,180.))
        np.testing.assert_allclose(s['pH'],c.bulk_ph,rtol=0,atol=1e-7)
        np.testing.assert_allclose(s['cl_molm3'],c.salt_molm3,rtol=0,atol=1e-8)

    def test_cathodic_ph_increases(self):
        s=self.sol
        self.assertGreater(s['pH'][-1,0],self.case.bulk_ph)
        self.assertGreater(s['pH'][-1,0],s['pH'][-1,-1])

    def test_imposed_charge_flux(self):
        q=self.sol['flux']['q']
        np.testing.assert_allclose(q,-self.case.cathodic_current_Am2/F,rtol=0,atol=1e-18)

    def test_no_cation_leakage_through_electrode(self):
        f=self.sol['flux']
        self.assertEqual(f['na'][0],0.)
        self.assertEqual(f['cl'][0],0.)
        self.assertAlmostEqual(f['a'][0],self.case.cathodic_current_Am2/F)

    def test_grid_convergence(self):
        vals=[]
        for n in (20,40,80):
            s=solve_np(replace(self.case,n_cells=n),times=(0.,180.))
            vals.append(float(s['pH'][-1,0]))
        self.assertLess(abs(vals[2]-vals[1]),abs(vals[1]-vals[0]))
        self.assertLess(abs(vals[2]-vals[1]),0.025)

    def test_integrated_alkalinity_balance(self):
        c=replace(self.case,n_cells=30)
        times=tuple(np.linspace(0,180,121))
        s=solve_np(c,times=times)
        dx=c.layer_m/c.n_cells
        # right boundary alkalinity flux can be reconstructed at each sampled time
        # from numerical conservation by using a dense derivative estimate of total moles;
        # here use end-points and compare to trapz to check finite-volume accounting.
        from stage2_np import _face_fluxes
        J=c.cathodic_current_Am2/F
        out=[]
        a0=initial_alkalinity(c.bulk_ph)
        for cl,a in zip(s['cl_molm3'],s['alkalinity_molm3']):
            r=_face_fluxes(np.array([cl[-1]]), np.array([a[-1]]),
                           np.array([c.salt_molm3]),np.array([a0]),np.array([dx/2]),-J)
            out.append(float(r['a'][0]))
        inventory=np.sum((s['alkalinity_molm3']-a0)*dx,axis=1)
        integral=np.trapezoid(J-np.array(out),times)
        self.assertLess(abs(inventory[-1]-integral)/(J*180),1e-4)

if __name__ == '__main__': unittest.main()
