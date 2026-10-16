import unittest
from stage5_thermo import Snapshot,davies_gamma,davies_reference,speciation_at_ph
class TestStage5(unittest.TestCase):
    def setUp(self):self.s=Snapshot('test',9.3,10,.5,2,7.5,.418,0,0)
    def test_positive_ionic_strength(self):
        a=davies_reference(self.s)
        self.assertGreater(a['ionic_strength_M'],0)
        self.assertLess(a['gamma2'],a['gamma1'])
    def test_electroneutrality(self):
        a=davies_reference(self.s)
        self.assertLess(abs(a['charge_mol_L']),1e-10)
    def test_dic_balance(self):
        a=davies_reference(self.s)
        self.assertAlmostEqual(a['CO2_M']+a['HCO3_M']+a['CO3_M'],self.s.dic_mM*1e-3,places=10)
    def test_gamma_zero(self):self.assertAlmostEqual(davies_gamma(2,0),1)
    def test_fixed_ph(self):
        a=davies_reference(self.s,False)
        self.assertAlmostEqual(a['pH_activity'],7.5)
    def test_invalid_high_ionic_strength(self):
        with self.assertRaises(ValueError):davies_gamma(1,0.51)
if __name__=='__main__':unittest.main()
