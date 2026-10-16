"""Run `python -m unittest discover -s tests -v` from repository folder."""
import unittest
from dataclasses import replace
from model import Case, solve_case, ph_from_alkalinity, alkalinity_from_ph

class Tests(unittest.TestCase):
    def test_inverse_ph(self):
        for dic in [0., 0.5, 2., 5.]:
            for ph in [6., 7.5, 9.]:
                self.assertAlmostEqual(ph_from_alkalinity(alkalinity_from_ph(ph,dic),dic),ph,places=7)
    def test_zero_current(self):
        s=solve_case(replace(Case(),current_Am2=0.),times=(0.,180.))
        self.assertTrue((abs(s['pH']-7.5)<1e-6).all())
    def test_cathode_increases_ph(self):
        s=solve_case(Case(),times=(0.,180.))
        self.assertGreater(s['pH'][-1,0],7.5)
        self.assertGreater(s['pH'][-1,0],s['pH'][-1,-1])
    def test_buffer_effect(self):
        a=solve_case(replace(Case(),dic_molm3=0.),times=(0.,180.))
        b=solve_case(replace(Case(),dic_molm3=2.),times=(0.,180.))
        self.assertGreater(a['pH'][-1,0],b['pH'][-1,0])
if __name__=='__main__':unittest.main()
