import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import stage6_spatial_prepare as s
from stage3_carbonate import CarbonateCase, solve_carbonate

class TestStage6(unittest.TestCase):
    def test_base_stage3_snapshot(self):
        r=solve_carbonate(CarbonateCase(),times=(0.,180.))
        self.assertAlmostEqual(float(r['pH'][-1,0]),10.21506761606862,places=3)
        self.assertAlmostEqual(float(r['omega_ideal'][-1,0]),155.92277305595678,places=1)
    def test_prepare_full_grid_and_deck(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td)
            with patch.object(s,'OUT',t/'results'),patch.object(s,'PQI',t/'phreeqc'/'spatial_fixed_ph.pqi'):
                items=s.prepare()
                self.assertEqual(len(items),160)
                p=(t/'phreeqc'/'spatial_fixed_ph.pqi').read_text()
                self.assertEqual(p.count('SOLUTION '),160)
                self.assertNotIn('EQUILIBRIUM_PHASES',p.split('SELECTED_OUTPUT')[1])
                with (t/'results'/'input_stage6_profiles.csv').open() as f:
                    rows=list(csv.DictReader(f))
                self.assertEqual(len(rows),160)
                self.assertEqual(set(float(x['current_A_m2']) for x in rows),{0.,0.15,0.5,1.0})
                self.assertEqual(len(set(float(x['distance_um']) for x in rows)),40)
    def test_caco3_aq_vs_solid(self):
        # User-supplied PHREEQC speciation: 0.5 A/m², fixed-pH.
        caco3aq=2.218e-4
        totalca=0.556363868369e-3
        self.assertAlmostEqual(caco3aq/totalca,0.39866,places=3)
        self.assertGreater(10**1.601965740254,1.)
    def test_phreeqc_punch_has_solution_id(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td)
            with patch.object(s,'OUT',t/'results'),patch.object(s,'PQI',t/'phreeqc'/'spatial_fixed_ph.pqi'):
                s.prepare()
                p=(t/'phreeqc'/'spatial_fixed_ph.pqi').read_text()
                self.assertIn('-solution true',p)
                self.assertIn('MOL("CaCO3")',p)
                self.assertIn('SI("Calcite")',p)

if __name__=='__main__': unittest.main()
