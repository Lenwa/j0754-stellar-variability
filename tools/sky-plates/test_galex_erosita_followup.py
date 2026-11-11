import unittest
import numpy as np
from astropy.table import Table
from analyze_galex_erosita_followup import flat_response,live_time,bin_edges,constant_test,OUT,rows

class FollowupSafeguards(unittest.TestCase):
    def test_flat_axis_orientation_and_sensitivity_epoch(self):
        flat=np.arange(1,17,dtype=float).reshape(4,4)
        t=np.array([822930000.]);scale=1.9946352-1.9679445e-9*t+9.3025231e-19*t*t
        self.assertAlmostEqual(flat_response(np.array([36000.]),np.array([0.]),t,flat,1)[0],flat[3,2]*scale[0])
    def test_shutter_gap_and_deadtime_reduce_exposure(self):
        live,dead,gap=live_time(0,.15,np.array([[0,1000],[2,1000]],float))
        self.assertAlmostEqual(gap,.05)
        self.assertAlmostEqual(live,.1*(1-dead))
        self.assertGreater(dead,0)
    def test_impossible_deadtime_rejected(self):
        with self.assertRaises(ValueError):live_time(0,.05,np.array([[0,10**9]],float))
    def test_short_tail_is_merged(self):
        self.assertEqual(bin_edges(0,82,20).tolist(),[0,20,40,60,82])
    def test_constant_and_variable_reference(self):
        args=([10]*4,[20]*4,[10]*4,[10]*4,.1)
        self.assertEqual(constant_test(*args,simulations=1000)['null_simulation_p'],1)
        self.assertLess(constant_test([0,0,40,0],*args[1:],simulations=1000)['null_simulation_p'],.01)
    def test_api_errors_are_not_empty_tables(self):
        with self.assertRaises(ValueError):rows('shutter_0_0.json')
    def test_released_counterpart_is_exact_gaia_source(self):
        t=Table.read(OUT/'erosita_gaia_counterparts.ecsv')
        self.assertEqual(len(t),1)
        self.assertEqual(int(t['UID'][0]),311707200169)
        self.assertEqual(int(t['GDR3_source_id'][0]),667401913089980544)
        self.assertEqual(int(t['NWAY_match_flag'][0]),1)

if __name__=='__main__':unittest.main()
