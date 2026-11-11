"""Search safeguards: reject impulses/gap joins and count photons only once."""
import unittest
import numpy as np
try:
    import search_tess_flares as s
except ImportError:
    s = None

class SearchTests(unittest.TestCase):
    def test_module_exists(self):
        self.assertIsNotNone(s, 'Flare search implementation is missing')

    @unittest.skipIf(s is None, 'Implementation missing')
    def test_runs_reject_single_impulse_and_split_at_gap(self):
        t = np.array([0., .01, .02, 1., 1.01, 1.02, 1.03])
        self.assertEqual(s.threshold_runs(t, np.array([0, 10, 0, 0, 0, 0, 0.]), 1.), [])
        self.assertEqual(s.threshold_runs(t, np.array([0, 0, 6, 6, 6, 0, 0.]), 1.), [])
        r = s.threshold_runs(t, np.array([0, 0, 0, 4, 6, 4, 0.]), 1.)
        self.assertEqual(len(r), 1)
        self.assertEqual(r[0]['peak_index'], 4)

    @unittest.skipIf(s is None, 'Implementation missing')
    def test_exposure_union_does_not_double_count_or_bridge_gaps(self):
        self.assertAlmostEqual(s.union_days([[(0., 1.), (3., 4.)], [(.5, 2.)]]), 3.)

    @unittest.skipIf(s is None, 'Implementation missing')
    def test_fred_injection_recovered_with_rotating_background(self):
        rng = np.random.default_rng(81)
        t = np.arange(0., 12., 1/144)
        y = .01*np.sin(t/3) + rng.normal(0, .001, len(t))
        added = .025*s.flare_shape(t, 6., 4.6)
        result = s.screen(t, y+added)
        self.assertTrue(any(abs(x['peak_btjd']-6.) < .1 for x in result['events']))

if __name__ == '__main__':
    unittest.main()
