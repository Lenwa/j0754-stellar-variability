"""Numerical contracts for the fixed candidate-period photometry experiment."""
import unittest
import numpy as np

try:
    import orbit_photometry as op
except ModuleNotFoundError:
    op = None


class FixedPeriodTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(op, "Fixed-period numerical implementation is absent")

    def test_joint_fit_recovers_signal_despite_camera_and_slow_trend(self):
        rng = np.random.default_rng(431)
        t = np.sort(rng.uniform(2457000, 2459000, 900))
        group = np.where(np.arange(len(t)) % 2, "a", "b")
        n = op.nuisance(t, group, "spline365")
        s = op.signal(t, 35.9204566, 1)
        beta = np.array([0.037, -0.021])
        y = n @ rng.normal(0, .2, n.shape[1]) + s @ beta
        result = op.fit(t, y, np.full(len(t), .03), group, "spline365")
        np.testing.assert_allclose(result["beta"], beta, atol=1e-11)

    def test_daily_groups_preserve_brightness_extreme_and_camera_identity(self):
        rows = [dict(t=2457000.1, y=10., err=.1, group="a", day="x"),
                dict(t=2457000.2, y=20., err=.1, group="a", day="x"),
                dict(t=2457000.3, y=99., err=.1, group="b", day="x")]
        daily = op.daily(rows)
        self.assertEqual(len(daily), 2)
        self.assertEqual(daily[0]["y"], 15.)
        self.assertEqual(daily[1]["y"], 99.)

    def test_frozen_wrong_phase_signal_hurts_heldout_prediction(self):
        t = np.linspace(2457000, 2458000, 650)
        y = op.signal(t, 35.9204566, 1) @ np.array([.03, .02])
        score = op.heldout(t, y, np.ones(len(t)), np.full(len(t), "a"),
                          "spline365", np.array([-.03, -.02]))
        self.assertLess(score["fraction_sse_improvement"], -2.9)

    def test_partial_sector_projection_exposes_signal_attenuation(self):
        t = np.linspace(2457000, 2457008, 180)
        group = np.full(len(t), "sector")
        f = op.fit(t, np.sin(2*np.pi*(t-2457000)/35.9204566),
                   np.ones(len(t)), group, "season")
        self.assertLess(min(f["signal_energy_retained"]), .05)


if __name__ == "__main__":
    unittest.main()
