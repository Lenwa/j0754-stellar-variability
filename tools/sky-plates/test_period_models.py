import unittest
import numpy as np

from period_models import PERIOD, EPOCH, bin_series, fit_model, replay_pairs


class PeriodTests(unittest.TestCase):
    def test_stable_cycle_is_recovered_with_camera_offsets(self):
        t = EPOCH+np.arange(0,3600,11.)
        camera = np.where(np.arange(len(t))%2,"a","b")
        y = .15*np.sin(2*np.pi*(t-EPOCH)/PERIOD)+np.where(camera=="a",12.,14.)
        model = fit_model(t,y,camera,"periodic",harmonics=1)
        self.assertLess(np.max(np.abs(y-model["prediction"])),1e-10)

    def test_slow_drift_is_not_forced_into_a_periodic_success(self):
        t = EPOCH+np.arange(0,3600,11.)
        camera = np.array(["a"]*len(t))
        y = .0001*(t-EPOCH)
        periodic = fit_model(t,y,camera,"periodic",harmonics=2)
        slow = fit_model(t,y,camera,"slow")
        total = np.sum((y-np.mean(y))**2)
        self.assertLess(1-np.sum((y-periodic["prediction"])**2)/total,.2)
        self.assertLess(np.max(np.abs(y-slow["prediction"])),1e-10)

    def test_future_brightness_cannot_change_training_coefficients(self):
        t = EPOCH+np.arange(0,3600,11.)
        camera = np.array(["a"]*len(t))
        y = np.sin(2*np.pi*(t-EPOCH)/PERIOD)
        train = t < EPOCH+1800
        first = fit_model(t,y,camera,"periodic",train=train)
        modified = y.copy();modified[~train] += 100
        second = fit_model(t,modified,camera,"periodic",train=train)
        np.testing.assert_array_equal(first["beta"],second["beta"])

    def test_unseen_forecast_camera_is_not_calibrated_on_holdout(self):
        t = EPOCH+np.arange(40.)*50
        with self.assertRaises(ValueError):
            fit_model(t,np.zeros(40),np.array(["a"]*20+["b"]*20),"null",train=np.arange(40)<20)

    def test_replay_uses_time_only_same_camera_and_unique_pairs(self):
        t = EPOCH+np.array([0,20,2*PERIOD+1,2*PERIOD+21]*2)
        camera = np.array(["a"]*4+["b"]*4)
        y = np.array([12,12.1,12,12.1,15,15.1,15,15.1])
        pairs = replay_pairs(t,y,camera)
        self.assertEqual(len(pairs),4)
        self.assertEqual(len({p["first_index"] for p in pairs}),4)
        self.assertTrue(all(abs(p["delta"])<1e-10 for p in pairs))

    def test_many_frames_on_one_date_do_not_outweigh_another_date(self):
        rows = [dict(jd=EPOCH+1,day="one",camera="a",band="g",mag=12.)]*10
        rows += [dict(jd=EPOCH+2,day="two",camera="a",band="g",mag=13.)]
        binned = bin_series(rows)
        self.assertEqual(len(binned),1)
        self.assertEqual(binned[0]["days"],2)
        self.assertEqual(binned[0]["mag"],12.5)
        with self.assertRaises(ValueError):
            bin_series(rows+[dict(jd=EPOCH+3,day="three",camera="a",band="V",mag=12.)])

    def test_descriptive_slow_reference_does_not_fabricate_a_forecast(self):
        from analyze_j0754_period import basic
        rows=[dict(jd=EPOCH+i,day=str(i),camera="a",band="g",mag=12+.0001*i)
              for i in range(400,3300,30)]
        result=basic(bin_series(rows))
        self.assertEqual(result["models"]["slow_reference"]["forward_forecast"]["status"],
                         "NOT_IDENTIFIABLE_FOR_EXTRAPOLATION")

    def test_good_single_wave_is_not_penalized_for_extra_harmonic_overfit(self):
        from analyze_j0754_period import decide
        good=dict(full=dict(r2=.8),blocked_cv=dict(rmse_gain=.5),forward_forecast=dict(rmse_gain=.5))
        overfit=dict(full=dict(r2=.82),blocked_cv=dict(rmse_gain=.4),forward_forecast=dict(rmse_gain=-.1))
        dataset=dict(models=dict(one_harmonic=good,two_harmonics=overfit),
                     nearby_training_selected=dict(forecast=dict(rmse_gain=-.1)),
                     two_cycle_replay=dict(rms=.03),cadence_simulation=dict(replay_rms_p95=.08))
        self.assertEqual(decide({"one":dataset,"two":dataset}),"SUPPORTED_DOMINANT_912_CYCLE")


if __name__ == "__main__":
    unittest.main()
