"""Numerical safeguards for extracted-lightcurve TESS diagnostics."""
import unittest
import json
import tempfile
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

try:
    import analyze_tess_variability as tv
except ImportError:
    tv = None


class TestModule(unittest.TestCase):
    def test_analysis_module_exists(self):
        self.assertIsNotNone(tv, "TESS diagnostic implementation is missing")


@unittest.skipIf(tv is None, "Implementation missing")
class TestTessNumerics(unittest.TestCase):
    def test_quality_keeps_real_large_brightening_and_dimming(self):
        t = np.arange(6.)
        f = np.array([1., 100., .01, -1., np.nan, 1.])
        mask = tv.quality_mask(t, f, np.array([0, 0, 0, 0, 0, 1]))
        np.testing.assert_array_equal(mask, [True, True, True, False, False, False])

    def test_native_errors_screen_only_when_present(self):
        t = np.arange(3.)
        self.assertEqual(tv.quality_mask(t, np.ones(3), np.zeros(3)).sum(), 3)
        self.assertEqual(tv.quality_mask(t, np.ones(3), np.zeros(3), [1., 0., np.nan]).sum(), 1)

    def test_binning_preserves_gaps_counts_and_injection_mean(self):
        t = np.array([0., .002, .2, .202])
        b = tv.bin_series(t, np.array([1., 3., 5., 7.]), width=.01)
        np.testing.assert_allclose(b['y'], [2., 6.])
        np.testing.assert_array_equal(b['count'], [2, 2])
        np.testing.assert_allclose(tv.bin_values(np.array([2., 4., 0., 2.]), b['inverse'], 2), [3., 1.])
        self.assertGreater(np.diff(b['t'])[0], .1)

    def test_equivalent_duration_does_not_bridge_gap(self):
        self.assertAlmostEqual(tv.equivalent_duration(np.array([0., 1., 10., 11.]), np.ones(4), max_gap=2.), 2*86400.)
        self.assertAlmostEqual(tv.equivalent_duration(np.array([0., 1.]), -np.ones(2), max_gap=2.), -86400.)

    def test_joint_harmonics_recover_despite_linear_baseline(self):
        t = np.linspace(0., 26., 1000)
        s = tv.periodic_design(t)
        beta = np.array([.01, -.02, .03, .005])
        y = s@beta + 3. + .05*t
        f = tv.template_fit(t, y, linear=True)
        np.testing.assert_allclose(f['beta'], beta, atol=1e-10)
        self.assertLess(f['residual_rms'], 1e-10)

    def test_box_screen_never_bridges_long_data_gap(self):
        t = np.array([0., .01, .5, .51])
        result = tv.box_scores(t, np.ones(4), noise=1., durations=[.6])
        self.assertFalse(np.any(np.isfinite(result['scores'])))

    def test_box_screen_sign_and_visible_injection(self):
        t = np.arange(0., 1., .01)
        y = np.zeros(len(t)); y[(t>=.4)&(t<.5)] = -.02
        result = tv.box_scores(t, y, noise=.001, durations=[.1])
        self.assertLess(np.nanmin(result['scores']), -20.)
        shape = tv.event_shape(t, .4, .1, 'box')
        self.assertEqual(int(np.sum(shape)), 10)
        self.assertTrue(np.all(tv.event_shape(t, 2., .1, 'box') == 0))

    def test_running_baseline_uses_time_and_preserves_gap(self):
        t = np.array([0., .01, 10., 10.01])
        y = np.array([1., 1., 3., 3.])
        np.testing.assert_allclose(tv.local_median(t, y, 1.), y)

    def test_candidates_are_native_json_values(self):
        t = np.arange(0., 2., 1/72)
        result = tv.candidates(t, .01*np.sin(t*10), .001, 10., {})
        json.dumps(result, allow_nan=False)

    def test_local_injection_response_matches_full_pipeline(self):
        t = np.r_[np.arange(0., 2., .01), np.arange(3., 5., .01)]
        y = .002*t + .001*np.sin(t*20)
        baseline = tv.local_median(t, y, 1.)
        added = .02*tv.event_shape(t, 1.5, .1, 'flare')
        actual = tv.injected_residual(t, y, baseline, added)
        expected = y+added-tv.local_median(t, y+added, 1.)
        np.testing.assert_allclose(actual, expected, atol=1e-14)

    def test_pixel_difference_localizes_stationary_source(self):
        t = np.arange(0., 10., .05)
        yy,xx = np.mgrid[:7,:7]
        star = np.exp(-((xx-3)**2+(yy-3)**2)/.8)
        cube = np.repeat((100+1000*star)[None],len(t),axis=0)
        cube[(t>=4)&(t<=4.2)] += 50*star
        result = tv.pixel_difference(t,cube,np.zeros(len(t)),[4,4.2],[3,3])
        np.testing.assert_allclose(result['difference_positive_centroid_xy'],[3,3],atol=.01)
        self.assertAlmostEqual(result['stationary_plus_motion_model']['fractional_stationary_scale'],.05,places=8)

    def test_pixel_difference_rejects_background_only_as_target_brightening(self):
        t = np.arange(0.,10.,.05);yy,xx=np.mgrid[:7,:7]
        star=np.exp(-((xx-3)**2+(yy-3)**2)/.8)
        cube=np.repeat((100+1000*star)[None],len(t),axis=0)
        cube[(t>=4)&(t<=4.2)] += 20+xx*.5
        result=tv.pixel_difference(t,cube,np.zeros(len(t)),[4,4.2],[3,3])
        self.assertAlmostEqual(result['stationary_plus_motion_model']['fractional_stationary_scale'],0.,places=8)

    def test_pixel_difference_distinguishes_displaced_positive_source(self):
        t=np.arange(0.,10.,.05);yy,xx=np.mgrid[:7,:7]
        star=np.exp(-((xx-3)**2+(yy-3)**2)/.8)
        other=np.exp(-((xx-5)**2+(yy-3)**2)/.3)
        cube=np.repeat((100+1000*star)[None],len(t),axis=0)
        cube[(t>=4)&(t<=4.2)] += 100*other
        result=tv.pixel_difference(t,cube,np.zeros(len(t)),[4,4.2],[3,3])
        self.assertGreater(result['difference_centroid_offset_pixels'],1.5)

    def test_native_pixel_table_wcs_and_result_export(self):
        t=np.arange(2571.5,2572.7,.01);yy,xx=np.mgrid[:7,:7]
        star=np.exp(-((xx-3)**2+(yy-3)**2)/.8)
        cube=np.repeat((100+1000*star)[None],len(t),axis=0)
        cube[(t>=2572.0260444438954)&(t<=2572.1371604654987)]+=50*star
        cols=[fits.Column(name='TIME',format='D',array=t),fits.Column(name='FLUX',format='49E',dim='(7,7)',unit='e-/s',array=cube),fits.Column(name='QUALITY',format='J',array=np.zeros(len(t),dtype=int)),fits.Column(name='FFI_FILE',format='20A',array=['synthetic']*len(t))]
        cols += [fits.Column(name=k,format='E',array=np.zeros(len(t))) for k in ['POS_CORR1','POS_CORR2']]
        primary=fits.PrimaryHDU();primary.header['SECTOR']=46;primary.header['CAMERA']=1;primary.header['CCD']=4
        w=WCS(naxis=2);w.wcs.crpix=[4,4];w.wcs.crval=[118.69117410685,16.694621231895];w.wcs.cdelt=[-.0058,.0058];w.wcs.ctype=['RA---TAN','DEC--TAN']
        with tempfile.TemporaryDirectory(dir=Path.cwd(),prefix='tess-test-') as temp:
            out=Path(temp);p=out/'pixels.fits';fits.HDUList([primary,fits.BinTableHDU.from_columns(cols),fits.ImageHDU(np.ones((7,7)),header=w.to_header())]).writeto(p)
            result=tv.analyze_pixels(p,out)
            self.assertAlmostEqual(result['stationary_only_fractional_scale'],.05,places=6)
            self.assertTrue((out/'sector46_pixel_assessment.json').exists())


if __name__ == '__main__':
    unittest.main()
