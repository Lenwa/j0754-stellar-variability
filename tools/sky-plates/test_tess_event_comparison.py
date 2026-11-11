"""Physical-unit, finite-cadence and source-position safeguards."""
import unittest
import numpy as np
import compare_j0754_tess_events as c

class EventComparisonTests(unittest.TestCase):
    def test_exposure_average_across_peak(self):
        value=c.exposure_profile(np.array([0.]),np.array([2.]),0.,.5,2.)[0]
        expected=(.5*(1-np.exp(-2))+2*(1-np.exp(-.5)))/2
        self.assertAlmostEqual(value,expected,places=10)

    def test_profile_integral_and_tail_order(self):
        t=np.linspace(-20,100,100001)
        p=c.exposure_profile(t,np.full_like(t,.001),0.,.5,3.)
        self.assertAlmostEqual(np.trapezoid(p,t),3.5,delta=.002)
        self.assertGreater(p[np.argmin(abs(t-2))],p[np.argmin(abs(t+2))])

    def test_fit_recovers_area_despite_gap_and_curved_baseline(self):
        t=np.arange(-24,30,1/6);t=t[(t<2)|(t>3)]
        y=.02*c.exposure_profile(t,np.full_like(t,1/6),.2,.6,4.)+.001+.00003*t+.000002*t*t
        fit=c.fit_event(t,y,np.full_like(t,1/6),2)
        self.assertAlmostEqual(fit['equivalent_duration_s'],.02*4.6*3600/(1+.001+.00003*.2+.000002*.2**2),delta=.2)
        self.assertAlmostEqual(fit['fwhm_h'],np.log(2)*4.6,delta=.003)

    def test_energy_equal_spectral_temperature(self):
        wave=np.linspace(500e-9,1100e-9,50);response=np.ones(50)
        self.assertAlmostEqual(c.energy_factor(wave,response,5000,5000),1.,places=12)

    def test_signed_area_keeps_negative_noise_and_gaps(self):
        self.assertAlmostEqual(c.observed_area(np.array([0.,1.,8.,9.]),np.array([1.,-1.,1.,-1.]),2.),0.)
        self.assertAlmostEqual(c.observed_area(np.array([0.,1.,8.,9.]),np.ones(4),2.),7200.)

    def test_offset_fit_finds_displaced_source_with_background_plane(self):
        yy,xx=np.indices((9,9));base=np.exp(-((xx-4.)**2/1.3+(yy-4.)**2/1.8)/2)
        # Analytic shifted image, independent of the interpolation used in fitting.
        data=50*np.exp(-((xx-4.45)**2/1.3+(yy-3.7)**2/1.8)/2)+2+.2*xx-.1*yy
        fit=c.fit_image_offset(data,base)
        np.testing.assert_allclose(fit['offset_xy'],[.45,-.3],atol=.035)

    def test_absent_positive_source_has_no_localization(self):
        yy,xx=np.indices((7,7));base=np.exp(-((xx-3)**2+(yy-3)**2)/2)
        self.assertFalse(c.fit_image_offset(np.ones((7,7)),base)['localized'])

    def test_nonzero_image_edges_do_not_pin_the_shift_at_zero(self):
        from scipy.ndimage import shift
        yy,xx=np.indices((7,7));base=1000*(np.exp(-((xx-3)**2+(yy-3)**2)/3)+.04*np.cos(xx+yy))
        image=.02*shift(base,[.22,-.31],order=3,mode='grid-constant')+.001*xx
        fit=c.fit_image_offset(image,base)
        np.testing.assert_allclose(fit['offset_xy'],[-.31,.22],atol=.03)

    def test_retained_detector_template_recovers_known_displacement(self):
        from scipy.ndimage import shift
        p=c.pixel_products()[46]
        maps=c.tv.pixel_difference(p['t'],p['cube'],p['quality'],c.EVENTS[2][2],p['xy'],guard=.15,side=.3)['_maps']
        template=c.target_template(p,maps)
        image=.02*shift(template,[.22,-.31],order=3,mode='grid-constant')
        np.testing.assert_allclose(c.fit_image_offset(image,template)['offset_xy'],[-.31,.22],atol=.03)

    def test_control_preserves_observed_sampling_and_rejects_flags(self):
        p=dict(t=np.arange(8,dtype=float),good=np.array([1,1,0,1,1,1,1,1],bool),cadence_days=1.)
        np.testing.assert_array_equal(c.matched_indices(p,[.01,3.01,5.01]),[0,3,5])
        with self.assertRaises(ValueError):c.matched_indices(p,[2.])
        with self.assertRaises(ValueError):c.matched_indices(p,[1.,1.01])
        with self.assertRaises(ValueError):c.matched_indices(p,[3.5])

if __name__=='__main__':unittest.main()
