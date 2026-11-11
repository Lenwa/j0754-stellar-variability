"""Forced radio flux and photon-counting safeguards."""
import unittest
import numpy as np
from astropy.coordinates import SkyCoord
try:
    import analyze_radio_uv as a
except ImportError:
    a=None

class FollowupTests(unittest.TestCase):
    def test_module_exists(self):self.assertIsNotNone(a)

    @unittest.skipIf(a is None,'Implementation missing')
    def test_beam_fwhm_orientation(self):
        b=a.beam(np.array([0.,0.,1.]),np.array([0.,2.,0.]),4.,2.,0.)
        np.testing.assert_allclose(b,[1.,.5,.5],atol=1e-12)

    @unittest.skipIf(a is None,'Implementation missing')
    def test_signed_flux_survives_background_plane(self):
        y,x=np.mgrid[-10:11,-10:11];b=np.exp(-4*np.log(2)*(x*x/2**2+y*y/4**2))
        for amplitude in [-.2,.7]:
            image=amplitude*b+1+.003*x-.002*y
            measured=a.forced_fit(x,y,image,0.,0.,4.,2.,0.)
            self.assertAlmostEqual(measured,amplitude,places=10)

    @unittest.skipIf(a is None,'Implementation missing')
    def test_photon_visit_test_accounts_for_different_exposures(self):
        self.assertAlmostEqual(a.poisson_deviance([100,200],[10,20]),0.,places=10)
        self.assertGreater(a.poisson_deviance([200,100],[10,20]),50.)

    @unittest.skipIf(a is None,'Implementation missing')
    def test_wcs_frame_attributes_do_not_move_source(self):
        world=a.TARGET.transform_to('fk5')
        east,north=a.tangent_offsets(world)
        self.assertAlmostEqual(float(east),0.,places=5)
        self.assertAlmostEqual(float(north),0.,places=5)

if __name__=='__main__':unittest.main()
