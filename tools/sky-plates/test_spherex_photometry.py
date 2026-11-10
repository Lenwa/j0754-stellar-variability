import unittest
import numpy as np
from astropy.io.fits import Header

from spherex_photometry import render_epsf,fit_linear_source,surface_to_ujy,calibration_value,native_xy,pixel_mask,fit_centroid_sensitivity


class SpherexTests(unittest.TestCase):
    def epsf(self):
        yy,xx=np.indices((33,33));p=np.exp(-((xx-16)**2+(yy-16)**2)/(2*3.8**2));return p/p.sum()

    def test_pixel_convolved_sampling_has_unit_flux_and_correct_phase(self):
        p=render_epsf(self.epsf(),(21,21),10.27,9.81,(5,5))
        yy,xx=np.indices(p.shape)
        self.assertAlmostEqual(p.sum(),1)
        self.assertAlmostEqual(float((p*xx).sum()),10.27,places=2)
        self.assertAlmostEqual(float((p*yy).sum()),9.81,places=2)

    def test_weighted_plane_and_mask_recover_injected_flux(self):
        p=render_epsf(self.epsf(),(21,21),10.2,10.1,(5,5));yy,xx=np.indices(p.shape)
        image=47*p+3+.02*(xx-10)-.01*(yy-10);variance=np.full_like(p,.01)
        bad=np.zeros_like(p,dtype=bool);bad[10,11]=True;image[10,11]+=1e5
        r=fit_linear_source(image,variance,p,bad,10.2,10.1)
        self.assertAlmostEqual(r['amplitude'],47,places=9)
        self.assertGreater(r['amplitude_error'],0)

    def test_negative_source_is_retained(self):
        p=render_epsf(self.epsf(),(21,21),10.,10.,(5,5))
        r=fit_linear_source(2-3*p,np.ones_like(p),p,np.zeros_like(p,dtype=bool),10.,10.)
        self.assertAlmostEqual(r['amplitude'],-3,places=10)

    def test_edge_crop_does_not_redistribute_missing_source_wings(self):
        full=render_epsf(self.epsf(),(21,21),10.27,10.81,(5,5))
        cropped=render_epsf(self.epsf(),(12,21),10.27,10.81,(5,5))
        np.testing.assert_allclose(cropped,full[:12],rtol=1e-12,atol=1e-14)
        self.assertLess(cropped.sum(),.95)

    def test_fine_centroid_check_recovers_subpixel_offset(self):
        p=render_epsf(self.epsf(),(21,21),10.24,10.04,(5,5))
        yy,xx=np.indices(p.shape);image=100*p+2+.01*xx
        r=fit_centroid_sensitivity(self.epsf(),image,np.full_like(p,.001),np.zeros_like(p,dtype=bool),10.2,10.1,(5,5))
        self.assertTrue(r['success'])
        self.assertAlmostEqual(r['dx'],.04,places=4)
        self.assertAlmostEqual(r['dy'],-.06,places=4)
        self.assertAlmostEqual(r['amplitude'],100,places=3)

    def test_surface_brightness_uses_pixel_solid_angle(self):
        self.assertAlmostEqual(surface_to_ujy(1,37.7645606995),887.6344938,places=5)
        self.assertAlmostEqual(surface_to_ujy(2,37.7645606995),2*surface_to_ujy(1,37.7645606995))

    def test_calibration_is_bilinear_in_native_pixel_coordinates(self):
        yy,xx=np.indices((10,10));a=2*xx+3*yy+1
        self.assertAlmostEqual(calibration_value(a,3.2,4.7),2*3.2+3*4.7+1)
        with self.assertRaises(ValueError):calibration_value(a,-1,3)

    def test_cutout_mapping_retains_fits_one_based_reference(self):
        h=Header({'CTYPE1A':'LINEAR','CTYPE2A':'LINEAR','CRPIX1A':-1127.,'CRPIX2A':-412.,'CRVAL1A':0.,'CRVAL2A':0.,'CDELT1A':1.,'CDELT2A':1.})
        x,y=native_xy(h,7.464454,7.347760)
        self.assertAlmostEqual(x,1135.464454);self.assertAlmostEqual(y,420.347760)

    def test_source_flag_is_not_an_artifact_mask(self):
        flags=np.zeros((5,5),np.int64);flags[2,2]=1<<21;flags[1,1]=1<<15
        mask=pixel_mask(flags)
        self.assertFalse(mask[2,2]);self.assertTrue(mask[1,1])

    def test_native_alternate_mapping_never_applies_primary_sip_distortion(self):
        h=Header({'CTYPE1A':'LINEAR','CTYPE2A':'LINEAR','CRPIX1A':-1127.,'CRPIX2A':-412.,'CRVAL1A':0.,'CRVAL2A':0.,'CDELT1A':1.,'CDELT2A':1.,'A_ORDER':2,'B_ORDER':2,'A_2_0':.001,'B_0_2':.001})
        x,y=native_xy(h,7.464454,7.347760)
        self.assertAlmostEqual(x,1135.464454);self.assertAlmostEqual(y,420.347760)


if __name__=='__main__':unittest.main()
