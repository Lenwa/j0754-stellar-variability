"""Contracts for source attribution, shared trends and motion alternatives."""
import unittest
import numpy as np


class J0830PixelTests(unittest.TestCase):
    def test_source_masked_plane_ignores_a_variable_star(self):
        from j0830_pixel_core import background_planes
        yy,xx=np.indices((13,13));plane=20+2*xx-3*yy
        cube=np.stack([plane.copy(),plane.copy()+100]);cube[0,6,6]+=1000;cube[1,6,6]+=3000
        fitted=background_planes(cube,[(6,6)])
        np.testing.assert_allclose(fitted,np.stack([plane,plane+100]),atol=1e-9)

    def test_aperture_preserves_signed_flux_and_rejects_missing_pixels(self):
        from j0830_pixel_core import aperture_flux
        cube=np.ones((3,9,9));cube[1,4,4]=-10;cube[2,4,4]=np.nan
        f=aperture_flux(cube,(4,4),1.)
        np.testing.assert_allclose(f[:2],[5,-6]);self.assertTrue(np.isnan(f[2]))

    def test_ensemble_removes_shared_throughput_without_removing_target_wave(self):
        from j0830_pixel_core import differential_flux
        t=np.linspace(0,80,401);common=1+.04*np.sin(t/13);wave=1+.025*np.cos(2*np.pi*t/38.4)
        result,ensemble=differential_flux(100*common*wave,np.column_stack([50*common,200*common]))
        np.testing.assert_allclose(result/np.median(result),wave/np.median(wave),atol=1e-12)
        self.assertGreater(np.ptp(ensemble),.05)

    def test_signed_map_finds_the_injected_source_not_constant_neighbour(self):
        from j0830_pixel_core import periodic_image
        t=np.linspace(0,80,401)+2457000;yy,xx=np.indices((11,11));psf=np.exp(-((xx-6.2)**2+(yy-4.3)**2)/1.2)
        wave=np.cos(2*np.pi*(t-2457000)/38.4);cube=1000*psf[None,:,:]*(1+.025*wave[:,None,None])+30
        cube[:,8,2]+=2000
        r=periodic_image(t,cube,np.zeros(len(t)),38.4,[1,0],'global_linear')
        np.testing.assert_allclose(r['signed'],25*psf,atol=1e-9)
        self.assertEqual(np.unravel_index(np.argmax(r['signed']),psf.shape),(4,6))

    def test_joint_spatial_fit_separates_a_stationary_source_from_motion(self):
        from j0830_pixel_core import spatial_components
        yy,xx=np.indices((15,15));template=1000*np.exp(-((xx-7.1)**2+(yy-6.8)**2)/1.5)
        dy,dx=np.gradient(template);signed=.03*template+.008*dx-.004*dy+2+.1*xx
        r=spatial_components(template,signed,(7.1,6.8))
        self.assertAlmostEqual(r['stationary_fraction'],.03,places=10)
        np.testing.assert_allclose(r['motion_coefficients'],[.008,-.004],atol=1e-10)

    def test_pure_motion_does_not_masquerade_as_stationary_brightness(self):
        from j0830_pixel_core import spatial_components
        yy,xx=np.indices((15,15));template=1000*np.exp(-((xx-7.1)**2+(yy-6.8)**2)/1.5)
        _,dx=np.gradient(template);r=spatial_components(template,.02*dx,(7.1,6.8))
        self.assertAlmostEqual(r['stationary_fraction'],0.,places=10)

    def test_comparison_centroids_align_motion_without_using_target_variability(self):
        from scipy.ndimage import shift
        from j0830_pixel_core import register_on_control
        yy,xx=np.indices((25,25));control=np.exp(-((xx-5.2)**2+(yy-5.1)**2)/1.3);target=np.exp(-((xx-15.1)**2+(yy-15.3)**2)/1.3)
        amplitudes=np.array([.98,1.,1.02]);movement=[(-.15,.1),(0,0),(.15,-.1)]
        cube=np.stack([shift(1000*control+800*a*target,m,order=3) for a,m in zip(amplitudes,movement)])
        aligned,offsets=register_on_control(cube,(5.2,5.1))
        # Compare with the known reverse translation. Undersampling means that
        # interpolation itself changes individual peak pixels, even for exact shifts.
        expected=np.stack([shift(frame,(-m[0],-m[1]),order=3,mode='constant',cval=np.nan) for frame,m in zip(cube,movement)])
        np.testing.assert_allclose(aligned[:,3:8,3:8],expected[:,3:8,3:8],atol=1.5)
        np.testing.assert_allclose(offsets[2]-offsets[0],[-.2,.3],atol=.025)
        measured=np.sum(aligned[:,13:18,13:18],axis=(1,2));np.testing.assert_allclose(measured/measured[1],amplitudes,atol=.004)


if __name__=='__main__':unittest.main()
