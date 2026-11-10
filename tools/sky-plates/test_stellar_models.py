import unittest
import numpy as np

from stellar_models import Bandpass,band_flux,magnitudes,fit_temperature_scale,fit_partial_cover,fit_two_component,fit_spots,read_bandpass,covering_predictions


class StellarModelTests(unittest.TestCase):
    def test_detector_weighting_and_union_grid(self):
        wave=np.array([1.,2.,3.]);flux=wave.copy()
        energy=Bandpass('x',np.array([1.,1.5,2.,2.5,3.]),np.ones(5),0)
        photon=Bandpass('x',energy.wave,energy.transmission,1)
        self.assertAlmostEqual(band_flux(wave,flux,energy),2.)
        self.assertGreater(band_flux(wave,flux,photon),2.)
        self.assertAlmostEqual(band_flux(wave,flux*3,photon)/band_flux(wave,flux,photon),3.)

    def test_incomplete_wavelength_coverage_fails(self):
        b=Bandpass('x',np.array([1.,2.,3.]),np.ones(3),0)
        with self.assertRaises(ValueError):band_flux(np.array([1.5,2.,2.5]),np.ones(3),b)

    def test_two_band_temperature_and_scale_is_recovered(self):
        temps=np.array([3500.,4000.,4500.,5000.])
        flux=np.stack([(temps/5000)**5,(temps/5000)**2,(temps/5000)**1.5],axis=1)
        bright=flux[-1];faint=flux[1]*.95
        observed=magnitudes(faint,bright)
        result=fit_temperature_scale(temps,flux,5000,observed)
        self.assertAlmostEqual(result['temperature_faint'],4000,places=4)
        self.assertAlmostEqual(result['area_ratio'],.95,places=6)
        np.testing.assert_allclose(result['prediction'],observed,atol=1e-7)
        changed=observed.copy();changed[2]+=3
        other=fit_temperature_scale(temps,flux,5000,changed)
        np.testing.assert_allclose(other['prediction'],result['prediction'])

    def test_partial_cover_recovers_physical_screen(self):
        tau=np.array([.1,1,10,100]);k=np.array([1.,.05,.025])
        deficit=1-np.exp(-tau[:,None]*k)
        observed=-2.5*np.log10(1-.25*deficit[2])
        result=fit_partial_cover(tau,deficit,observed)
        self.assertAlmostEqual(result['covering_fraction'],.25)
        self.assertAlmostEqual(result['tau'],10)
        np.testing.assert_allclose(result['prediction'],observed,atol=1e-10)

    def test_two_component_occultation_recovery_and_physical_bounds(self):
        hot=np.array([10.,3.,2.]);cool=np.array([1.,2.,2.]);area=2.;cover=.4
        observed=magnitudes(hot*(1-cover)+area*cool,hot+area*cool)
        result=fit_two_component(hot,cool,observed)
        self.assertAlmostEqual(result['companion_area_ratio'],area)
        self.assertAlmostEqual(result['covered_fraction'],cover)
        np.testing.assert_allclose(result['prediction'],observed,atol=1e-10)
        self.assertIsNone(fit_two_component(hot,hot,observed))

    def test_spot_covering_preserves_starting_spots(self):
        hot=np.array([10.,3.,2.]);cool=np.array([[1.,1.5,1.2],[2.,2.,1.5]])
        start=.1;end=.4
        observed=magnitudes((1-end)*hot+end*cool[0],(1-start)*hot+start*cool[0])
        result=fit_spots(hot,cool,np.array([3500.,4000.]),observed,start)
        self.assertAlmostEqual(result['spot_fraction_faint'],end)
        self.assertEqual(result['spot_temperature'],3500)
        np.testing.assert_allclose(result['prediction'],observed,atol=1e-10)

    def test_filter_identity_is_checked_despite_case_insensitive_filenames(self):
        from io import BytesIO
        xml='''<VOTABLE version="1.3"><RESOURCE><PARAM name="filterID" datatype="char" arraysize="*" value="Generic/Cousins.I"/><PARAM name="DetectorType" datatype="int" value="0"/><PARAM name="WavelengthUnit" datatype="char" arraysize="*" value="Angstrom"/><TABLE><FIELD name="Wavelength" datatype="double" unit="Angstrom"/><FIELD name="Transmission" datatype="double"/><DATA><TABLEDATA><TR><TD>4000</TD><TD>1</TD></TR><TR><TD>5000</TD><TD>1</TD></TR></TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>'''
        with self.assertRaises(ValueError):read_bandpass(BytesIO(xml.encode()),'SLOAN/SDSS.i')
        b=read_bandpass(BytesIO(xml.encode()),'Generic/Cousins.I');self.assertEqual(b.detector,0)

    def test_envelope_keeps_compatible_nonoptimal_covering_and_ignores_check_band(self):
        deficit=np.array([[.99,.4,.2],[1.,.5,.3]])
        observed=-2.5*np.log10(1-.25*deficit[0])
        result=covering_predictions(deficit,observed,.02)
        self.assertGreater(len(result),2)
        self.assertTrue(any(abs(r['fraction']-.25)>.005 for r in result))
        for r in result:self.assertLessEqual(max(abs(r['prediction'][:2]-observed[:2])),.020000001)
        changed=observed.copy();changed[2]=99
        other=covering_predictions(deficit,changed,.02)
        np.testing.assert_allclose([r['prediction'] for r in result],[r['prediction'] for r in other])


if __name__=='__main__':unittest.main()
