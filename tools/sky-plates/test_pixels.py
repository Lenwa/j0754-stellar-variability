"""Checks for scientific failure modes in independent plate measurements."""
import io
import unittest
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from plate_pixels import decode_band, aperture, fit_calibration, predict_magnitude, blank_positions, fit_isotonic, crop_origin


class PlatePixelsTest(unittest.TestCase):
    def test_unsigned_fits_band_and_crop_preserve_pixels_and_sky(self):
        data = (np.arange(120).reshape(10, 12) * 487).astype(np.uint16)
        h = fits.Header(dict(CTYPE1="RA---TAN", CTYPE2="DEC--TAN",
                            CRPIX1=6., CRPIX2=5., CRVAL1=118.69,
                            CRVAL2=16.69, CDELT1=-.001, CDELT2=.001))
        stream = io.BytesIO()
        fits.PrimaryHDU(data, header=h).writeto(stream)
        stream.seek(0)
        original = fits.Header.fromfile(stream)
        offset = stream.tell()
        band = stream.getvalue()[offset+2*12*2:offset+6*12*2]
        crop, ch = decode_band(band, original, 3, 2, 5, 4)
        np.testing.assert_array_equal(crop, data[2:6, 3:8])
        np.testing.assert_allclose(WCS(ch).pixel_to_world_values(2., 1.),
                                   WCS(original).pixel_to_world_values(5., 3.))

    def test_forced_aperture_preserves_negative_noise(self):
        image = np.full((101, 101), 30000.)
        image[49:52, 49:52] -= 50
        value = aperture(image, 50, 50, radius=6, polarity=1)
        self.assertLess(value["signal"], 0)

    def test_background_and_scan_polarity_do_not_change_star_signal(self):
        yy, xx = np.indices((101, 101))
        star = 1200 * np.exp(-((xx-50.2)**2+(yy-49.8)**2)/8.)
        a = aperture(30000+star, 50.2, 49.8, radius=6, polarity=1)
        b = aperture(51000-star, 50.2, 49.8, radius=6, polarity=-1)
        self.assertAlmostEqual(a["signal"], b["signal"], places=5)

    def test_calibration_does_not_extrapolate_to_fainter_target(self):
        mag = np.linspace(11, 14, 30)
        signal = 10**((18-mag)/2.)
        model = fit_calibration(signal, mag)
        self.assertAlmostEqual(predict_magnitude(model, 10**((18-13.2)/2.)), 13.2, places=7)
        self.assertIsNone(predict_magnitude(model, 10**((18-15)/2.)))

    def test_calibration_requires_enough_independent_stars(self):
        self.assertIsNone(fit_calibration(np.array([100., 200.]), np.array([13., 12.])) )

    def test_missing_catalogue_coordinates_do_not_erase_empty_sky(self):
        positions = blank_positions(np.array([50., np.nan]), np.array([50., np.nan]), 101, 101, count=20)
        self.assertEqual(len(positions), 20)
        self.assertTrue(all(np.hypot(x-50, y-50) > 25 for x, y in positions))

    def test_nonparametric_calibration_stays_monotonic_with_noisy_controls(self):
        signals = np.logspace(3, 6, 20)
        mags = 23-2*np.log10(signals)
        mags[8] += .4
        model = fit_isotonic(signals, mags)
        predictions = [predict_magnitude(model, s) for s in np.logspace(3,6,100)]
        self.assertTrue(np.all(np.diff(predictions) <= 0))
        self.assertIsNone(predict_magnitude(model, 100.))

    def test_equal_density_controls_keep_their_weight_in_monotonic_fit(self):
        signals = np.array([1., 1., *range(2,15)])
        mags = np.array([17., 13., *np.linspace(14.5, 8.5, 13)])
        model = fit_isotonic(signals, mags)
        self.assertAlmostEqual(predict_magnitude(model, 1.), 15.)

    def test_enlarged_edge_crop_preserves_target_and_never_reads_outside_scan(self):
        self.assertEqual(crop_origin(1079.,1361.,28616,22930,2200,clamp=True), (0,261))
        with self.assertRaises(ValueError):
            crop_origin(1079.,1361.,28616,22930,2200)


if __name__ == "__main__":
    unittest.main()
