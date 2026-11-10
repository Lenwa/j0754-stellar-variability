import unittest

from historical_extension import (applause_reasons, asas_aperture_status, dasch_reasons,
                                  parse_asas, spherical_polygon_contains)


class HistoricalExtensionTests(unittest.TestCase):
    def test_asas_keeps_datasets_apertures_grade_and_hjd(self):
        text = '''#ndata= 1
#dataset= 1 ; 1 F0800+16_438
#ra= 7.912720 07:54:45.8
#dec= 16.694115 16:41:38.8
2635.71095 12.139 12.188 29.999 99.999 12.287 0.032 0.033 0.031 0.036 0.040 A 32257
#ndata= 1
#dataset= 2 ; 1 F0800+24_483
2635.71095 12.239 12.288 12.296 12.335 12.387 0.042 0.043 0.041 0.046 0.050 D 32257
'''
        rows, datasets = parse_asas(text)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(datasets), 2)
        self.assertEqual(rows[0]['hjd'], 2452635.71095)
        self.assertEqual(rows[0]['mag_4'], 12.287)
        self.assertEqual(rows[0]['grade'], 'A')
        self.assertEqual(rows[1]['dataset'], '2')
        self.assertEqual(datasets[0]['field'], 'F0800+16_438')
        self.assertEqual(datasets[0]['ndata'], '1')
        self.assertAlmostEqual(float(datasets[0]['ra_hours']), 7.912720)

    def test_missing_and_negative_flux_sentinels_are_not_magnitudes(self):
        self.assertEqual(asas_aperture_status(29.999), 'below_threshold_no_numeric_limit')
        self.assertEqual(asas_aperture_status(99.999), 'negative_flux')
        self.assertEqual(asas_aperture_status(12.1), 'measured')
        self.assertEqual(asas_aperture_status(float('nan')), 'invalid')

    def test_asas_rejects_incomplete_data_instead_of_silent_skip(self):
        with self.assertRaises(ValueError):
            parse_asas('#dataset= 1 ; 1 F\n2635.7 12.0 A 1')

    def test_dasch_reasons_preserve_simultaneous_failures(self):
        row = dict(magcal_magdep='13', magcal_local_rms='.4', limiting_mag_local='13.1',
                   ra_deg='1', dec_deg='2', ra_cat_corrected='1', dec_cat_corrected='2',
                   aflags=str((1 << 6) | (1 << 11)))
        reasons = dasch_reasons(row)
        self.assertIn('error_above_0.3', reasons)
        self.assertIn('within_0.5mag_of_limit', reasons)
        self.assertIn('aflag_bit_6', reasons)
        self.assertIn('aflag_bit_11', reasons)

    def test_dasch_valid_row_and_missing_detection(self):
        row = dict(magcal_magdep='13', magcal_local_rms='.1', limiting_mag_local='14',
                   ra_deg='1', dec_deg='2', ra_cat_corrected='1', dec_cat_corrected='2', aflags='0')
        self.assertEqual(dasch_reasons(row), [])
        row['magcal_magdep'] = ''
        self.assertIn('no_photometric_detection', dasch_reasons(row))

    def test_applause_keeps_calibration_outlier_but_rejects_depth_artifact_and_rim(self):
        row = dict(natmag='13', cat_natmag='13.2', natmag_error='.2',
                   model_prediction='.99', phot_range_flags='0', annular_bin='4',
                   phot_calib_flags='2')
        self.assertEqual(applause_reasons(row), [])
        row.update(phot_range_flags='2', model_prediction='.1', annular_bin='9')
        self.assertEqual(set(applause_reasons(row)),
                         {'outside_calibrated_depth_range', 'artifact_probability_gate', 'plate_rim'})

    def test_spherical_footprint_handles_ra_wrap_and_rejects_antipode(self):
        polygon = 'Polygon ICRS 359 -1 1 -1 1 1 359 1'
        self.assertTrue(spherical_polygon_contains(polygon, 0, 0))
        self.assertTrue(spherical_polygon_contains(polygon, 359.5, .5))
        self.assertFalse(spherical_polygon_contains(polygon, 180, 0))
        self.assertFalse(spherical_polygon_contains(polygon, 2, 0))

    def test_spherical_footprint_reversed_vertices_and_invalid_shape(self):
        polygon = 'Polygon ICRS 359 1 1 1 1 -1 359 -1'
        self.assertTrue(spherical_polygon_contains(polygon, 0, 0))
        with self.assertRaises(ValueError):
            spherical_polygon_contains('Box ICRS 0 0 1 1', 0, 0)


if __name__ == '__main__':
    unittest.main()
