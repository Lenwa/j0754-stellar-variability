"""Adversarial synthetic checks; these do not validate astronomical data."""
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import analysis_core as core


class PlateAnalysisTests(unittest.TestCase):
    def test_applause_matches_catalogue_synthetic_magnitudes_when_g_missing(self):
        self.assertTrue(callable(getattr(core, "eligible_applause", None)))
        target={"cat_natmag":"13.3", "gaiaedr3_bp_rp":"1.25", "gaiaedr3_gmag":""}
        row={"cat_natmag":"13.1", "gaiaedr3_bp_rp":"1.35", "gaiaedr3_gmag":""}
        self.assertTrue(core.eligible_applause(row,target))
        self.assertFalse(core.eligible_applause(dict(row,cat_natmag=""),target))
        self.assertFalse(core.eligible_applause(dict(row,cat_natmag="11.0"),target))
        self.assertFalse(core.eligible_applause(dict(row,gaiaedr3_bp_rp="0.5"),target))

    def test_empty_catalogue_count_excluded_without_changing_matching_rules(self):
        self.assertTrue(callable(getattr(core, "eligible_catalog", None)))
        target = {"ref_number":"1", "stdmag":"13.13", "color":"1.01"}
        candidate = {"ref_number":"2", "stdmag":"13.0", "color":"1.1", "num_matches":"600"}
        self.assertTrue(core.eligible_catalog(candidate, target))
        self.assertFalse(core.eligible_catalog(dict(candidate, num_matches=""), target))
        self.assertFalse(core.eligible_catalog(dict(candidate, color=""), target))
        self.assertFalse(core.eligible_catalog(dict(candidate, stdmag="14.5"), target))

    def test_shared_plate_step_removed_without_erasing_target_dip(self):
        self.assertTrue(callable(getattr(core, "differential_curve", None)))
        years = [1911, 1912, 1913, 1934, 1935, 1936, 1951, 1952, 1953]
        controls = {}
        target = []
        for star in range(6):
            controls[str(star)] = [
                {"key": str(y), "year": float(y), "mag": 12 + star / 10 + (0.8 if y >= 1930 else 0)}
                for y in years
            ]
        for y in years:
            target.append({"key": str(y), "year": float(y), "mag": 13 + (0.8 if y >= 1930 else 0) + (1.2 if 1930 <= y < 1940 else 0)})
        result = core.differential_curve(target, controls, minimum_controls=5)
        baseline = [r["corrected_mag"] for r in result if r["year"] < 1930]
        event = [r["corrected_mag"] for r in result if 1930 <= r["year"] < 1940]
        self.assertAlmostEqual(sum(event) / len(event) - sum(baseline) / len(baseline), 1.2, places=10)

    def test_missing_controls_cannot_manufacture_a_dip(self):
        self.assertTrue(callable(getattr(core, "differential_curve", None)))
        years = [1911, 1912, 1913, 1935, 1951, 1952, 1953]
        controls = {str(s): [{"key": str(y), "year": float(y), "mag": 10 + s} for y in years if not (y == 1935 and s >= 4)] for s in range(6)}
        target = [{"key": str(y), "year": float(y), "mag": 13.0} for y in years]
        result = core.differential_curve(target, controls, minimum_controls=5)
        self.assertNotIn(1935.0, [r["year"] for r in result])
        self.assertEqual(len(result), 6)
        self.assertTrue(all(abs(r["corrected_mag"] - 13) < 1e-10 for r in result))

    def test_quality_masks_preserve_variable_catalogue_flag(self):
        self.assertTrue(callable(getattr(core, "quality_dasch", None)))
        row = {"magcal_magdep": "14.0", "magcal_local_rms": "0.12", "limiting_mag_local": "15.2", "aflags": str(1 << 19), "ra_deg": "118.0", "dec_deg": "16.0", "ra_cat_corrected": "118.0", "dec_cat_corrected": "16.0"}
        self.assertTrue(core.quality_dasch(row))
        self.assertFalse(core.quality_dasch(dict(row, aflags=str(1 << 13))))
        self.assertFalse(core.quality_dasch(dict(row, magcal_magdep="")))
        self.assertFalse(core.quality_dasch(dict(row, limiting_mag_local="14.2")))
        self.assertFalse(core.quality_dasch(dict(row, ra_deg="118.1")))

    def test_year_blocks_prevent_one_dense_year_dominating_contrast(self):
        self.assertTrue(callable(getattr(core, "annual_contrast", None)))
        rows = []
        for year, mag, count in [(1911, 0, 3), (1912, 0, 3), (1951, 0, 3), (1952, 0, 3), (1934, 1, 3), (1935, 1, 3), (1936, 9, 300)]:
            rows.extend({"year": float(year), "corrected_mag": float(mag)} for _ in range(count))
        result = core.annual_contrast(rows, replicates=200, minimum_per_year=3)
        self.assertEqual(result["event_years"], 3)
        self.assertEqual(result["baseline_years"], 4)
        self.assertAlmostEqual(result["delta_mag"], 1.0)
        self.assertEqual(result, core.annual_contrast(rows, replicates=200, minimum_per_year=3))


if __name__ == "__main__":
    unittest.main()
