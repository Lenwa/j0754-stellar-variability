"""Guard the distinction between missing archive data and a null result."""
import unittest
import argparse
import io
from pathlib import Path
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq

from summarize_asassn import camera_filter, coverage, annual
from fetch_asassn import save_response
from fetch_asassn_job import bounded_minutes


class AsassnTests(unittest.TestCase):
    def test_polling_limit_rejects_nonfinite_and_unbounded_inputs(self):
        for value in ("nan", "inf", "-inf", "-1", "1e308"):
            with self.subTest(value=value), self.assertRaises(argparse.ArgumentTypeError):
                bounded_minutes(value)
        self.assertEqual(bounded_minutes("12"), 12)

    def test_comparisons_cannot_substitute_for_absent_target(self):
        index = [dict(asas_sn_id=1, gaia_id=10), dict(asas_sn_id=2, gaia_id=20)]
        rows = [dict(asas_sn_id=2, camera="bi", jd=2458200, mag=12,
                     mag_err=.02, flux=10, flux_err=.1, quality="G")]
        result = coverage(index, rows, target_id=1)
        self.assertEqual(result["target_rows"], 0)
        self.assertEqual(result["target_status"], "NO_RETURNED_MEASUREMENTS")
        self.assertEqual(result["matched_catalogue_sources"], 2)
        self.assertEqual(result["sources_with_measurements"], 1)

    def test_old_years_cannot_supply_the_missing_2022_comparison(self):
        result = annual([dict(jd=2458200+i, mag=12+i*.001) for i in range(20)])
        self.assertEqual(result[0]["year"], 2018)
        self.assertEqual(result[0]["n"], 20)
        self.assertNotIn(2022, [r["year"] for r in result])

    def test_camera_case_changes_filter_and_unknown_stays_unknown(self):
        self.assertEqual(camera_filter("ba"), "V")
        self.assertEqual(camera_filter("bA"), "g")
        self.assertIsNone(camera_filter("unknown"))

    def test_malformed_first_host_cannot_poison_next_hosts_cache(self):
        with tempfile.TemporaryDirectory(prefix="sky-asassn-test-", dir=Path(__file__).resolve().parent) as directory:
            out = Path(directory)
            with self.assertRaises((ValueError, OSError)):
                save_response(out, "block.parquet", b"not parquet")
            self.assertFalse((out / "block.parquet").exists())
            self.assertEqual((out / "block.parquet.rejected.bin").read_bytes(), b"not parquet")
            buf = io.BytesIO()
            pq.write_table(pa.table({"asas_sn_id": [1]}), buf)
            save_response(out, "block.parquet", buf.getvalue())
            self.assertEqual(pq.read_table(out / "block.parquet")["asas_sn_id"].to_pylist(), [1])


if __name__ == "__main__":
    unittest.main()
