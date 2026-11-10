import unittest

from asassn_replication import parse_row, contrast


def sample(year, day, mag, camera="bi"):
    return dict(year=year, month=1, day=f"{year}-01-{day:02}", camera=camera,
                band="g", mag=mag, mag_err=.02, flux=10., flux_err=.1, fwhm=2.)


class ReplicationTests(unittest.TestCase):
    def test_limits_and_sentinels_are_not_detections(self):
        base = {"HJD":"2458000", "UT Date":"2018-01-01.5", "Camera":"bi", "Filter":"g",
                "FWHM":"2", "mag":"12", "mag_err":".02", "flux(mJy)":"10", "flux_err":".1"}
        self.assertTrue(parse_row(base)["accepted"])
        self.assertFalse(parse_row(dict(base, mag=">12"))["accepted"])
        self.assertFalse(parse_row(dict(base, mag="99.990", mag_err="99.990"))["accepted"])

    def test_positive_magnitude_change_means_fading_and_days_are_blocks(self):
        rows = [sample(y,d,m) for y,m in [(2018,12.),(2022,12.3)]
                for d in range(1,21) for _ in range(3)]
        r = contrast(rows,2018,2022,draws=100)
        self.assertEqual(r["days"],[20,20])
        self.assertEqual(r["observations"],[60,60])
        self.assertAlmostEqual(r["delta"],.3)
        self.assertGreater(r["ci95"][0],0)

    def test_camera_offsets_cannot_be_mixed(self):
        with self.assertRaises(ValueError):
            contrast([sample(2018,1,12.),sample(2022,1,15.,"bm")],2018,2022)

    def test_many_frames_on_few_days_are_insufficient(self):
        rows = [sample(y,d,12.) for y in (2018,2022) for d in range(1,6) for _ in range(3)]
        self.assertEqual(contrast(rows,2018,2022)["status"],"INSUFFICIENT")


if __name__ == "__main__":
    unittest.main()
