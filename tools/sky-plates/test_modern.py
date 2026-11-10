import unittest
from summarize_modern import summarize_ztf


class ModernPhotometryTest(unittest.TestCase):
    def row(self, oid, mag, band="zg", flags="0"):
        return dict(oid=oid,mag=str(mag),magerr="0.02",catflags=flags,
                    filtercode=band,mjd="59000",ra="118.0",dec="16.0")

    def test_object_and_filter_zero_points_are_never_combined(self):
        rows=[self.row("a",12),self.row("a",12.02),self.row("b",13),self.row("b",13.02),
              self.row("a",11,"zr"),self.row("a",11.02,"zr")]
        result=summarize_ztf(rows)
        self.assertEqual(len(result["groups"]),3)
        self.assertTrue(all(g["robust_scatter_mag"]<.02 for g in result["groups"]))

    def test_rare_fading_is_retained_and_bad_flags_are_counted(self):
        rows=[self.row("a",12)]*40+[self.row("a",13),self.row("a",15,flags="32768")]
        result=summarize_ztf(rows)
        self.assertEqual(result["raw_points"],42)
        self.assertEqual(result["good_points"],41)
        self.assertEqual(result["groups"][0]["fainter_than_median_by_0_3"],1)


if __name__=="__main__":unittest.main()
