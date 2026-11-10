import unittest
import numpy as np
from analyze_spherex_co import residual_index


class COTests(unittest.TestCase):
    def test_index_cancels_gray_normalization(self):
        w=np.array([4.05,4.15,4.6,4.7,4.8]);m=np.array([4.,3.9,3.2,3.,2.9])
        self.assertAlmostEqual(residual_index(w,7*m,m),1.)

    def test_localized_fillin_has_positive_index(self):
        w=np.array([4.05,4.15,4.6,4.7,4.8]);m=np.ones(5);f=np.array([2.,2.,2.2,2.2,2.2])
        self.assertAlmostEqual(residual_index(w,f,m),1.1)

    def test_missing_reference_is_rejected(self):
        with self.assertRaises(ValueError):residual_index(np.array([4.6,4.7]),np.ones(2),np.ones(2))


if __name__=='__main__':unittest.main()
