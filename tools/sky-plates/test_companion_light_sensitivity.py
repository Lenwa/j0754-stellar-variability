import unittest
import numpy as np
from companion_light_sensitivity import shape_score, mix

class TestCompanionLight(unittest.TestCase):
    def test_scale_invariance(self):
        y=np.array([1.,2.,3.]); self.assertLess(float(shape_score(7*y,y)),1e-12)
    def test_known_mixture(self):
        p=np.array([1.,2.,3.]); c=np.array([3.,2.,1.]); y=mix(p,c,.3)
        scores=[float(shape_score(y,mix(p,c,f))) for f in np.arange(0,1.01,.1)]
        self.assertEqual(int(np.argmin(scores)),3)
    def test_identical_degeneracy(self):
        p=np.array([1.,2.,3.]); self.assertTrue(np.allclose(mix(p,p,.9),p))
    def test_positive(self):
        with self.assertRaises(ValueError): shape_score(np.array([0.,1.]),np.ones(2))
if __name__=='__main__': unittest.main()
