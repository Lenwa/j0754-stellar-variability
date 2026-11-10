import unittest
import numpy as np
from check_existing_chromaticity import design,predictions

class Checks(unittest.TestCase):
    def test_endpoints(self):
        a,b=predictions([12.,12.3],12.,12.3,9.,9.1)
        np.testing.assert_allclose(a,[9.,9.1]);np.testing.assert_allclose(b,[9.,9.1])
    def test_identical_gray(self):
        g=np.linspace(12.,12.4,20);a,b=predictions(g,12.,12.3,9.,9.3)
        np.testing.assert_allclose(a,g-3);np.testing.assert_allclose(a,b)
    def test_joint_coefficients(self):
        t=np.linspace(2458000,2458200,200);c=np.array(['a','b']*100);band=np.where(c=='a','V','g')
        X=design(t,c,band);coef=np.arange(X.shape[1])*.013
        np.testing.assert_allclose(np.linalg.pinv(X)@(X@coef),coef,atol=1e-12)

if __name__=='__main__':unittest.main()
