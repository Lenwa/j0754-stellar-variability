import unittest
import numpy as np
from check_rotation_aliases import prepare,residual,basis,power,PERIODS
class AliasTests(unittest.TestCase):
    def test_injected_and_nuisance(self):
        t=np.linspace(2457000,2457220,180);e=np.ones(len(t))*.02;g=np.array(['one']*len(t))
        q=prepare(t,e,g);p=PERIODS['Gaia'];y=3+.002*(t-t.mean())+.03*np.sin(2*np.pi*(t-2457000)/p)
        z=residual(y,e,q)
        self.assertAlmostEqual(power(z,basis(t,e,q,p)),1,places=10)
        self.assertLess(power(z,basis(t,e,q,PERIODS['year_plus'])),.6)
    def test_offsets_irrelevant(self):
        t=np.arange(2457000,2457200,1.1);e=np.ones(len(t));g=np.array(['one']*len(t));q=prepare(t,e,g)
        y=np.cos(2*np.pi*t/35.9)
        np.testing.assert_allclose(residual(y,e,q),residual(y+20,e,q),atol=1e-12)
    def test_alias_frequency(self):
        self.assertAlmostEqual(1/PERIODS['year_plus']-1/PERIODS['Gaia'],1/365.25)
        self.assertAlmostEqual(1/PERIODS['Gaia']-1/PERIODS['year_minus'],1/365.25)
if __name__=='__main__':unittest.main()
