import unittest
import numpy as np
from check_tess_oscillations import residual,comb,spectrum,FREQ

class OscillationTests(unittest.TestCase):
    def test_trend_projection(self):
        t=np.linspace(0,25,2000);y=2+3*t
        self.assertLess(np.max(abs(residual(t,y))),1e-10)
    def test_comb_transfer(self):
        t=np.linspace(0,25,2000);s=comb(t)
        self.assertAlmostEqual(np.std(s),1e-4,places=12)
        self.assertGreater(np.std(residual(t,s))/np.std(s),.99)
    def test_frequency_recovery(self):
        t=np.arange(0,25,1/144);y=np.cos(2*np.pi*60e-6*t*86400)
        p=spectrum(t,residual(t,y))
        self.assertLess(abs(FREQ[np.argmax(p)]-60),.1)
    def test_regular_cadence_integrated_power_units(self):
        t=np.arange(0,25,1/144);y=1e-4*np.cos(2*np.pi*60e-6*t*86400)
        p=spectrum(t,residual(t,y))
        # The contained sinusoid has variance 5000 ppm squared.
        self.assertAlmostEqual(np.trapezoid(p,FREQ)/5000,1,delta=.02)

if __name__=='__main__':unittest.main()
