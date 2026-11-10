import unittest
import numpy as np
from fast_relative_rv import shift_fit,C,sky_center
class Tests(unittest.TestCase):
    def test_shift(self):
        w=np.linspace(4900,5500,1500);y=-np.exp(-.5*((w-5200)/3)**2)-.5*np.exp(-.5*((w-5300)/4)**2)
        z=np.interp(w/(1+10/C),w,y)
        self.assertAlmostEqual(shift_fit(w,z,w,y,5000,5450)['shift_kms'],10,places=3)
    def test_amplitude_offset(self):
        w=np.linspace(4900,5500,1500);y=-np.exp(-.5*((w-5200)/3)**2)
        self.assertAlmostEqual(shift_fit(w,y*1.7+.2,w,y,5000,5450)['shift_kms'],0,places=3)
    def test_sky(self):
        w=np.linspace(5550,5600,100);y=3+.2*(w-5577)+10*np.exp(-.5*((w-5577.8)/2)**2)
        self.assertAlmostEqual(sky_center(w,y)['center_A'],5577.8,places=4)
if __name__=='__main__':unittest.main()
