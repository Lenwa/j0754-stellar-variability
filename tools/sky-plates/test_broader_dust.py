import unittest
import numpy as np
from broader_dust import cover_ratio, gray_selective, select_rows,final_cover,fixed_cover_from_columns
from stellar_models import Bandpass
from compare_stellar_models import extinction_deficits
from check_dust_baseline_color import fnu

class Tests(unittest.TestCase):
    def test_cover_limits(self):
        d=np.array([0.9,0.2,0.1])
        np.testing.assert_allclose(cover_ratio(d,0.2,0.2),1)
        np.testing.assert_allclose(cover_ratio(d,0,1),1-d)
    def test_gray_selective(self):
        self.assertAlmostEqual(gray_selective(0.2,0.4,np.ones(3))[0],0.6)
    def test_selection_ignores_check(self):
        p=np.array([[0.312,0.1175,99],[0.4,0.1175,0.08875]])
        self.assertEqual(select_rows(p,np.array([0.312,0.1175,0.08875])).tolist(),[True,False])
    def test_exact_anchor(self):
        d=np.array([0.4,0.8,0.999])
        for c0 in [0,0.3,0.7]:
            c1=final_cover(d,c0,0.25)
            np.testing.assert_allclose(cover_ratio(d,c0,c1),0.75)
    def test_integrated_gray(self):
        wave=np.linspace(3000,50000,1000);flux=np.ones_like(wave)
        b=[Bandpass('a',np.array([4000.,5000.,6000.]),np.array([0,1,0]),1),
           Bandpass('b',np.array([30000.,35000.,40000.]),np.array([0,1,0]),1)]
        grain=np.zeros((100,6));grain[:,0]=np.geomspace(.1,100,100);grain[:,3]=1
        d=extinction_deficits(wave,flux,b,grain,np.array([0,0.3]))
        np.testing.assert_allclose(d[0],0,atol=1e-15)
        np.testing.assert_allclose(d[1],1-np.exp(-0.3))
    def test_fullcover_column_ratio(self):
        t0=.2;t1=.7;d0=1-np.exp(-t0);d1=1-np.exp(-t1)
        self.assertAlmostEqual((1-d1)/(1-d0),np.exp(-(t1-t0)))
    def test_column_anchor(self):
        d0=.4;d1=.9;r=.75
        c=fixed_cover_from_columns(d0,d1,r)
        self.assertAlmostEqual((1-c*d1)/(1-c*d0),r)
    def test_flux_density_convention(self):
        wave=np.linspace(3000,50000,1000)
        b=Bandpass('test',np.array([4000.,5000.,6000.]),np.array([0.,1.,0.]),1)
        self.assertAlmostEqual(fnu(wave,2.99792458e18/wave**2*3.7,b),3.7)

if __name__=='__main__':unittest.main()
