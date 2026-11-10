import unittest
import numpy as np
from decadal_patterns import annual, contrast, template_fit

class DecadalTests(unittest.TestCase):
    def test_annual_actual_date(self):
        rows=[dict(year=2001.1,corrected_mag=1,series='a'),dict(year=2001.3,corrected_mag=3,series='a')]
        self.assertAlmostEqual(annual(rows)[0]['time'],2001.2)
        self.assertEqual(annual(rows)[0]['mag'],2)
    def test_linear_baseline(self):
        rows=[dict(year=y+.1,corrected_mag=2*y,series='a') for y in range(10)]
        self.assertAlmostEqual(contrast(rows,[(0,3),(3,6),(6,10)])['peak_minus_linear_baseline'],0)
    def test_template_known_amplitude(self):
        t=np.linspace(0,10,101); z=np.exp(-(t-3)**2/5)
        rows=[dict(year=x+50,corrected_mag=12+0.6*v,series='a') for x,v in zip(t,z)]
        self.assertAlmostEqual(template_fit(rows,t,z,50)['amplitude'],.6)
    def test_outside_not_extrapolated(self):
        t=np.linspace(0,10,101);z=np.sin(t)
        rows=[dict(year=x,corrected_mag=12+v,series='a') for x,v in zip(t,z)]
        rows.append(dict(year=100,corrected_mag=999,series='a'))
        r=template_fit(rows,t,z,0)
        self.assertEqual(r['excluded_outside_template'],1)
        self.assertAlmostEqual(r['amplitude'],1)

if __name__=='__main__':unittest.main()
