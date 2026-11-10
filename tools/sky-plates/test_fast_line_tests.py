"""Numerical contracts for retrospective FAST diagnostics (no attachment measurements)."""
import unittest
import importlib.util
from pathlib import Path
import numpy as np

_spec = importlib.util.find_spec("fast_line_tests")
m = None if _spec is None else __import__("fast_line_tests")

class NumericalTests(unittest.TestCase):
    def api(self):
        self.assertIsNotNone(m, "FAST numerical API is not implemented")
        return m

    def profile(self):
        x=np.linspace(-40,40,401)
        y=1-1.6/np.sqrt(2*np.pi)/2*np.exp(-.5*(x/2)**2)
        return x,y,np.eye(len(x))*1e-6,(-15,15),[(-35,-20),(20,35)]

    def test_broadband_offset_and_continuum_bias_differ(self):
        f=self.api();x,y,c,line,sides=self.profile()
        w=f.equivalent_width(x,y,c,line,sides)["ew"]
        shifted=f.equivalent_width(x,y+.015,c,line,sides)["ew"]
        biased=f.equivalent_width(x,y,c,line,sides,continuum_bias=.015)["ew"]
        self.assertAlmostEqual(shifted,w/1.015,places=9)
        self.assertAlmostEqual(biased,w+.015/1.015*(30-w),places=9)

    def test_local_emission_subtracts_its_equivalent_area(self):
        f=self.api();x,y,c,line,sides=self.profile()
        emission=.3/np.sqrt(2*np.pi)*np.exp(-.5*x*x)
        a=f.equivalent_width(x,y,c,line,sides)["ew"]
        b=f.equivalent_width(x,y+emission,c,line,sides)["ew"]
        self.assertAlmostEqual(a-b,.3,places=7)

    def test_gradient_includes_fitted_continuum(self):
        f=self.api();x,y,c,line,sides=self.profile()
        r=f.equivalent_width(x,y,c,line,sides)
        for j in [40,150,200,270,350]:
            delta=np.zeros_like(y);delta[j]=1e-5
            fd=(f.equivalent_width(x,y+delta,c,line,sides)["ew"]-
                f.equivalent_width(x,y-delta,c,line,sides)["ew"])/2e-5
            self.assertAlmostEqual(r["gradient"][j],fd,places=7)
        self.assertGreater(r["variance_continuum"],0)
        self.assertGreater(r["sigma"],np.sqrt(r["variance_line"]))

    def test_common_resolution_recovers_equal_injected_lines(self):
        f=self.api();x=np.arange(-50,50,.25);out=np.arange(-40,40,.5)
        values=[]
        for native in [2.3,7]:
            sig=np.hypot(.7,native/2.354820045)
            y=1-1.6/np.sqrt(2*np.pi)/sig*np.exp(-.5*(x/sig)**2)
            A=f.resolution_operator(x,out,native,10)
            z,c=f.apply_operator(A,y,np.full(len(x),1e-6))
            self.assertLess(np.max(abs(A@np.ones(len(x))-1)),1e-12)
            values.append(f.equivalent_width(out,z,c,(-20,20),[(-38,-25),(25,38)])["ew"])
        self.assertLess(abs(values[0]-values[1]),.002)
        self.assertLess(abs(values[0]-1.6),.002)
        with self.assertRaises(ValueError):f.resolution_operator(x,out,10,7)

    def test_covariance_matches_native_pixel_monte_carlo(self):
        f=self.api();rng=np.random.default_rng(1296)
        x=np.arange(-50,50,.5);out=np.arange(-40,40,.8)
        y=1-.2*np.exp(-.5*(x/4)**2);v=np.full(len(x),.003**2)
        A=f.resolution_operator(x,out,6,10);z,c=f.apply_operator(A,y,v)
        r=f.equivalent_width(out,z,c,(-15,15),[(-35,-22),(22,35)])
        draws=(y+rng.normal(size=(3000,len(x)))*np.sqrt(v))@A.T
        empirical=np.std([f.equivalent_width(out,d,c,(-15,15),[(-35,-22),(22,35)])["ew"] for d in draws],ddof=1)
        self.assertLess(abs(empirical/r["sigma"]-1),.06)
        self.assertGreater(np.max(abs(c-np.diag(np.diag(c)))),0)

    def test_repeat_exposures_are_one_night(self):
        f=self.api()
        out=f.night_summary([{"night":"A","ew":1,"sigma":.1},{"night":"A","ew":2,"sigma":.1},{"night":"B","ew":3,"sigma":.2}])
        self.assertEqual(len(out),2)
        self.assertEqual(out["A"]["n_exposures"],2)
        self.assertAlmostEqual(out["A"]["ew"],1.5)
        self.assertAlmostEqual(out["A"]["sigma_formal"],.1/np.sqrt(2))

    def test_independent_night_sensitivity_preserves_repeat_grouping(self):
        f=self.api()
        self.assertTrue(hasattr(f,"independent_night_bounds"),"Missing sensitivity-envelope helper")
        rows=[]
        for night,file,values in [("A","a",[1,1.2]),("B","b",[1.2,1.4]),
                                  ("E","e1",[.8,1]),("E","e2",[.9,1.1])]:
            rows += [{"night":night,"file":file,"ew":v} for v in values]
        r=f.independent_night_bounds(rows,baseline=("A","B"),event="E")
        self.assertAlmostEqual(r["delta_min_A"],.05)
        self.assertAlmostEqual(r["delta_max_A"],.45)

if __name__=="__main__":unittest.main()
