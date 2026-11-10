import unittest
import numpy as np
from multiday_eclipse_check import trapezoid,projection,residualize,run,P,REF

class Check(unittest.TestCase):
    def test_shape(self):
        np.testing.assert_allclose(trapezoid(np.array([REF-1,REF-.6,REF,REF+.6,REF+1]),0,2),[0,1,1,1,0],atol=1e-11)
        np.testing.assert_allclose(trapezoid(np.arange(10.),2,4),trapezoid(np.arange(10.)+P,2,4),atol=1e-10)
    def test_joint_recovery(self):
        t=REF+np.arange(0,70,.1);g=(t>REF+35).astype(int);q=projection(t,g,2)
        s=-trapezoid(t,5,4);sr=residualize(s,q)
        y=.03*s+.2+.001*(t-REF)+.00001*(t-REF)**2
        self.assertAlmostEqual(sr@residualize(y,q)/(sr@sr),.03,12)
    def test_increment(self):
        t=REF+np.arange(0,60,.1);q=projection(t,np.zeros(len(t)),1);s=-trapezoid(t,10,2);sr=residualize(s,q)
        y=np.sin(t)*.02;beta=sr@residualize(y,q)/(sr@sr)
        self.assertAlmostEqual(sr@residualize(y+.01*s,q)/(sr@sr)-beta,.01,12)
    def test_blind_nondefault_period(self):
        period=P+.03;t=REF+np.arange(0,80,.25)
        y=-.03*trapezoid(t,period/4,2,period)+.1+.0001*(t-REF)
        result=run([dict(key='synthetic',native_t=t,native_y=y)],1,period)
        self.assertAlmostEqual(result['observed_best']['depth'],.03,10)
        self.assertAlmostEqual(result['observed_best']['phase_days'],period/4,10)
        self.assertEqual(result['observed_best']['duration_days'],2)
