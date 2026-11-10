import unittest
import numpy as np

from dust_models import quality, differential, visits, extinction_ratio, cooling_prediction, sampling_interval


class DustTests(unittest.TestCase):
    def good(self):
        return dict(separation=.1,qual_frame='10',qi_fact='1',saa_sep='10',nb='1',
                    moon_masked='00',cc_flags='0000',ph_qual='AA',w1snr='50',
                    w2snr='50',w1sigmpro='.02',w2sigmpro='.02',w1sat='0',w2sat='0',
                    w1rchi2='1',w2rchi2='1',w1mpro='10',w2mpro='10',
                    w1frtr='0000',w2frtr='0000')

    def test_flags_are_band_specific_and_errors_required(self):
        r=self.good();r['cc_flags']='0D00'
        self.assertEqual(quality(r,1),[])
        self.assertIn('cc_flags',quality(r,2))
        r['w1sigmpro']=''
        self.assertIn('error',quality(r,1))

    def test_cuts_do_not_reject_brightness_and_sensitivity_is_explicit(self):
        r=self.good();r['w1mpro']='20';r['w1rchi2']='4';r['w1frtr']='1100'
        self.assertEqual(quality(r,1,'loose'),[])
        self.assertIn('chi2',quality(r,1))
        self.assertIn('transient',quality(r,1))

    def test_ensemble_removes_common_shift_preserves_target_fade(self):
        rows=[]
        for sid,base in [('target',10),('a',11),('b',12),('c',13),('d',14)]:
            for k in range(4):
                rows.append(dict(sid=sid,key=('scan',str(k)),mjd=k,
                                 mag=base+[0,.1,-.1,.2][k]+(.3 if sid=='target' and k>1 else 0)))
        corrected=differential(rows,'target')
        values=[r['corrected'] for r in corrected]
        np.testing.assert_allclose(np.diff(values),[0,.3,0],atol=1e-12)
        self.assertTrue(all(r['n_controls']==4 for r in corrected))

    def test_different_frame_is_not_a_control(self):
        rows=[dict(sid=s,key=('s',s),mjd=1,mag=10) for s in ['target','a','b','c']]
        self.assertEqual(differential(rows,'target'),[])

    def test_visits_require_five_and_gap_greater_than_thirty(self):
        rows=[dict(mjd=float(d),mag=10+d/100) for d in [0,1,2,3,4,34,65,66,67,68]]
        v=visits(rows,'mag')
        self.assertEqual([r['n'] for r in v],[6])
        self.assertAlmostEqual(v[0]['mag'],10.025)

    def test_extinction_interpolation_accepts_descending_grid(self):
        wave=np.array([10.,1.,.1]);cross=2*wave**-1
        self.assertAlmostEqual(extinction_ratio(wave,cross,2,1),.5)

    def test_fixed_area_cooling_is_redder_and_matches_optical(self):
        result=cooling_prediction(4500,.30)
        self.assertAlmostEqual(result['g'],.30,places=10)
        self.assertLess(result['W1'],.10)
        self.assertGreater(result['W1'],.03)
        self.assertLess(result['temperature_faint'],4500)

    def test_sampling_resamples_dates_not_individual_exposures(self):
        rows=[dict(mjd=day+.1+k/100,mag=value) for day,value in [(0,10),(1,10.02),(200,10.1),(201,10.12)] for k in range(3)]
        a=sampling_interval(rows,'mag',(0,100),(200,300),draws=300)
        b=sampling_interval(rows*5,'mag',(0,100),(200,300),draws=300)
        np.testing.assert_allclose(a['interval_95'],b['interval_95'])
        self.assertEqual(a['dates_per_visit'],[2,2])
        self.assertLessEqual(a['interval_95'][0],.1)
        self.assertGreaterEqual(a['interval_95'][1],.1)


if __name__=='__main__':
    unittest.main()
