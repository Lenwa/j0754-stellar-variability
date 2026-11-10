import unittest

from historical_extension_controls import frame_ensemble


class ControlTests(unittest.TestCase):
    def test_exact_frames_distinct_stars_and_dataset_baselines(self):
        controls = []
        for star in ('a', 'b', 'c'):
            for field, offset in (('f1', 0), ('f2', 2)):
                controls += [dict(star=star, field=field, dataset='1', year=2004,
                                  hjd=i, frame=str(i), mag=12+offset) for i in (1,2,3)]
                controls.append(dict(star=star, field=field, dataset='1', year=2007,
                                     hjd=4, frame='4', mag=12.1+offset))
        target = [dict(hjd=4, frame='4', mag=13.2, year=2007)]
        result = frame_ensemble(target, controls)
        self.assertEqual(result[0]['controls'], 3)
        self.assertAlmostEqual(result[0]['control_offset'], .1)
        self.assertAlmostEqual(result[0]['corrected_mag'], 13.1)

    def test_no_fake_frame_matching_or_duplicate_star_votes(self):
        controls = [dict(star='a', field='f', dataset='1', year=2004,
                         hjd=i, frame=str(i), mag=12) for i in (1,2,3)]
        controls += [dict(star='a', field='f', dataset='1', year=2007,
                          hjd=4, frame='4', mag=12.1)] * 5
        self.assertEqual(frame_ensemble([dict(hjd=4, frame='4', mag=13, year=2007)], controls), [])
        self.assertEqual(frame_ensemble([dict(hjd=4.0001, frame='4', mag=13, year=2007)], controls), [])

    def test_shared_actual_asas_apertures_are_one_vote_across_gaia_ids(self):
        controls=[]
        for star, aperture in (('a','one'),('b','one'),('c','two'),('d','two'),('e','three'),('f','three')):
            controls += [dict(star=star,field='f',dataset='1',aperture_identity=aperture,
                              year=2004,hjd=i,frame=str(i),mag=12) for i in (1,2,3)]
            controls.append(dict(star=star,field='f',dataset='1',aperture_identity=aperture,
                                 year=2007,hjd=4,frame='4',mag=12.1))
        result=frame_ensemble([dict(hjd=4,frame='4',mag=13,year=2007)],controls)
        self.assertEqual(result[0]['controls'],3)
        self.assertEqual(frame_ensemble([dict(hjd=4,frame='4',mag=13,year=2007)],
                                       [r for r in controls if r['star'] not in ('e','f')]),[])


if __name__ == '__main__':
    unittest.main()
