"""Scientific selection safeguards for the exploratory comparison screen."""
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np


class GiantScreeningTests(unittest.TestCase):
    def test_missing_stellar_parameters_cannot_be_best_analog(self):
        from screen_giant_analogs import analog_score, REFERENCE
        self.assertIsNone(analog_score(dict(Teff=5000, Rad='', Lbol=40)))
        self.assertIsNone(analog_score(dict(Teff=5000, Rad=-1, Lbol=40)))
        self.assertAlmostEqual(analog_score(dict(Teff=REFERENCE['teff'], Rad=REFERENCE['radius'], Lbol=REFERENCE['luminosity'])), 0)

    def test_dwarf_labels_are_excluded_even_with_identical_parameters(self):
        from screen_giant_analogs import rank_catalog, REFERENCE
        row=dict(KIC='1', Teff=REFERENCE['teff'], Rad=REFERENCE['radius'], Lbol=REFERENCE['luminosity'])
        result=rank_catalog([dict(row,LumClass='dwarf'), dict(row,KIC='2',LumClass='giant')])
        self.assertIsNone(result[0]['physical_rank'])
        self.assertEqual(result[1]['physical_rank'],1)

    def test_service_error_is_not_an_empty_catalogue(self):
        from screen_giant_analogs import parse_vizier_tsv
        with self.assertRaises(ValueError):parse_vizier_tsv('#INFO\tError=Table not found\n')

    def test_units_and_separator_are_not_data_and_blank_stays_blank(self):
        from screen_giant_analogs import parse_vizier_tsv
        raw='#Name: test\nKIC\tTeff\tRad\n \tK\tRsun\n---\t----\t----\n 11\t5000\t \n'
        self.assertEqual(parse_vizier_tsv(raw),[dict(KIC='11',Teff='5000',Rad='')])

    def test_year_summary_requires_distinct_nights(self):
        from screen_giant_analogs import annual_nights
        summary=annual_nights(np.repeat(np.arange(59000,59005),20),np.ones(100))
        self.assertEqual(summary,[])

    def test_duplicate_points_do_not_change_nightly_or_yearly_median(self):
        from screen_giant_analogs import annual_nights
        t=np.arange(59000,59020); y=np.arange(20)*.01
        a=annual_nights(t,y)
        b=annual_nights(np.repeat(t,7),np.repeat(y,7))
        self.assertEqual(a,b)
        self.assertEqual(a[0]['nights'],20)

    def test_html_table_ignores_duplicate_tex_annotations(self):
        from fetch_giant_screening import scientific_tables
        tables=scientific_tables('<table><tr><td>2852961</td><td><math><semantics><mn>35.5</mn><annotation>35.5</annotation></semantics></math></td></tr></table>')
        self.assertEqual(tables,[[['2852961','35.5']]])

    def test_rotation_screen_profiles_out_the_same_linear_trend(self):
        from screen_giant_analogs import rotation_screen
        t=np.arange(200,dtype=float)
        y=.2*np.sin(2*np.pi*t/39.4)+.005*t
        result=rotation_screen(t,y,np.full(len(t),.02))
        self.assertLess(abs(result['best_period_days']-39.4),.1)
        self.assertGreater(result['power'],.999)

    def test_short_segments_do_not_get_a_rotation_measurement(self):
        from screen_giant_analogs import rotation_screen
        self.assertIsNone(rotation_screen(np.arange(15.),np.ones(15),np.ones(15)))

    def test_plot_lines_leave_unobserved_gaps_open(self):
        from screen_giant_analogs import gap_series
        t,y=gap_series([0,1,10],[2,3,4],3)
        self.assertTrue(np.isnan(t[2]) and np.isnan(y[2]))
        self.assertEqual(list(t[[0,1,3]]),[0,1,10])

    def test_shortlist_does_not_overwrite_enriched_results_on_cache_reads(self):
        import screen_giant_analogs as screen
        raw='KIC\tTeff\tRad\tLbol\tLumClass\n \tK\tRsun\tLsun\t \n---\t----\t----\t----\t---\n1\t5000\t9\t50\tgiant\n'
        with tempfile.TemporaryDirectory() as folder:
            out=Path(folder);(out/'olah2022_stars.tsv').write_text(raw)
            (out/'screened_catalog.csv').write_text('retained literature enrichment')
            with patch.object(screen,'OUT',out):rows,selected=screen.shortlist()
            self.assertEqual(selected[0]['KIC'],'1')
            self.assertEqual((out/'screened_catalog.csv').read_text(),'retained literature enrichment')
            self.assertFalse((out/'physical_shortlist.csv').exists())

    def test_published_period_pairs_are_preserved_without_becoming_one_period(self):
        from screen_giant_analogs import period_components
        self.assertEqual(period_components('23.458, 29.136'),[23.458,29.136])
        self.assertEqual(period_components('4.304, 4.187'),[4.304,4.187])
        self.assertEqual(period_components('—'),[])
        self.assertEqual(period_components('<5'),[])


if __name__=='__main__':unittest.main()
