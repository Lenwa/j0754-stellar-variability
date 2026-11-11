"""Scientific contracts for the two historical fading siblings."""
import unittest
import numpy as np


class TangSiblingTests(unittest.TestCase):
    def test_time_zero_points_do_not_collapse_2015_into_2021(self):
        from tang_sibling_core import absolute_time
        np.testing.assert_allclose(absolute_time([2306.6],dict(BJDREFI=2454833,BJDREFF=.0001,TIMEZERO=.02)),[2457139.6201],rtol=0,atol=1e-8)
        np.testing.assert_allclose(absolute_time([2500.2],dict(BJDREFI=2457000)),[2459500.2],rtol=0,atol=1e-8)

    def test_quality_default_does_not_mean_all_flags_are_good(self):
        from tang_sibling_core import quality_mask
        q=np.array([0,393216,1048576,4])
        self.assertEqual(quality_mask(q,'strict').tolist(),[True,False,False,False])
        self.assertEqual(quality_mask(q,'kepler_default').tolist(),[True,True,False,False])

    def test_joint_fit_recovers_a_known_wave_despite_sector_offsets_and_slopes(self):
        from tang_sibling_core import periodic_fit,signal_matrix
        t=np.linspace(0,80,401);g=np.floor(t/27).astype(int)
        beta=np.array([.013,-.007]);y=signal_matrix(t,39.42)@beta+g*.4+.002*t
        fit=periodic_fit(t,y,g,39.42,'sector_linear')
        np.testing.assert_allclose(fit['beta'],beta,atol=1e-11)
        self.assertAlmostEqual(fit['power'],1.)

    def test_short_arcs_disclose_lost_periodic_information(self):
        from tang_sibling_core import periodic_fit
        t=np.linspace(0,8,100);y=np.sin(2*np.pi*t/39)
        f=periodic_fit(t,y,np.zeros(100),39,'sector_linear')
        self.assertLess(min(f['phase_energy_retention']),.01)

    def test_two_harmonics_can_represent_a_half_period_wave_without_relabeling_it(self):
        from tang_sibling_core import periodic_fit
        t=np.linspace(0,80,401);y=.03*np.cos(2*np.pi*t/(38.4/2))
        f=periodic_fit(t,y,np.zeros(len(t)),38.4,'global_linear',harmonics=2)
        self.assertLess(f['harmonic_semiamplitudes'][0],1e-10)
        self.assertAlmostEqual(f['harmonic_semiamplitudes'][1],.03)

    def test_bin_groups_and_exact_duplicates_do_not_change_the_data(self):
        from tang_sibling_core import bin_series
        t=np.array([0,.1,.2,.3]);y=np.array([1,3,5,7]);g=np.array([0,0,1,1])
        a=bin_series(t,y,g,width=.4)
        b=bin_series(np.repeat(t,3),np.repeat(y,3),np.repeat(g,3),width=.4)
        np.testing.assert_allclose(a['y'],[2,6]);np.testing.assert_allclose(a['y'],b['y'])
        self.assertEqual(a['count'].tolist(),[2,2])

    def test_flat_data_has_no_fractional_variance_improvement(self):
        from tang_sibling_core import periodic_fit
        t=np.linspace(0,80,400)
        f=periodic_fit(t,np.ones(400),np.zeros(400),39,'global_linear')
        self.assertEqual(f['power'],0.)

    def test_flat_scan_reports_no_identifiable_period(self):
        from tang_sibling_core import scan_periods
        t=np.linspace(0,80,400)
        result,power=scan_periods(t,np.ones(400),np.zeros(400),np.array([10.,39.,80.]))
        np.testing.assert_array_equal(power,np.zeros(3))
        self.assertIsNone(result['best_period_days'])

    def test_arc_labels_use_time_gaps_and_sector_identity_in_original_order(self):
        from tang_sibling_core import arc_labels
        t=np.array([4.,0.,.2,4.2,4.3,4.4]);g=np.array(['A','A','A','A','B','B'])
        self.assertEqual(arc_labels(t,g).tolist(),['A_arc1','A_arc0','A_arc0','A_arc1','B_arc0','B_arc0'])

    def test_arc_offsets_and_sector_slopes_recover_injected_modulation(self):
        from tang_sibling_core import periodic_fit,signal_matrix,arc_labels
        t=np.concatenate([np.linspace(start,start+11,60) for start in [0,14,27,41,54,68]])
        g=np.repeat(['A','B','C'],120);arcs=arc_labels(t,g)
        offsets={name:i*.05 for i,name in enumerate(sorted(set(arcs)))}
        beta=np.array([.021,-.013])
        y=signal_matrix(t,38.4)@beta+np.array([offsets[a] for a in arcs])+.001*t
        fit=periodic_fit(t,y,g,38.4,'arc_sector_linear')
        np.testing.assert_allclose(fit['beta'],beta,atol=1e-10)
        self.assertAlmostEqual(fit['power'],1.)

    def test_flagged_time_gaps_cannot_be_joined_into_an_excursion(self):
        from tang_sibling_core import excursion_runs
        self.assertEqual(excursion_runs([0,.01,1],[True,True,True],cadence=.01),[])
        self.assertEqual(excursion_runs([0,.01,.02],[True,True,True],cadence=.01),[[0,1,2]])

    def test_frozen_correct_phase_predicts_heldout_data_and_opposite_phase_fails(self):
        from tang_sibling_core import heldout_score,signal_matrix
        t=np.linspace(55,80,200);g=np.zeros(len(t));beta=np.array([.02,.03]);y=signal_matrix(t,38.4)@beta+.01*t
        good=heldout_score(t,y,g,38.4,'global_linear',beta)
        bad=heldout_score(t,y,g,38.4,'global_linear',-beta)
        self.assertAlmostEqual(good['fraction_improvement'],1.)
        self.assertLess(bad['fraction_improvement'],-2.9)

    def test_k2sc_time_component_is_restored_for_variability_studies(self):
        from tang_sibling_core import restore_k2sc_variability
        t=np.linspace(0,80,401);wave=3*np.sin(2*np.pi*t/38.4)
        corrected=restore_k2sc_variability(np.full(401,100.),100+wave)
        np.testing.assert_allclose(corrected-np.median(corrected),wave-np.median(wave),atol=1e-12)

    def test_one_missing_native_cadence_interrupts_an_excursion(self):
        from tang_sibling_core import excursion_runs
        runs=excursion_runs([0,.02,.03],[True]*3,.01,cadence_numbers=[0,10,15],cadence_step=5)
        self.assertEqual(runs,[])

    def test_noise_threshold_is_local_for_heterogeneous_noise(self):
        from tang_sibling_core import local_residual_scale
        t=np.arange(80)*.01;residual=np.r_[np.tile([-.001,.001],20),np.tile([-.01,.01],20)]
        scales=local_residual_scale(t,residual,half_window=.15,minimum=5)
        self.assertGreater(scales[60]/scales[20],8)
        self.assertTrue(.01>5*scales[20] and .01<5*scales[60])

    def test_floating_point_ties_count_as_equal_null_outcomes(self):
        from tang_sibling_core import count_exceedances
        self.assertEqual(count_exceedances([.96333-2e-16,.96333,.96333-1e-4],.96333),2)

    def test_identical_block_signs_reproduce_the_observed_peak(self):
        from tang_sibling_core import block_sign_diagnostic
        t=np.arange(0,78,.2);y=np.sin(2*np.pi*t/38.4)
        r=block_sign_diagnostic(t,y,np.zeros(len(t)),np.array([20.,38.4,60.]),draws=8,block_days=(14,))['14']
        self.assertEqual(r['exact_patterns'],64)
        self.assertEqual(r['exact_exceedances'],2)


if __name__=='__main__':unittest.main()
