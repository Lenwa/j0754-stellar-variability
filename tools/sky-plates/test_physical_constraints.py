"""Physical limits and numerical anchors for the conditional Argon audit."""
import unittest

try:
    import physical_constraints as pc
except ModuleNotFoundError:
    pc = None


class PhysicalConstraintsTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(pc, 'The new conditional physical-model module is required')

    def test_cgs_gravity_is_converted_before_scale_height(self):
        self.assertAlmostEqual(pc.pressure_scale_height(4500., 1.7)/1e7, 5.7426, delta=.002)

    def test_gravity_radius_obeys_mass_and_round_trip(self):
        radius=pc.radius_from_logg(1.2, 1.7)
        self.assertAlmostEqual(radius, 25.62, delta=.03)
        self.assertAlmostEqual(pc.logg_from_radius(1.2, radius), 1.7, places=12)
        self.assertAlmostEqual(radius/pc.radius_from_logg(.8,1.7), 1.5**.5, places=12)

    def test_ellipsoidal_coefficient_matches_primary_literature(self):
        self.assertAlmostEqual(pc.ellipsoidal_coefficient(.8,.32), 1.422, places=12)
        self.assertAlmostEqual(pc.ellipsoidal_coefficient(.8,.32)/(15*(15+.8)*(1+.32)/(32*(3-.8))), .32, places=12)

    def test_mass_function_uses_actual_candidate_K(self):
        result=pc.mass_function(35.9204566326,10.7244716399,.02297498158)
        # Astropy nominal solar-mass constant; the earlier hand result was rounded.
        self.assertAlmostEqual(result, .00458720, delta=2e-8)
        self.assertGreater(result, 2*pc.mass_function(35.9204566326,8.,.02297498158))

    def test_minimum_companion_mass_reproduces_mass_function(self):
        mass=pc.minimum_companion_mass(1., .00458706)
        self.assertAlmostEqual(mass,.186186,delta=2e-6)
        self.assertAlmostEqual(mass**3/(1+mass)**2,.00458706,places=11)

    def test_primary_roche_lobe_uses_primary_to_companion_mass_ratio(self):
        self.assertAlmostEqual(pc.roche_fraction(1.),.3789205,places=6)
        self.assertGreater(pc.roche_fraction(2.4),pc.roche_fraction(1/2.4))

    def test_unphysical_RV_inclination_is_not_clipped_into_a_solution(self):
        self.assertIsNone(pc.inclination_from_K(30.,20.))
        self.assertAlmostEqual(pc.inclination_from_K(10.,20.),30.,places=12)

    def test_edge_on_roundoff_is_not_a_false_inconsistent_orbit(self):
        self.assertEqual(pc.inclination_from_K(20.*(1+1e-14),20.),90.)

    def test_planck_contrast_uses_wavelength_in_microns(self):
        self.assertAlmostEqual(pc.planck_contrast(.65,3800.,4500.),.4023,delta=.0001)
        self.assertAlmostEqual(pc.planck_contrast(.55,4500.,4500.),1.,places=12)

    def test_spot_fraction_anchors_actual_g_not_invented_V(self):
        contrast=[pc.planck_contrast(w,3800.,4500.) for w in [.48,.65,.8,3.4,4.6]]
        result=pc.spot_prediction(contrast,.312)
        self.assertAlmostEqual(result['fade_mag'][0],.312,places=12)
        self.assertGreater(result['fade_mag'][1],result['fade_mag'][2])
        self.assertGreater(result['fade_mag'][2],result['fade_mag'][3])

    def test_partial_screen_inversion_reproduces_two_anchors(self):
        model=pc.partial_powerlaw(.312,.1175,.48,3.4,.5)
        values=pc.partial_prediction([.48,3.4,4.6],model)
        self.assertAlmostEqual(values[0],.312,places=12)
        self.assertAlmostEqual(values[1],.1175,places=12)
        self.assertLess(values[2],.1175)
        with self.assertRaises(ValueError):pc.partial_powerlaw(.312,.1175,.48,3.4,.1)

    def test_photometric_radius_stays_explicitly_conditional(self):
        r=pc.photometric_radius(12.35,1668.,.15,-.5,4500.)
        self.assertAlmostEqual(r['luminosity_solar'],45.7,delta=.15)
        self.assertAlmostEqual(r['radius_solar'],11.12,delta=.03)


if __name__=='__main__':
    unittest.main()
