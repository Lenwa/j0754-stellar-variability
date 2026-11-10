import unittest
import numpy as np
from compare_spherex import rectangular_fnu, source_filename


class ComparisonTests(unittest.TestCase):
    def test_flat_fnu_is_preserved_in_every_rectangular_band(self):
        wave=np.linspace(5000.,60000.,5000)
        # F_lambda = F_nu c/lambda^2; same arbitrary F_nu units returned.
        flux=7.2*2.99792458e18/wave**2
        result=rectangular_fnu(wave,flux,np.array([.8,2.,4.8]),np.array([.02,.06,.04]))
        np.testing.assert_allclose(result,7.2,rtol=2e-6)

    def test_archive_id_maps_to_exact_detector_file(self):
        self.assertEqual(source_filename('ivo://irsa.ipac/spherex_qr2?2025W42_2A_0454_1/D1'),
                         '2025W42_2A_0454_1D1.fits')


if __name__=='__main__':unittest.main()
