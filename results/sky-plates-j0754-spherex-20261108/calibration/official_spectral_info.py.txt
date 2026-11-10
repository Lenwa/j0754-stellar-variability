__all__ = ["DetectorSpectralInfo"]

import numpy as np
import astropy.units as u
from functools import lru_cache
from typing import Optional
from numpy.typing import ArrayLike
from .lvf import get_wavelength_mesh
from .spx_param import get_lvf_values


class DetectorSpectralInfo:

    def __init__(self, detector: int,
                 central_wl: Optional[np.ndarray] = None, central_bw: Optional[np.ndarray] = None,
                 wl_unit=u.micron, bp_unit=u.micron):
        """Per-pixel spectral info for a detector.

        Parameters
        ----------
        detector : `int`
            Detector id 1-6.
        central_wl : `numpy.ndarray`, optional
            Central wavelength for each active pixel of the detector.
            If `None`, central wavelength is based on simulator's Tabular
            Bandpass class.
        central_bw : `numpy.ndarray`, optional
            Central bandwidth for each active pixel of the detector.
            If `None` central bandwidth is derived from central wavelength as
            central wavelength divided by per-detector spectral resolution.
        wl_unit : `str`
            Wavelength unit, default - micron.
        bp_unit : `str`
            Bandpass unit, default - micron.
        """

        self.array_id = detector

        if central_wl is None:
            # use wavelength model based on simulator's Tabular Bandpass
            self.central_wl = get_wavelength_mesh(detector)
        else:
            self.central_wl = central_wl
        self.wl_unit = wl_unit

        lvf = get_lvf_values(detector)
        self.spectral_r = lvf["R"]
        if central_bw is None:
            self.central_bw = None
        else:
            self.central_bw = central_bw
        self.bp_unit = bp_unit

    @classmethod
    @lru_cache(maxsize=6)
    def create(cls, detector: int):
        """Factory method to create and cache DetectorSpectralInfo instances.

        Parameters
        ----------
        detector : `int`

        Returns
        -------
        DetectorSpectralInfo
            A cached or newly created instance of DetectorSpectralInfo
            for the given detector.
        """
        return cls(detector)

    def get_wl_mid_fwhm(self, x: ArrayLike, y: ArrayLike) -> tuple[ArrayLike, ArrayLike]:
        """Get wavelength and bandwidth at the center of given pixels.

        Parameters
        ----------
        x : `numpy.ndarray` or `list` [`int`]
            Zero-based indexes in x or columns direction.
        y : `numpy.ndarray` or `list` [`int`]
            Zero-based indexes in y or rows direction.
        Returns
        -------
        result : `tuple` [`~numpy.typing.ArrayLike`, `~numpy.typing.ArrayLike`]
            Central wavelength and bandwidth.
        """
        wl = self.central_wl[y, x]
        if self.central_bw is None:
            # approximate bandwidth using spectral resolution for each array
            bw = wl / self.spectral_r
        else:
            bw = self.central_bw[y, x]
        return wl, bw
