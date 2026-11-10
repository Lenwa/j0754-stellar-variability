# Bring the Update Photometric Calibration code in
# locally to modify and measure the gain function
__all__ = ["UpdatePhotometricCalibrationConfig", "UpdatePhotometricCalibrationTask"]

# Standard library imports
from dataclasses import dataclass, field
import importlib.resources
import logging
import tempfile

# Related third party imports
from astropy import constants as c
from astropy.coordinates import ICRS, SkyCoord, get_body, CartesianRepresentation
from astropy.time import Time
import astropy.units as u
from astropy.table import Table, Column
from astropy.wcs import Sip
from lsst.pex.config import Field, DictField
import numpy as np
from scipy.interpolate import interp1d
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
from scipy.interpolate import UnivariateSpline
from joblib import Parallel, delayed, parallel_backend
from pathlib import Path

# Local application/library specific imports
from spx_data import SPHERExImage
from spx_utils.utils import get_data_file
from spx_utils.catalog_info import binned_catalog_descriptions
from spx_utils.primary_calibrator_info import primary_calibrator_names
from spherex_parameters import load_spherex_parameters
from spx_utils.coarse_astro_utils import get_initial_linear_wcs, get_initial_wcs_with_sip
from spx_utils.update_photo_calib_utils import CalspecPredict

# Simulator dependencies
try:
    from SPHEREx_Simulator_Tools import data_filename
    import SPHEREx_InstrumentSimulator as SPinst
    import SPHEREx_SkySimulator as SPsky
    from SPHEREx_SkySimulator import Catalog_to_Simulate
    hasSPHERExSimulator = True
except ImportError:
    hasSPHERExSimulator = False
# ==============================================================================


_param_dict: dict = load_spherex_parameters(append_quantity=True)


@dataclass
class ArrayCalValues:
    cal_values: dict = field(default_factory=lambda: {
        1: 8.286700087641215,
        2: 8.38949569816677,
        3: 7.896646949035173,
        4: 9.867400104358904,
        5: 2.0516901757300006,
        6: 1.682176698588182
    })


def process_single_observation(row: dict,
                               coord: SkyCoord,
                               time: Time,
                               obspos: CartesianRepresentation,
                               sunpos: CartesianRepresentation,
                               spectab: Table,
                               zodi_model: SPsky.zodicalc.ModifiedKelsallModelWithHPFT,
                               bandpass_model: SPinst.Tabular_Bandpass,
                               wl_unit: u) -> tuple:
    """Process a single observation to compute zodiacal light, true flux,
      mean wavelength, and color correction for a given source and detector.

    Parameters
    ----------
    row : dict
        A dictionary representing a single row of observational data with keys:
        - "WAVELENGTH" : Observed wavelength in `wl_unit`.
        - "X", "Y" : Pixel coordinates of the observation.
        - "DETECTOR" : Detector ID used in the observation.
        - "SOURCE_ID" : Unique identifier for the source.
    coord : `astropy.coordinates.SkyCoord`
        Sky coordinates of the observation in the ICRS frame.
    time : `astropy.time.Time`
        Time of the observation.
    obspos : `astropy.coordinates.CartesianRepresentation`
        Observer's position in Cartesian coordinates.
    sunpos : `astropy.coordinates.CartesianRepresentation` or None
        Sun's position in Cartesian coordinates.
    spectab : Tab
        Spectral table containing source spectra. Must provide:
        - `interp`: Dictionary mapping source IDs to spectral interpolators.
        - `flux_unit`: Dictionary mapping source IDs to flux units.
        - `Spectra`: Dictionary with spectral data for each source.
          Each entry must include:
          - "lambda": Wavelength array for the source spectrum.
          - "FLUX": Flux values at each wavelength.
    zodi_model : `SPsky.zodicalc.ModifiedKelsallModelWithHPFT`
        Zodiacal light model function.
    bandpass_model : `SPinst.Tabular_Bandpass`
        Bandpass model function.
    wl_unit : `astropy.units.Unit`
        Unit of the observed wavelength, typically `u.um`.

    Returns
    -------
    Z : `astropy.units.Quantity`
        Zodiacal light for the observation in MJy/sr.
    true_flux : `float`
        Integrated flux of the source in mJy.
    wl_mean : `float`
        Mean wavelength in the bandpass in microns.
    color_correction : `float`
        Color correction factor due to zodiacal light.
    """
    wl = np.array([row["WAVELENGTH"]]) << wl_unit

    # Calculate the zodi background
    Z = zodi_model(time, coord, wl, obspos=obspos, sunpos=sunpos)

    # Get the bandpass of the pixel
    wls, trans = bandpass_model(row["X"], row["Y"], array=row["DETECTOR"], sparse=True)

    # Integrate the spectrum over the bandpass
    _int = np.sum(trans * spectab.interp[str(row["SOURCE_ID"])](wls)) / np.sum(trans)
    true_flux = (_int * spectab.flux_unit[row["SOURCE_ID"]]).to_value(u.mJy)

    # Calculate mean wavelength
    wl_mean = np.mean(wls)

    # Normalize zodiacal light
    calspec_wl = spectab.Spectra[row["SOURCE_ID"]]["lambda"]
    spx_wl = row["WAVELENGTH"]
    n_i = np.argmin(np.abs(calspec_wl - spx_wl))

    # Ensure calspec_wl has the same shape as obscrd
    reshaped_wl = np.broadcast_to(calspec_wl[:, np.newaxis], (len(calspec_wl), 1)) << u.um
    Z_wavelength = zodi_model(time, coord, reshaped_wl)
    s_zodi = Z_wavelength[:, 0] / Z_wavelength[n_i, 0]
    s_src = spectab.Spectra[row["SOURCE_ID"]]["FLUX"] / spectab.Spectra[row["SOURCE_ID"]]["FLUX"][n_i]

    zodi_interp = interp1d(calspec_wl, s_zodi, fill_value=0.0, bounds_error=False)
    src_interp = interp1d(calspec_wl, s_src, fill_value=0.0, bounds_error=False)

    int_zodi = np.sum(trans * zodi_interp(wls)) / np.sum(trans)
    int_src = np.sum(trans * src_interp(wls)) / np.sum(trans)
    color_correction = int_zodi / int_src

    return Z, true_flux, wl_mean, color_correction


def calc_pred_flux_and_color_corr_parallel(spectab: Table,
                                           obstab: Table,
                                           bandpass_model: SPinst.Tabular_Bandpass,
                                           zodi_model: SPsky.zodicalc.ModifiedKelsallModelWithHPFT,
                                           wl_unit: u,
                                           njobs: int = 8):
    """Compute predicted fluxes, zodiacal light, mean wavelengths,
       and color corrections for a set of observations using parallel processing.

    Parameters
    ----------
    spectab : `astropy.table.Table`
        Spectral table containing source spectra. Must provide:
        - `interp`: Dictionary mapping source IDs to spectral interpolators.
        - `flux_unit`: Dictionary mapping source IDs to flux units.
        - `Spectra`: Dictionary with spectral data for each source.
          Each entry must include:
          - "lambda": Wavelength array for the source spectrum.
          - "FLUX": Flux values at each wavelength.
    obstab : `astropy.table.Table`
        Table of observations containing the following columns:
        - "RA", "DEC": Sky coordinates of the observation in degrees.
        - "JD": Julian Date of the observation.
        - "X", "Y": Pixel coordinates of the observation.
        - "DETECTOR": Detector ID.
        - "WAVELENGTH": Observed wavelength in `wl_unit`.
        - "SOURCE_ID": Unique identifier for the source.
    bandpass_model : `SPinst.Tabular_Bandpass``
        Bandpass model function with signature:
        `bandpass_model(x, y, array=None, sparse=True)`.
        Returns:
        - `wls`: Wavelengths in the bandpass.
        - `trans`: Transmission values at the corresponding wavelengths.
    zodi_model : `SPsky.zodicalc.ModifiedKelsallModelWithHPFT`
        Zodiacal light model function with signature:
        `zodi_model(time, coord, wavelength, obspos=None, sunpos=None)`.
    wl_unit : `astropy.units.Unit`
        Unit of the observed wavelength, typically `u.um`.
    njobs : `int`, optional
        Number of parallel jobs to use. Default is 8.

    Returns
    -------
    obstab : `astropy.table.Table`
        Input observation table with the following new columns:
        - "zodi": Zodiacal light contamination (in MJy/sr).
        - "predicted_fluxes": Predicted fluxes for each observation (in mJy).
        - "wavelength_mean": Mean wavelengths for each observation (in microns).
        - "zodi_cc": Color correction factors for each observation.
    """
    coords = SkyCoord(obstab["RA"], obstab["DEC"], frame=ICRS()).geocentrictrueecliptic
    times = Time(obstab["JD"], format="jd")

    zodiname = zodi_model.__class__.__name__
    if zodiname == "ModifiedKelsallModelWithHPFT":
        suns = get_body("sun", times).geocentrictrueecliptic
        obss = iter(int, None)
    elif zodiname == "ModifiedKelsallModel":
        suns = iter(int, None)
        obss = get_body("earth", times).heliocentrictrueecliptic
    else:
        raise ValueError(f"Unknown zodi model: {zodiname}")

    # Parallelize the row-level processing
    with parallel_backend('multiprocessing'):
        results = Parallel(n_jobs=njobs, verbose=5)(
            delayed(process_single_observation)(
                row, coord, time, obspos, sunpos, spectab, zodi_model, bandpass_model, wl_unit
            )
            for row, coord, time, obspos, sunpos in zip(obstab, coords, times, obss, suns)
        )

    # Unpack results and update the table
    zodi, true_fluxes, wl_mean, zodi_cc = zip(*results)

    obstab["zodi"] = Column(zodi, unit=u.MJy / u.sr)
    obstab["predicted_fluxes"] = Column(true_fluxes, unit=u.mJy)
    obstab["wavlength_mean"] = Column(wl_mean, unit=u.um)
    obstab["zodi_cc"] = Column(zodi_cc, unit=u.dimensionless_unscaled)

    return obstab


def calc_pred_flux_and_color_corr_local(spectab: Catalog_to_Simulate,
                                        obstab: Table,
                                        bandpass_model: SPinst.Tabular_Bandpass,
                                        zodi_model: SPsky.zodicalc.ModifiedKelsallModelWithHPFT,
                                        wl_unit: u) -> Table:
    """Calculate predicted fluxes and color corrections for observations.

    Parameters
    ----------
    spectab : `Catalog_to_Simulate`
        A spectral data object that contains interpolated spectral functions
        for each source and the associated flux unit for each source.

    obstab : `astropy.table.Table`
        The observation table that contains the following columns:
        - "RA" : Right Ascension (degrees).
        - "DEC" : Declination (degrees).
        - "JD" : Julian Date of observation.
        - "X" : X-coordinate of the pixel.
        - "Y" : Y-coordinate of the pixel.
        - "DETECTOR" : Detector ID.
        - "WAVELENGTH" : Wavelength of observation (microns).
        - "SOURCE_ID" : ID of the observed source.

    bandpass_model : `SPinst.Tabular_Bandpass`
        A model that returns the wavelength and transmission for the pixel
        given its X, Y coordinates, and the detector ID. This is used to
        calculate the bandpass for the observation.

    zodi_model : `SPsky.zodicalc.ModifiedKelsallModelWithHPFT`
        A zodiacal light model that returns the zodiacal light contamination
        at the specified position, time, and wavelength. It can be either
        `ModifiedKelsallModelWithHPFT` or `ModifiedKelsallModel`.

    wl_unit : `astropy.units.Unit`
        The unit for wavelength, typically `u.um` for microns.

    Returns
    -------
    obstab : `astropy.table.Table`
        The input observation table with the following new columns:
        - "zodi" : Zodiacal light contamination (in MJy/sr).
        - "predicted_fluxes" : The predicted flux for each observation (in mJy).
        - "wavlength_mean" : The mean wavelength for the observation (in microns).
        - "zodi_cc" : The color correction factor due to zodiacal light.

    Notes
    -----
    Keeping this function during validation, will be removed after assessment.
    """
    coords = SkyCoord(obstab["RA"], obstab["DEC"], unit="deg", frame=ICRS()).geocentrictrueecliptic
    times = Time(obstab["JD"], format="jd")

    zodiname = zodi_model.__class__.__name__
    if zodiname == "ModifiedKelsallModelWithHPFT":
        suns = get_body("sun", times).geocentrictrueecliptic
        obss = iter(int, None)  # dummy iterator to give None
    elif zodiname == "ModifiedKelsallModel":
        suns = iter(int, None)  # dummy iterator to give None
        obss = get_body("earth", times).heliocentrictrueecliptic
    else:
        # Modify if new model is introduced in the future...
        raise ValueError(f"Unknown zodi model: {zodiname}")

    zodi = []
    true_fluxes = []
    zodi_cc = []
    wl_mean = []

    for row, coord, time, obspos, sunpos in zip(obstab, coords, times, obss, suns):
        wl = row["WAVELENGTH"] << wl_unit

        # Calculate the zodi background at this position, time, and wavelength
        Z = zodi_model(time, coord, wl, obspos=obspos, sunpos=sunpos)
        zodi.append(Z)
        # Get the bandpass of the pixel where this landed
        wls, trans = bandpass_model(row["X"], row["Y"], array=row["DETECTOR"], sparse=True)
        # integrate the spectrum of the source over the SPHEREx bandpass
        _int = np.sum(trans * spectab.interp[str(row["SOURCE_ID"])](wls)) / np.sum(trans)
        # Let's make sure we convert everything to milliJansky
        true_fluxes.append((_int * spectab.flux_unit[row["SOURCE_ID"]]).to_value(u.mJy))
        # Wavelength mean
        wl_mean.append(np.mean(wls))

        # Normalize to the wavelength where the observation occurs
        calspec_wl = spectab.Spectra[row["SOURCE_ID"]]["lambda"]
        spx_wl = row["WAVELENGTH"]
        n_i = np.argmin(np.abs(calspec_wl - spx_wl))

        # Zodiacal light normalization
        nwl = len(calspec_wl)
        boresight_icrs = SkyCoord(row['RA']*np.ones(nwl), row['DEC']*np.ones(nwl), frame='icrs', unit='deg')
        Z_wavelength = zodi_model(time, boresight_icrs, calspec_wl << u.um)

        s_zodi = Z_wavelength / Z_wavelength[n_i]  # Normalize zodi to the observation wavelength
        # Normalize source flux
        s_src = spectab.Spectra[row["SOURCE_ID"]]["FLUX"] / spectab.Spectra[row["SOURCE_ID"]]["FLUX"][n_i]

        zodi_interp = interp1d(calspec_wl, s_zodi, fill_value=0.0, bounds_error=False)
        src_interp = interp1d(calspec_wl, s_src, fill_value=0.0, bounds_error=False)

        int_zodi = np.sum(trans * zodi_interp(wls)) / np.sum(trans)
        int_src = np.sum(trans * src_interp(wls)) / np.sum(trans)

        color_correction = int_zodi/int_src
        zodi_cc.append(color_correction)

    # Convert lists to numpy arrays
    zodi = np.array(zodi)
    true_fluxes = np.array(true_fluxes)
    wl_mean = np.array(wl_mean)
    zodi_cc = np.array(zodi_cc)

    # Add new columns to obstab
    obstab["zodi"] = Column(zodi, unit=u.MJy / u.sr)
    obstab["predicted_fluxes"] = Column(true_fluxes, unit=u.mJy)
    obstab["wavlength_mean"] = Column(wl_mean, unit=u.um)
    obstab["zodi_cc"] = Column(zodi_cc, unit=u.dimensionless_unscaled)

    return obstab


def flux_engineering_units_local(obstab: Table,
                                 flat_fields: np.ndarray,
                                 flat_field_spline: UnivariateSpline,
                                 detector_id: int,
                                 initial_distortion_table: Table = None,
                                 sip: Sip = None) -> Table:
    """Calculate and append flux in engineering units for each entry
       in the observation table.

    Parameters
    ----------
    obstab : `astropy.table.Table`
        A table containing observational data. The table must have the
        following columns:
        - "FLUX" : Flux values in the original units.
        - "FLUX_ERR" : Flux error values.
        - "X" : Pixel X-coordinates.
        - "Y" : Pixel Y-coordinates.
        - "DETECTOR" : Detector IDs.
    flat_fields : `np.ndarray`
        An array of flat field data for each detector.
    detector_id : `int`
        Detector id number.
    initial_distortion_table : `astropy.table.Table`
        Initial distortion model from ray-tracing model.
    sip : `Sip`
        Most up-to-date SIP distortion solution.

    Returns
    -------
    obstab : `astropy.table.Table`
        The input table with the following new or updated columns:
        - "flux_eng" : Flux values converted to engineering units (electron/s).
        - "flux_eng_err" : Flux error values in engineering units (electron/s).
        - "omega_pix" : The solid angle of each pixel in steradians (sr).

    Notes
    -----
    This function computes the flux values in engineering units (electron/s)
    for the input observation table (`obstab`), along with the associated error
    and the pixel solid angle (`omega_pix`). The calculation uses
    instrument parameters including aperture size, pixel scale, and
    filter properties. Flat-field corrections are applied to the flux measurements.
    """
    xx, yy = np.meshgrid(np.arange(2049) - 0.5, np.arange(2049) - 0.5)
    obstab["flux_eng"] = np.zeros(len(obstab)) << (u.electron / u.s)
    obstab["flux_eng_err"] = np.zeros(len(obstab)) << (u.electron / u.s)
    obstab["omega_pix"] = np.zeros(len(obstab)) << u.sr

    _ap = _param_dict["nominal_instrument_models"]["aperture"]["quantity"]
    _area = np.pi * (_ap.to_value(u.m) / 2.0) ** 2
    _ncm = _param_dict["nominal_camera_model"]
    _pix_sr = (_ncm["pix_scale"]["quantity"] ** 2).to_value(u.sr)
    _det = _param_dict["detector"]
    eta_total = np.array(_det["eta_total"]["value"])
    _fil = _param_dict["filter"]
    R = np.array(_fil["R"]["value"])
    _factor = 1.0e-20 / c.h.si.value * (u.electron / u.second) / (u.MJy / u.sr)

    cal_values = ArrayCalValues()

    # TODO 6.1 include as input most up to date SIP solution
    # to test distortion -> solid angle
    # Linear WCS
    wcs_linear = get_initial_linear_wcs(0.0, 0.0, array=detector_id)
    # WCS -SIP with initial distortion model from ray tracing.
    wcs = get_initial_wcs_with_sip(wcs_linear, initial_distortion_table, array=detector_id)
    rr, dd = wcs.all_pix2world(xx, yy, 0)
    rr[rr > 180.0] -= 360.0
    diff1 = (rr[:-1, :-1] - rr[:-1, 1:]) << u.deg
    diff2 = (dd[:-1, :-1] - dd[1:, :-1]) << u.deg
    Omega = (diff1 * diff2).to(u.sr)

    ii = obstab["DETECTOR"] == detector_id

    obstab["omega_pix"][ii] = Omega[(obstab["X"][ii]).astype("int16"), (obstab["Y"][ii]).astype("int16")]

    _ap = _param_dict["nominal_instrument_models"]["aperture"]["quantity"]
    _area = np.pi * (_ap.to_value(u.m) / 2.0) ** 2
    _pix_sr = (_ncm["pix_scale"]["quantity"] ** 2).to_value(u.sr)
    # TODO CAL from image header (temporary short cut)
    cal_value = cal_values.cal_values[detector_id]
    logging.info(f"Array {detector_id} CAL value: {cal_value}")

    # TODO 6.1 use Gain matrix from module L2.16
    Mjysr_to_e_per_s = (_area * _pix_sr * eta_total / R) * _factor
    Mjysr_to_e_per_s_value = Mjysr_to_e_per_s[ii-1]
    logging.info(f"Mjysr_to_e_per_s = {Mjysr_to_e_per_s[ii-1] }")

    # Assert if all elements in the array are nearly equal to cal_value
    differences = np.abs(Mjysr_to_e_per_s_value.value - cal_value)
    relative_diff = np.abs(differences / cal_value)
    logging.info(f"Absolute differences: {differences}")
    logging.info(f"Relative differences: {relative_diff}")

    Omega_ref = (6.15*u.arcsec)**2
    conv = Mjysr_to_e_per_s[ii-1][0] / obstab["omega_pix"][ii] * (Omega_ref.to(u.sr))
    logging.info(f"Array {detector_id} CONV value: {conv}")

    # Here we do something different, DCM.
    # Use the spectral component of instrument gain.
    ff_spline_res = flat_field_spline(obstab["WAVELENGTH"])
    ff_spline_res = flat_field_spline(obstab["WAVELENGTH"][obstab["DETECTOR"] == 1])
    conv *= ff_spline_res

    obstab["flux_eng"][ii] = obstab["FLUX"][ii] * conv
    obstab["flux_eng_err"][ii] = obstab["FLUX_ERR"][ii] * conv

    return obstab


def compute_corr_gain_and_weight_local(obstab: Table) -> Table:
    """Compute the gain factor and weight for each observation and
       append them to the observation table.

    Parameters
    ----------
    obstab : `astropy.table.Table`
        The observation table that contains the following columns:
        - "predicted_fluxes" : The predicted flux values (mJy).
        - "flux_eng" : The flux in engineering units (electron/s).
        - "FLUX_ERR" : The flux error values.

    Returns
    -------
    obstab : `astropy.table.Table`
        The input table with two new columns added:
        - "gain_factor" : The gain factor, which is the ratio of
          predicted flux to flux in  engineering units (dimensionless).
        - "weight" : The weight, computed as the inverse square of the
          flux error (dimensionless).
    """
    # Compute gain_factor and weight
    gain_factor = obstab["predicted_fluxes"] / obstab['flux_eng']
    weight = 1.0 / (obstab['FLUX_ERR']**2)

    # Append gain_factor and weight to the obstab table with appropriate units
    obstab["gain_factor"] = Column(gain_factor, unit=u.dimensionless_unscaled)
    obstab["weight"] = Column(weight, unit=u.dimensionless_unscaled)

    # Return the updated obstab table
    return obstab


def generate_binned_calibration_data_local(obstab: Table,
                                           bandpass_model: SPinst.Tabular_Bandpass,
                                           detector_id: int,
                                           nbins: int = 51,
                                           version: str = "R2",
                                           snr_threshold: float = 10.0) -> list[Table]:
    """Generate binned calibration data for observations based on gain,
    zodiacal correction, and pixel solid angle.

    Parameters
    ----------
    obstab : `astropy.table.Table`
        The observation table containing the following columns:
        - "FLUX" : Flux values.
        - "FLUX_ERR" : Flux error values.
        - "DETECTOR" : Detector ID.
        - "WAVELENGTH" : Wavelength of observation (microns).
        - "weight" : Weight of each observation (inverse of flux error squared).
        - "gain_factor" : The gain factor for the observation (dimensionless).
        - "zodi_cc" : Zodiacal light color correction factor.
        - "omega_pix" : The pixel solid angle in steradians (sr).

    bandpass_model : `SPinst.Tabular_Bandpass`
        A function or model that returns the wavelength and transmission for each pixel
        given its X, Y coordinates and the detector array. Used to determine the
        central wavelength and transmission for binning.

    detector_id : `int`
        Detector id number.

    nbins : `int`, optional
        The number of wavelength bins to divide the data into. Default is 51.

    version : `str`, optional
        The version string to include in the output table's metadata. Default is "R2".

    snr_threshold : `float`, optional
        The minimum signal-to-noise ratio (SNR) required for an observation to be
        included in the binned data. Default is 10.0.

    Returns
    -------
    Tbin : `list[astropy.table.Table]`
        A table containing the binned calibration data with the following columns:
        - "wl_min" : The minimum wavelength of each bin (microns).
        - "wl_max" : The maximum wavelength of each bin (microns).
        - "a_bin" : The weighted average gain factor in each bin.
        - "z_bin" : The weighted average zodiacal light color
          correction factor in each bin.
        - "omega_bin" : The weighted average pixel solid angle in each bin (steradians).
    """
    xx, yy = np.meshgrid(np.arange(2049) - 0.5, np.arange(2049) - 0.5)

    wl, t = bandpass_model(xx, yy, array=detector_id, central_bandpass_only=True)
    snr = obstab["FLUX"] / obstab["FLUX_ERR"]
    T = obstab[(obstab["DETECTOR"] == detector_id) & (snr > snr_threshold)]
    bin_edge = np.linspace(wl.min(), wl.max(), num=nbins + 1)
    w = T["weight"]
    wl = T["WAVELENGTH"]
    cal = T["gain_factor"]
    z = T["zodi_cc"]
    omega = T["omega_pix"]
    vals_gain = (
        np.histogram(wl, weights=(w * cal), bins=bin_edge)[0]
        / np.histogram(wl, weights=(w), bins=bin_edge)[0]
    )
    vals_zodi = (
        np.histogram(wl, weights=(w * z), bins=bin_edge)[0]
        / np.histogram(wl, weights=(w), bins=bin_edge)[0]
    )
    vals_omega = (
        np.histogram(wl, weights=(w * omega), bins=bin_edge)[0]
        / np.histogram(wl, weights=(w), bins=bin_edge)[0]
    )
    Tbin = Table(
        data=[bin_edge[:-1], bin_edge[1:], vals_gain, vals_zodi, vals_omega],
        names=["wl_min", "wl_max", "a_bin", "z_bin", "omega_bin"],
        units=["micron", "micron", "", "", "sr"],
    )
    Tbin.meta["VER"] = (version, "Calibration version")
    Tbin.meta["DETECTOR"] = (detector_id, "SPHEREx detector")
    Tbin.meta["NBINS"] = (nbins, "Number of bins")
    for key in binned_catalog_descriptions.keys():
        if key in Tbin.colnames:
            Tbin[key].description = binned_catalog_descriptions[key]
    # Binned Calibration data table
    return Tbin


def load_spline_from_file(filename: str) -> UnivariateSpline:
    """Load a UnivariateSpline object from a file.

    Parameters
    ----------
    filename : `str`
        The filename where the spline is saved.

    Returns
    -------
    spline : `UnivariateSpline`
        The reconstructed spline.
    """
    # Load the knots, coefficients, and degree
    data = np.load(filename)
    knots = data['knots']
    coeffs = data['coeffs']
    degree = int(data['degree'])
    # Reconstruct the spline
    return UnivariateSpline._from_tck((knots, coeffs, degree))


class UpdatePhotometricCalibrationConnections(pipeBase.PipelineTaskConnections,
                                              dimensions=("instrument", "detector"),
                                              defaultTemplates={}):
    ff_prev = cT.Input(
        name='flat',
        doc="Relative flat field estimate.",
        storageClass="SPHERExImage",
        dimensions=("instrument", "detector"),
        isCalibration=True
    )
    initial_distortion_table = cT.Input(
        name="initial_distortion_table",
        doc="Initial data to fit distortion from ray-tracing.",
        storageClass="InitDistTable",
        dimensions=["instrument"],
        isCalibration=True
    )
    sip = cT.Input(
        name="sip",
        doc="Output distortion as SIP (Simple Imaging Polynomial).",
        storageClass="AstropySIP",
        dimensions=["instrument", "detector"],
        isCalibration=True
    )
    photo_calib = cT.Output(
        name='photo_calib',
        doc="Input Binned Photometric Calibration data.",
        storageClass="SpxPhotoCalib",
        dimensions=("instrument", "detector"),
        isCalibration=True
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class UpdatePhotometricCalibrationConfig(pipeBase.PipelineTaskConfig,
                                         pipelineConnections=UpdatePhotometricCalibrationConnections):
    param_int = Field(
        doc="An example of int configuration parameter."
        "An example of a two-line doc.",
        dtype=int,
        default=42)

    param_float = Field(
        "An example of float configuration parameter.",
        dtype=float,
        default=4.2)

    param_dict = DictField(
        doc="An example of dictionary configuration parameter.",
        keytype=int,
        itemtype=float,
        default={1: 12.5, 2: 12.5, 3: 12.5, 4: 10.5, 5: 10.5, 6: 10.5})

    calspec_spx_path = Field(
        doc=("Path to the CalSpec SPHEREx observation catalog FITS table. "),
        dtype=str,
        default=""
    )


@dataclass
class Output:
    image: SPHERExImage


class UpdatePhotometricCalibrationTask(pipeBase.PipelineTask):
    """Implementation of L3 Update Photometric Calibration.

    Attributes
    ----------
    ConfigClass : `lsst.pex.config.Config`
        Task configuration class.
    full_name : `str`
        Full module name as it appears in SSDC document.
    ssdc_pipe_level: `int`
        Module pipeline level.
    ssdc_doc_tag: `str`
        SPHEREx Science Data Center (SSDC) document ID.
        Document GitHub page: https://github.com/SPHEREx/SSDC-MS-037

    Notes
    -----
    - Use config to customize algorithm parameters.
    - Pass module inputs as arguments of the ``run`` method.
    - Module outputs are returned by the ``run`` method as
      attributes of a `~dataclasses.dataclass` object.
    """

    _DefaultName = "update-photometric-calibration"
    ConfigClass = UpdatePhotometricCalibrationConfig
    full_name = "Update Photometric Calibration"
    ssdc_pipe_level = 3
    ssdc_doc_tag = "SSDC-MS-037, v1.5"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self,
            calspec_sources: Table = None,
            calspec_spx: Table = None,
            ff_prev: np.array = None,
            ff_prev_spline: UnivariateSpline = None,
            initial_distortion_table: Table = None,
            sip: Sip = None) -> pipeBase.Struct:
        """
        Notes
        ----
        - calspec_sources provides metadata about the sources
        - calspec_spectra provides the high-resolution spectra for those sources.
        """
        # Reload the spline if not provided
        if ff_prev_spline is None:
            spline_file = get_data_file("data/splfit_data.npz", package="spx_utils")
            ff_prev_spline = load_spline_from_file(spline_file)

        detector_id = ff_prev.getArrayId()
        # Load defaults if not provided
        if calspec_sources is None:
            calspec_sources = Table.read(data_filename("data/templates/calspec/calspec_catalog.fits",
                                                       "SPHEREx_SkySimulator"))
        # Primary catalog
        if calspec_spx is None:
            if not self.config.calspec_spx_path:
                raise ValueError(
                    "config.calspec_spx_path must be set to an absolute path "
                    "when calspec_spx is not provided."
                )

            path = Path(self.config.calspec_spx_path).expanduser()

            if not path.is_file():
                raise FileNotFoundError(f"CalSpec SPX catalog not found: {path}")

            self.log.info(f"Loading CalSpec SPX catalog from {path}")
            calspec_spx = Table.read(path)

        logging.info("Loading CalSpec sources to simulate")
        # Each file in calspec_templates is a FITS file that
        # contains the high-resolution spectrum of a source.
        # TODO use `data_filename`, see: SP-1153
        calspec_spectra = sorted([
            str(file)
            for file in importlib.resources.files("spx_utils.data.templates.calspec").iterdir()
            if file.name.endswith(".fits") and len(file.name.split('_')[0]) == 3])

        logging.info(f'calspec_template length = {len(calspec_spectra)}')

        if not calspec_spectra:
            raise FileNotFoundError("No CalSpec template files found. Please ensure the "
                                    "data/templates/calspec/ directory contains valid files.")

        calspec_sources_pred = CalspecPredict()

        # Loop through templates and simulate sources
        for L in calspec_spectra:
            TT = Table.read(L)
            mag = (TT["fnu"]).to(u.ABmag)
            for row in calspec_sources:
                if row["name"] in L:
                    name = row["name"]
                    if name in primary_calibrator_names:
                        logging.info(f"{name} {np.min(mag):0.2f}")
                        TT.rename_column("wavelength", "lambda")
                        TT.rename_column("fnu", "FLUX")

                        # Use tempfile to create a temporary file
                        with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as tmp_file:
                            TT.write(tmp_file.name, overwrite=True)
                            calspec_sources_pred.load_single(name,
                                                             row["ra"] * u.deg,
                                                             row["dec"] * u.deg,
                                                             tmp_file.name)

        # Use in place of All-Sky-Spectral Catalog (primary version)
        # Use a subset of the large dataset for testing

        # Initialize models
        bandpass_model = SPinst.Tabular_Bandpass(max_wl=None, min_wl=None)
        zodi_model = SPsky.zodicalc.ModifiedKelsallModelWithHPFT()

        # Run photometric calibration
        # NOTE currently keeping this during validation.
        # obstab = calc_pred_flux_and_color_corr_local(calspec_sources_pred,
        #                                       calspec_spx,
        #                                       bandpass_model,
        #                                       zodi_model,
        #                                       wl_unit=u.um)

        # Run photometric calibration
        obstab = calc_pred_flux_and_color_corr_parallel(calspec_sources_pred,
                                                        calspec_spx,
                                                        bandpass_model,
                                                        zodi_model,
                                                        wl_unit=u.um)

        # Apply flux engineering units
        obstab_with_engineeing_units = flux_engineering_units_local(obstab,
                                                                    ff_prev,
                                                                    ff_prev_spline,
                                                                    detector_id,
                                                                    initial_distortion_table,
                                                                    sip)

        # Compute correction gain and weight
        tab = compute_corr_gain_and_weight_local(obstab_with_engineeing_units)

        # Generate binned calibration data
        binned_calib_data = generate_binned_calibration_data_local(tab, bandpass_model, detector_id)

        # Create an Output object
        output = pipeBase.Struct(tab=tab, photo_calib=binned_calib_data)

        return output
