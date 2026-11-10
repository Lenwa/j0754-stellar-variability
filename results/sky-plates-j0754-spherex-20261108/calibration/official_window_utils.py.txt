"""Utility functions for performing photometry on image subregions.
"""

__all__ = ["midpoint", "slide_window", "estimate_loc_bkg", "run_tractor_window", "combine_measurements",
           "combine_fits"]

# Standard library imports
import itertools
from collections import namedtuple
from typing import NamedTuple

# Related third party imports
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D, VarianceUncertainty
from astropy.stats import SigmaClip
from astropy.table import Table
from astropy.units import CompositeUnit, Quantity
from astropy.wcs import WCS
from lsst.pex.config import Config
from photutils.background import Background2D
import numpy as np
from numba import njit

# Local application/ library specific imports
from spx_data import PSFZoneData, SPHERExImage
from spx_utils.photometry.tractor_utils import tractor_fit_spherex
from spx_utils import SourcePipelineClass


# Memory-efficient tuple subclass to store window metadata
#   cut_<axis><N>:           central cutout boundaries
#   window_<axis><N>:        window (cutout + frame edge) boundaries
#   bkgd_<axis><N>:          local background estimation cutout boundaries
#   x1, x2:                  left/right edge
#   y1, y2:                  lower/upper edge
WindowInfo = namedtuple(
    "WindowInfo",
    ['row', 'col', 'cut_x1', 'cut_x2', 'cut_y1', 'cut_y2', 'window_x1', 'window_x2', 'window_y1',
     'window_y2', 'bkgd_x1', 'bkgd_x2', 'bkgd_y1', 'bkgd_y2'])


@njit(cache=True)
def _compute_mask(src_ximg, src_yimg, x1, x2, y1, y2):
    """
    Numba JIT-compiled mask function (lazy compile on first call)
    to return window mask
    src_ximg, src_yimg : `numpy.ndarray`
        Source positions in X, Y pixel coordinates on the detector, calculated
        assuming zero-based indexes.
    x1,x2,y1,y2: `float`
        Window parameter offset positions.
    """

    n = src_ximg.shape[0]
    mask = np.empty(n, dtype=np.bool_)
    for i in range(n):
        mask[i] = (src_ximg[i] > x1) & (src_ximg[i] < x2) & \
                  (src_yimg[i] > y1) & (src_yimg[i] < y2)
    return mask


def window_masking_np(src_ximg: np.ndarray, src_yimg: np.ndarray, src_cat_np: np.ndarray,
                      window_params: NamedTuple):
    """
    src_ximg, src_yimg : `numpy.ndarray`
        Source positions in X, Y pixel coordinates on the detector, calculated
        assuming zero-based indexes.
    src_cat_np : `np.ndarray`
        Numpy ndarray representation of reference catalog sources within
        or close to the image. Includes HEALPix index for each source.
    window_params : `collections.namedtuple`
        Named tuple with field names corresponding to different window
        parameters. All X, Y values are pixel coordinates on the detector.
    """
    x1 = window_params.window_x1 - 0.5
    x2 = window_params.window_x2 - 0.5
    y1 = window_params.window_y1 - 0.5
    y2 = window_params.window_y2 - 0.5

    # Mask is computed at C speed after first compile
    window_mask = _compute_mask(src_ximg, src_yimg, x1, x2, y1, y2)

    # Slice NumPy structured array (fast)
    window_cat_np = src_cat_np[window_mask]

    return window_mask, window_cat_np


def midpoint(x1: float, x2: float, y1: float, y2: float) -> tuple[float, float]:
    """Calculates the midpoint of a window.

    Parameters
    ----------
    x1 : `float`
        Initial X position/left side of the window.
    x2 : `float`
        Final X position/right side of the window.
    y1 : `float`
        Initial Y position/bottom of the window.
    y2 : `float`
        Final Y position/top of the window.

    Returns
    -------
    `tuple`[`float`, `float`]
        The X, Y position of the midpoint. All values are returned in the same
        coordinate system assumed for the inputs.
    """
    return ((x1+x2)/2.0, (y1+y2)/2.0)


def slide_window(spximg: SPHERExImage, config: Config, use_local: bool) -> np.ndarray:
    """Slides a window across an image to identify boundaries for photometry.

    Parameters
    ----------
    spximg : `spx_data.SPHERExImage`
        LVF image.
    config : `lsst.pex.config.Config`
        Contains fields with sliding window parameters.
    use_local : `bool`
        Whether or not to determine bounds for local background estimation.

    Returns
    -------
    img_windows : `np.ndarray`
        Array of `collections.namedtuple` elements. The field names of each
        named tuple correspond to different window parameters. All X, Y values
        are pixel coordinates on the detector.
    """
    # NOTE: FITS NAXIS1 (NAXIS2) is numpy array axes 1 (0). See
    # https://docs.astropy.org/en/stable/io/fits/usage/image.html for more.
    cut_x1_vals = np.arange(0, spximg.data.shape[1], config.cut_side_length)
    cut_y1_vals = np.arange(0, spximg.data.shape[0], config.cut_side_length)
    num_x_windows = len(cut_x1_vals)
    num_y_windows = len(cut_y1_vals)
    bkg_border = config.bkg_side_length - config.cut_side_length

    window_idx = 0
    img_windows = np.empty(shape=(num_x_windows * num_y_windows), dtype=WindowInfo)
    # Slide window along Y-axis
    for cut_y1, y_window in zip(cut_y1_vals, range(num_y_windows)):
        cut_y2 = cut_y1 + config.cut_side_length
        if y_window == 0:  # along bottom edge
            window_y1 = cut_y1  # no frame below
            window_y2 = cut_y2 + config.frame_width
            bkgd_y1 = cut_y1 if use_local else np.nan
            bkgd_y2 = cut_y2 + bkg_border if use_local else np.nan
        elif y_window == num_y_windows - 1:  # along top edge
            window_y1 = cut_y1 - config.frame_width
            window_y2 = cut_y2  # no frame above
            bkgd_y1 = cut_y1 - bkg_border if use_local else np.nan
            bkgd_y2 = cut_y2 if use_local else np.nan
        else:
            window_y1 = cut_y1 - config.frame_width
            window_y2 = cut_y2 + config.frame_width
            bkgd_y1 = cut_y1 - bkg_border/2.0 if use_local else np.nan
            bkgd_y2 = cut_y2 + bkg_border/2.0 if use_local else np.nan

        # Slide window along X-axis
        for cut_x1, x_window in zip(cut_x1_vals, range(num_x_windows)):
            cut_x2 = cut_x1 + config.cut_side_length
            if x_window == 0:  # along left edge
                window_x1 = cut_x1  # no frame to left
                window_x2 = cut_x2 + config.frame_width
                bkgd_x1 = cut_x1 if use_local else np.nan
                bkgd_x2 = cut_x2 + bkg_border if use_local else np.nan
            elif x_window == num_x_windows - 1:  # along right edge
                window_x1 = cut_x1 - config.frame_width
                window_x2 = cut_x2  # no frame to right
                bkgd_x1 = cut_x1 - bkg_border if use_local else np.nan
                bkgd_x2 = cut_x2 if use_local else np.nan
            else:
                window_x1 = cut_x1 - config.frame_width
                window_x2 = cut_x2 + config.frame_width
                bkgd_x1 = cut_x1 - bkg_border/2.0 if use_local else np.nan
                bkgd_x2 = cut_x2 + bkg_border/2.0 if use_local else np.nan

            # Populate array with corresponding window parameters
            img_windows[window_idx] = WindowInfo(
                y_window, x_window, cut_x1, cut_x2, cut_y1, cut_y2, window_x1, window_x2, window_y1,
                window_y2, bkgd_x1, bkgd_x2, bkgd_y1, bkgd_y2
            )
            window_idx += 1
    return img_windows


def estimate_loc_bkg(window_params: NamedTuple, img_data: np.ndarray, img_units: CompositeUnit,
                     img_wcs: WCS, config: Config) -> Cutout2D:
    """Estimates the local background for a window within a crowded image.

    Parameters
    ----------
    window_params : `collections.namedtuple`
        Named tuple with field names corresponding to different window
        parameters. All X, Y values are pixel coordinates on the detector.
    img_data : `numpy.ndarray`
         Crowded SPHEREx image CCD data with no global background solution.
    img_units : `astropy.units.core.CompositeUnit`
        SPHEREx image CCD data units.
    img_wcs : `astropy.wcs.WCS`
        Object containing world coordinate system keywords.
    config : `lsst.pex.config.Config`
        Contains fields with local background estimation parameters.

    Returns
    -------
    bkgd_cut : `astropy.nddata.Cutout2D`
        Cutout of the local background estimation result with the same
        origin/size as the window cutout.
    """
    # Create larger cutout for background estimation region
    est_region_ctr = midpoint(window_params.bkgd_x1, window_params.bkgd_x2,
                              window_params.bkgd_y1, window_params.bkgd_y2)
    est_region_cut = Cutout2D(Quantity(value=img_data, unit=img_units, copy=False), position=est_region_ctr,
                              size=(config.bkg_side_length, config.bkg_side_length), wcs=img_wcs)

    # Estimate local background
    sigma_clip = SigmaClip(sigma_lower=config.bkg_sc_sigma_lower, sigma_upper=config.bkg_sc_sigma_upper,
                           maxiters=config.bkg_sc_maxiters)
    bkgd = Background2D(est_region_cut.data, box_size=config.bkg_box_size, sigma_clip=sigma_clip,
                        exclude_percentile=config.bkg_exclude_pct)

    # Create local background cutout overlapping with window cutout
    win_relative_bkgd = [window_params.window_x1-window_params.bkgd_x1,
                         window_params.window_x2-window_params.bkgd_x1,
                         window_params.window_y1-window_params.bkgd_y1,
                         window_params.window_y2-window_params.bkgd_y1]
    bkgd_ctr = midpoint(*win_relative_bkgd)
    size = ((window_params.window_y2 - window_params.window_y1),
            (window_params.window_x2 - window_params.window_x1))
    bkgd_cut = Cutout2D(Quantity(value=bkgd.background, unit=img_units, copy=False), position=bkgd_ctr,
                        size=size, wcs=est_region_cut.wcs)

    return bkgd_cut


def run_tractor_window(img_data: np.ndarray, img_units: CompositeUnit, img_wcs: WCS,
                       img_var: VarianceUncertainty, img_flags: np.ndarray, window_params: NamedTuple,
                       src_ximg: np.ndarray, src_yimg: np.ndarray, src_cat: np.ndarray, cat_coords: SkyCoord,
                       use_local: bool, config: Config, psf: PSFZoneData, detector: int,
                       keep_frame_sources: bool = False):
    """Sets up input parameters and runs Tractor fitting for a window.

    Parameters
    ----------
    img_data : `numpy.ndarray`
        SPHEREx image CCD data.
    img_units : `astropy.units.core.CompositeUnit`
        SPHEREx image CCD data units.
    img_wcs : `astropy.wcs.WCS`
        Object containing world coordinate system keywords.
    img_var : `astropy.nddata.nduncertainty.VarianceUncertainty`
        SPHEREx image variance.
    img_flags : `numpy.ndarray`
        Flags giving information about each pixel.
    window_params : `collections.namedtuple`
        Named tuple with field names corresponding to different window
        parameters. All X, Y values are pixel coordinates on the detector.
    src_ximg, src_yimg : `numpy.ndarray`
        Source positions in X, Y pixel coordinates on the detector, calculated
        assuming zero-based indexes.
    src_cat : `np.ndarray`
        Numpy ndarray representation of reference catalog sources within
        or close to the image. Includes HEALPix index for each source.
    cat_coords : `astropy.coordinates.SkyCoord`
        ICRS RA/DEC of nearby reference catalog sources, in degrees.
    use_local : `bool`
        Whether or not to perform local background subtraction.
    config : `lsst.pex.config.Config`
        Contains fields with local background estimation parameters.
    psf : `spx_data.PSFZoneData`
        SPHEREx PSF.
    detector : `int`
        Detector number (1-6).
    keep_frame_sources : `bool`
        If `True`, returns photometry for sources in the outer frame of the
        cutout. Default is `False`.

    Returns
    -------
    tractor_results : `numpy.ndarray` or `dict`
        Returns array of photometry measurements by default. If
        `config`.``return_fit`` is ``True``, returns dictionary with flux
         results and best-fit model/residual image.

    """
    # Identify sources to photometer
    # NOTE: Subtract 0.5 from edges since input source positions assume zero-
    # based indexing.
    window_mask, window_cat_np = window_masking_np(src_ximg, src_yimg, src_cat, window_params)
    # window_cat = Table(window_cat_np, copy=False)
    if len(window_cat_np) > 0:
        window_cat = Table(window_cat_np, copy=False)
        # Add columns for pixel positions
        window_cat.add_column(col=src_ximg[window_mask], name="x_image", copy=False)
        window_cat.add_column(col=src_yimg[window_mask], name="y_image", copy=False)

        # Flag sources in central cutout
        central_mask = ((window_cat["x_image"] >= window_params.cut_x1 - 0.5)
                        & (window_cat["x_image"] <= window_params.cut_x2 - 0.5)
                        & (window_cat["y_image"] >= window_params.cut_y1 - 0.5)
                        & (window_cat["y_image"] <= window_params.cut_y2 - 0.5))
        window_cat.add_column(col=central_mask, name="central_source", copy=False)

        # Filter source classes if necessary
        if config.class_filter != SourcePipelineClass.ALL.value:
            window_cls_flgs = set(itertools.chain.from_iterable(
                SourcePipelineClass(int(val)).get_names_as_list() for val in window_cat["spherex_class"]))
            filter_names = SourcePipelineClass(config.class_filter).get_names_as_list()
            class_in_window = any(spx_class in window_cls_flgs for spx_class in filter_names)
        else:
            class_in_window = True

        if (len(window_cat[central_mask]) > 0) and class_in_window:
            # Create image cutout
            window_center = midpoint(window_params.window_x1, window_params.window_x2,
                                     window_params.window_y1, window_params.window_y2)
            size = ((window_params.window_y2 - window_params.window_y1),
                    (window_params.window_x2 - window_params.window_x1))
            # NOTE: position in (x, y) order, but size in (ny, nx) order
            cutout = Cutout2D(Quantity(value=img_data, unit=img_units, copy=False), position=window_center,
                              size=size, wcs=img_wcs)

            # Subtract local background for relevant pointings
            if use_local:
                bkgd_cutout = estimate_loc_bkg(window_params, img_data, img_units, img_wcs, config)
                cutout.data = cutout.data - bkgd_cutout.data

            # Create variance and flag subregions
            slice_y, slice_x = cutout.slices_original
            var_slice = img_var.array[slice_y, slice_x]
            flags_slice = img_flags[slice_y, slice_x]

            # Convert positions to X, Y pixel coordinates in cutout
            window_coords = cat_coords[window_mask]
            src_xwin, src_ywin = (cutout.wcs).world_to_pixel(window_coords)
            # Check expected sources are actually within cutout
            # NOTE: Subtract/add 0.5 because <axis min/max>_cutout are indices
            cut_mask = ((src_xwin > cutout.xmin_cutout - 0.5)
                        & (src_xwin < cutout.xmax_cutout + 0.5)
                        & (src_ywin > cutout.ymin_cutout - 0.5)
                        & (src_ywin < cutout.ymax_cutout + 0.5))
            window_cat = window_cat[cut_mask]
            window_cat.add_column(col=src_xwin[cut_mask], name="x_window", copy=False)
            window_cat.add_column(col=src_ywin[cut_mask], name="y_window", copy=False)

            # Create generative model using Tractor
            window_origin = (window_params.window_x1, window_params.window_y1)
            tractor_results = tractor_fit_spherex(
                spx_cutout=cutout, spx_variance=var_slice, spx_flags=flags_slice,
                spx_bkgd=bkgd_cutout if use_local else None,
                source_list=window_cat, detector=detector, cutout_origin=window_origin, psf=psf,
                hips_name=f"hips{config.hips_order}",
                keep_frame_sources=keep_frame_sources, return_fit=config.return_fit,
                fit_shapes=config.fit_shapes)
            return tractor_results


def combine_measurements(combined_arr: np.ndarray, result: dict | np.ndarray | None, start: int,
                         end: int) -> tuple[int, int]:
    """Combine Tractor photometry measurements from a single window into an
    array.

    Parameters
    ----------
    combined_arr : `numpy.ndarray`
        Array to hold photometry measurements for a whole image. It can be
        empty to start, or it may include existing indices for sources in other
        windows.
    result : `dict`, `numpy.ndarray`, or `None`
        Result from running Tractor on a single window. If `dict`, the
        `flux_result` key contains the measurement array.
    start : `int`
        Starting index for array slice.
    end : `int`
        Ending index for array slice.

    Returns
    -------
    start, end : `tuple` of `int`
        Modified start/end indices.
    """
    # Filter out windows with no sources
    if result is not None:
        if isinstance(result, dict):
            end = start + len(result["flux_result"])
            combined_arr[start:end] = result["flux_result"]
        else:
            end = start + len(result)
            combined_arr[start:end] = result
        start = end
    return start, end


def combine_fits(tractor_fits: dict, result: dict | None) -> dict:
    """Organize Tractor fit results from a single window into a dictionary.

    Parameters
    ----------
    tractor_fits : `dict`
        Dictionary to hold fitting results for a whole image. It can be empty
        to start, or it may include existing key-value pairs for other windows.
    result : `dict` or `None`
        Result from running Tractor on a single window. If `dict`, the `fit`
        key contains a list, in which the first index holds a position string
        and the second index holds an `astropy.io.fits.HDUList`.

    Returns
    -------
    tractor_fits : `dict`
        Modified dictionary.
    """
    # Filter out windows with no sources:
    if result is not None:
        tractor_fits[result["fit"][0]] = result["fit"][1]
    return tractor_fits
