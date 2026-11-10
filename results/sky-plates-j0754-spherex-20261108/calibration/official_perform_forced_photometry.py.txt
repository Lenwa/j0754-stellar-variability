__all__ = ["PerformForcedPhotometryConfig", "PerformForcedPhotometryTask"]

# Standard library imports
from itertools import repeat
import copy
import gc
import multiprocessing as mp
import time

# Related third party imports
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header, HDUList
from astropy.table import Table, Column
from lsst.daf.butler import DeferredDatasetHandle
from lsst.pex.config import Field, ListField
import numpy as np
import healpy as hp
import astropy.units as u
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT

# Local application/library specific imports
from spx_data import (
    DetectorSpectralInfo, Flag, ObservationId, PSFZoneData, SPHERExImage, SPSFlag)
from spx_data.spx_param import get_npix_x_active, get_npix_y_active
from spx_data.spx_version import SSDC_PIPELINES_VERSION
from spx_utils.photometry.window_utils import (
    slide_window, run_tractor_window, combine_measurements, combine_fits)
from spx_utils.photometry.tractor_utils import MJysr_to_uJy
from spx_utils import SourcePipelineClass, get_current_rss
from spx_utils.utils import get_refcat_provenance


# ==============================================================================
# TODO: Incorporate x/y pixel position info
def calculate_mjd(hdr: Header):
    mjd = np.average([hdr["MJD-BEG"], hdr["MJD-END"]])
    return mjd


def get_deep_flag(hdr: Header):
    if "SFLAG" in hdr:
        # one-year simulated images
        sps_flag_int = SPSFlag.from_string(hdr["SFLAG"])
    elif "SPS_FLAGS" in hdr:
        # SSDC level1 and level2 images
        sps_flag_int = SPSFlag.from_string(hdr["SPS_FLAGS"])
    else:
        sps_flag_int = 0

    return SPSFlag.is_deep_field(sps_flag_int)


def check_cut_side_length(value: int) -> bool:
    """Check that cutout side length (in pixels) is a factor of the number of
       active pixels along a side.
    """
    img_side_length = get_npix_x_active()
    if (img_side_length % value == 0):
        return True
    else:
        return False


def add_srccat_pass_thru(idx_map: np.array, input_cat: Table, col_name: str, output_tab: Table):
    """Add a pass-through column from the input catalog to the output
    photometry measurement table.

    Parameters
    ----------
    idx_map : `numpy.ndarray`
        Reorders the input catalog to match the order of the output table.
    input_cat : `astropy.table.Table`
        Filtered, image-level input source catalog.
    col_name : `str`
        Column name in input source catalog.
    output_tab : `astropy.table.Table`
        Output photometry measurement table.
    """
    # Input catalog columns that can be added to the output table
    pass_thru_lookup = {
        "name": {
            "description": "User-specified name of the source ('Target NNN' if no name provided)",
            "dtype": "U",
            "unit": None,
            "idx": 1
        },
        "spherex_class": {
            "description": "Source origin in reference catalog",
            "dtype": "int32",
            "unit": None,
            "idx": 1
        }
    }

    # NOTE: Whether adding or replacing columns, we avoid copying to save
    # memory. If future module changes require a column's input data to be
    # modified *after* the data has been added to the table, switch to use
    # `copy=True`.
    if col_name in pass_thru_lookup.keys():
        output_tab.add_column(
            col=Column(data=input_cat[col_name][idx_map],
                       name=col_name,
                       unit=pass_thru_lookup[col_name]["unit"],
                       description=pass_thru_lookup[col_name]["description"],
                       dtype=pass_thru_lookup[col_name]["dtype"]),
            index=pass_thru_lookup[col_name]["idx"], copy=False
        )


class PerformForcedPhotometryConnections(pipeBase.PipelineTaskConnections,
                                         dimensions=("instrument", "exposure"),
                                         defaultTemplates={}):
    spximg = cT.Input(
        name="spximg",  # input is output from previous module
        doc="Input image.",
        storageClass="SPHERExImage",
        dimensions=("instrument", "exposure"),
    )
    global_bkgd = cT.Input(
        name="global_bkgd",  # input is output from previous module
        doc="Estimated global background.",
        storageClass="NumpyArray",
        dimensions=("instrument", "exposure"),
    )
    srccat = cT.Input(
        name="catalog",
        doc="Input catalog.",
        storageClass="AstropyParquet",
        dimensions=("instrument", "dichroicpair", "spsplan", "slewgroup"),
        deferLoad=True,
    )
    solid_angle_pixel_map = cT.PrerequisiteInput(
        name="solid_angle_pixel_map",
        doc="Input solid angle pixel map.",
        storageClass="SPHERExImage",
        dimensions=["instrument", "detector"],
        isCalibration=True
    )
    spectral_wcs = cT.PrerequisiteInput(
        name="spectral_wcs",
        doc="Spectral WCS with lookup table.",
        storageClass="AstropyImage",
        dimensions=["instrument", "detector"],
        isCalibration=True,
    )
    phot_meas = cT.Output(
        name="phot_meas",
        doc="Output photometry measurement table (per HEALPix).",
        storageClass="AstropyParquet",
        # dimensions should include a table with spatial dimension (exposure)
        # skypix dimension should match commonSkyPix in dimensions.yaml
        dimensions=("instrument", "exposure", "hips6"),
        multiple=True,
    )
    phot_meas_combined = cT.Output(
        name="phot_meas_combined",
        doc="Output photometry measurement table (full image).",
        storageClass="AstropyParquet",
        dimensions=("instrument", "exposure"),
    )
    tractor_fits = cT.Output(
        name="tractor_fits",
        doc="Tractor results",
        storageClass="StructuredDataDict",
        dimensions=("instrument", "exposure"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if config.combinedOutput:
            del self.phot_meas
        else:
            del self.phot_meas_combined

        if not config.return_fit:
            del self.tractor_fits

        if config.use_simulator_bandpass:
            del self.spectral_wcs


class PerformForcedPhotometryConfig(pipeBase.PipelineTaskConfig,
                                    pipelineConnections=PerformForcedPhotometryConnections):
    hips_order = Field(
        doc="HEALPix order for grouping the measurements",
        dtype=int,
        default=6)

    combinedOutput = Field(
        doc="If True, output the photometry measurements for all sources in an image in one table."
            " Default is False.",
        dtype=bool,
        default=False)

    pass_thru_cols = ListField(
        doc="List of columns to pass through to the output catalog, in order of relative placement in the"
            "output catalog.",
        dtype=str,
        default=[])

    gc_interval = Field(
        doc="Triggers manual garbage collection very N iterations during window loop. Default is 0 (no"
            " garbage collection).",
        dtype=int,
        default=0)

    cat_phot_cols = ListField(
        doc="Minimal list of input source catalog columns required to perform photometry.",
        dtype=str,
        default=["SPHERExRefID", "ra", "dec", "ra_error", "dec_error", "ref_ra", "ref_dec", "sersic",
                 "shape_r", "ab_ratio", "phi", "for_l3phot"]
    )

    merge_srccat = Field(
        doc="If True, adds `pass_thru_cols` from the input source catalog to the output photometry table."
            " Default is False.",
        dtype=bool,
        default=False)

    # Sliding window parameters
    cut_side_length = Field(
        doc="Side length (in pixels) of the NxN square cutouts used to step across the detector. Sources "
            "inside the cutout region are included in the Tractor model, and their photometry results are "
            "stored in the output of this module.",
        dtype=int,
        check=check_cut_side_length,
        default=5)

    frame_width = Field(
        doc="Number of extra pixels around a cutout to include in Tractor fitting. Sources in the frame "
            "region are included in the Tractor model, but their photometry results are not recorded by "
            "default. The frame is added along 2 (corner), 3 (edge), or 4 (inner) sides. For example, the "
            "window size for a corner with `cut_side_length`=5 and `frame_width`=3 would be 8x8 pixels.",
        dtype=int,
        default=3)

    parallel_fitting = Field(
        doc="Whether to use parallel window fitting. Only suggested for local testing with single images "
            "using Python scripts/notebooks. Default is False.",
        dtype=bool,
        default=False)

    class_filter = Field(
        doc="Filters target windows to photometer based on `spherex_class` from SPHEREx Reference Catalog. "
            "Default is `SourcePipelineClass.ALL`.",
        dtype=int,
        default=SourcePipelineClass.ALL)

    # Local background estimation parameters
    bkg_side_length = Field(
        doc="Side length (in pixels) of the NxN square cutouts used to estimate the"
            " local background.",
        dtype=int,
        default=15)

    bkg_sc_sigma_lower = Field(
        doc="Number of standard deviations to use as the lower bound for the background clipping limit.",
        dtype=float,
        default=4)

    bkg_sc_sigma_upper = Field(
        doc="Number of standard deviations to use as the upper bound for the background clipping limit.",
        dtype=float,
        default=4)

    bkg_sc_maxiters = Field(
        doc="Maximum number of background sigma-clipping iterations to perform.",
        dtype=int,
        default=10)

    bkg_filter_threshold = Field(
        doc="The threshold value used for selective median filtering of"
            " the low-resolution 2D background map."
            " The median filter will be applied to only the background boxes"
            " with values larger than filter_size."
            " A value of 1 or (1, 1) disables filtering.",
        dtype=int,
        optional=True,
        default=None)

    bkg_box_size = Field(
        doc="The size of the square box used along each axis. The data should"
            " be covered by an integer number of boxes.",
        dtype=int,
        default=3)

    bkg_filter_size = ListField(
        doc="The window size of the 2D median filter to apply to the low-resolution"
            " background map. A value of 1 or (1, 1) disables filtering.",
        dtype=int,
        default=(3, 3))

    bkg_exclude_pct = Field(
        doc="Percentage of masked pixels allowed in a box. If the percent of"
            " masked pixels exceeds this threshold, the box will be excluded"
            " from the low-resolution map.",
        dtype=float,
        default=50.0)

    # Fitting settings
    fit_shapes = Field(
        doc="If True, treat galaxies as shapes during modeling. If False, "
            "treat galaxies as point sources. Default is True.",
        dtype=bool,
        default=True)

    return_fit = Field(
        doc="If True, fit model and residual image are returned. Default is"
            " False.",
        dtype=bool,
        default=False)

    use_simulator_bandpass = Field(
        doc="Use simulator wavelength model.",
        dtype=bool,
        default=False)

    def update(self, **kw):
        """Override default update method to validate upon field changes.
        """
        super().update(**kw)
        self.validate()

    def validate(self):
        """Override default validate method to include inter-field checks.
        """
        super().validate()

        # Frame width is smaller than cutout
        frame_size_check = (self.frame_width < self.cut_side_length)
        if not frame_size_check:
            raise ValueError("frame_width must be less than cut_side_length.")

        # Cutout for local BG estimation is larger than photometry cutout
        bkg_check_1 = (self.bkg_side_length > self.cut_side_length)
        if not bkg_check_1:
            raise ValueError("bkg_side_length must be greater than cut_side_length. Received:"
                             f" bkg_side_length={self.bkg_side_length},"
                             f" cut_side_length={self.cut_side_length}")

        # Difference in BG estimation size and cutout size is an even number
        # (for non-corner cutouts, proper Background2D and Cutout2D behavior
        #  requires the background "frame" on each side to be an integer)
        bkg_check_2 = ((self.bkg_side_length - self.cut_side_length) % 2.0 == 0)
        if not bkg_check_2:
            raise ValueError("Difference between bkg_side_length and cut_side_length must be an"
                             f" even number. Received: bkg_side_length={self.bkg_side_length},"
                             f" cut_side_length={self.cut_side_length}")

        # Background box size must be a factor of background side length
        bkg_check_3 = (self.bkg_side_length % self.bkg_box_size == 0)
        if not bkg_check_3:
            raise ValueError("bkg_box_size must be a factor of bkg_side_length. Received:"
                             f" bkg_side_length={self.bkg_side_length},"
                             f" bkg_box_size={self.bkg_box_size}")

        # Non-default filter class requires spherex_class column
        class_missing = "spherex_class" not in self.pass_thru_cols
        filtering_enabled = (self.class_filter != SourcePipelineClass.ALL.value)
        if (filtering_enabled and class_missing):
            raise ValueError("Pass-through column 'spherex_class' must be included to perform class "
                             "filtering.")

        # If merging srccat data, 1+ pass-through columns must be specified
        if self.merge_srccat and (len(self.pass_thru_cols) == 0):
            raise ValueError("To merge data from the srccat, 1+ pass-through columns must be specified.")


class PerformForcedPhotometryTask(pipeBase.PipelineTask):
    """Implementation of L3 Perform Forced Photometry.

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
        Document GitHub page: https://github.com/SPHEREx/SSDC-MS-033

    Notes
    -----
    - Use config to customize algorithm parameters.
    - Pass module inputs as arguments of the ``run`` method.
    - Module outputs are returned by the ``run`` method as
      attributes of a `~dataclasses.dataclass` object.
    """

    _DefaultName = "perform-forced-photometry"
    ConfigClass = PerformForcedPhotometryConfig
    full_name = "Perform Forced Photometry"
    ssdc_pipe_level = 3
    ssdc_doc_tag = "SSDC-MS-033 3.1"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        if self.config.combinedOutput is True:
            super().runQuantum(butlerQC, inputRefs, outputRefs)
        else:
            inputs = butlerQC.get(inputRefs)
            result = self.run(**inputs)
            if result.phot_meas is not None:
                healpix_IDs = result.phot_meas.keys()
                for name, ref_or_lst in outputRefs:
                    if name == "phot_meas":
                        for ref in ref_or_lst:
                            healpix_id = ref.dataId.get("hips6")
                            if healpix_id in healpix_IDs:
                                table = result.phot_meas[healpix_id]
                                butlerQC.put(table, ref)
                    else:
                        # TODO handle other types of output
                        print(f"ERROR: output for connection {name} is not saved")

    def run(self, spximg: SPHERExImage, srccat: DeferredDatasetHandle | Table, global_bkgd: np.ndarray,
            solid_angle_pixel_map: SPHERExImage, spectral_wcs: HDUList = None) -> pipeBase.Struct:
        """Perform Forced Photometry.

        Parameters
        ----------
        spximg : `spx_data.SPHERExImage`
            LVF image with catalog of sources to be measured
        srccat : `astropy.table.Table`|`lsst.daf.butler.DeferredDatasetHandle`
            Catalog of sources within or close to the image.
        global_bkgd : `numpy.ndarray`
            Global background estimate.
        solid_angle_pixel_map : `spx_data.SPHERExImage`
            Solid angle pixel map derived by Update Field Distortion module.
        spectral_wcs : `~astropy.io.fits.HDUList`, optional
            Spectral WCS calibration with central wavelength and bandwidth.
            If not specified, Simulator's wavelength model is used.

        Returns
        -------
        return_obj : `pipeBase.Struct`
            Container including output connection fields. By default, it only
            includes `phot_meas_combined` (`astropy.table.Table` of photometry
            measurements for whole image). If `self.config.return_fit` is
            True, it includes `tractor_fits` (dictionary of
            `astropy.fits.io.HDUList`). If `self.config.combinedOutput` is
            False, it includes `phot_meas` (dictionary of `astropy.table.Table`
            grouped by HEALPix).
        """
        # =====================================================================
        # Prepare inputs for photometry fitting
        task_t0 = time.time()
        self.log.info("Starting forced photometry...")

        # Only read necessary columns from input catalog
        if isinstance(srccat, DeferredDatasetHandle):
            # NOTE: Need to convert `lsst.pex.config.listField.List` to
            # `list` for table item access
            read_cols = list(self.config.cat_phot_cols) + list(self.config.pass_thru_cols) if \
                self.config.merge_srccat else list(self.config.cat_phot_cols)
            try:
                srccat = srccat.get(
                    parameters={
                        "columns": read_cols
                    }
                )
            except KeyError as e:
                raise ValueError(f"Requested column {e} is not in the input catalog.")

        # Extract refcat provenance (version) from input srccat once
        ref_ver = None
        try:
            _, ref_ver = get_refcat_provenance(srccat)
        except Exception:
            # Optional: try the meta dict if your pipeline stores it there
            try:
                _, ref_ver = get_refcat_provenance(getattr(srccat, "meta", {}))
            except Exception:
                ref_ver = None

        if ref_ver:
            # Record in task metadata (shows up in Butler task metadata)
            self.metadata.add("refcat_version", str(ref_ver))

        # Convert positions/their errors to X, Y pixel coordinates on detector
        srctab = copy.copy(srccat)  # avoid directly modifying input cat
        cat_coords = SkyCoord(ra=srctab["ra"], dec=srctab["dec"], unit=u.deg, frame="icrs")
        src_xdet, src_ydet = (spximg.wcs).world_to_pixel(cat_coords)

        # Select sources whose positions fall within image boundaries
        # NOTE: Subtract 0.5 from edges since WCS.world_to_pixel returns
        # positions assuming zero-based indices.
        img_mask = ((src_xdet > 0 - 0.5)
                    & (src_xdet < get_npix_x_active() - 0.5)
                    & (src_ydet > 0 - 0.5)
                    & (src_ydet < get_npix_y_active() - 0.5))
        src_ximg = src_xdet[img_mask]
        src_yimg = src_ydet[img_mask]
        cat_coords = cat_coords[img_mask]
        srctab = srctab[img_mask]

        # If "for_l3phot" column exists, filter out sources not flagged for L3
        if "for_l3phot" in srctab.colnames:
            l3phot_mask = (srctab["for_l3phot"] == 1)
            src_ximg = src_ximg[l3phot_mask]
            src_yimg = src_yimg[l3phot_mask]
            cat_coords = cat_coords[l3phot_mask]
            srctab = srctab[l3phot_mask]
        self.log.debug(get_current_rss("Image-filtered catalog created"))

        # Add HEALPix index for each source, using refcat RA/DEC if available
        # to ensure consistency over time
        hips_order = self.config.hips_order
        hips_col_name = f"hips{hips_order}"
        hips_ra = "ref_ra" if "ref_ra" in self.config.cat_phot_cols and srctab.colnames else "ra"
        hips_dec = "ref_dec" if "ref_dec" in self.config.cat_phot_cols and srctab.colnames else "dec"
        col = Column(data=hp.ang2pix(2 ** hips_order, srctab[hips_ra], srctab[hips_dec], nest=True,
                                     lonlat=True),
                     name=hips_col_name, unit="")  # leaving unit None, makes "deg" unit
        # NOTE: We avoid copying columns to save memory. If future module
        # changes require a column's input data to be modified after the data
        # has been added to the table, switch to use `copy=True`.
        srctab.add_column(col, copy=False)

        # Modify pixels that should be ignored by Tractor
        # NOTE: Tractor doesn't support masking, so resetting these values
        #       effectively forces it to discount pixels in the fit.
        var_nanmax_init = np.nanmax(spximg.uncertainty.array)
        var_nanmean_init = np.nanmean(spximg.uncertainty.array)

        # Flagged pixels (in image and variance)
        ignore_pixel_flags = (Flag.SUR_ERROR | Flag.NONFUNC | Flag.MISSING_DATA | Flag.HOT | Flag.COLD
                              | Flag.NONLINEAR | Flag.PERSIST | Flag.BLOOM | Flag.SNOWBALL | Flag.HALO
                              | Flag.SATELLITE_HALO)
        ignore_mask = spximg.getFlagsAsMask(ignore_pixel_flags)
        spximg.data[ignore_mask] = 0
        spximg.uncertainty.array[ignore_mask] = 1e3 * var_nanmax_init

        # Remaining non-finite pixels (in image and variance)
        n_nonfin_img = np.count_nonzero(~np.isfinite(spximg.data))
        if n_nonfin_img > 0:
            spximg.data[~np.isfinite(spximg.data)] = 0
            self.log.info(f"Assigned {n_nonfin_img} non-finite-valued image pixels to 0.")
        n_nonfin_var = np.count_nonzero(~np.isfinite(spximg.uncertainty.array))
        if n_nonfin_var > 0:
            spximg.uncertainty.array[~np.isfinite(spximg.uncertainty.array)] = 1e3 * var_nanmax_init
            self.log.info(f"Assigned {n_nonfin_var} non-finite-valued variance pixels to initial max "
                          "variance value.")

        # Zeros (in variance)
        n_zeros = np.count_nonzero(spximg.uncertainty.array == 0)
        if n_zeros > 0:
            spximg.uncertainty.array[spximg.uncertainty.array == 0] = var_nanmean_init
            self.log.info(f"Assigned {n_zeros} zero-valued variance pixels to initial mean variance value.")

        # Subtract global background if suitable; if not, use local background
        use_local = np.all(np.isnan(global_bkgd))
        if not use_local:
            spximg.data -= global_bkgd
            self.log.info("Subtracting global background estimate...")
        else:
            self.log.info("Local background estimation and subtraction will occur during window fitting...")

        # Save arguments for data, unit, variance, etc. in memory
        data = spximg.data
        units = spximg.unit
        variance = spximg.uncertainty
        flags = spximg.flags
        wcs = spximg.wcs
        psf = PSFZoneData.from_hdu(spximg.layers["PSF"])
        detector = spximg.getArrayId()
        hdr = spximg.header

        # Only pass the columns relevant for photometry fitting
        relevant_cols = ["SPHERExRefID", "ra", "dec", "shape_r", "ab_ratio", "phi", "sersic",
                         hips_col_name]
        if (
            (self.config.class_filter != SourcePipelineClass.ALL.value)
            and ("spherex_class" in self.config.pass_thru_cols)
        ):
            relevant_cols += ["spherex_class"]
        fit_tab = srctab[relevant_cols]

        # =====================================================================
        # Perform fitting
        # Step across the detector, storing bounds for "windows" or subregions
        img_windows = slide_window(spximg, self.config, use_local)
        self.log.debug(get_current_rss("Array of window objects created"))

        # Set up photometry measurement array/fitting results dictionary
        combined_arr = np.empty(len(srctab), dtype=[("SPHERExRefID", "int64"),
                                                    ("ra", "float64"),
                                                    ("dec", "float64"),
                                                    ("x_image", "float32"),
                                                    ("y_image", "float32"),
                                                    ("mjd", "float64"),
                                                    ("lambda", "float32"),
                                                    ("lambda_width", "float32"),
                                                    ("flux", "float32"),
                                                    ("flux_err", "float32"),
                                                    ("flux_bkg", "float32"),
                                                    ("local_bkg_flg", "bool"),
                                                    ("ext_flg", "bool"),
                                                    ("flags", "int64"),
                                                    ("fit_ql", "float32"),
                                                    ("deep_flg", "bool"),
                                                    ("lvf_id", "int64"),
                                                    ("det_id", "uint8"),
                                                    (hips_col_name, "int32")])
        tractor_fits = {} if self.config.return_fit else None
        start = 0
        end = 0

        self.log.info("Starting window fitting with Tractor...")
        fitting_t0 = time.time()
        fit_tab_np = fit_tab.as_array()
        self.log.debug(get_current_rss("Empty fitting array created"))
        if not self.config.parallel_fitting:
            # Fit windows sequentially
            for i, window in enumerate(img_windows):
                result = run_tractor_window(
                    img_data=data, img_units=units, img_wcs=wcs, img_var=variance, img_flags=flags,
                    window_params=window, src_ximg=src_ximg, src_yimg=src_yimg, src_cat=fit_tab_np,
                    cat_coords=cat_coords, use_local=use_local, config=self.config, psf=psf,
                    detector=detector, keep_frame_sources=False)
                start, end = combine_measurements(combined_arr, result, start, end)
                if self.config.return_fit:
                    tractor_fits = combine_fits(tractor_fits, result)
                if (self.config.gc_interval > 0) and (i % self.config.gc_interval == 0):
                    n_objs = gc.collect()
                    self.log.info(f"{n_objs} objects garbage collected")
        else:
            # Fit windows in parallel (for small datasets/quick local tests)
            iterator = zip(repeat(data), repeat(units), repeat(wcs), repeat(variance), repeat(flags),
                           img_windows, repeat(src_ximg), repeat(src_yimg), repeat(fit_tab_np),
                           repeat(cat_coords), repeat(use_local), repeat(self.config), repeat(psf),
                           repeat(detector), repeat(False))
            mp.freeze_support()
            with mp.pool.Pool() as pool:
                for result in pool.starmap(run_tractor_window, iterator):
                    start, end = combine_measurements(combined_arr, result, start, end)
                    if self.config.return_fit:
                        tractor_fits = combine_fits(tractor_fits, result)
        fitting_t1 = time.time()
        # Truncate unfilled elements
        combined_arr = combined_arr[:end]

        # Convert specific intensity to spectral flux density
        pix_map = solid_angle_pixel_map.data << solid_angle_pixel_map.unit  # use << to avoid copying
        ypos_idx, xpos_idx = (spximg.wcs).world_to_array_index(
            SkyCoord(combined_arr["ra"], combined_arr["dec"], unit=u.deg, frame="icrs"))  # keep array order
        # NOTE: We pass an oversample factor of 1 because the Tractor
        # result flux was concentrated in a single 6.2" x 6.2" pixel.
        combined_arr["flux"] = MJysr_to_uJy(combined_arr["flux"]*(u.MJy/u.sr), 1, pix_map, xpos_idx, ypos_idx)
        combined_arr["flux_err"] = MJysr_to_uJy(combined_arr["flux_err"]*(u.MJy/u.sr), 1, pix_map,
                                                xpos_idx, ypos_idx)

        # To save memory and avoid creating more cutouts, populate the flux_bkg
        # column in the global background case all at once
        if not use_local:
            combined_arr["flux_bkg"] = global_bkgd[ypos_idx, xpos_idx]  # replaces 0s

        # =====================================================================
        # Add back in values for non-fitting parameters
        # Get central wavelengths and bandwidths
        # Use position indices initialized in order of combined_array
        if spectral_wcs is None:
            # Use Simulator bandpass model
            spec_info = DetectorSpectralInfo.create(detector)
            central_wls, central_bws = spec_info.get_wl_mid_fwhm(xpos_idx, ypos_idx)
        else:
            # x corresponds to zero-based indexes in x or columns direction
            # y corresponds to zero-based indexes in y or row direction
            central_wls = spectral_wcs["CWAVE"].data[ypos_idx, xpos_idx]
            central_bws = spectral_wcs["CBAND"].data[ypos_idx, xpos_idx]

        # Construct exposure ID
        exp_id_int = ObservationId(hdr["OBSID"]).pointing_id + hdr["DETECTOR"]

        # Quantities based on image/detector data
        combined_arr["mjd"] = calculate_mjd(hdr)
        combined_arr["lambda"] = central_wls
        combined_arr["lambda_width"] = central_bws
        combined_arr["deep_flg"] = get_deep_flag(hdr)
        combined_arr["lvf_id"] = exp_id_int
        combined_arr["det_id"] = detector

        # Construct output table
        out_fmt_arr = np.array([
            ("source_id", "Unique SPHEREx ID", "", "int64"),
            ("ra", "Right ascension in input source catalog", u.deg, "float64"),
            ("dec", "Declination in input source catalog", u.deg, "float64"),
            ("x_image", "x pixel position on LVF image used for fitting", u.pixel, "float32"),
            ("y_image", "y pixel position on LVF image used for fitting", u.pixel, "float32"),
            ("mjd", "Modified Julian Date", u.day, "float64"),
            ("lambda", "SPHEREx wavelengths", u.um, "float32"),
            ("lambda_width", "SPHEREx wavelength widths", u.um, "float32"),
            ("flux", "SPHEREx observed flux", u.uJy, "float32"),
            ("flux_err", "SPHEREx observed flux error", u.uJy, "float32"),
            ("flux_bkg", "SPHEREx background surface brightness at source position", u.MJy/u.sr, "float32"),
            ("local_bkg_flg", "True if background estimation was local", "", "bool"),
            ("ext_flg", "True if treated as extended source", "", "bool"),
            ("flags", "Photometry flags", "", "int64"),
            ("fit_ql", "Photometry metric", "", "float32"),
            ("deep_flg", "True if in Deep Field", "", "bool"),
            ("lvf_id", "ID of LVF exposure used for photometry", "", "int64"),
            ("det_id", "Detector ID", "", "uint8"),
            (hips_col_name, f"HEALPix tile (order {self.config.hips_order})", "", "int32")],
            dtype=[("name", "U20"), ("parameter", "U60"), ("unit", "U10"), ("format", "U10")]
        )
        combined_tab = Table(data=combined_arr,
                             names=out_fmt_arr["name"].tolist(),
                             descriptions=out_fmt_arr["parameter"].tolist(),
                             units=out_fmt_arr["unit"].tolist(),
                             copy=False)

        # Add requested pass-through columns from srccat
        if self.config.merge_srccat and len(self.config.pass_thru_cols) > 0:
            srctab.sort(keys="SPHERExRefID")
            combined_tab.sort(keys="source_id")
            if (
                (self.config.class_filter == SourcePipelineClass.ALL.value)
                and (len(srctab) == len(combined_tab))
            ):
                # Use contiguous indexing to avoid copying
                index_map = slice(None)
            else:
                # Manually find corresponding indices (table lengths unequal)
                index_map = [np.nonzero(srctab["SPHERExRefID"] == combined_tab["source_id"][idx])[0][0] for
                             idx in range(len(combined_tab))]

            for col_name in self.config.pass_thru_cols:
                add_srccat_pass_thru(index_map, srctab, col_name, combined_tab)

        # =====================================================================
        # Add metadata and create return object based on config parameters
        combined_tab.meta["version_level2"] = spximg.header["VERSION"]
        combined_tab.meta["version_phot_meas"] = SSDC_PIPELINES_VERSION

        if ref_ver:
            combined_tab.meta["refcat_version"] = str(ref_ver)

        # Photometry measurement table keyword args
        if self.config.combinedOutput:
            # Group photometry measurements by image
            kwargs = {"phot_meas_combined": combined_tab}
        else:
            # Group photometry measurements by HEALPix
            cat_dict = {}
            grouped_table = combined_tab.group_by(hips_col_name)
            for key, group in zip(grouped_table.groups.keys[hips_col_name], grouped_table.groups):
                group.meta = dict(combined_tab.meta)
                cat_dict[key] = group
            kwargs = {"phot_meas": cat_dict}

        # Fitting results keyword args
        if self.config.return_fit:
            kwargs["tractor_fits"] = tractor_fits

        # Create return object
        return_obj = pipeBase.Struct(**kwargs)

        task_t1 = time.time()
        self.log.info("Fitting time in seconds (minutes): %s (%s)", f"{(fitting_t1-fitting_t0):.2f}",
                      f"{((fitting_t1-fitting_t0)/60.):.2f}")
        self.log.info("Task run time in seconds (minutes): %s (%s)", f"{(task_t1-task_t0):.2f}",
                      f"{((task_t1-task_t0)/60.):.2f}")
        return return_obj
