"""Utility classes and functions for performing Tractor fitting. Adapted
   from SPHEREx-Sky-Simulator/py/SPHEREx_SkySimulator/tractor_utils.py and
   SPHEREx-Sky-Simulator/py/SPHEREx_SkySimulator/quickcatalog.py.
"""

__all__ = ["tractor_fit_spherex", "MJysr_to_uJy"]

# Related third party imports
from astropy.coordinates import SkyCoord
from astropy.io.fits import HDUList, PrimaryHDU, ImageHDU, BinTableHDU, Header
from astropy.nddata import Cutout2D
from astropy.table import Table
from cv2 import INTER_AREA, resize
from regions import CirclePixelRegion, PixCoord
from tractor import (
    ConstantSky, Flux, GalaxyShape, LinearPhotoCal, NullWCS, PixPos, PointSource)
from tractor.psf import lanczos_shift_image
from tractor.sersic import SersicGalaxy, SersicIndex
from tractor.spherex_utils import getUnitFluxModelPatchSPHEREx
import astropy.units as u
import numpy as np
import tractor

# Local application/ library specific imports
from spx_data import Flag, PhotometryFlag, PSFZoneData, spx_param

tractor.galaxy.ProfileGalaxy.getUnitFluxModelPatch = \
    getUnitFluxModelPatchSPHEREx  # override to support arbitrary PSF oversampling


# =============================================================================
# Abstracted Tractor PSF class
def _subimg_center(height, width, subimg):
    subimg_h, subimg_w = subimg.shape

    if height >= subimg_h:
        pad_y0 = (height - subimg_h) // 2
        subimg_cy = pad_y0 + np.floor((subimg_h - 1.0) / 2)
    else:
        pad_y0 = 0
        cut = (subimg_h - height) // 2
        subimg = subimg[cut: cut + height, :]
        subimg_cy = np.floor((subimg_h - 1.0) / 2) - cut

    if width >= subimg_w:
        pad_x0 = (width - subimg_w) // 2
        subimg_cx = pad_x0 + np.floor((subimg_w - 1.0) / 2)
    else:
        pad_x0 = 0
        cut = (subimg_w - width) // 2
        subimg = subimg[:, cut: cut + width]
        subimg_cx = np.floor((subimg_w - 1.0) / 2) - cut

    subimg_h, subimg_w = subimg.shape

    pad = np.zeros((height, width), subimg.dtype)
    pad[pad_y0: pad_y0 + subimg_h, pad_x0: pad_x0 + subimg_w] = subimg

    return pad, subimg_cx, subimg_cy


class SPHERExTractorPSF(tractor.ducks.Params, tractor.ducks.ImageCalibration):
    """Tractor PSF object based on SPHEREx instrument.

    Parameters
    ----------
    psf : `spx_data.psf.PSFZoneData`
        Instantiated SPHEREx PSF object.
    detector : {1, 2, 3, 4, 5, 6}
        ID of the SPHEREx detector. Should be one of 1, 2, 3, 4, 5, and
        6.
    xshift, yshift : `int`, optional
        x- and y-axis shift in the SPHEREx detector pixel coordinates. If
        initializing the PSF for a window cutout, these are the positions of
        the `tractor.Image` origin relative to the SPHEREx detector origin.
        Default is 0.

    Notes
    -----
    Note that this currently does not work with extended sources.

    """
    _psf = None      # instantiated SPHEREx PSF object
    detector = None  # ID of the SPHEREx detector
    xshift = 0       # x-axis shift in the SPHEREx detector pixel coordinates
    yshift = 0       # y-axis shift in the SPHEREx detector pixel coordinates
    radius = 1       # radius of PSF in the unit of pixels

    def __init__(self, spherex_psf, detector, xshift=0, yshift=0):
        if not isinstance(spherex_psf, PSFZoneData):
            raise TypeError(
                "input parameter `spherex_psf` not a "
                "spx_data.psf.PSFZoneData object"
            )
        self._psf = spherex_psf
        self.detector = int(detector)
        self.xshift = xshift
        self.yshift = yshift

        # set the shape of PSF taking into account the oversampling
        self.psf_shape = [
            int(x / self._psf.oversampling) for x in self._psf.get_blended_psf(0, 0).shape
        ]
        # set the radius of PSF taking into account the oversampling
        self.radius = max(self.psf_shape) * 0.5 / self._psf.oversampling

        return

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return (
            f"{self.__class__.__name__}: detector={self.detector}, "
            f"xshift={self.xshift}, yshift={self.yshift}"
        )

    @property
    def shape(self):
        return self.psf_shape

    def copy(self):
        return self.__class__(
            self._psf, self.detector, xshift=self.xshift, yshift=self.yshift
        )

    def getRadius(self):
        return self.radius

    def getPointSourcePatch(
        self, px: float, py: float, minval: float = 0, radius: float = None, modelMask=None, **kwargs
    ):
        """Return a patch of PSF of a point source at the given pixel
        coordinates.

        Notes
        -----
        This obtains high-resolution, oversampled PSF from SPHEREx
        instrument simulator, resamples it for sub-pixel shift, and downscales
        to get the detector-pixel-scaled PSF.

        Parameters
        ----------
        px, py : `float`
            x- and y-axis position of the point source, respectively, in
            `tractor.Image`/SPHEREx window cutout pixel coordinates.
        minval : `float`, optional
        radius : `float`, optional
        modelMask : optional

        Returns
        -------
        patch : `tractor.patch.Patch`
            Patch of PSF
        """

        if radius is not None:
            raise NotImplementedError("input parameter `radius` not " "supported yet")
        if modelMask is not None:
            raise NotImplementedError(
                "input parameter `modelMask` not " "supported yet"
            )

        # get high-resolution PSF
        hires_psf = self._psf.get_blended_psf(px + self.xshift, py + self.yshift)

        # pixel coordinates of the left-lower corner of "first" pixel of
        # full-scale PSF in the scene with resolution same to that of
        # full-scale PSF
        # For simplicity, assume that the "left-lower" corner of a pixel has
        # integer index in this scene.
        # We use a complex equation for displacement from the given pixel
        # coordinates (`px`, `py`), because the centroid of PSF is at the
        # center of left-lower pixel of the exact center of image, if the image
        # is even-shaped.
        x0_hires = (px + 0.5) * self._psf.oversampling - (hires_psf.shape[1] - 1) // 2 - 0.5
        y0_hires = (py + 0.5) * self._psf.oversampling - (hires_psf.shape[0] - 1) // 2 - 0.5

        # integer parts of `x0_hires` and `y0_hires`
        ix0_hires = round(x0_hires)
        iy0_hires = round(y0_hires)

        # float part of `x0_hires` and `y0_hires`
        dx_hires = x0_hires - ix0_hires
        dy_hires = y0_hires - iy0_hires

        # subpixel shift of the hires-scale PSF
        hires_psf = lanczos_shift_image(hires_psf, dx_hires, dy_hires)

        # pad zero-value pixels around the hires-scale PSF to align it to the
        # real resolution scene
        npad_left = ix0_hires % self._psf.oversampling
        hires_psf = np.c_[np.zeros((hires_psf.shape[0], npad_left)), hires_psf]

        npad_bottom = iy0_hires % self._psf.oversampling
        hires_psf = np.r_[np.zeros((npad_bottom, hires_psf.shape[1])), hires_psf]

        if nover := hires_psf.shape[1] % self._psf.oversampling > 0:
            npad_right = self._psf.oversampling - nover
            hires_psf = np.c_[hires_psf, np.zeros((hires_psf.shape[0], npad_right))]

        if nover := hires_psf.shape[0] % self._psf.oversampling > 0:
            npad_top = self._psf.oversampling - nover
            hires_psf = np.r_[hires_psf, np.zeros((npad_top, hires_psf.shape[1]))]

        # pixel coordinates of the "first" pixel of padded hires-scale PSF
        # in the scene with resolution same to the hires-scale PSF
        x0_hires = ix0_hires - npad_left
        y0_hires = iy0_hires - npad_bottom

        # They should be multiple of `self.downsampling_factor`.
        assert x0_hires % self._psf.oversampling == 0
        assert y0_hires % self._psf.oversampling == 0

        # downsampling
        outimg = resize(hires_psf, dsize=None, fx=1/self._psf.oversampling, fy=1/self._psf.oversampling,
                        interpolation=INTER_AREA)
        outimg /= outimg.sum()
        ix0 = x0_hires // self._psf.oversampling
        iy0 = y0_hires // self._psf.oversampling

        return tractor.patch.Patch(ix0, iy0, outimg)

    def getFourierTransformSize(self, radius):
        """Calculate the power-of-two size next to the given `radius`."""
        size = 2 ** int(np.ceil(np.log2(2 * radius)))
        return size

    def getFourierTransform(self, px, py, radius):
        """Fourier transform of PSF, with the next-power-of-2 size up from
        `radius`.

        Parameters
        ----------
        px : float
        py : float
        radius : float

        Returns
        -------
        p : array_like
            Fourier transform of PSF.
        (x0, y0) : tuple of floats
            Pixel coordinates of the PSF center in the PSF subimage.
        (imh, imw) : tuple of ints
            Shape of the padded PSF image.
        (v, w) : tuple of floats
            v=np.fft.rfftfreq(imw), w=np.fft.fftfreq(imh)

        """

        size = self.getFourierTransformSize(radius)

        pad, cx, cy = self._padInImage(size, size, px, py)

        p = np.fft.rfft2(pad)
        p = p.astype(np.complex64)

        pad_height, pad_width = pad.shape

        v = np.fft.rfftfreq(pad_width)
        w = np.fft.fftfreq(pad_height)

        return p, (cx, cy), (pad_height, pad_width), (v, w)

    def _padInImage(self, height, width, px, py, img=None):
        """Embed PSF image into shape (`height`, `width`).

        Parameters
        ----------
        height, width : int
            Embed PSF image into a larger or smaller image of shape, (`height`,
            `width`).
        px, py : float
        img : array_like, optional

        Returns
        -------
        pad : array_like
        cx : int
            *x*-axis coordinate of the PSF center in `pad`.
        cy : int
            *y*-axis coordinate of the PSF center in `pad`.

        Notes
        -----
        Return pixel coordinates follow the FITS standard with 0-based indices.

        """

        if img is None:
            point_source_patch = self.getPointSourcePatch(px, py)
            subimg = point_source_patch.patch
        else:
            subimg = img.copy()

        # Center of PSF is at the "exact" center of the image, `pad`, if the
        # image is odd-shaped.  It is at the center of lower-left pixel (with
        # origin = lower) of the exact center of the image, if the image is
        # even-shaped.

        return _subimg_center(height, width, subimg)

    def getFourierTransformHiRes(self, px, py, radius):

        oversample = self._psf.oversampling
        size = self.getFourierTransformSize(oversample * radius)

        pad, cx, cy = self._padInImageHiRes(size, size, px, py)

        p = np.fft.rfft2(pad)
        p = p.astype(np.complex64)

        pad_height, pad_width = pad.shape

        v = np.fft.rfftfreq(pad_width)
        w = np.fft.fftfreq(pad_height)

        return p, (cx, cy), (pad_height, pad_width), (v, w)

    def _padInImageHiRes(self, height, width, px, py, img=None):
        if img is None:
            # get high-resolution PSF
            psfimg = self._psf.get_blended_psf(px, py)
            subimg = psfimg / psfimg.sum()
        else:
            subimg = img.copy()

        subimg_h, subimg_w = subimg.shape

        # Center of PSF is at the "exact" center of the image, `pad`, if the
        # image is odd-shaped.  It is at the center of lower-left pixel (with
        # origin = lower) of the exact center of the image, if the image is
        # even-shaped.
        if height >= subimg_h:
            pad_y0 = (height - subimg_h) // 2
            subimg_cy = pad_y0 + np.floor((subimg_h - 1.0) / 2)
        else:
            pad_y0 = 0
            cut = (subimg_h - height) // 2
            subimg = subimg[cut:cut + height, :]
            subimg_cy = np.floor((subimg_h - 1.0) / 2) - cut

        if width >= subimg_w:
            pad_x0 = (width - subimg_w) // 2
            subimg_cx = pad_x0 + np.floor((subimg_w - 1.0) / 2)
        else:
            pad_x0 = 0
            cut = (subimg_w - width) // 2
            subimg = subimg[:, cut:cut + width]
            subimg_cx = np.floor((subimg_w - 1.0) / 2) - cut

        subimg_h, subimg_w = subimg.shape

        pad = np.zeros((height, width), subimg.dtype)
        pad[pad_y0:pad_y0 + subimg_h, pad_x0:pad_x0 + subimg_w] = subimg

        return pad, subimg_cx, subimg_cy


# =============================================================================
# Fitting functions
def MJysr_to_uJy(flux_MJysr: u.Quantity, oversample: int, pix_map: np.ndarray,
                 xpos_idx: np.ndarray | list[int], ypos_idx: np.ndarray | list[int]) -> u.Quantity:
    """Convert specific intensity to spectral flux density.

    Parameters
    ----------
    flux_MJysr : `astropy.units.Quantity`
        Input source flux (MJy/sr).
    oversample : `int`
        Oversampling factor.
    pix_map : `astropy.units.Quantity`
        Solid angle pixel map (arcsec^2).
    xpos_idx : `numpy.ndarray` or `list`[`int`]
        Zero-based indexes closest to source X positions on the detector.
    ypos_idx : `numpy.ndarray` or `list`[`int`]
        Zero-based indexes closest to source Y positions on the detector.

    Returns
    -------
    `astropy.units.Quantity`
        Converted source flux (uJy).
    """
    assert (flux_MJysr.unit == "MJy/sr"), "Input flux is not in MJy/sr"
    assert (pix_map.unit == u.arcsec**2), "Input pixel map is not in arcsec^2"
    flux_uJy = flux_MJysr.to(u.uJy / u.arcsec**2) * pix_map[ypos_idx, xpos_idx] / oversample**2
    return flux_uJy


def tractor_fit_spherex(
    spx_cutout: Cutout2D,
    spx_variance: np.ndarray,
    spx_flags: np.ndarray,
    spx_bkgd: Cutout2D | None,
    source_list: Table,
    detector: int,
    cutout_origin: tuple,
    psf: PSFZoneData,
    hips_name: str,
    keep_frame_sources: bool = False,
    return_fit: bool = False,
    fit_shapes: bool = True,
) -> Table | dict[str, Table | HDUList]:
    """
    Computes forced photometry with Tractor of sources in a SPHEREx cutout.

    Parameters
    ----------
    spx_cutout : `astropy.nddata.Cutout2D`
        A small subregion of a SPHEREx image. WCS info should be included.
    spx_variance : `numpy.ndarray`
        Variance for the relevant subregion.
    spx_flags : `numpy.ndarray`
        Flags for the relevant subregion.
    spx_bkgd : `astropy.nddata.Cutout2D` or `None`
        Local background estimate for the relevant subregion.
    source_list : `astropy.table.Table`
        Sources to fit. Table includes: ID, RA, DEC, SersicRadius (arcsec),
        SersicAB, SersicNu, SersicPhi (deg), X, Y. For sources treated as point
        sources, set Sersic Radius to 0 or a negative number. If X, Y are not
        NaN, they are treated as the source coordinates and the WCS is
        bypassed.
    detector : `int`
        Detector number (1-6).
    cutout_origin : `tuple`
        Position (x,y) in pixels of the cutout origin on full SPHEREx detector.
        Used to query the PSF.
    psf : `spx_data.psf.PSFZoneData`
        Contains the PSF, estimated per zone.
    hips_name : `str`
        Name of the column in `source_list` storing HEALPix indices.
    keep_frame_sources : `bool`
        If `True`, returns photometry for sources in the outer frame of the
        cutout. Default is `False`.
    return_fit : `bool`, optional
        If `True`, fit model and residual image are returned. Default is
        `False`.
    fit_shapes : bool, optional
        If `True`, treat galaxies as shapes during modeling. If False, treat
        galaxies as point sources. Default is `True`.

    Returns
    -------
    phot_meas_arr: `numpy.ndarray`
        Fluxes/errors of the sources provided in `source_list`. Includes
        photometry for frame sources if `keep_frame_sources` is ``True``.
    fluxfit: `dict`[`str`, `numpy.ndarray` | `list`]
        If `return_fit` is ``True``, returns a dictionary including the flux
        measurement array and a list of information for plotting residuals.
    """
    spherex_psf = SPHERExTractorPSF(
        psf, detector, xshift=cutout_origin[0], yshift=cutout_origin[1]
    )

    tractor_source_objs = np.empty(shape=len(source_list), dtype=object)
    for i, (_x, _y, _ra, _dec, _rser, _ab, _phi, _nu) in enumerate(zip(
        source_list["x_window"],
        source_list["y_window"],
        source_list["ra"],
        source_list["dec"],
        source_list["shape_r"],  # LS half-light radius in arcsec
        source_list["ab_ratio"],  # derived axis ratio a/b
        source_list["phi"],  # derived position angle
        source_list["sersic"]  # LS sersic index
    )):
        # NOTE: We assign the flux assuming the same units as the input cutout,
        # MJy/sr (since model is initialized with LinearPhotCal).
        _flux = Flux(np.random.uniform(high=2))
        if (not fit_shapes) or (_rser < 1) or (_rser is None):
            # Start off with point source
            _src = PointSource(PixPos(_x, _y), _flux)
        else:
            # Start off with extended source
            ba = 1. / _ab
            cut_angle = np.arccos(spx_cutout.wcs.wcs.pc[0, 0]) * 180. / np.pi
            # Get position angle sign following "East of North" convention
            pc_det = np.linalg.det(spx_cutout.wcs.wcs.pc)
            if pc_det < 0:
                # E is CCW of N (right-handed Cartesian)
                pos_ang = _phi + cut_angle
            elif pc_det > 0:
                # E is CW of N (left-handed Cartesian)
                pos_ang = - (_phi + cut_angle)
            _src = SersicGalaxy(
                PixPos(_x, _y),
                _flux,
                GalaxyShape(_rser, ba, pos_ang),
                SersicIndex(_nu)
            )

        # Freeze position for all sources (pure forced photometry)
        _src.freezeAllRecursive()
        _src.thawParam("brightness")
        tractor_source_objs[i] = _src

    tim = tractor.Image(
        data=spx_cutout.data.value,
        invvar=1 / spx_variance,
        psf=spherex_psf,
        wcs=NullWCS(pixscale=spx_param.get_pix_scale().value),
        photocal=LinearPhotoCal(1.0),
        sky=ConstantSky(0.0),
    )

    tim.freezeAllRecursive()  # Image parameters fixed

    # Build Tractor object and double-check source types
    if not fit_shapes:
        source_types = [type(source) for source in tractor_source_objs]
        assert SersicGalaxy not in source_types
    trac_spherex = tractor.Tractor([tim], tractor_source_objs)

    # Optimize the fit with Tractor, positions frozen
    optres = trac_spherex.optimize_forced_photometry(variance=True)
    optres_inv_var = optres.IV

    # Build the output table with photometry measurements
    if keep_frame_sources is False:  # filter frame sources out
        tractor_source_objs = tractor_source_objs[source_list["central_source"]]
        optres_inv_var = optres_inv_var[source_list["central_source"]]
        source_list = source_list[source_list["central_source"]]

    phot_meas_arr = np.empty(shape=len(tractor_source_objs),
                             dtype=[("SPHERExRefID", "int64"),
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
                                    (hips_name, "int32")])
    for i in range(0, len(source_list)):
        # Calculate source flux/flux error
        src_flux = tractor_source_objs[i].getParams()[0]
        src_flux_err = 1 / np.sqrt(optres_inv_var[i])

        # Get background from local estimate cutout if available, otherwise set
        # to zero (will be populated with global background estimate later in
        # Perform Forced Photometry run method)
        position_yidx, position_xidx = (spx_cutout.wcs).world_to_array_index(
            SkyCoord(source_list["ra"][i], source_list["dec"][i], unit=u.deg))

        if spx_bkgd is not None:
            src_bkgd = spx_bkgd.data[position_yidx, position_xidx].value
            local_bkg_flg = True
        else:
            src_bkgd = 0
            local_bkg_flg = False

        # Get source region pixels
        source_center = PixCoord(x=position_xidx, y=position_yidx)
        reg_mask = CirclePixelRegion(center=source_center, radius=2.5).to_mask(mode="center")

        # To account for edge/corner sources: get the bounding box slices
        # overlapping the window (with respect to the cutout/bbox origin),
        # then extract masked pixels
        slices_cutout, slices_mask = reg_mask.get_overlap_slices(shape=spx_cutout.data.shape)
        submask = reg_mask.data[slices_mask].astype(bool)

        # Calculate fit quality
        fit_ql = np.nanmean(
            np.abs((spx_cutout.data.value - trac_spherex.getModelImage(0))[slices_cutout][submask])
            / np.sqrt(spx_variance[slices_cutout][submask]))

        # Set photometry flag
        bad_pixel_flags = (Flag.TRANSIENT | Flag.OVERFLOW | Flag.SUR_ERROR | Flag.NONFUNC | Flag.DICHROIC
                           | Flag.MISSING_DATA | Flag.HOT | Flag.COLD | Flag.NONLINEAR | Flag.CROSSTALK
                           | Flag.GHOST | Flag.GHOST_FPA | Flag.GHOST_EXT | Flag.STREAK | Flag.PERSIST
                           | Flag.BLOOM | Flag.SNOWBALL | Flag.HALO | Flag.SATELLITE_HALO)
        src_reg_flags = spx_flags[slices_cutout][submask]
        phot_flag = np.bitwise_or.reduce(src_reg_flags)
        bad_phot = np.any((src_reg_flags & bad_pixel_flags) != 0)
        if bad_phot:
            phot_flag |= PhotometryFlag.CONTAINS_BAD_PIXEL

        # If fit failed, set fit quality to NaN and add fit error flag
        if (np.isnan(src_flux) or np.isnan(src_flux_err)):
            fit_ql = np.nan
            phot_flag |= PhotometryFlag.FIT_ERROR

        # Populate array element with column values
        # NOTE: To save memory, we populate flux_bkg (in global case), mjd,
        # lambda, lambda_width, deep_flg, lvf_id, and det_id with their true
        # values in the Perform Forced Photometry run method.
        phot_meas_arr[i] = (
            source_list["SPHERExRefID"][i],                                      # SPHERExRefID
            source_list["ra"][i],                                                # ra
            source_list["dec"][i],                                               # dec
            source_list["x_image"][i],                                           # x_image
            source_list["y_image"][i],                                           # y_image
            0,                                                                   # mjd
            0,                                                                   # lambda
            0,                                                                   # lambda_width
            src_flux,                                                            # flux
            src_flux_err,                                                        # flux_err
            src_bkgd,                                                            # flux_bkg
            local_bkg_flg,                                                       # local_bkg_flg
            True if type(tractor_source_objs[i]) is SersicGalaxy else False,     # ext_flg
            phot_flag,                                                           # flags
            fit_ql,                                                              # fit_ql
            0,                                                                   # deep_flg
            0,                                                                   # lvf_id
            0,                                                                   # det_id
            source_list[hips_name][i]                                            # <hips_name>
        )

    if return_fit:
        # Return measurement table AND relevant info to examine residuals
        pos_str = f"{cutout_origin[0]},{cutout_origin[1]}"
        hdu_model = PrimaryHDU(trac_spherex.getModelImage(0)*spx_cutout.data.unit,
                               Header(cards={"BUNIT": spx_cutout.data.unit.to_string()}))
        hdu_cutout = ImageHDU(spx_cutout.data, name="CUTOUT",
                              header=Header(cards={"BUNIT": spx_cutout.data.unit.to_string()}))
        hdu_src = BinTableHDU(data=source_list, name="WINDOW CAT")
        hdu_var = ImageHDU(spx_variance, name="VARIANCE",
                           header=Header(cards={"BUNIT": (spx_cutout.data.unit**2).to_string()}))
        hdu_chi = ImageHDU(trac_spherex.getChiImage(0), name="CHI")
        hdul = HDUList([hdu_model, hdu_cutout, hdu_src, hdu_var, hdu_chi])
        fluxfit = {"flux_result": phot_meas_arr, "fit": [pos_str, hdul]}
        return fluxfit
    else:
        return phot_meas_arr
