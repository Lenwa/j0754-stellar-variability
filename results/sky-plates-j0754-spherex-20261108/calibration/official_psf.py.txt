__all__ = ["PSFZoneData", "spx_psf_reader", "spx_psf_writer"]

# Standard library imports
import copy
import logging
from typing import TypeVar, Union, Any
from pathlib import Path

# Related third party imports
import astropy.units as u
from astropy.io import fits, registry
from astropy import wcs
from astropy.nddata import NDIOMixin
from numpy.typing import ArrayLike
import numba
import numpy as np

from .spx_param import get_pix_scale

# Local application/library specific imports

TPSFZoneData = TypeVar("TPSFZoneData", bound="PSFZoneData")


@numba.njit
def manual_dot(psfs: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Manually computes the tensor dot product, summing the product of
    the weight and the PSF over all corners of an inner PSF rectangle.

    Parameters
    ----------
    psfs : `numpy.ndarray`
        Array of corner PSFs with shape (4, N, N), in order lower left,
        upper left, lower right, upper right.
    weights : `numpy.ndarray`
        Array of corner weights with shape (4, ), in order lower left,
        upper left, lower right, upper right.

    Returns
    -------
    blended_psf : `numpy.ndarray`
        Blended/interpolated PSF with shape (N, N).
    """
    j, k = psfs.shape[1], psfs.shape[2]
    blended_psf = np.zeros(psfs[0].shape, psfs[0].dtype)
    for i in range(weights.shape[0]):
        for jj in range(j):
            for kk in range(k):
                blended_psf[jj, kk] += weights[i] * psfs[i, jj, kk]
    return blended_psf


class PSFZoneData(NDIOMixin):
    """PSF normalized to unity integral and centered to create PSF grid.

    This class is used to represent the exposure-averaged effective PSF
    (optical+jitter) calculated from the multiple exposures, and the
    per-exposure effective PSF (optical+jitter).

    PSF is calculated per zone, defined by the x center, y center, x width, and
    y width, where coordinates are in FITS.

    Parameters
    ----------
    data : `~numpy.typing.ArrayLike`
        The PSF grid.
    x_centers : `~numpy.ndarray`
        Array of x coordinates of zone centers in the FITS convention.
    y_centers : `~numpy.ndarray`
        Array of y coordinates of zone centers in the FITS convention.
    oversampling : `float`, optional
        PSF oversampling value, by default 10.

    Notes
    -----
    Using FITS coords in the original image.
    We maybe have to apply the offsets when the coordinates are used elsewhere.
    """

    def __init__(self, data: ArrayLike, x_centers: np.ndarray, y_centers: np.ndarray,
                 oversampling: float = 10) -> None:
        self.data: ArrayLike = data
        # Explicitly initialize `<axis>_centers` as arrays in order to use
        # `repeat` method and convert any input lists
        self.x_centers = np.asarray(x_centers)
        self.y_centers = np.asarray(y_centers)
        self.oversampling = oversampling
        # Calculate zone IDs using centers
        self.zone_ids = range(1, (len(x_centers) * len(y_centers)))
        self.x_zones = len(x_centers)
        self.y_zones = len(y_centers)
        # Determine the widths using the centers
        self.x_widths = np.diff(np.asarray([-x_centers[0]] + list(x_centers)))
        self.y_widths = np.diff(np.asarray([-y_centers[0]] + list(y_centers)))
        # Store precomputed zone center coords & quads to get interpolated PSF
        self.header = self._get_header()
        self.zone_xy = np.column_stack((self.x_centers.repeat(self.y_zones),
                                        np.tile(self.y_centers, self.x_zones)))
        self.quads = self._get_psf_quads()
        self.zone_idx_map = {
            (x, y): idx for idx, (x, y) in enumerate(self.zone_xy)
        }
        # Checks for the PSF cube.
        if len(np.asarray(self.data).shape) != 3:
            raise ValueError(f"PSF data is not a cube, it has {len(np.asarray(self.data).shape)}"
                             " dimension(s).")
        if len(self.x_centers)*len(self.y_centers) != np.asarray(self.data).shape[0]:
            raise ValueError(f"Third dimension of PSF cube, {np.asarray(self.data).shape[0]},"
                             " does not equal x and y centers multiplied:"
                             f" {len(self.x_centers)}, {len(self.y_centers)}.")

    def _get_header(self) -> fits.Header:
        """Builds the header and WCS for PSF grid.

        Returns
        -------
        `fits.Header`
            Header for PSF grid.
        """
        # Build WCS
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [(self.data[0].shape[0]//2)+1, (self.data[0].shape[1]//2)+1]
        w.wcs.cdelt = np.repeat(get_pix_scale().to(u.arcsec).value / self.oversampling, 2)
        w.wcs.ctype = ["XANG", "YANG"]
        w.wcs.cname = ["Relative offset in detector X from center of PSF, in arcseconds",
                       "Relative offset in detector Y from center of PSF, in arcseconds"]
        w.wcs.cunit = ["arcsec", "arcsec"]
        header = w.to_header()

        # Remove headers not needed
        del header["LON*"]
        del header["LAT*"]

        header["OVERSAMP"] = self.oversampling
        # Build header using zone IDs
        zone_id = 1
        for ii in range(self.x_zones):
            for jj in range(self.y_zones):
                location = f"({ii}, {jj})"
                header[f"XCTR_{zone_id}"] = self.x_centers[ii], f"Center of x zone {location}"
                header[f"YCTR_{zone_id}"] = self.y_centers[jj], f"Center of y zone {location}"
                header[f"XWID_{zone_id}"] = self.x_widths[ii], f"Width of x zone {location}"
                header[f"YWID_{zone_id}"] = self.y_widths[jj], f"Width of y zone {location}"
                zone_id += 1

        # PSF is normalized to unity integral
        header["BUNIT"] = "", "Unit of pixel value"
        return header

    def _get_psf_quads(self):
        """Get the coordinate boundaries (along zone centers or edges)
        for each quadrilateral enclosing where 4 PSF zones meet.

        Notes
        -----
        1) Definition of a 'quad'
           Example grid with zone IDs
            [3][6][9]
            [2][5][8]
            [1][4][7]

           Example quad
            [2][5]
            [1][4]

        2) Distinguishing between centers and edges
           A quad dictionary has 2 keys. "zone_ctrs" ("zone_edges") outline the
           boundaries of the inner (outer) rectangle within the PSF quad:

           edge ___________ edge
               |     |     |
               | ctr | ctr |
               |_____|_____|
               |     |     |
               | ctr | ctr |
               |_____|_____|
           edge             edge

        3) Using coordinate boundaries
           From the 4-tuple value associated with each quad dictionary key,
           we can determine the corners of a rectangle:

           (xmin, ymax)     (xmax, ymax)
                UL ------------ UR
                 |              |
                 |              |
                LL ------------ LR
           (xmin, ymin)     (xmax, ymin)

        Returns
        -------
        quads : `numpy.ndarray` of `dict`
            An array where each element is a dictionary of quadrilaterals. Keys
            "zone_ctrs" and "zone_edges" map to different coordinate
            boundary values. The coordinates are in order xmin,xmax,ymin,ymax.
        """
        def find_value_coords(grid_xy, position, axis):
            # Find all coordinate sets matching an input position along 1 axis
            if axis == "x":
                match_idx = np.nonzero(grid_xy[:, 0] == position)
            if axis == "y":
                match_idx = np.nonzero(grid_xy[:, 1] == position)

            return grid_xy[match_idx]

        quads = np.empty(shape=((self.x_zones - 1) * (self.y_zones - 1)), dtype=dict)
        quad_idx = 0
        for (ctr_xmin, ctr_ymin) in zip(self.zone_xy[:, 0], self.zone_xy[:, 1]):
            i_xmin = np.nonzero(self.x_centers == ctr_xmin)[0][0]  # find first occurence where true
            i_ymin = np.nonzero(self.y_centers == ctr_ymin)[0][0]

            # If min coordinates belong to last index along an axis
            # (right/top edge of grid), there's no quad
            if (i_xmin == self.x_zones - 1) or (i_ymin == self.y_zones - 1):
                continue

            # Search for a complete quad
            d_ix = 1
            d_iy = 1
            i_xtest = i_xmin + d_ix
            i_ytest = i_ymin + d_iy
            found_quad = False
            while (i_xtest < self.x_zones) and (i_ytest < self.y_zones):
                coords_xtest = find_value_coords(self.zone_xy, self.x_centers[i_xtest], 'x')
                coords_ytest = find_value_coords(self.zone_xy, self.y_centers[i_ytest], 'y')
                # Intersection of the two lists of coordinates
                int_idx = (coords_xtest[:, None] == coords_ytest).all(-1).any(1)  # broadcast to compare rows
                coords_xytest = coords_xtest[int_idx]
                if len(coords_xtest) > 0 and len(coords_ytest) > 0 and len(coords_xytest) > 0:
                    ctr_xmax, ctr_ymax = coords_xytest[0]
                    # Verified presence of LL and UR corners in list of grid
                    # coordinates; check for UL and LR corners
                    if (ctr_xmax, ctr_ymin) in self.zone_xy and (ctr_xmin, ctr_ymax) in self.zone_xy:
                        found_quad = True
                        break

                ''' If none of the above conditions were met, search for
                satisfactory grid points at incrementally larger distances.
                The following sequence stretches the block in each dimension
                while keeping the other dimension as thin as possible, for
                example, (1,1), (2,1), (1,2), (3,1), (1,3), etc. '''
                if d_ix == d_iy:
                    d_ix += 1
                elif d_ix > d_iy:
                    delta = d_ix - d_iy
                    d_ix = copy.deepcopy(d_iy)
                    d_iy += delta
                elif d_iy > d_ix:
                    delta = d_iy - d_ix
                    d_iy = copy.deepcopy(d_ix)
                    d_ix += delta + 1
                i_xtest = i_xmin + d_ix
                i_ytest = i_ymin + d_iy

            if found_quad:
                # Calculate boundaries along edge; use -0.5/2039.5 at margins
                # due to source positions assuming zero-based indices
                edge_xmin = -0.5 if (i_xmin == 0) else ctr_xmin
                edge_ymin = -0.5 if (i_ymin == 0) else ctr_ymin
                edge_xmax = 2039.5 if (i_xtest == self.x_zones - 1) else ctr_xmax
                edge_ymax = 2039.5 if (i_ytest == self.y_zones - 1) else ctr_ymax

                quad = {}
                quad["zone_ctrs"] = (ctr_xmin, ctr_xmax, ctr_ymin, ctr_ymax)
                quad["zone_edges"] = (edge_xmin, edge_xmax, edge_ymin, edge_ymax)
                quads[quad_idx] = quad
                quad_idx += 1

        return quads

    def _get_inner_rectangle(self, x: float, y: float) -> tuple[float]:
        """Get the coordinate boundaries for the inner PSF rectangle
        enclosing the point (x, y).

        Parameters
        ----------
        x, y : `float`
            Source position (in detector coordinates).

        Returns
        -------
        `tuple` of `float`
            Four-tuple of coordinates for the inner PSF rectangle,
            in order (xmin, xmax, ymin, ymax).

        Raises
        ------
        ValueError
            Raised if input positions are outside of all PSF zone edges.
        """
        inner_rect = None
        for quad in self.quads:
            # Use zone edges to avoid excluding query points near image margins
            x1, x2, y1, y2 = quad["zone_edges"]

            if (x1 <= x <= x2) and (y1 <= y <= y2):
                inner_rect = quad["zone_ctrs"]
                return inner_rect
        if inner_rect is None:
            raise ValueError(f"Input positions do not fall within the boundaries of any PSF quads: X={x},"
                             f" Y={y}.")

    def _get_corner_psf(self, rectangle: tuple[float],
                        corner_type: str) -> np.ndarray:
        """Helper to extract one of the four PSFs in an inner PSF rectangle.

        Parameters
        ----------
        rectangle : `tuple` of `float`
            Four-tuple of coordinates for the inner PSF rectangle,
            in order (xmin, xmax, ymin, ymax).
        corner_type : `str`
            The corner for which to extract the PSF estimate, specified as one
            of: {'lower left', 'upper left', 'lower right', 'upper right'}.

        Returns
        -------
        `numpy.ndarray`
            The PSF for the queried corner.
        """
        x1, x2, y1, y2 = rectangle

        if corner_type == "lower left":
            zone_coords = (x1, y1)
        elif corner_type == "upper left":
            zone_coords = (x1, y2)
        elif corner_type == "lower right":
            zone_coords = (x2, y1)
        elif corner_type == "upper right":
            zone_coords = (x2, y2)

        if zone_coords in self.zone_idx_map:
            idx = self.zone_idx_map[zone_coords]
            return self.data[idx]
        else:
            raise ValueError(f"Coordinates {zone_coords} not found in `PSFZoneData.zone_idx_map`. Valid"
                             f"coordinates must correspond to a PSF zone center.")

    def _get_dist_weights(self, x_in: float, y_in: float, rectangle: tuple[float]) -> np.ndarray[float]:
        """Calculate the distance-based weights for PSF interpolation.

        Notes
        -----
        These weights are identical to those of a bilinear interpolation.

        Since we may not have PSF estimates on the boundaries of the detector,
        it is necessary to extrapolate to the margins using PSFs from the
        interior. To avoid numerical pathologies in the extrapolation, use the
        nearest part of the PSF rectangle to determine the PSF for a point on
        the boundary.

        Parameters
        ----------
        x_in, y_in : `float`
            Source position (in detector coordinates).
        rectangle : `tuple` of `float`
            Four-tuple of coordinates for the inner PSF rectangle,
            in order (xmin, xmax, ymin, ymax).

        Returns
        -------
        weights : `np.ndarray` of `float`
            The distance-based weight for each corner, in order (lower left,
            upper left, lower right, upper right).
        """
        x1, x2, y1, y2 = rectangle

        x_eff = max(min(x_in, x2), x1)
        y_eff = max(min(y_in, y2), y1)

        rectangle_area = (x2 - x1) * (y2 - y1)

        weights = np.array([
            (x2 - x_eff) * (y2 - y_eff) / rectangle_area,  # lower left
            (x2 - x_eff) * (y_eff - y1) / rectangle_area,  # upper left
            (x_eff - x1) * (y2 - y_eff) / rectangle_area,  # lower right
            (x_eff - x1) * (y_eff - y1) / rectangle_area   # upper right
        ])

        return weights

    def get_blended_psf(self, x: float, y: float) -> np.ndarray:
        """Retrieves the interpolated PSF for a given source position on the
        detector.

        Notes
        -----
        This method and its helpers are adapted from Ari Cukierman's `get_psf`
        method used in the SPHEREx Sky Simulator. It smoothly varies across the
        detector.

        Parameters
        ----------
        x, y : `float`
            Source position (in detector coordinates).

        Returns
        -------
        psf : `numpy.ndarray`
            The blended PSF.
        """
        psf_rectangle = self._get_inner_rectangle(x, y)
        corner_types = ["lower left", "upper left", "lower right", "upper right"]
        corner_psfs = np.array([self._get_corner_psf(psf_rectangle, c) for c in corner_types])
        corner_weights = self._get_dist_weights(x, y, psf_rectangle)
        psf = manual_dot(corner_psfs, corner_weights)
        weight_sum = np.sum(corner_weights)
        if not (abs(weight_sum - 1.0) <= 1e-3):
            logging.warning(f"Weights sum to {weight_sum:.6f}, expected ~1.0.")
        return psf

    def to_hdu(self, name: str = "PSF-DATA-CUBE") -> fits.ImageHDU:
        """Returns the PSF zone data and header as an image HDU.

        Parameters
        ----------
        name : `str`, optional
            Name of the image HDU, by default "PSF-DATA-CUBE".

        Returns
        -------
        `~astropy.io.fits.ImageHDU`
            The image HDU of the PSF zone data and header.
        """
        return fits.ImageHDU(data=self.data.astype(np.float32, copy=False), header=self.header, name=name)

    @classmethod
    def from_hdu(cls, hdu: fits.ImageHDU) -> TPSFZoneData:
        """Builds a class instance from an image HDU.

        Parameters
        ----------
        hdu : `~astropy.io.fits.ImageHDU`
            HDU to build class instance from.

        Returns
        -------
        `PSFZoneData`
            The class instance built from the image HDU.

        Raises
        ------
        ValueError
            Raised if HDU type is not an image HDU.
        """
        if isinstance(hdu, fits.ImageHDU):
            data = np.array(hdu.data)
            header = hdu.header
            x_centers = np.array(sorted(set(dict(**header["XCTR_*"]).values())))
            y_centers = np.array(sorted(sorted(set(dict(**header["YCTR_*"]).values()))))
            oversampling = header["OVERSAMP"]
            psf = PSFZoneData(data, x_centers, y_centers, oversampling)
            return psf
        else:
            raise ValueError("Unexpected HDU type for PSF layer.")

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, PSFZoneData) and \
            np.allclose(self.x_centers, other.x_centers) and \
            np.allclose(self.y_centers, other.y_centers) and \
            np.allclose(self.x_widths, other.x_widths) and \
            np.allclose(self.y_widths, other.y_widths) and \
            self.zone_ids == other.zone_ids and \
            np.allclose(self.data, other.data)

    def __repr__(self) -> str:
        return f"SPHEREx PSF Data Object - {len(self.data)} kernels"


def spx_psf_writer(psf: PSFZoneData, filename: Union[str, Path], **kwargs) -> None:
    """Writes PSF zone data object to a file.

    Parameters
    ----------
    psf : `PSFZoneData`
        The instance to write to a file.
    filename : `Union[str, Path]`
        Filename and path to write to.
    """
    primary_header = fits.Header()
    primary_header["COMMENT"] = "PSF zone data defined by x and y centers and x and y zones"
    primary_hdu = fits.PrimaryHDU(data=None, header=primary_header)
    psf_hdu = psf.to_hdu()
    hdul = fits.HDUList([primary_hdu, psf_hdu])
    hdul.writeto(filename, **kwargs)


def spx_psf_reader(filename: Union[str, Path], **kwargs) -> PSFZoneData:
    """Reads in FITS file and creates PSF instance.

    Parameters
    ----------
    filename : `Union[str, Path]`
        Path to FITS file to read in.

    Returns
    -------
    `PSFZoneData`
        Instance of PSF built from FITS file.

    Notes
    -----
    FITS files that contain scaled data (e.g. unsigned integer images) will
    be scaled and the keywords used to manage scaled data in
    :mod:`astropy.io.fits` are disabled.
    """
    with fits.open(filename, **kwargs) as hdul:
        return PSFZoneData.from_hdu(hdul[1])


# Register read/write methods for PSFZoneData
with registry.delay_doc_updates(PSFZoneData):
    registry.register_reader(data_format="fits", data_class=PSFZoneData,
                             function=spx_psf_reader)
    registry.register_writer(data_format="fits", data_class=PSFZoneData,
                             function=spx_psf_writer)
    registry.register_identifier(data_format="fits", data_class=PSFZoneData,
                                 identifier=fits.connect.is_fits)
