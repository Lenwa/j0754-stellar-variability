__all__ = ["Flag", "flag2bit", "flag_description", "PhotometryFlag",
           "phot_flag_description"]

from enum import IntFlag
from math import log2
from typing import Union


class Flag(IntFlag):
    """Bit flags"""
    TRANSIENT = 1 << 0
    OVERFLOW = 1 << 1
    SUR_ERROR = 1 << 2
    PHANTOM = 1 << 4
    REFERENCE = 1 << 5
    NONFUNC = 1 << 6
    DICHROIC = 1 << 7
    MISSING_DATA = 1 << 9
    HOT = 1 << 10
    COLD = 1 << 11
    FULLSAMPLE = 1 << 12
    SUR_MISMATCH = 1 << 13
    PHANMISS = 1 << 14
    NONLINEAR = 1 << 15
    VAR_UNDER = 1 << 16
    PERSIST = 1 << 17
    PERSIST_UNK = 1 << 18
    OUTLIER = 1 << 19
    CROSSTALK = 1 << 20
    SOURCE = 1 << 21
    GHOST = 1 << 22
    GHOST_FPA = 1 << 23
    GHOST_EXT = 1 << 24
    STREAK = 1 << 25
    BLOOM = 1 << 26
    SNOWBALL = 1 << 27
    HALO = 1 << 28
    SATELLITE_HALO = 1 << 29


class PhotometryFlag(IntFlag):
    CONTAINS_BAD_PIXEL = 1 << 32
    FIT_ERROR = 1 << 33


def flag2bit(flag: Union[Flag, PhotometryFlag, int]) -> int:
    """Convert flag value to bit position.

    Parameters
    ----------
    flag : `Flag` or `PhotometryFlag` or `int`
        Flag value
    Returns
    -------
    bitpos : `int`
        Bit position of a flag, the least significant bit is 0
    """
    return int(log2(flag))


def flag_description(flag: Flag) -> str:
    """Return flag description.

    Parameters
    ----------
    flag : `Flag`
        A bit flag.
    Returns
    -------
    description: `str`
        Short description of the flag, which would be good for FITS comment.
    """
    desc = None
    match flag:
        case Flag.TRANSIENT:
            desc = "Transient detected during SUR"
        case Flag.OVERFLOW:
            desc = "Overflow detected during SUR"
        case Flag.SUR_ERROR:
            desc = "Error in onboard processing"
        case Flag.PHANTOM:
            desc = "Phantom pixel"
        case Flag.REFERENCE:
            desc = "Reference pixel"
        case Flag.NONFUNC:
            desc = "Permanently unusable"
        case Flag.DICHROIC:
            desc = "Low efficiency due to dichroic"
        case Flag.MISSING_DATA:
            desc = "Onboard data lost"
        case Flag.HOT:
            desc = "Hot pixel"
        case Flag.COLD:
            desc = "Anomalously low signal"
        case Flag.FULLSAMPLE:
            desc = "Pixel full sample history is available"
        case Flag.SUR_MISMATCH:
            desc = "Reconstructed and onboard slopes do not match"
        case Flag.PHANMISS:
            desc = "Phantom correction was not applied"
        case Flag.NONLINEAR:
            desc = "Nonlinearity correction cannot be applied reliably"
        case Flag.VAR_UNDER:
            desc = "Unknown integration time, underestimated variance"
        case Flag.PERSIST:
            desc = "Persistent charge above threshold"
        case Flag.PERSIST_UNK:
            desc = "Persistence model is known to be inaccurate"
        case Flag.OUTLIER:
            desc = "Pixel flagged by Detect Outliers"
        case Flag.CROSSTALK:
            desc = "Pixel affected by crosstalk."
        case Flag.SOURCE:
            desc = "Pixel mapped to a known source"
        case Flag.GHOST:
            desc = "Affected by a bright source inside the exposure"
        case Flag.GHOST_FPA:
            desc = "Affected by a bright source on another detector"
        case Flag.GHOST_EXT:
            desc = "Pixel affected by a bright source outside FOV"
        case Flag.STREAK:
            desc = "Pixel affected by a streak"
        case Flag.BLOOM:
            desc = "Pixel part of a bloom artifact"
        case Flag.SNOWBALL:
            desc = "Pixel part of a snowball artifact"
        case Flag.HALO:
            desc = "Pixel part of a transient halo"
        case Flag.SATELLITE_HALO:
            desc = "Pixel part of a satellite streak halo"

    if desc is None:
        raise ValueError(f"Unknown bit flag {flag}")

    return desc


def phot_flag_description(flag: PhotometryFlag) -> str:
    """Return flag description.

    Parameters
    ----------
    flag : `PhotometryFlag`
        A 64-bit flag.
    Returns
    -------
    description: `str`
        Short description of the flag.
    """
    desc = None
    match flag:
        case PhotometryFlag.CONTAINS_BAD_PIXEL:
            desc = "1 or more of the source pixels have a harmful flag set"
        case PhotometryFlag.FIT_ERROR:
            desc = "Tractor fit failed to return a valid flux or flux error"

    if desc is None:
        raise ValueError(f"Unknown bit flag {flag}")

    return desc
