"""
Module for parsing astronomical object names to extract embedded coordinates
eg: '2MASS J06495091-0737408'.
"""

import re

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

RA_REGEX = r"()([0-2]\d)([0-5]\d)([0-5]\d)\.?(\d{0,3})"
DEC_REGEX = r"([+-])(\d{1,2})([0-5]\d)([0-5]\d)\.?(\d{0,3})"
JCOORD_REGEX = "(.*?J)" + RA_REGEX + DEC_REGEX
JPARSER = re.compile(JCOORD_REGEX)


def _sexagesimal(g):
    # convert matched regex groups to sexigesimal array
    sign, h, m, s, frac = g
    sign = -1 if (sign == "-") else 1
    s = ".".join((s, frac))
    return sign * np.array([h, m, s], float)


def search(name, raise_=False):
    """Regex match for coordinates in name."""
    # extract the coordinate data from name
    match = JPARSER.search(name)
    if match is None and raise_:
        raise ValueError("No coordinate match found!")
    return match


def to_ra_dec_angles(name):
    """get RA in hourangle and DEC in degrees by parsing name."""
    groups = search(name, True).groups()
    prefix, hms, dms = np.split(groups, [1, 6])
    ra = (_sexagesimal(hms) / (1, 60, 60 * 60) * u.hourangle).sum()
    dec = (_sexagesimal(dms) * (u.deg, u.arcmin, u.arcsec)).sum()
    return ra, dec


def to_skycoord(name, frame="icrs"):
    """Convert to `name` to `SkyCoords` object."""
    return SkyCoord(*to_ra_dec_angles(name), frame=frame)


def shorten(name):
    """Produce a shortened version of the full object name.

    The shortened name is built from the prefix (usually the survey name) and RA (hour,
    minute), DEC (deg, arcmin) parts.
    e.g.: '2MASS J06495091-0737408' --> '2MASS J0649-0737'

    Parameters
    ----------
    name : str
        Full object name with J-coords embedded.

    Returns
    -------
    shortName: str
    """
    match = search(name)
    return "".join(match.group(1, 3, 4, 7, 8, 9))
