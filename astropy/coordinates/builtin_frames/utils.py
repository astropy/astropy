# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains functions/values used repeatedly in different modules of
the ``builtin_frames`` package.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import warnings

import numpy as np

from ... import units as u
from ... import _erfa as erfa
from ...time import Time
from ...utils import iers
from ...utils.exceptions import AstropyWarning

from ...extern.six.moves import range

# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
EQUINOX_J2000 = Time('J2000', scale='utc')
EQUINOX_B1950 = Time('B1950', scale='tai')

# This is a time object that is the default "obstime" when such an attribute is
# necessary.  Currently, we use J2000.
DEFAULT_OBSTIME = Time('J2000', scale='utc')

PIOVER2 = np.pi / 2.

# comes from the mean of the 1962-2014 IERS B data
_DEFAULT_PM = (0.035, 0.29)*u.arcsec


def get_polar_motion(time):
    """
    gets the two polar motion components in radians for use with apio13
    """
    # Get the polar motion from the IERS table
    xp, yp, status = iers.IERS_Auto.open().pm_xy(time, return_status=True)

    wmsg = None
    if np.any(status == iers.TIME_BEFORE_IERS_RANGE):
        wmsg = ('Tried to get polar motions for times before IERS data is '
                'valid. Defaulting to polar motion from the 50-yr mean for those. '
                'This may affect precision at the 10s of arcsec level')
        xp.ravel()[status.ravel() == iers.TIME_BEFORE_IERS_RANGE] = _DEFAULT_PM[0]
        yp.ravel()[status.ravel() == iers.TIME_BEFORE_IERS_RANGE] = _DEFAULT_PM[1]

        warnings.warn(wmsg, AstropyWarning)

    if np.any(status == iers.TIME_BEYOND_IERS_RANGE):
        wmsg = ('Tried to get polar motions for times after IERS data is '
                'valid. Defaulting to polar motion from the 50-yr mean for those. '
                'This may affect precision at the 10s of arcsec level')

        xp.ravel()[status.ravel() == iers.TIME_BEYOND_IERS_RANGE] = _DEFAULT_PM[0]
        yp.ravel()[status.ravel() == iers.TIME_BEYOND_IERS_RANGE] = _DEFAULT_PM[1]

        warnings.warn(wmsg, AstropyWarning)

    return xp.to_value(u.radian), yp.to_value(u.radian)


def _warn_iers(ierserr):
    """
    Generate a warning for an IERSRangeerror

    Parameters
    ----------
    ierserr : An `~astropy.utils.iers.IERSRangeError`
    """
    msg = '{0} Assuming UT1-UTC=0 for coordinate transformations.'
    warnings.warn(msg.format(ierserr.args[0]), AstropyWarning)


def get_dut1utc(time):
    """
    This function is used to get UT1-UTC in coordinates because normally it
    gives an error outside the IERS range, but in coordinates we want to allow
    it to go through but with a warning.
    """
    try:
        return time.delta_ut1_utc
    except iers.IERSRangeError as e:
        _warn_iers(e)
        return np.zeros(time.shape)


def get_jd12(time, scale):
    """
    Gets ``jd1`` and ``jd2`` from a time object in a particular scale.

    Parameters
    ----------
    time : `~astropy.time.Time`
        The time to get the jds for
    scale : str
        The time scale to get the jds for

    Returns
    -------
    jd1 : float
    jd2 : float
    """
    if time.scale == scale:
        newtime = time
    else:
        try:
            newtime = getattr(time, scale)
        except iers.IERSRangeError as e:
            _warn_iers(e)
            newtime = time

    return newtime.jd1, newtime.jd2


def norm(p):
    """
    Normalise a p-vector.
    """
    if np.__version__ == '1.14.0':
        # there is a bug in numpy v1.14.0 (fixed in 1.14.1) that causes
        # this einsum call to break with the default of optimize=True
        # see https://github.com/astropy/astropy/issues/7051
        return p / np.sqrt(np.einsum('...i,...i', p, p, optimize=False))[..., np.newaxis]
    else:
        return p / np.sqrt(np.einsum('...i,...i', p, p))[..., np.newaxis]


def get_cip(jd1, jd2):
    """
    Find the X, Y coordinates of the CIP and the CIO locator, s.

    Parameters
    ----------
    jd1 : float or `np.ndarray`
        First part of two part Julian date (TDB)
    jd2 : float or `np.ndarray`
        Second part of two part Julian date (TDB)

    Returns
    --------
    x : float or `np.ndarray`
        x coordinate of the CIP
    y : float or `np.ndarray`
        y coordinate of the CIP
    s : float or `np.ndarray`
        CIO locator, s
    """
    # classical NPB matrix, IAU 2006/2000A
    rpnb = erfa.pnm06a(jd1, jd2)
    # CIP X, Y coordinates from array
    x, y = erfa.bpn2xy(rpnb)
    # CIO locator, s
    s = erfa.s06(jd1, jd2, x, y)
    return x, y, s


def aticq(ri, di, astrom):
    """
    A slightly modified version of the ERFA function ``eraAticq``.

    ``eraAticq`` performs the transformations between two coordinate systems,
    with the details of the transformation being encoded into the ``astrom`` array.

    The companion function ``eraAtciqz`` is meant to be its inverse. However, this
    is not true for directions close to the Solar centre, since the light deflection
    calculations are numerically unstable and therefore not reversible.

    This version sidesteps that problem by artificially reducing the light deflection
    for directions which are within 90 arcseconds of the Sun's position. This is the
    same approach used by the ERFA functions above, except that they use a threshold of
    9 arcseconds.

    Parameters
    ----------
    ri : float or `~numpy.ndarray`
        right ascension, radians
    di : float or `~numpy.ndarray`
        declination, radians
    astrom : eraASTROM array
        ERFA astrometry context, as produced by, e.g. ``eraApci13`` or ``eraApcs13``

    Returns
    --------
    rc : float or `~numpy.ndarray`
    dc : float or `~numpy.ndarray`
    """
    # RA, Dec to cartesian unit vectors
    pos = erfa.s2c(ri, di)

    # Bias-precession-nutation, giving GCRS proper direction.
    ppr = erfa.trxp(astrom['bpn'], pos)

    # Aberration, giving GCRS natural direction
    d = np.zeros_like(ppr)
    for j in range(2):
        before = norm(ppr-d)
        after = erfa.ab(before, astrom['v'], astrom['em'], astrom['bm1'])
        d = after - before
    pnat = norm(ppr-d)

    # Light deflection by the Sun, giving BCRS coordinate direction
    d = np.zeros_like(pnat)
    for j in range(5):
        before = norm(pnat-d)
        after = erfa.ld(1.0, before, before, astrom['eh'], astrom['em'], 5e-8)
        d = after - before
    pco = norm(pnat-d)

    # ICRS astrometric RA, Dec
    rc, dc = erfa.c2s(pco)
    return erfa.anp(rc), dc


def atciqz(rc, dc, astrom):
    """
    A slightly modified version of the ERFA function ``eraAtciqz``.

    ``eraAtciqz`` performs the transformations between two coordinate systems,
    with the details of the transformation being encoded into the ``astrom`` array.

    The companion function ``eraAticq`` is meant to be its inverse. However, this
    is not true for directions close to the Solar centre, since the light deflection
    calculations are numerically unstable and therefore not reversible.

    This version sidesteps that problem by artificially reducing the light deflection
    for directions which are within 90 arcseconds of the Sun's position. This is the
    same approach used by the ERFA functions above, except that they use a threshold of
    9 arcseconds.

    Parameters
    ----------
    rc : float or `~numpy.ndarray`
        right ascension, radians
    dc : float or `~numpy.ndarray`
        declination, radians
    astrom : eraASTROM array
        ERFA astrometry context, as produced by, e.g. ``eraApci13`` or ``eraApcs13``

    Returns
    --------
    ri : float or `~numpy.ndarray`
    di : float or `~numpy.ndarray`
    """
    # BCRS coordinate direction (unit vector).
    pco = erfa.s2c(rc, dc)

    # Light deflection by the Sun, giving BCRS natural direction.
    pnat = erfa.ld(1.0, pco, pco, astrom['eh'], astrom['em'], 5e-8)

    # Aberration, giving GCRS proper direction.
    ppr = erfa.ab(pnat, astrom['v'], astrom['em'], astrom['bm1'])

    # Bias-precession-nutation, giving CIRS proper direction.
    # Has no effect if matrix is identity matrix, in which case gives GCRS ppr.
    pi = erfa.rxp(astrom['bpn'], ppr)

    # CIRS (GCRS) RA, Dec
    ri, di = erfa.c2s(pi)
    return erfa.anp(ri), di


def prepare_earth_position_vel(time):
    """
    Get barycentric position and velocity, and heliocentric position of Earth

    Parameters
    -----------
    time : `~astropy.time.Time`
        time at which to calculate position and velocity of Earth

    Returns
    --------
    earth_pv : `np.ndarray`
        Barycentric position and velocity of Earth, in au and au/day
    earth_helio : `np.ndarray`
        Heliocentric position of Earth in au
    """
    # this goes here to avoid circular import errors
    from ..solar_system import (get_body_barycentric, get_body_barycentric_posvel)
    # get barycentric position and velocity of earth
    earth_pv = get_body_barycentric_posvel('earth', time)

    # get heliocentric position of earth, preparing it for passing to erfa.
    sun = get_body_barycentric('sun', time)
    earth_heliocentric = (earth_pv[0] -
                          sun).get_xyz(xyz_axis=-1).to_value(u.au)

    # Also prepare earth_pv for passing to erfa, which wants xyz in last
    # dimension, and pos/vel in one-but-last.
    # (Note could use np.stack once our minimum numpy version is >=1.10.)
    earth_pv = np.concatenate((earth_pv[0].get_xyz(xyz_axis=-1).to(u.au)
                               [..., np.newaxis, :].value,
                               earth_pv[1].get_xyz(xyz_axis=-1).to(u.au/u.d)
                               [..., np.newaxis, :].value), axis=-2)
    return earth_pv, earth_heliocentric
