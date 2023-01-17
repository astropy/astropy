# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains functions/values used repeatedly in different modules of
the ``builtin_frames`` package.
"""

import warnings

import erfa
import numpy as np

from astropy import units as u
from astropy.coordinates.earth import EarthLocation
from astropy.coordinates.representation import CartesianDifferential
from astropy.time import Time
from astropy.utils import iers
from astropy.utils.exceptions import AstropyWarning

# We use tt as the time scale for this equinoxes, primarily because it is the
# convention for J2000 (it is unclear if there is any "right answer" for B1950)
# while #8600 makes this the default behavior, we show it here to ensure it's
# clear which is used here
EQUINOX_J2000 = Time("J2000", scale="tt")
EQUINOX_B1950 = Time("B1950", scale="tt")

# This is a time object that is the default "obstime" when such an attribute is
# necessary.  Currently, we use J2000.
DEFAULT_OBSTIME = Time("J2000", scale="tt")

# This is an EarthLocation that is the default "location" when such an attribute is
# necessary. It is the centre of the Earth.
EARTH_CENTER = EarthLocation(0 * u.km, 0 * u.km, 0 * u.km)

PIOVER2 = np.pi / 2.0

# comes from the mean of the 1962-2014 IERS B data
_DEFAULT_PM = (0.035, 0.29) * u.arcsec


def get_polar_motion(time):
    """
    gets the two polar motion components in radians for use with apio.
    """
    # Get the polar motion from the IERS table
    iers_table = iers.earth_orientation_table.get()
    xp, yp, status = iers_table.pm_xy(time, return_status=True)

    wmsg = (
        "Tried to get polar motions for times {} IERS data is "
        "valid. Defaulting to polar motion from the 50-yr mean for those. "
        "This may affect precision at the arcsec level. Please check your "
        "astropy.utils.iers.conf.iers_auto_url and point it to a newer "
        "version if necessary."
    )
    if np.any(status == iers.TIME_BEFORE_IERS_RANGE):
        xp[status == iers.TIME_BEFORE_IERS_RANGE] = _DEFAULT_PM[0]
        yp[status == iers.TIME_BEFORE_IERS_RANGE] = _DEFAULT_PM[1]

        warnings.warn(wmsg.format("before"), AstropyWarning)

    if np.any(status == iers.TIME_BEYOND_IERS_RANGE):
        xp[status == iers.TIME_BEYOND_IERS_RANGE] = _DEFAULT_PM[0]
        yp[status == iers.TIME_BEYOND_IERS_RANGE] = _DEFAULT_PM[1]

        warnings.warn(wmsg.format("after"), AstropyWarning)

    return xp.to_value(u.radian), yp.to_value(u.radian)


def _warn_iers(ierserr):
    """
    Generate a warning for an IERSRangeerror.

    Parameters
    ----------
    ierserr : An `~astropy.utils.iers.IERSRangeError`
    """
    msg = "{0} Assuming UT1-UTC=0 for coordinate transformations."
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
    return p / np.sqrt(np.einsum("...i,...i", p, p))[..., np.newaxis]


def pav2pv(p, v):
    """
    Combine p- and v- vectors into a pv-vector.
    """
    pv = np.empty(np.broadcast(p, v).shape[:-1], erfa.dt_pv)
    pv["p"] = p
    pv["v"] = v
    return pv


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
    -------
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


def aticq(srepr, astrom):
    """
    A slightly modified version of the ERFA function ``eraAticq``.

    ``eraAticq`` performs the transformations between two coordinate systems,
    with the details of the transformation being encoded into the ``astrom`` array.

    There are two issues with the version of aticq in ERFA. Both are associated
    with the handling of light deflection.

    The companion function ``eraAtciqz`` is meant to be its inverse. However, this
    is not true for directions close to the Solar centre, since the light deflection
    calculations are numerically unstable and therefore not reversible.

    This version sidesteps that problem by artificially reducing the light deflection
    for directions which are within 90 arcseconds of the Sun's position. This is the
    same approach used by the ERFA functions above, except that they use a threshold of
    9 arcseconds.

    In addition, ERFA's aticq assumes a distant source, so there is no difference between
    the object-Sun vector and the observer-Sun vector. This can lead to errors of up to a
    few arcseconds in the worst case (e.g a Venus transit).

    Parameters
    ----------
    srepr : `~astropy.coordinates.SphericalRepresentation`
        Astrometric GCRS or CIRS position of object from observer
    astrom : eraASTROM array
        ERFA astrometry context, as produced by, e.g. ``eraApci13`` or ``eraApcs13``

    Returns
    -------
    rc : float or `~numpy.ndarray`
        Right Ascension in radians
    dc : float or `~numpy.ndarray`
        Declination in radians
    """
    # ignore parallax effects if no distance, or far away
    srepr_distance = srepr.distance
    ignore_distance = srepr_distance.unit == u.one

    # RA, Dec to cartesian unit vectors
    pos = erfa.s2c(srepr.lon.radian, srepr.lat.radian)

    # Bias-precession-nutation, giving GCRS proper direction.
    ppr = erfa.trxp(astrom["bpn"], pos)

    # Aberration, giving GCRS natural direction
    d = np.zeros_like(ppr)
    for j in range(2):
        before = norm(ppr - d)
        after = erfa.ab(before, astrom["v"], astrom["em"], astrom["bm1"])
        d = after - before
    pnat = norm(ppr - d)

    # Light deflection by the Sun, giving BCRS coordinate direction
    d = np.zeros_like(pnat)
    for j in range(5):
        before = norm(pnat - d)
        if ignore_distance:
            # No distance to object, assume a long way away
            q = before
        else:
            # Find BCRS direction of Sun to object.
            # astrom['eh'] and astrom['em'] contain Sun to observer unit vector,
            # and distance, respectively.
            eh = astrom["em"][..., np.newaxis] * astrom["eh"]
            # unit vector from Sun to object
            q = eh + srepr_distance[..., np.newaxis].to_value(u.au) * before
            sundist, q = erfa.pn(q)
            sundist = sundist[..., np.newaxis]
            # calculation above is extremely unstable very close to the sun
            # in these situations, default back to ldsun-style behaviour,
            # since this is reversible and drops to zero within stellar limb
            q = np.where(sundist > 1.0e-10, q, before)

        after = erfa.ld(1.0, before, q, astrom["eh"], astrom["em"], 1e-6)
        d = after - before
    pco = norm(pnat - d)

    # ICRS astrometric RA, Dec
    rc, dc = erfa.c2s(pco)
    return erfa.anp(rc), dc


def atciqz(srepr, astrom):
    """
    A slightly modified version of the ERFA function ``eraAtciqz``.

    ``eraAtciqz`` performs the transformations between two coordinate systems,
    with the details of the transformation being encoded into the ``astrom`` array.

    There are two issues with the version of atciqz in ERFA. Both are associated
    with the handling of light deflection.

    The companion function ``eraAticq`` is meant to be its inverse. However, this
    is not true for directions close to the Solar centre, since the light deflection
    calculations are numerically unstable and therefore not reversible.

    This version sidesteps that problem by artificially reducing the light deflection
    for directions which are within 90 arcseconds of the Sun's position. This is the
    same approach used by the ERFA functions above, except that they use a threshold of
    9 arcseconds.

    In addition, ERFA's atciqz assumes a distant source, so there is no difference between
    the object-Sun vector and the observer-Sun vector. This can lead to errors of up to a
    few arcseconds in the worst case (e.g a Venus transit).

    Parameters
    ----------
    srepr : `~astropy.coordinates.SphericalRepresentation`
        Astrometric ICRS position of object from observer
    astrom : eraASTROM array
        ERFA astrometry context, as produced by, e.g. ``eraApci13`` or ``eraApcs13``

    Returns
    -------
    ri : float or `~numpy.ndarray`
        Right Ascension in radians
    di : float or `~numpy.ndarray`
        Declination in radians
    """
    # ignore parallax effects if no distance, or far away
    srepr_distance = srepr.distance
    ignore_distance = srepr_distance.unit == u.one

    # BCRS coordinate direction (unit vector).
    pco = erfa.s2c(srepr.lon.radian, srepr.lat.radian)

    # Find BCRS direction of Sun to object
    if ignore_distance:
        # No distance to object, assume a long way away
        q = pco
    else:
        # Find BCRS direction of Sun to object.
        # astrom['eh'] and astrom['em'] contain Sun to observer unit vector,
        # and distance, respectively.
        eh = astrom["em"][..., np.newaxis] * astrom["eh"]
        # unit vector from Sun to object
        q = eh + srepr_distance[..., np.newaxis].to_value(u.au) * pco
        sundist, q = erfa.pn(q)
        sundist = sundist[..., np.newaxis]
        # calculation above is extremely unstable very close to the sun
        # in these situations, default back to ldsun-style behaviour,
        # since this is reversible and drops to zero within stellar limb
        q = np.where(sundist > 1.0e-10, q, pco)

    # Light deflection by the Sun, giving BCRS natural direction.
    pnat = erfa.ld(1.0, pco, q, astrom["eh"], astrom["em"], 1e-6)

    # Aberration, giving GCRS proper direction.
    ppr = erfa.ab(pnat, astrom["v"], astrom["em"], astrom["bm1"])

    # Bias-precession-nutation, giving CIRS proper direction.
    # Has no effect if matrix is identity matrix, in which case gives GCRS ppr.
    pi = erfa.rxp(astrom["bpn"], ppr)

    # CIRS (GCRS) RA, Dec
    ri, di = erfa.c2s(pi)
    return erfa.anp(ri), di


def prepare_earth_position_vel(time):
    """
    Get barycentric position and velocity, and heliocentric position of Earth.

    Parameters
    ----------
    time : `~astropy.time.Time`
        time at which to calculate position and velocity of Earth

    Returns
    -------
    earth_pv : `np.ndarray`
        Barycentric position and velocity of Earth, in au and au/day
    earth_helio : `np.ndarray`
        Heliocentric position of Earth in au
    """
    # this goes here to avoid circular import errors
    from astropy.coordinates.solar_system import (
        get_body_barycentric,
        get_body_barycentric_posvel,
        solar_system_ephemeris,
    )

    # get barycentric position and velocity of earth

    ephemeris = solar_system_ephemeris.get()

    # if we are using the builtin erfa based ephemeris,
    # we can use the fact that epv00 already provides all we need.
    # This avoids calling epv00 twice, once
    # in get_body_barycentric_posvel('earth') and once in
    # get_body_barycentric('sun')
    if ephemeris == "builtin":
        jd1, jd2 = get_jd12(time, "tdb")
        earth_pv_heliocentric, earth_pv = erfa.epv00(jd1, jd2)
        earth_heliocentric = earth_pv_heliocentric["p"]

    # all other ephemeris providers probably don't have a shortcut like this
    else:
        earth_p, earth_v = get_body_barycentric_posvel("earth", time)

        # get heliocentric position of earth, preparing it for passing to erfa.
        sun = get_body_barycentric("sun", time)
        earth_heliocentric = (earth_p - sun).get_xyz(xyz_axis=-1).to_value(u.au)

        # Also prepare earth_pv for passing to erfa, which wants it as
        # a structured dtype.
        earth_pv = pav2pv(
            earth_p.get_xyz(xyz_axis=-1).to_value(u.au),
            earth_v.get_xyz(xyz_axis=-1).to_value(u.au / u.d),
        )

    return earth_pv, earth_heliocentric


def get_offset_sun_from_barycenter(time, include_velocity=False, reverse=False):
    """
    Returns the offset of the Sun center from the solar-system barycenter (SSB).

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time at which to calculate the offset
    include_velocity : `bool`
        If ``True``, attach the velocity as a differential.  Defaults to ``False``.
    reverse : `bool`
        If ``True``, return the offset of the barycenter from the Sun.  Defaults to ``False``.

    Returns
    -------
    `~astropy.coordinates.CartesianRepresentation`
        The offset
    """
    if include_velocity:
        # Import here to avoid a circular import
        from astropy.coordinates.solar_system import get_body_barycentric_posvel

        offset_pos, offset_vel = get_body_barycentric_posvel("sun", time)
        if reverse:
            offset_pos, offset_vel = -offset_pos, -offset_vel
        offset_vel = offset_vel.represent_as(CartesianDifferential)
        offset_pos = offset_pos.with_differentials(offset_vel)

    else:
        # Import here to avoid a circular import
        from astropy.coordinates.solar_system import get_body_barycentric

        offset_pos = get_body_barycentric("sun", time)
        if reverse:
            offset_pos = -offset_pos

    return offset_pos
