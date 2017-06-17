# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Functions for computing velocity corrections.

Note that this is *not* meant to be a public module - it is an internal
implementation convenience, and the public API is on SkyCoord.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .solar_system import get_body_barycentric_posvel
from .representation import CartesianDifferential, UnitSphericalRepresentation
from .builtin_frames import GCRS
from .sky_coordinate import SkyCoord


__all__ = ['radial_velocity_correction', 'helio_vector', 'bary_vector']


def helio_vector(t, loc, ephemeris=None):
    """
    Computes the heliocentric velocity correction at a given time and place.

    Paramters
    ---------
    t : astropy.time.Time
        Time of the observation. Can be a Time array.
    loc : astropy.coordinates.EarthLocation
        The observer location at which to compute the correction.

    Returns
    -------
    vh : CartesianDifferential
        The heliocentric velocity vector
    """
    vsun = get_body_barycentric_posvel('sun', t, ephemeris=ephemeris)[1]
    vearth = get_body_barycentric_posvel('earth', t, ephemeris=ephemeris)[1]

    vsunearth = vearth - vsun

    gcrs_p, gcrs_v = loc.get_gcrs_posvel(t)

    return (vsunearth + gcrs_v).represent_as(CartesianDifferential)


def bary_vector(t, loc, ephemeris=None):
    """
    Computes the Solar System barycenter velocity correction at a given time and
    place.

    Paramters
    ---------
    t : astropy.time.Time
        Time of the observation. Can be a Time array.
    loc : astropy.coordinates.EarthLocation
        The observer location at which to compute the correction.

    Returns
    -------
    vh : CartesianDifferential
        The barycentric velocity vector
    """
    vearth = get_body_barycentric_posvel('earth', t, ephemeris=ephemeris)[1]

    gcrs_p, gcrs_v = loc.get_gcrs_posvel(t)

    return (vearth + gcrs_v).represent_as(CartesianDifferential)


def radial_velocity_correction(t, loc, target, bary=True, ephemeris=None):
    """
    This is the private API for the SkyCoord.radial_velocity_correction method.
    See the SkyCoord docstrings for a more detailed set of algorithm notes
    (remember that users generally will *not* see this docstring.)

    Parameters
    ----------
    t : astropy.time.Time
        Time of the observation. Can be a Time array.
    loc : astropy.coordinates.EarthLocation
        The observer location at which to compute the correction.
    target : astropy.coordinates.SkyCoord
        The on-sky location at which to compute the correction.
    bary : bool
        If True, do barycentric correction.  If False, do heliocentric.
    ephemeris, optional
        The ephemeris to use for the caclulation.  Must be something that
        `bary_vector`/`helio_vector` accept.

    Returns
    -------
    vcorr : astropy.units.Quantity with velocity units
        The  correction with a positive sign.  I.e., *add* this
        to an observed radial velocity to get the heliocentric velocity.
    """
    if bary:
        vcorr_cart = bary_vector(t, loc, ephemeris=ephemeris)
    else:
        vcorr_cart = helio_vector(t, loc, ephemeris=ephemeris)

    gcrs_p, _ = loc.get_gcrs_posvel(t)

    gtarg = target.transform_to(GCRS(obstime=t, obsgeoloc=gcrs_p))
    targcart = gtarg.represent_as(UnitSphericalRepresentation).to_cartesian()
    return targcart.dot(vcorr_cart)
