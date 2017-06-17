# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools for computing velocities and velocity corrections using Astropy
coordinates and related machinery.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .. import units as u

from .solar_system import get_body_barycentric_posvel
from .matrix_utilities import matrix_product
from .representation import (CartesianRepresentation, CartesianDifferential,
                             UnitSphericalRepresentation)
from .builtin_frames import GCRS


__all__ = ['helio_corr', 'helio_vector']

KPS = u.km/u.s

def helio_vector(t, loc):
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
    vsun = get_body_barycentric_posvel('sun', t)[1]
    vearth = get_body_barycentric_posvel('earth', t)[1]

    vsunearth = vearth - vsun

    gcrs_p, gcrs_v = loc.get_gcrs_posvel(t)

    return (vsunearth + gcrs_v).represent_as(CartesianDifferential)


def bary_vector(t, loc):
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
    vearth = get_body_barycentric_posvel('earth', t)[1]

    gcrs_p, gcrs_v = loc.get_gcrs_posvel(t)

    return (vearth + gcrs_v).represent_as(CartesianDifferential)


def radial_velocity_correction(t, loc, target, kind='barycentric'):
    """
    Compute the correction required to convert a radial velocity at a given
    time and place to a barycentric or heliocentric velocity.

    Parameters
    ----------
    t : astropy.time.Time
        Time of the observation. Can be a Time array.
    loc : astropy.coordinates.EarthLocation
        The observer location at which to compute the correction.
    target : astropy.coordinates.SkyCoord
        The on-sky location at which to compute the correction.
    kind : str
        The kind of correction. Only 'barycentric' and 'heliocentric' are
        allowed values.

    Returns
    -------
    vcorr : astropy.units.Quantity with velocity units
        The  correction with a positive sign.  I.e., *add* this
        to an observed radial velocity to get the heliocentric velocity.

    Notes
    -----
    The algorithm here is sufficient to perform corrections at the ~1 to
    10 m/s level, but has not been validated at higher prevision.  Future
    versions of Astropy will likely aim to improve this.

    Additionally, this function may be deprecated in the future in favor of a
    more complete representation of velocities built into the coordinate frame
    classes. There will be ample warning if this occurs, however (following the
    standard Astropy deprecation rules).

    """
    if kind == 'barycentric':
        vsuntarg_cart = bary_vector(t, loc)
    elif kind == 'heliocentric':
        vsuntarg_cart = helio_vector(t, loc)
    else:
        raise ValueError('Invalid "kind" in radial_velocity_correction: "{}"'.format(kind))

    gcrs_p, _ = loc.get_gcrs_posvel(t)

    gtarg = target.transform_to(GCRS(obstime=t, obsgeoloc=gcrs_p))
    targxyz = gtarg.represent_as(UnitSphericalRepresentation).to_cartesian().xyz

    res = matrix_product(vsuntarg_cart.xyz, targxyz)

    if hasattr(res, 'unit'):
        return res
    else:  # for unclear reasons, matrix_product here drops the unit for scalars
        return res*KPS
