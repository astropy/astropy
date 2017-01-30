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
from .representation import CartesianRepresentation, UnitSphericalRepresentation
from .builtin_frames import GCRS


__all__ = ['helio_corr', 'helio_vector']

KPS = u.km/u.s

def helio_vector(t, loc):
    """
    Compute the heliocentric velocity correction at a given time and place.

    Paramters
    ---------
    t : astropy.time.Time
        Time of the observation. Can be a Time array.
    loc : astropy.coordinates.EarthLocation
        The observer location at which to compute the correction.
    """
    vsun = get_body_barycentric_posvel('sun', t)[1]
    vearth = get_body_barycentric_posvel('earth', t)[1]

    vsunearth = vearth - vsun

    gcrs_p, gcrs_v = loc.get_gcrs_posvel(t)

    vsuntargxyz = vsunearth.xyz + gcrs_v.xyz
    return CartesianRepresentation(vsuntargxyz.to(KPS))


def helio_corr(t, loc, target):
    """
    Compute the correction required to convert a radial velocity at a given
    time and place to a heliocentric velocity.

    Paramters
    ---------
    t : astropy.time.Time
        Time of the observation. Can be a Time array.
    loc : astropy.coordinates.EarthLocation
        The observer location at which to compute the correction.
    target : astropy.coordinates.SkyCoord
        The on-sky location at which to compute the correction.

    Returns
    -------
    vcorr : astropy.units.Quantity with velocity units
        The heliocentric correction with a positive sign.  I.e., *add* this
        to an observed radial velocity to get the heliocentric velocity.
    """
    vsuntarg_cartrepr = helio_vector(t, loc)
    gcrs_p, _ = loc.get_gcrs_posvel(t)

    gtarg = target.transform_to(GCRS(obstime=t, obsgeoloc=gcrs_p))
    targxyz = gtarg.represent_as(UnitSphericalRepresentation).to_cartesian().xyz

    res = matrix_product(vsuntarg_cartrepr.xyz, targxyz)

    if hasattr(res, 'unit'):
        return res
    else:  # for unclear reasons, matrix_product here drops the unit for scalars
        return res*KPS
