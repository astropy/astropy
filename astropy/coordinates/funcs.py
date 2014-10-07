# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains convinience functions for coordinate-related functionality.

This is generally just wrapping around the object-oriented coordinates
framework, but it is useful for some users who are used to more functional
interfaces.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from .. import units as u

__all__ = ['cartesian_to_spherical', 'spherical_to_cartesian']


def cartesian_to_spherical(x, y, z):
    """
    Converts 3D rectangular cartesian coordinates to spherical polar
    coordinates.

    Note that the resulting angles are latitude/longitude or
    elevation/azimuthal form.  I.e., the origin is along the equator
    rather than at the north pole.

    .. note::
        This function simply wraps functionality provided by the
        `~astropy.coordinates.CartesianRepresentation` and
        `~astropy.coordinates.SphericalRepresentation` classes.  In general,
        for both performance and readability, we suggest using these classes
        directly.  But for situations where a quick one-off conversion makes
        sense, this function is provided.

    Parameters
    ----------
    x : scalar, array-like, or `~astropy.units.Quantity`
        The first cartesian coordinate.
    y : scalar, array-like, or `~astropy.units.Quantity`
        The second cartesian coordinate.
    z : scalar, array-like, or `~astropy.units.Quantity`
        The third cartesian coordinate.

    Returns
    -------
    r : `~astropy.units.Quantity`
        The radial coordinate (in the same units as the inputs).
    lat : `~astropy.units.Quantity`
        The latitude in radians
    lon : `~astropy.units.Quantity`
        The longitude in radians
    """
    from .representation import SphericalRepresentation, CartesianRepresentation

    if not hasattr(x, 'unit'):
        x = x * u.dimensionless_unscaled
    if not hasattr(y, 'unit'):
        y = y * u.dimensionless_unscaled
    if not hasattr(z, 'unit'):
        z = z * u.dimensionless_unscaled

    cart = CartesianRepresentation(x, y, z)
    sph = cart.represent_as(SphericalRepresentation)

    return sph.distance, sph.lat, sph.lon


def spherical_to_cartesian(r, lat, lon):
    """
    Converts spherical polar coordinates to rectangular cartesian
    coordinates.

    Note that the input angles should be in latitude/longitude or
    elevation/azimuthal form.  I.e., the origin is along the equator
    rather than at the north pole.

    .. note::
        This is a low-level function used internally in
        `astropy.coordinates`.  It is provided for users if they really
        want to use it, but it is recommended that you use the
        `astropy.coordinates` coordinate systems.

    Parameters
    ----------
    r : scalar, array-like, or `~astropy.units.Quantity`
        The radial coordinate (in the same units as the inputs).
    lat : scalar, array-like, or `~astropy.units.Quantity`
        The latitude (in radians if array or scalar)
    lon : scalar, array-like, or `~astropy.units.Quantity`
        The longitude (in radians if array or scalar)

    Returns
    -------
    x : float or array
        The first cartesian coordinate.
    y : float or array
        The second cartesian coordinate.
    z : float or array
        The third cartesian coordinate.


    """
    from .representation import SphericalRepresentation, CartesianRepresentation

    if not hasattr(r, 'unit'):
        r = r * u.dimensionless_unscaled
    if not hasattr(lat, 'unit'):
        lat = lat * u.radian
    if not hasattr(lon, 'unit'):
        lon = lon * u.radian

    sph = SphericalRepresentation(distance=r, lat=lat, lon=lon)
    cart = sph.represent_as(CartesianRepresentation)

    return cart.x, cart.y, cart.z
