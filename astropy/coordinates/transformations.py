# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the functionality for transforming points from
one coordinate system to another (e.g. equatorial to galactic). The public
interface for this functionality is in the file core.py.
"""

import math

import numpy as np

pi = math.pi

__all__ = ['cartesian_to_spherical', 'spherical_to_cartesian']


def cartesian_to_spherical(x, y, z):
    """
    Converts 3D rectangular cartesian coordinates to spherical polar
    coordinates.

    Note that the resulting angles are latitude/longitude or
    elevation/azimuthal form.  I.e., the origin is along the equator
    rather than at the north pole.

    .. note::
        This is a low-level function used internally in
        `astropy.coordinates`.  It is provided for users if they really
        want to use it, but it is recommended that you use the
        `astropy.coordinates` coordinate systems.

    Parameters
    ----------
    x : scalar or array-like
        The first cartesian coordinate.
    y : scalar or array-like
        The second cartesian coordinate.
    z : scalar or array-like
        The third cartesian coordinate.

    Returns
    -------
    r : float or array
        The radial coordinate (in the same units as the inputs).
    lat : float or array
        The latitude in radians
    lng : float or array
        The longitude in radians
    """

    xsq = x ** 2
    ysq = y ** 2
    zsq = z ** 2

    r = (xsq + ysq + zsq) ** 0.5
    s = (xsq + ysq) ** 0.5

    if np.isscalar(x) and np.isscalar(y) and np.isscalar(z):
        lng = math.atan2(y, x)
        lat = math.atan2(z, s)
    else:
        lng = np.arctan2(y, x)
        lat = np.arctan2(z, s)

    return r, lat, lng


def spherical_to_cartesian(r, lat, lng):
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
    r : scalar or array-like
        The radial coordinate (in the same units as the inputs).
    lat : scalar or array-like
        The latitude in radians
    lng : scalar or array-like
        The longitude in radians

    Returns
    -------
    x : float or array
        The first cartesian coordinate.
    y : float or array
        The second cartesian coordinate.
    z : float or array
        The third cartesian coordinate.


    """

    if np.isscalar(r) and np.isscalar(lat) and np.isscalar(lng):
        x = r * math.cos(lat) * math.cos(lng)
        y = r * math.cos(lat) * math.sin(lng)
        z = r * math.sin(lat)
    else:
        x = r * np.cos(lat) * np.cos(lng)
        y = r * np.cos(lat) * np.sin(lng)
        z = r * np.sin(lat)

    return x, y, z
