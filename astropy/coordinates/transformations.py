# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the functionality for transforming points from
one coordinate system to another (e.g. equatorial to galactic). The public
interface for this functionality is in the file core.py.
"""

import math

import numpy as np

pi = math.pi
piotwo = pi / 2.

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
        lat = piotwo - math.atan2(s, z)
    else:
        lng = np.arctan2(y, x)
        lat = piotwo - np.arctan2(s, z)


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

    polarang = piotwo - lng

    if np.isscalar(r) and np.isscalar(lat) and np.isscalar(lng):
        x = r * math.sin(polarang) * math.cos(lng)
        y = r * math.sin(polarang) * math.sin(lng)
        z = r * math.cos(polarang)
    else:
        x = r * np.sin(polarang) * np.cos(lng)
        y = r * np.sin(polarang) * np.sin(lng)
        z = r * np.cos(polarang)

    return x,y,z