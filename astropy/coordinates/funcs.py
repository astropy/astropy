# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains convenience functions for coordinate-related functionality.

This is generally just wrapping around the object-oriented coordinates
framework, but it is useful for some users who are used to more functional
interfaces.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .. import units as u
from .. import _erfa as erfa
from ..io import ascii
from ..utils import isiterable, data
from .sky_coordinate import SkyCoord
from .builtin_frames import GCRS, FK5
from .representation import SphericalRepresentation, CartesianRepresentation

__all__ = ['cartesian_to_spherical', 'spherical_to_cartesian', 'get_sun',
           'concatenate', 'get_constellation']


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
    if not hasattr(r, 'unit'):
        r = r * u.dimensionless_unscaled
    if not hasattr(lat, 'unit'):
        lat = lat * u.radian
    if not hasattr(lon, 'unit'):
        lon = lon * u.radian

    sph = SphericalRepresentation(distance=r, lat=lat, lon=lon)
    cart = sph.represent_as(CartesianRepresentation)

    return cart.x, cart.y, cart.z


def get_sun(time):
    """
    Determines the location of the sun at a given time, in
    geocentric coordinates.

    Parameters
    ----------
    table : `~astropy.time.Time`
        The time at which to compute the location of the sun.

    Returns
    -------
    newsc : `~astropy.coordinates.SkyCoord`
        The location of the sun as a `~astropy.coordinates.SkyCoord` in the
        `~astropy.coordinates.GCRS` frame.


    Notes
    -----
    The algorithm for determining the sun/earth relative position is based
    on the simplified version of VSOP2000 that is part of ERFA. Compared to
    JPL's ephemeris, it should be good to about 4 km (in the Sun-Earth
    vector) from 1900-2100 C.E., 8 km for the 1800-2200 span, and perhaps
    250 km over the 1000-3000.

    """
    earth_pv_helio, earth_pv_bary = erfa.epv00(time.jd1, time.jd2)
    x = -earth_pv_helio[..., 0, 0] * u.AU
    y = -earth_pv_helio[..., 0, 1] * u.AU
    z = -earth_pv_helio[..., 0, 2] * u.AU
    cartrep = CartesianRepresentation(x=x, y=y, z=z)
    return SkyCoord(cartrep, frame=GCRS)


def concatenate(coords):
    """
    Combine multiple coordinate objects into a single
    `~astropy.coordinates.SkyCoord`.

    "Coordinate objects" here mean frame objects with data,
    `~astropy.coordinates.SkyCoord`, or representation objects.  Currently,
    they must all be in the same frame, but in a future version this may be
    relaxed to allow inhomogenous sequences of objects.

    Parameters
    ----------
    coords : sequence of coordinate objects
        The objects to concatenate

    Returns
    -------
    cskycoord : SkyCoord
        A single sky coordinate with its data set to the concatenation of all
        the elements in ``coords``
    """
    if getattr(coords, 'isscalar', False) or not isiterable(coords):
        raise TypeError('The argument to concatenate must be iterable')
    return SkyCoord(coords)

_constellation_frame = FK5(equinox='B1875')
_constellation_data = {}

def get_constellation(coord, short=False):
    """
    Determines the constellation(s) of the coordinates this `SkyCoord`
    contains.

    Parameters
    ----------
    coords : coordinate object
        The object to determine the constellation
    short : bool
        If True, the returned names are the IAU-sanctioned abbreviated
        names.  Otherwise, full names for the constellations are used.

    Returns
    -------
    constellation : str or string array
        If ``coords`` contains a scalar coordinate, returns the name of the
        constellation.  If it is an array `SkyCoord`, it returns an array of
        names.
    """
    if not _constellation_data:
        cdata = data.get_pkg_data_contents('data/constellation_data_roman87.dat')
        ctable = ascii.read(cdata, names=['ral', 'rau', 'decl', 'name'])
        cnames = data.get_pkg_data_contents('data/constellation_names.dat')
        cnames_short_to_long = dict([(l[:3], l[4:])
                                     for l in cnames.split('\n')
                                     if not l.startswith('#')])

        _constellation_data['ctable'] = ctable
        _constellation_data['cnames_short_to_long'] = cnames_short_to_long
    else:
        ctable = _constellation_data['ctable']
        cnames_short_to_long = _constellation_data['cnames_short_to_long']

    constel_coord = coord.transform_to(_constellation_frame)
    if constel_coord.isscalar:
        constel_coord = [constel_coord]
        scalar = True
    else:
        scalar = False

    names = []
    for coo in constel_coord:
        rah = coo.ra.hour
        decd = coo.dec.deg
        for row in ctable:
            if row['ral'] < rah < row['rau'] and decd < row['decl']:
                if short:
                    names.append(row['name'])
                else:
                    names.append(cnames_short_to_long[row['name']])
                break

    if scalar:
        return names[0]
    else:
        return names

