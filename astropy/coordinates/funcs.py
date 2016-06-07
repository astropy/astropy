# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains convenience functions for coordinate-related functionality.

This is generally just wrapping around the object-oriented coordinates
framework, but it is useful for some users who are used to more functional
interfaces.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from .. import units as u
from ..constants import c
from .. import _erfa as erfa
from ..io import ascii
from ..utils import isiterable, data
from .sky_coordinate import SkyCoord
from .builtin_frames import GCRS, PrecessedGeocentric
from .representation import SphericalRepresentation, CartesianRepresentation
from .builtin_frames.utils import get_jd12

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
    Determines the location of the sun at a given time (or times, if the input
    is an array `~astropy.time.Time` object), in geocentric coordinates.

    Parameters
    ----------
    time : `~astropy.time.Time`
        The time(s) at which to compute the location of the sun.

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
    earth_pv_helio, earth_pv_bary = erfa.epv00(*get_jd12(time, 'tdb'))

    # We have to manually do aberration because we're outputting directly into
    # GCRS
    earth_p = earth_pv_helio[..., 0, :]
    earth_v = earth_pv_bary[..., 1, :]

    # convert barycentric velocity to units of c, but keep as array for passing in to erfa
    earth_v /= c.to(u.au/u.d).value

    dsun = np.sqrt(np.sum(earth_p**2, axis=-1))
    invlorentz = (1-np.sum(earth_v**2, axis=-1))**0.5
    properdir = erfa.ab(earth_p/dsun.reshape(dsun.shape + (1,)),
                        -earth_v, dsun, invlorentz)

    cartrep = CartesianRepresentation(x=-dsun*properdir[..., 0] * u.AU,
                                      y=-dsun*properdir[..., 1] * u.AU,
                                      z=-dsun*properdir[..., 2] * u.AU)
    return SkyCoord(cartrep, frame=GCRS(obstime=time))


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


# global dictionary that caches repeatedly-needed info for get_constellation
_constellation_data = {}
def get_constellation(coord, short_name=False, constellation_list='iau'):
    """
    Determines the constellation(s) a given coordinate object contains.

    Parameters
    ----------
    coords : coordinate object
        The object to determine the constellation of.
    short_name : bool
        If True, the returned names are the IAU-sanctioned abbreviated
        names.  Otherwise, full names for the constellations are used.
    constellation_list : str
        The set of constellations to use.  Currently only ``'iau'`` is
        supported, meaning the 88 "modern" constellations endorsed by the IAU.

    Returns
    -------
    constellation : str or string array
        If ``coords`` contains a scalar coordinate, returns the name of the
        constellation.  If it is an array coordinate object, it returns an array
        of names.

    Notes
    -----
    To determine which constellation a point on the sky is in, this precesses
    to B1875, and then uses the Delporte boundaries of the 88 modern
    constellations, as tabulated by
    `Roman 1987 <http://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/42>`_.
    """
    if constellation_list != 'iau':
        raise ValueError("only 'iau' us currently supported for constellation_list")

    # read the data files and cache them if they haven't been already
    if not _constellation_data:
        cdata = data.get_pkg_data_contents('data/constellation_data_roman87.dat')
        ctable = ascii.read(cdata, names=['ral', 'rau', 'decl', 'name'])
        cnames = data.get_pkg_data_contents('data/constellation_names.dat', encoding='UTF8')
        cnames_short_to_long = dict([(l[:3], l[4:])
                                     for l in cnames.split('\n')
                                     if not l.startswith('#')])
        cnames_long = np.array([cnames_short_to_long[nm] for nm in ctable['name']])

        _constellation_data['ctable'] = ctable
        _constellation_data['cnames_long'] = cnames_long
    else:
        ctable = _constellation_data['ctable']
        cnames_long = _constellation_data['cnames_long']

    isscalar = coord.isscalar

    # if it is geocentric, we reproduce the frame but with the 1875 equinox,
    # which is where the constellations are defined
    constel_coord = coord.transform_to(PrecessedGeocentric(equinox='B1875'))
    if isscalar:
        rah = constel_coord.ra.ravel().hour
        decd = constel_coord.dec.ravel().deg
    else:
        rah = constel_coord.ra.hour
        decd = constel_coord.dec.deg

    constellidx = -np.ones(len(rah), dtype=int)

    notided = constellidx == -1  # should be all
    for i, row in enumerate(ctable):
        msk = (row['ral'] < rah) & (rah < row['rau']) & (decd > row['decl'])
        constellidx[notided & msk] = i
        notided = constellidx == -1
        if np.sum(notided) == 0:
            break
    else:
        raise ValueError('Could not find constellation for coordinates {0}'.format(constel_coord[notided]))

    if short_name:
        names = ctable['name'][constellidx]
    else:
        names = cnames_long[constellidx]

    if isscalar:
        return names[0]
    else:
        return names
