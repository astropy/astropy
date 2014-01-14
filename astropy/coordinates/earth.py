# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from .. import units as u
from . import CartesianPoints, Longitude, Latitude

try:
    # Not guaranteed available at setup time
    from ..time import erfa_time
except ImportError:
    if not _ASTROPY_SETUP_:
        raise

__all__ = ['EarthLocation']

ELLIPSOIDS = {'WGS84': 1, 'GRS80': 2, 'WGS72': 3}


def _check_ellipsoid(ellipsoid=None, default='WGS84'):
    if ellipsoid is None:
        ellipsoid = default
    if ellipsoid not in ELLIPSOIDS:
        raise ValueError('Ellipsoid {0} not among known ones ({1})'
                         .format(ellipsoid, ELLIPSOIDS.keys()))
    return ellipsoid


class EarthLocation(CartesianPoints):
    """
    Location on Earth

    Initialization is first attempted assuming geocentric (x, y, z) coordinates
    are given; if that fails, another attempt is made assuming geodetic
    coordinates (longitude, latitude, height above a reference ellipsoid).
    Internally, the coordinates are stored as geocentric.

    To ensure a specific type of coordinates is used, either use the
    corresponding class methods (`from_geocentric` and `from_geodetic`) or
    initialize the arguments with names (`x`, `y`, `z` for geocentric;
    `lon`, `lat`, `height` for geodetic).  See the class methods for details.

    """

    ellipsoid = 'WGS84'

    def __new__(cls, *args, **kwargs):
        try:
            self = cls.from_geocentric(*args, **kwargs)
        except (u.UnitsError, TypeError) as exc_geocentric:
            try:
                self = cls.from_geodetic(*args, **kwargs)
            except Exception as exc_geodetic:
                raise TypeError('Coordinates could not be parsed as either '
                                'geocentric or geodetic, with respective '
                                'exceptions "{0}" and "{1}"'
                                .format(exc_geocentric, exc_geodetic))
        return self

    @classmethod
    def from_geocentric(cls, x, y=None, z=None, unit=None, dtype=None,
                        copy=True):
        """
        Location on Earth, initialized from geocentric coordinates

        Parameters
        ----------
        x, y, z : `~astropy.units.Quantity` or array-like
            Cartesian coordinates.  Have to have length units (m if not given).
            If `x` has length 3, `y` and `z` can be omitted.
        unit : `~astropy.units.UnitBase` object or None
            Physical unit of the coordinate values.  If `x`, `y`, and/or `z`
            are quantities, they will be converted to this unit.
        dtype : ~numpy.dtype, optional
            See `~astropy.units.Quantity`
        copy : bool, optional
            See `~astropy.units.Quantity`

        Raises
        ------
        astropy.units.UnitsError
            If the units on `x`, `y`, and `z` do not match or an invalid unit
            is given
        ValueError
            If `y` and `z` don't match `x`'s shape or `x` is not length-3
        TypeError
            If incompatible array types are passed into `x`, `y`, or `z`

        """
        return super(EarthLocation, cls).__new__(cls, x, y, z, unit,
                                                 dtype, copy)

    @classmethod
    def from_geodetic(cls, lon, lat, height=0., ellipsoid=None, dtype=None):
        """
        Location on Earth, initialized from geodetic coordinates

        Parameters
        ----------
        lon, lat : Angle or float
            Earth East longitude and latitude of observer.  Can be anything
            that initialises an `Angle` object (if float, should be decimal
            degrees).
        height : Quantity or float, optional
            Height above reference ellipsoid (if float, in meters; default: 0)
        ellipsoid : str, optional
            Name of the reference ellipsoid to use (default: 'WGS84')
            Available ellipoids are:  'WGS84', 'GRS80', 'WGS72'
        dtype : ~numpy.dtype, optional
            See `~astropy.units.Quantity`

        Raises
        ------
        astropy.units.UnitsError
            If the units on `lon` and `lat` are inconsistent with angular ones,
            or that on `height` with a length.
        ValueError
            If `lon`, `lat`, and `height` do not have the same shape, or
            if the ellipsoid is not recognized as among the ones implemented

        """
        ellipsoid = _check_ellipsoid(ellipsoid, default=cls.ellipsoid)
        lon = Longitude(lon, u.degree, wrap_angle=180*u.degree, copy=False)
        lat = Latitude(lat, u.degree, copy=False)
        height = u.Quantity(height, u.meter, copy=False)
        if not (lon.shape == lat.shape == height.shape):
            raise ValueError('`lon`, `lat`, and `height` do not have the '
                             'same length')
        xyz = erfa_time.era_gd2gc(ELLIPSOIDS[ellipsoid],
                                  np.atleast_1d(lon.to(u.radian).value),
                                  np.atleast_1d(lat.to(u.radian).value),
                                  np.atleast_1d(height.value))
        self = super(EarthLocation, cls).__new__(cls, xyz.squeeze(),
                                                 unit=u.meter, dtype=dtype,
                                                 copy=False)
        self.ellipsoid = ellipsoid
        return self

    @property
    def geodetic(self):
        """Convert to geodetic coordinates for the default ellipsoid."""
        return self.to_geodetic()

    def to_geodetic(self, ellipsoid=None):
        """Convert to geodetic coordinates

        Parameters
        ----------
        ellipsoid : str, optional
            Reference ellipsoid to use.  Default is the one the coordinates
            were initialized with.  Available are: 'WGS84', 'GRS80', 'WGS72'

        Returns
        -------
        (lon, lat, height) : tuple
            The tuple contains instances of `~astropy.coordinates.Longitude`,
            `~astropy.coordinates.Latitude`, and `~astropy.units.Quantity`
        """

        ellipsoid = _check_ellipsoid(ellipsoid, default=self.ellipsoid)
        self_value = self.to(u.meter).value
        if self_value.ndim == 1:
            self_value = self_value.reshape(-1, 1)
        lon, lat, height = erfa_time.era_gc2gd(ELLIPSOIDS[ellipsoid],
                                               self_value)
        return (Longitude(lon.squeeze() * u.radian, u.degree,
                          wrap_angle=180.*u.degree),
                Latitude(lat.squeeze() * u.radian, u.degree),
                u.Quantity(height.squeeze(), u.meter))

    @property
    def longitude(self):
        return self.geodetic[0]

    @property
    def latitude(self):
        return self.geodetic[1]

    @property
    def height(self):
        return self.geodetic[2]

    # mostly for symmetry with geodedic and to_geodetic
    @property
    def geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return self.to_xyz()

    def to_geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return (self.x, self.y, self.z)
