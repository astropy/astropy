# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from .. import units as u
from . import Longitude, Latitude

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


class EarthLocation(u.Quantity):
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
    _location_dtype = np.dtype({'names': ['x', 'y', 'z'],
                                'formats': [np.float64]*3})
    _array_dtype = np.dtype((np.float64, (3,)))

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
    def from_geocentric(cls, x, y=None, z=None, unit=None):
        """
        Location on Earth, initialized from geocentric coordinates

        Parameters
        ----------
        x, y, z : `~astropy.units.Quantity` or array-like
            Cartesian coordinates.  If not quantities, `unit` should be given.
        unit : `~astropy.units.UnitBase` object or None
            Physical unit of the coordinate values.  If `x`, `y`, and/or `z`
            are quantities, they will be converted to this unit.

        Raises
        ------
        astropy.units.UnitsError
            If the units on `x`, `y`, and `z` do not match or an invalid unit
            is given
        ValueError
            If the shapes of 'x', 'y', and 'z' do not match
        TypeError
            If 'x' is not a `~astropy.units.Quantity` and no unit is given

        """
        if unit is None:
            try:
                unit = x.unit
            except AttributeError:
                raise TypeError("Geocentric coordinates should be Quantities "
                                "unless an explicit unit is given.")
        else:
            unit = u.Unit(unit)

        if unit.physical_type != 'length':
            raise u.UnitsError("Geocentric coordinates should be in "
                               "units of length.")

        try:
            x = u.Quantity(x, unit, copy=False)
            y = u.Quantity(y, unit, copy=False)
            z = u.Quantity(z, unit, copy=False)
        except u.UnitsError:
            raise u.UnitsError("Geocentric coordinate units should all be "
                               "consistent.")

        x, y, z = np.broadcast_arrays(x, y, z)
        struc = np.empty(x.shape, cls._location_dtype)
        struc['x'], struc['y'], struc['z'] = x, y, z
        return super(EarthLocation, cls).__new__(cls, struc, unit, copy=False)

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
        # convert to float in units required for erfa routine, and ensure
        # all broadcast to same shape, and are at least 1-dimensional
        _lon, _lat, _height = np.broadcast_arrays(lon.to(u.radian).value,
                                                  lat.to(u.radian).value,
                                                  height.value)
        # get geocentric coordinates. Have to give one-dimensional array.
        xyz = erfa_time.era_gd2gc(ELLIPSOIDS[ellipsoid], _lon.ravel(),
                                  _lat.ravel(), _height.ravel())
        self = xyz.view(cls._location_dtype, cls).reshape(lon.shape)
        self._unit = u.meter
        self.ellipsoid = ellipsoid
        return self.to(height.unit)

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

        self_array = self.to(u.meter).view(self._array_dtype, np.ndarray)
        lon, lat, height = erfa_time.era_gc2gd(ELLIPSOIDS[ellipsoid],
                                               np.atleast_2d(self_array))
        return (Longitude(lon.squeeze() * u.radian, u.degree,
                          wrap_angle=180.*u.degree),
                Latitude(lat.squeeze() * u.radian, u.degree),
                u.Quantity(height.squeeze() * u.meter, self.unit))

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
        return self.to_geocentric()

    def to_geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return (self.x, self.y, self.z)

    @property
    def x(self):
        return self['x']

    @property
    def y(self):
        return self['y']

    @property
    def z(self):
        return self['z']

    def __getitem__(self, item):
        result = super(EarthLocation, self).__getitem__(item)
        if result.dtype is self.dtype:
            return result.view(self.__class__)
        else:
            return result.view(u.Quantity)

    def __eq__(self, other):
        # pass by Quantity
        return np.ndarray.__eq__(self, other)

    def __ne__(self, other):
        # pass by Quantity
        return np.ndarray.__ne__(self, other)

    def __len__(self):
        if self.shape == ():
            raise IndexError('0-d EarthLocation arrays cannot be indexed')
        else:
            return super(EarthLocation, self).__len__()

    def to(self, unit, equivalencies=[]):
        array_view = self.view(self._array_dtype, u.Quantity)
        converted = array_view.to(unit, equivalencies)
        return converted.view(self.dtype,
                              self.__class__).reshape(self.shape)
    to.__doc__ = u.Quantity.to.__doc__
