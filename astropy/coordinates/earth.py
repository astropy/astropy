# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from urllib2 import URLError
from warnings import warn

import numpy as np
from .. import units as u
from ..utils.exceptions import AstropyUserWarning
from . import Longitude, Latitude

try:
    # Not guaranteed available at setup time.
    from .. import _erfa as erfa
except ImportError:
    if not _ASTROPY_SETUP_:
        raise

__all__ = ['EarthLocation']

# Available ellipsoids (defined in erfam.h, with numbers exposed in erfa).
ELLIPSOIDS = ('WGS84', 'GRS80', 'WGS72')

_NO_ONLINE_SITES_WARNING_MSG = ('Could not access the online site list. Falling'
                                'back on the builtin version.')


def _check_ellipsoid(ellipsoid=None, default='WGS84'):
    if ellipsoid is None:
        ellipsoid = default
    if ellipsoid not in ELLIPSOIDS:
        raise ValueError('Ellipsoid {0} not among known ones ({1})'
                         .format(ellipsoid, ELLIPSOIDS))
    return ellipsoid


class EarthLocation(u.Quantity):
    """
    Location on the Earth.

    Initialization is first attempted assuming geocentric (x, y, z) coordinates
    are given; if that fails, another attempt is made assuming geodetic
    coordinates (longitude, latitude, height above a reference ellipsoid).
    When using the geodetic forms, Longitudes are measured increasing to the
    east, so west longitudes are negative. Internally, the coordinates are
    stored as geocentric.

    To ensure a specific type of coordinates is used, use the corresponding
    class methods (`from_geocentric` and `from_geodetic`) or initialize the
    arguments with names (``x``, ``y``, ``z`` for geocentric; ``lon``, ``lat``,
    ``height`` for geodetic).  See the class methods for details.


    Notes
    -----
    This class fits into the coordinates transformation framework in that it
    encodes a position on the `~astropy.coordinates.ITRS` frame.  To get a
    proper `~astropy.coordinates.ITRS` object from this object, use the ``itrs``
    property.
    """

    _ellipsoid = 'WGS84'
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
    def from_geocentric(cls, x, y, z, unit=None):
        """
        Location on Earth, initialized from geocentric coordinates.

        Parameters
        ----------
        x, y, z : `~astropy.units.Quantity` or array-like
            Cartesian coordinates.  If not quantities, ``unit`` should be given.
        unit : `~astropy.units.UnitBase` object or None
            Physical unit of the coordinate values.  If ``x``, ``y``, and/or
            ``z`` are quantities, they will be converted to this unit.

        Raises
        ------
        astropy.units.UnitsError
            If the units on ``x``, ``y``, and ``z`` do not match or an invalid
            unit is given.
        ValueError
            If the shapes of ``x``, ``y``, and ``z`` do not match.
        TypeError
            If ``x`` is not a `~astropy.units.Quantity` and no unit is given.
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
    def from_geodetic(cls, lon, lat, height=0., ellipsoid=None):
        """
        Location on Earth, initialized from geodetic coordinates.

        Parameters
        ----------
        lon : `~astropy.coordinates.Longitude` or float
            Earth East longitude.  Can be anything that initialises an
            `~astropy.coordinates.Angle` object (if float, in degrees).
        lat : `~astropy.coordinates.Latitude` or float
            Earth latitude.  Can be anything that initialises an
            `~astropy.coordinates.Latitude` object (if float, in degrees).
        height : `~astropy.units.Quantity` or float, optional
            Height above reference ellipsoid (if float, in meters; default: 0).
        ellipsoid : str, optional
            Name of the reference ellipsoid to use (default: 'WGS84').
            Available ellipsoids are:  'WGS84', 'GRS80', 'WGS72'.

        Raises
        ------
        astropy.units.UnitsError
            If the units on ``lon`` and ``lat`` are inconsistent with angular
            ones, or that on ``height`` with a length.
        ValueError
            If ``lon``, ``lat``, and ``height`` do not have the same shape, or
            if ``ellipsoid`` is not recognized as among the ones implemented.

        Notes
        -----
        For the conversion to geocentric coordinates, the ERFA routine
        ``gd2gc`` is used.  See https://github.com/liberfa/erfa
        """
        ellipsoid = _check_ellipsoid(ellipsoid, default=cls._ellipsoid)
        lon = Longitude(lon, u.degree, wrap_angle=180*u.degree, copy=False)
        lat = Latitude(lat, u.degree, copy=False)
        # don't convert to m by default, so we can use the height unit below.
        if not isinstance(height, u.Quantity):
            height = u.Quantity(height, u.m, copy=False)
        # convert to float in units required for erfa routine, and ensure
        # all broadcast to same shape, and are at least 1-dimensional.
        _lon, _lat, _height = np.broadcast_arrays(lon.to(u.radian).value,
                                                  lat.to(u.radian).value,
                                                  height.to(u.m).value)
        # get geocentric coordinates. Have to give one-dimensional array.
        xyz = erfa.gd2gc(getattr(erfa, ellipsoid), _lon.ravel(),
                                  _lat.ravel(), _height.ravel())
        self = xyz.view(cls._location_dtype, cls).reshape(lon.shape)
        self._unit = u.meter
        self._ellipsoid = ellipsoid
        return self.to(height.unit)

    @classmethod
    def of_site(cls, site_name):
        """
        Return an object of this class for a known observatory/site by name.

        This is intended as a quick convenience function to get basic site
        information, not a fully-featured exhaustive registry of observatories
        and all their properties.

        Note that when this function is called, it will first attempt to
        download site information from the data.astropy.org server.  If it
        cannot (i.e., an internet connection is not available), it will fall
        back on the list included with astropy (which is a limited and dated set
        of sites).

        Parameters
        ----------
        site_name : str
            Name of the observatory (case-insensitive).

        Returns
        -------
        site : This class (a `~astropy.coordinates.EarthLocation` or subclass)
            The location of the observatory.

        See Also
        --------
        get_site_names : the list of sites that this function can access
        """
        # need to import inside function to avoid circular dependencies
        from .sites import get_site

        try:
            el = get_site(site_name, online=True)  # this is always an EarthLocation
        except URLError:
            warn(AstropyUserWarning(_NO_ONLINE_SITES_WARNING_MSG))
            el = get_site(site_name, online=False)  # this is always an EarthLocation

        if cls is el.__class__:
            return el
        else:
            return cls.from_geodetic(*el.to_geodetic())

    @classmethod
    def get_site_names(cls):
        """
        Get list of names of observatories for use with
        `~astropy.coordinates.EarthLocation.of_site`.

        Note that when this function is called, it will first attempt to
        download site information from the data.astropy.org server.  If it
        cannot (i.e., an internet connection is not available), it will fall
        back on the list included with astropy (which is a limited and dated set
        of sites).


        Returns
        -------
        names : list of str
            List of valid observatory names

        See Also
        --------
        of_site : Gets the actual location object for one of the sites names
                  this returns.
        """
        # need to import inside function to avoid circular dependencies
        from .sites import get_site_names

        try:
            return get_site_names(online=True)  # this is always an EarthLocation
        except URLError:
            warn(AstropyUserWarning(_NO_ONLINE_SITES_WARNING_MSG))
            return get_site_names(online=False)  # this is always an EarthLocation

    @property
    def ellipsoid(self):
        """The default ellipsoid used to convert to geodetic coordinates."""
        return self._ellipsoid

    @ellipsoid.setter
    def ellipsoid(self, ellipsoid):
        self._ellipsoid = _check_ellipsoid(ellipsoid)

    @property
    def geodetic(self):
        """Convert to geodetic coordinates for the default ellipsoid."""
        return self.to_geodetic()

    def to_geodetic(self, ellipsoid=None):
        """Convert to geodetic coordinates.

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

        Raises
        ------
        ValueError
            if ``ellipsoid`` is not recognized as among the ones implemented.

        Notes
        -----
        For the conversion to geodetic coordinates, the ERFA routine
        ``gc2gd`` is used.  See https://github.com/liberfa/erfa
        """
        ellipsoid = _check_ellipsoid(ellipsoid, default=self.ellipsoid)
        self_array = self.to(u.meter).view(self._array_dtype, np.ndarray)
        lon, lat, height = erfa.gc2gd(getattr(erfa, ellipsoid), self_array)
        return (Longitude(lon * u.radian, u.degree,
                          wrap_angle=180.*u.degree),
                Latitude(lat * u.radian, u.degree),
                u.Quantity(height * u.meter, self.unit))

    @property
    def longitude(self):
        """Longitude of the location, for the default ellipsoid."""
        return self.geodetic[0]

    @property
    def latitude(self):
        """Latitude of the location, for the default ellipsoid."""
        return self.geodetic[1]

    @property
    def height(self):
        """Height of the location, for the default ellipsoid."""
        return self.geodetic[2]

    # mostly for symmetry with geodetic and to_geodetic.
    @property
    def geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return self.to_geocentric()

    def to_geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return (self.x, self.y, self.z)

    @property
    def itrs(self):
        """
        Generates an `~astropy.coordinates.ITRS` object with the coordinates of
        this object.
        """
        #potential circular imports prevent this from being up top
        from .builtin_frames import ITRS

        return ITRS(x=self.x, y=self.y, z=self.z)

    @property
    def x(self):
        """The X component of the geocentric coordinates."""
        return self['x']

    @property
    def y(self):
        """The Y component of the geocentric coordinates."""
        return self['y']

    @property
    def z(self):
        """The Z component of the geocentric coordinates."""
        return self['z']

    def __getitem__(self, item):
        result = super(EarthLocation, self).__getitem__(item)
        if result.dtype is self.dtype:
            return result.view(self.__class__)
        else:
            return result.view(u.Quantity)

    def __array_finalize__(self, obj):
        super(EarthLocation, self).__array_finalize__(obj)
        if hasattr(obj, '_ellipsoid'):
            self._ellipsoid = obj._ellipsoid

    def __len__(self):
        if self.shape == ():
            raise IndexError('0-d EarthLocation arrays cannot be indexed')
        else:
            return super(EarthLocation, self).__len__()

    def to(self, unit, equivalencies=[]):
        array_view = self.view(self._array_dtype, u.Quantity)
        converted = array_view.to(unit, equivalencies)
        return self._new_view(converted.view(self.dtype).reshape(self.shape),
                              unit)
    to.__doc__ = u.Quantity.to.__doc__
