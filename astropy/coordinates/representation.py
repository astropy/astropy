"""
In this module, we define the coordinate representation classes, which are
used to represent low-level cartesian, spherical, cylindrical, and other
coordinates.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import abc
import functools

import numpy as np
import astropy.units as u

from .angles import Angle, Longitude, Latitude
from .distances import Distance
from ..extern import six
from ..utils import OrderedDict
from ..utils.compat.numpy import broadcast_arrays

__all__ = ["BaseRepresentation", "CartesianRepresentation",
           "SphericalRepresentation", "UnitSphericalRepresentation",
           "PhysicsSphericalRepresentation", "CylindricalRepresentation"]

NUMPY_LT_1P7 = [int(x) for x in np.__version__.split('.')[:2]] < [1, 7]

# Module-level dict mapping representation string alias names to class.
# This is populated by the metaclass init so all representation classes
# get registered automatically.
REPRESENTATION_CLASSES = {}

class MetaBaseRepresentation(type):
    def __init__(cls, name, bases, dct):
        super(MetaBaseRepresentation, cls).__init__(name, bases, dct)

        if name != 'BaseRepresentation' and 'attr_classes' not in dct:
            raise NotImplementedError('Representations must have an '
                                      '"attr_classes" class attribute.')

        # Register representation name (except for BaseRepresentation)
        if cls.__name__ == 'BaseRepresentation':
            return

        REPRESENTATION_CLASSES[cls.get_name()] = cls


def _fstyle(precision, x):
    fmt_str = '{0:.{precision}f}'
    s = fmt_str.format(x, precision=precision)
    s_trunc = s.rstrip('0')
    if s_trunc[-1] == '.':
        # Ensure there is one trailing 0 after a bare decimal point
        return s_trunc + '0'
    else:
        return s_trunc


@six.add_metaclass(MetaBaseRepresentation)
class BaseRepresentation(object):
    """
    Base Representation object, for representing a point in a 3D coordinate
    system.

    Notes
    -----
    All representation classes should subclass this base representation
    class. All subclasses should then define a ``to_cartesian`` method and a
    ``from_cartesian`` class method. By default, transformations are done via
    the cartesian system, but classes that want to define a smarter
    transformation path can overload the ``represent_as`` method.
    Furthermore, all classes must define an ``attr_classes`` attribute, an
    `~collections.OrderedDict` which maps component names to the class that
    creates them.  They can also define a `recommended_units` dictionary, which
    maps component names to the units they are best presented to users in.  Note
    that frame classes may override this with their own preferred units.
    """

    recommended_units = {}  # subclasses can override

    def represent_as(self, other_class):
        if other_class == self.__class__:
            return self
        else:
            # The default is to convert via cartesian coordinates
            return other_class.from_cartesian(self.to_cartesian())

    @classmethod
    def from_representation(cls, representation):
        return representation.represent_as(cls)

    # Should be replaced by abstractclassmethod once we support only Python 3
    @abc.abstractmethod
    def from_cartesian(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def to_cartesian(self):
        raise NotImplementedError()

    @property
    def components(self):
        """A tuple with the in-order names of the coordinate components"""
        return tuple(self.attr_classes)

    @classmethod
    def get_name(cls):
        name = cls.__name__.lower()
        if name.endswith('representation'):
            name = name[:-14]
        return name

    def __getitem__(self, view):
        return self.__class__(*[getattr(self, component)[view]
                                for component in self.components])

    def __len__(self):
        if self.isscalar:
            raise TypeError("'{cls}' object with scalar values have no "
                            "len()".format(cls=self.__class__.__name__))
        else:
            return len(getattr(self, self.components[0]))

    def __nonzero__(self):  # Py 2.x
        return self.isscalar or len(self) != 0

    def __bool__(self):  # Py 3.x
        return self.isscalar or len(self) != 0

    @property
    def shape(self):
        return getattr(self, self.components[0]).shape

    @property
    def isscalar(self):
        return getattr(self, self.components[0]).isscalar

    @property
    def _values(self):
        """Turn the coordinates into a record array with the coordinate values.

        The record array fields will have the component names.
        """
        allcomp = np.array([getattr(self, component).value
                            for component in self.components])
        dtype = np.dtype([(str(component), getattr(self, component).dtype)
                          for component in self.components])
        return (np.rollaxis(allcomp, 0, len(allcomp.shape))
                .copy().view(dtype).squeeze())

    @property
    def _units(self):
        """Return a dictionary with the units of the coordinate components."""
        return dict([(component, getattr(self, component).unit)
                     for component in self.components])

    @property
    def _unitstr(self):
        units_set = set(self._units.values())
        if len(units_set) == 1:
            unitstr = units_set.pop().to_string()
        else:
            unitstr = '({0})'.format(
                ', '.join([self._units[component].to_string()
                           for component in self.components]))
        return unitstr

    def __str__(self):
        return '{0} {1:s}'.format(self._values, self._unitstr)

    def __repr__(self):
        prefixstr = '    '

        if self._values.shape == ():
            v = [tuple([self._values[nm] for nm in self._values.dtype.names])]
            v = np.array(v, dtype=self._values.dtype)
        else:
            v = self._values

        names = self._values.dtype.names
        precision = np.get_printoptions()['precision']
        fstyle = functools.partial(_fstyle, precision)
        format_val = lambda val: np.array2string(val, style=fstyle)
        formatter = {
            'numpystr': lambda x: '({0})'.format(
                ', '.join(format_val(x[name]) for name in names))
        }

        if NUMPY_LT_1P7:
            arrstr = np.array2string(v, separator=', ',
                                     prefix=prefixstr)

        else:
            arrstr = np.array2string(v, formatter=formatter,
                                     separator=', ',
                                     prefix=prefixstr)

        if self._values.shape == ():
            arrstr = arrstr[1:-1]

        unitstr = ('in ' + self._unitstr) if self._unitstr else '[dimensionless]'
        return '<{0} ({1}) {2:s}\n{3}{4}>'.format(
            self.__class__.__name__, ', '.join(self.components),
            unitstr, prefixstr, arrstr)


class CartesianRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cartesian coordinates.

    Parameters
    ----------
    x, y, z : `~astropy.units.Quantity`
        The x, y, and z coordinates of the point(s). If ``x``, ``y``, and
        ``z`` have different shapes, they should be broadcastable.

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('x', u.Quantity),
                                ('y', u.Quantity),
                                ('z', u.Quantity)])

    def __init__(self, x, y=None, z=None, copy=True):

        if y is None and z is None:
            x, y, z = x
        elif (y is None and z is not None) or (y is not None and z is None):
            raise ValueError("x, y, and z are required to instantiate CartesianRepresentation")

        if not isinstance(x, self.attr_classes['x']):
            raise TypeError('x should be a {0}'.format(self.attr_classes['x'].__name__))

        if not isinstance(y, self.attr_classes['x']):
            raise TypeError('y should be a {0}'.format(self.attr_classes['y'].__name__))

        if not isinstance(z, self.attr_classes['x']):
            raise TypeError('z should be a {0}'.format(self.attr_classes['z'].__name__))

        x = self.attr_classes['x'](x, copy=copy)
        y = self.attr_classes['y'](y, copy=copy)
        z = self.attr_classes['z'](z, copy=copy)

        if not (x.unit.physical_type == y.unit.physical_type == z.unit.physical_type):
            raise u.UnitsError("x, y, and z should have matching physical types")

        try:
            x, y, z = broadcast_arrays(x, y, z, subok=True)
        except ValueError:
            raise ValueError("Input parameters x, y, and z cannot be broadcast")

        self._x = x
        self._y = y
        self._z = z

    @property
    def x(self):
        """
        The x component of the point(s).
        """
        return self._x

    @property
    def y(self):
        """
        The y component of the point(s).
        """
        return self._y

    @property
    def z(self):
        """
        The z component of the point(s).
        """
        return self._z

    @property
    def xyz(self):
        return u.Quantity((self._x, self._y, self._z))

    @classmethod
    def from_cartesian(cls, other):
        return other

    def to_cartesian(self):
        return self


class UnitSphericalRepresentation(BaseRepresentation):
    """
    Representation of points on a unit sphere.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity` or str
        The longitude and latitude of the point(s), in angular units. The
        latitude should be between -90 and 90 degrees, and the longitude will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}

    def __init__(self, lon, lat, copy=True):

        if not isinstance(lon, u.Quantity) or isinstance(lon, Latitude):
            raise TypeError('lon should be a Quantity, Angle, or Longitude')

        if not isinstance(lat, u.Quantity) or isinstance(lat, Longitude):
            raise TypeError('lat should be a Quantity, Angle, or Latitude')
        # Let the Longitude and Latitude classes deal with e.g. parsing
        lon = self.attr_classes['lon'](lon, copy=copy)
        lat = self.attr_classes['lat'](lat, copy=copy)

        try:
            lon, lat = broadcast_arrays(lon, lat, subok=True)
        except ValueError:
            raise ValueError("Input parameters lon and lat cannot be broadcast")

        self._lon = lon
        self._lat = lat

    @property
    def lon(self):
        """
        The longitude of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The latitude of the point(s).
        """
        return self._lat

    # TODO: implement represent_as for efficient transformations

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        x = u.one * np.cos(self.lat) * np.cos(self.lon)
        y = u.one * np.cos(self.lat) * np.sin(self.lon)
        z = u.one * np.sin(self.lat)

        return CartesianRepresentation(x=x, y=y, z=z)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        s = np.hypot(cart.x, cart.y)

        lon = np.arctan2(cart.y, cart.x)
        lat = np.arctan2(cart.z, s)

        return cls(lon=lon, lat=lat)

    def represent_as(self, other_class):
        # Take a short cut if the other clsss is a spherical representation
        if issubclass(other_class, PhysicsSphericalRepresentation):
            return other_class(phi=self.lon,theta=90 * u.deg - self.lat, r=1.0)
        elif issubclass(other_class, SphericalRepresentation):
            return other_class(lon=self.lon, lat=self.lat, distance=1.0)
        else:
            return super(UnitSphericalRepresentation, self).represent_as(other_class)


class SphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity`
        The longitude and latitude of the point(s), in angular units. The
        latitude should be between -90 and 90 degrees, and the longitude will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    distance : `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, otherwise
        it is passed to the :class:`~astropy.units.Quantity` class.

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude),
                                ('distance', u.Quantity)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}

    _unitrep = UnitSphericalRepresentation

    def __init__(self, lon, lat, distance, copy=True):

        if not isinstance(lon, u.Quantity) or isinstance(lon, Latitude):
            raise TypeError('lon should be a Quantity, Angle, or Longitude')

        if not isinstance(lat, u.Quantity) or isinstance(lat, Longitude):
            raise TypeError('lat should be a Quantity, Angle, or Latitude')

        # Let the Longitude and Latitude classes deal with e.g. parsing
        lon = self.attr_classes['lon'](lon, copy=copy)
        lat = self.attr_classes['lat'](lat, copy=copy)

        distance = self.attr_classes['distance'](distance, copy=copy)
        if distance.unit.physical_type == 'length':
            distance = distance.view(Distance)

        try:
            lon, lat, distance = broadcast_arrays(lon, lat, distance,
                                                  subok=True)
        except ValueError:
            raise ValueError("Input parameters lon, lat, and distance cannot be broadcast")

        self._lon = lon
        self._lat = lat
        self._distance = distance

    @property
    def lon(self):
        """
        The longitude of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The latitude of the point(s).
        """
        return self._lat

    @property
    def distance(self):
        """
        The distance from the origin to the point(s).
        """
        return self._distance

    def represent_as(self, other_class):
        # Take a short cut if the other class is a spherical representation
        if issubclass(other_class, PhysicsSphericalRepresentation):
            return other_class(phi=self.lon, theta=90 * u.deg - self.lat,
                               r=self.distance)
        elif issubclass(other_class, UnitSphericalRepresentation):
            return other_class(lon=self.lon, lat=self.lat)
        else:
            return super(SphericalRepresentation, self).represent_as(other_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        # We need to convert Distance to Quantity to allow negative values.
        if isinstance(self.distance, Distance):
            d = self.distance.view(u.Quantity)
        else:
            d = self.distance

        x = d * np.cos(self.lat) * np.cos(self.lon)
        y = d * np.cos(self.lat) * np.sin(self.lon)
        z = d * np.sin(self.lat)

        return CartesianRepresentation(x=x, y=y, z=z)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        s = np.hypot(cart.x, cart.y)
        r = np.hypot(s, cart.z)

        lon = np.arctan2(cart.y, cart.x)
        lat = np.arctan2(cart.z, s)

        return cls(lon=lon, lat=lat, distance=r)


class PhysicsSphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates (using the physics
    convention of using ``phi`` and ``theta`` for azimuth and inclination
    from the pole).

    Parameters
    ----------
    phi, theta : `~astropy.units.Quantity` or str
        The azimuth and inclination of the point(s), in angular units. The
        inclination should be between 0 and 180 degrees, and the azimuth will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`.  If ``copy`` is False, `phi`
        will be changed inplace if it is not between 0 and 360 degrees.

    r : `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, otherwise
        it is passed to the :class:`~astropy.units.Quantity` class.

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('phi', Angle),
                                ('theta', Angle),
                                ('r', u.Quantity)])
    recommended_units = {'phi': u.deg, 'theta': u.deg}

    def __init__(self, phi, theta, r, copy=True):

        if not isinstance(phi, u.Quantity) or isinstance(phi, Latitude):
            raise TypeError('phi should be a Quantity or Angle')

        if not isinstance(theta, u.Quantity) or isinstance(theta, Longitude):
            raise TypeError('phi should be a Quantity or Angle')

        # Let the Longitude and Latitude classes deal with e.g. parsing
        phi = self.attr_classes['phi'](phi, copy=copy)
        theta = self.attr_classes['theta'](theta, copy=copy)

        # Wrap/validate phi/theta
        if copy:
            phi = phi.wrap_at(360 * u.deg)
        else:
            # necessary because the above version of `wrap_at` has to be a copy
            phi.wrap_at(360 * u.deg, inplace=True)
        if np.any(theta.value < 0.) or np.any(theta.value > 180.):
            raise ValueError('Inclination angle(s) must be within 0 deg <= angle <= 180 deg, '
                             'got {0}'.format(theta.to(u.degree)))

        r = self.attr_classes['r'](r, copy=copy)
        if r.unit.physical_type == 'length':
            r = r.view(Distance)

        try:
            phi, theta, r = broadcast_arrays(phi, theta, r, subok=True)
        except ValueError:
            raise ValueError("Input parameters phi, theta, and r cannot be broadcast")

        self._phi = phi
        self._theta = theta
        self._distance = r

    @property
    def phi(self):
        """
        The azimuth of the point(s).
        """
        return self._phi

    @property
    def theta(self):
        """
        The elevation of the point(s).
        """
        return self._theta

    @property
    def r(self):
        """
        The distance from the origin to the point(s).
        """
        return self._distance

    def represent_as(self, other_class):
        # Take a short cut if the other clsss is a spherical representation
        if issubclass(other_class, SphericalRepresentation):
            return other_class(lon=self.phi, lat=90 * u.deg - self.theta,
                               distance=self.r)
        elif issubclass(other_class, UnitSphericalRepresentation):
            return other_class(lon=self.phi, lat=90 * u.deg - self.theta)
        else:
            return super(PhysicsSphericalRepresentation, self).represent_as(other_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        # We need to convert Distance to Quantity to allow negative values.
        if isinstance(self.r, Distance):
            d = self.r.view(u.Quantity)
        else:
            d = self.r

        x = d * np.sin(self.theta) * np.cos(self.phi)
        y = d * np.sin(self.theta) * np.sin(self.phi)
        z = d * np.cos(self.theta)

        return CartesianRepresentation(x=x, y=y, z=z)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        s = np.hypot(cart.x, cart.y)
        r = np.hypot(s, cart.z)

        phi = np.arctan2(cart.y, cart.x)
        theta = np.arctan2(s, cart.z)

        return cls(phi=phi, theta=theta, r=r)


class CylindricalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cylindrical coordinates.

    Parameters
    ----------
    rho : `~astropy.units.Quantity`
        The distance from the z axis to the point(s).

    phi : `~astropy.units.Quantity`
        The azimuth of the point(s), in angular units, which will be wrapped
        to an angle between 0 and 360 degrees. This can also be instances of
        `~astropy.coordinates.Angle`,

    z : `~astropy.units.Quantity`
        The z coordinate(s) of the point(s)

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('rho', u.Quantity),
                                ('phi', Angle),
                                ('z', u.Quantity)])
    recommended_units = {'phi': u.deg}

    def __init__(self, rho, phi, z, copy=True):

        if not isinstance(phi, u.Quantity) or isinstance(phi, Latitude):
            raise TypeError('phi should be a Quantity or Angle')

        rho = self.attr_classes['rho'](rho, copy=copy)
        phi = self.attr_classes['phi'](phi, copy=copy)
        z = self.attr_classes['z'](z, copy=copy)

        if not (rho.unit.physical_type == z.unit.physical_type):
            raise u.UnitsError("rho and z should have matching physical types")

        try:
            rho, phi, z = broadcast_arrays(rho, phi, z, subok=True)
        except ValueError:
            raise ValueError("Input parameters rho, phi, and z cannot be broadcast")

        self._rho = rho
        self._phi = phi
        self._z = z

    @property
    def rho(self):
        """
        The distance of the point(s) from the z-axis.
        """
        return self._rho

    @property
    def phi(self):
        """
        The azimuth of the point(s).
        """
        return self._phi

    @property
    def z(self):
        """
        The height of the point(s).
        """
        return self._z

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to cylindrical polar
        coordinates.
        """

        rho = np.hypot(cart.x, cart.y)
        phi = np.arctan2(cart.y, cart.x)
        z = cart.z

        return cls(rho=rho, phi=phi, z=z)

    def to_cartesian(self):
        """
        Converts cylindrical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        x = self.rho * np.cos(self.phi)
        y = self.rho * np.sin(self.phi)
        z = self.z

        return CartesianRepresentation(x=x, y=y, z=z)
