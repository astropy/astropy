# In this file, we define the coordinate representation classes, which are
# used to represent low-level cartesian, spherical, cylindrical, and other #
# coordinate. All classes should define a to_cartesian method and a
# from_cartesian class method. By default, transformations are done via the
# cartesian system, but classes that want to define a smarter transformation
# path can overload the ``represent_as`` method.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import abc

import numpy as np
import astropy.units as u

from .angles import Angle, Longitude, Latitude
from .distances import Distance
from ..extern import six

# Suggestions to improve API
#
# - switch from allowing a representation to be passed as an argument to
# having a ``from_representation`` class method. Alternatively, we get rid of
# the representation argument and just allow the first argument to be a
# representation object, since this would already simplify __new__ somewhat.
#
# - change PhysicistSphericalRepresentation to PhysicsSphericalRepresentation
# (implemented below).

__all__ = ["CartesianRepresentation", "SphericalRepresentation",
           "PhysicsSphericalRepresentation", "CylindricalRepresentation"]


def broadcast_quantity(*args):
    """
    A Quantity-aware version of np.broadcast_arrays
    """
    new_arrays = np.broadcast_arrays(*args)
    new_quantities = []
    for i in range(len(new_arrays)):
        new_quantities.append(args[i].__class__(new_arrays[i] * args[i].unit))
    return tuple(new_quantities)


@six.add_metaclass(abc.ABCMeta)
class BaseRepresentation(object):
    """
    Base Representation object, for representing a point in a 3D coordinate system
    """

    def __new__(cls, *args, **kwargs):

        representation = kwargs.pop('representation', None)

        if representation is None:
            return super(BaseRepresentation, cls).__new__(cls)
        else:
            if any([x is not None for x in kwargs.values()]):
                raise ValueError("If representation is passed, no other arguments can be passed")
            return cls.from_representation(representation)

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


class CartesianRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cartesian coordinates.

    Parameters
    ----------
    x, y, z : `~astropy.units.Quantity` or float or `~numpy.ndarray`, optional
        The x, y, and z coordinates of the point(s), which should either be
        `~astropy.units.Quantity` instances, or can be passed as
        float or `numpy.ndarray` provided that the ``unit`` parameter is
        specified. If ``x``, ``y``, and ``z`` have different shapes, they
        should be broadcastable.

    unit : `~astropy.units.Unit`, optional
        If ``x``, ``y``, or ``z`` are specified as float or
        ``numpy.ndarray``, then ``unit`` should be specified to indicate the
        units for these parameters. If ``x``, ``y``, or ``z`` are
        `~astropy.units.Quantity` instances, and ``unit`` is specified, they
        are converted to ``unit``.

    representation : BaseRepresentation, optional
        A pre-existing Representation object to convert to cartesian
        coordinates.

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    def __init__(self, x=None, y=None, z=None, unit=None, representation=None,
                 copy=True):

        if representation is not None:
            return

        if x is None or y is None or z is None:
            raise ValueError('x, y, and z are required to instantiate CartesianRepresentation')

        if unit is not None:
            unit = u.Unit(unit)

        if isinstance(x, u.Quantity) and x.unit.physical_type != 'length':
            raise u.UnitsError("x should have units of length")

        if isinstance(y, u.Quantity) and y.unit.physical_type != 'length':
            raise u.UnitsError("y should have units of length")

        if isinstance(z, u.Quantity) and z.unit.physical_type != 'length':
            raise u.UnitsError("z should have units of length")

        if unit is not None and unit.physical_type != 'length':
            raise u.UnitsError("unit should be a unit of length")

        x = u.Quantity(x, unit=unit, copy=copy)
        y = u.Quantity(y, unit=unit, copy=copy)
        z = u.Quantity(z, unit=unit, copy=copy)

        try:
            x, y, z = broadcast_quantity(x, y, z)
        except ValueError:
            raise ValueError("Input parameters x, y, and z cannot be broadcast")

        self._x = x
        self._y = y
        self._z = z

    @property
    def x(self):
        """
        The x position of the point(s).
        """
        return self._x

    @property
    def y(self):
        """
        The y position of the point(s).
        """
        return self._y

    @property
    def z(self):
        """
        The z position of the point(s).
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


class SphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity` or str, optional
        The longitude and latitude of the point(s). The input values are
        passed to the `~astropy.coordinates.Longitude` and
        `~astropy.coordinates.Latitude` class respectively, so any valid
        input for these classes is acceptable. This includes
        `~astropy.units.Quantity` instances, strings, lists of strings, and
        so on. `~astropy.coordinates.Longitude` instances can only be passed
        to ``lon``, and `~astropy.coordinates.Latitude` instances can only be
        passed to ``lat``.

    distance : `~astropy.units.Quantity`, optional
        The distance to the point(s). The input value is passed to the
        `~astropy.coordinates.Distance` class, so any valid input to that
        class is acceptable.

    representation : BaseRepresentation, optional
        A pre-existing Representation object to convert to spherical
        coordinates.

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    def __init__(self, lon=None, lat=None, distance=None, representation=None, copy=True):

        if representation is not None:
            return

        if lon is None or lat is None:
            raise ValueError('lon and lat are required to instantiate SphericalRepresentation')

        # Let the Longitude and Latitude classes deal with e.g. parsing
        lon = Longitude(lon, copy=copy)
        lat = Latitude(lat, copy=copy)

        if distance is not None:
            distance = Distance(distance, copy=copy)
        else:
            distance = None

        try:
            if distance is None:
                lon, lat = broadcast_quantity(lon, lat)
            else:
                lon, lat, distance = broadcast_quantity(lon, lat, distance)
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
        # Take a short cut if the other clsss is a spherical representation
        if other_class is PhysicsSphericalRepresentation:
            return PhysicsSphericalRepresentation(phi=self.lon, theta=self.lat, distance=self.distance)
        else:
            return super(SphericalRepresentation, self).represent_as(other_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        if self.distance is None:
            raise ValueError("can only convert to cartesian coordinates if distance is set")

        # We need to convert Distance to Quantity to allow negative values.
        # At the moment, there is no easy way to convert Distance objects to
        # Quantity objects (https://github.com/astropy/astropy/issues/2259)
        d = self.distance.value * self.distance.unit

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

        xsq = cart.x ** 2
        ysq = cart.y ** 2
        zsq = cart.z ** 2

        r = (xsq + ysq + zsq) ** 0.5
        s = (xsq + ysq) ** 0.5

        lon = np.arctan2(cart.y, cart.x)
        lat = np.arctan2(cart.z, s)

        return SphericalRepresentation(lon=lon, lat=lat, distance=r)


class PhysicsSphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates (using the physics
    convention of using ``phi`` and ``theta`` for longitude and latitude).

    Parameters
    ----------
    phi, theta : `~astropy.units.Quantity` or str, optional
        The longitude and latitude of the point(s). The input values are
        passed to the `~astropy.coordinates.Longitude` and
        `~astropy.coordinates.Latitude` class respectively, so any valid
        input for these classes is acceptable. This includes
        `~astropy.units.Quantity` instances, strings, lists of strings, and
        so on. `~astropy.coordinates.Longitude` instances can only be passed
        to ``phi``, and `~astropy.coordinates.Latitude` instances can only be
        passed to ``theta``.

    distance : `~astropy.units.Quantity`, optional
        The distance to the point(s). The input value is passed to the
        `~astropy.coordinates.Distance` class, so any valid input to that
        class is acceptable.

    representation : BaseRepresentation, optional
        A pre-existing Representation object to convert to spherical
        coordinates.

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    def __init__(self, phi=None, theta=None, distance=None, representation=None, copy=True):

        if representation is not None:
            return

        if phi is None or theta is None:
            raise ValueError('phi and theta are required to instantiate PhysicsSphericalRepresentation')

        # Let the Longitude and Latitude classes deal with e.g. parsing
        phi = Longitude(phi, copy=copy)
        theta = Latitude(theta, copy=copy)

        if distance is not None:
            distance = Distance(distance, copy=copy)
        else:
            distance = None

        try:
            if distance is None:
                phi, theta = broadcast_quantity(phi, theta)
            else:
                phi, theta, distance = broadcast_quantity(phi, theta, distance)
        except ValueError:
            raise ValueError("Input parameters phi, theta, and distance cannot be broadcast")

        self._phi = phi
        self._theta = theta
        self._distance = distance

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
    def distance(self):
        """
        The distance from the origin to the point(s).
        """
        return self._distance

    def represent_as(self, other_class):
        # Take a short cut if the other clsss is a spherical representation
        if other_class is SphericalRepresentation:
            return SphericalRepresentation(lon=self.phi, lat=self.theta, distance=self.distance)
        else:
            return super(PhysicsSphericalRepresentation, self).represent_as(other_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        if self.distance is None:
            raise ValueError("can only convert to cartesian coordinates if distance is set")

        # We need to convert Distance to Quantity to allow negative values.
        # At the moment, there is no easy way to convert Distance objects to
        # Quantity objects (https://github.com/astropy/astropy/issues/2259)
        d = self.distance.value * self.distance.unit

        x = d * np.cos(self.theta) * np.cos(self.phi)
        y = d * np.cos(self.theta) * np.sin(self.phi)
        z = d * np.sin(self.theta)

        return CartesianRepresentation(x=x, y=y, z=z)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        xsq = cart.x ** 2
        ysq = cart.y ** 2
        zsq = cart.z ** 2

        r = (xsq + ysq + zsq) ** 0.5
        s = (xsq + ysq) ** 0.5

        phi = np.arctan2(cart.y, cart.x)
        theta = np.arctan2(cart.z, s)

        return PhysicsSphericalRepresentation(phi=phi, theta=theta, distance=r)


class CylindricalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cylindrical coordinates.

    Parameters
    ----------
    rho : `~astropy.units.Quantity`, optional
        The distance from the z axis to the point(s).

    phi : `~astropy.units.Quantity` or str, optional
        The azimuth of the point(s). The input is passed to the
        `~astropy.coordinates.Angle` class, so any valid input for that class
        is acceptable

    z : `~astropy.units.Quantity`, optional
        The z coordinate(s) of the point(s)

    representation : BaseRepresentation, optional
        A pre-existing Representation object to convert to cylindrical
        coordinates.

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    def __init__(self, rho=None, phi=None, z=None, representation=None, copy=True):

        if representation is not None:
            return

        if rho is None or phi is None or z is None:
            raise ValueError('rho, phi, and z are required to instantiate CylindricalRepresentation')

        rho = u.Quantity(rho, copy=copy)
        phi = Angle(phi, copy=copy)

        if isinstance(z, u.Quantity) and z.unit.physical_type != 'length':
            raise u.UnitsError("z should have units of length")

        z = u.Quantity(z, copy=copy)

        try:
            rho, phi, z = broadcast_quantity(rho, phi, z)
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

        rho = np.sqrt(cart.x ** 2 + cart.y ** 2)

        phi = np.zeros(cart.x.shape) * u.deg
        phi[rho > 0] = np.arctan2(cart.y, cart.x)

        z = cart.z

        return CylindricalRepresentation(rho=rho, phi=phi, z=z)

    def to_cartesian(self):
        """
        Converts cylindrical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        x = self.rho * np.cos(self.phi)
        y = self.rho * np.sin(self.phi)
        z = self.z

        return CartesianRepresentation(x=x, y=y, z=z)
