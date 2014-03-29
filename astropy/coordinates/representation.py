import abc

import numpy as np
import astropy.units as u

from .angles import Angle, Longitude, Latitude
from .distances import Distance
from ..extern import six


@six.add_metaclass(abc.ABCMeta)
class BaseRepresentation(object):
    """
    Base Representation object, for representing a point in a 3D system
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
        # The default is to convert via cartesian coordinates        
        return other_class.from_cartesian(self.to_cartesian())

    @classmethod
    def from_representation(cls, representation):
        return representation.represent_as(cls)


class CartesianRepresentation(BaseRepresentation):
    """
    Representation of a point on three cartesian axes, x, y and z
    
    Parameters
    ----------
    x: Quantity or float
        The x value, either a quanity or a value to be passed to Quantity with unit.

    y: Quantity or float
        The y value, either a quanity or a value to be passed to Quantity with unit.

    z: Quantity or float
        The z value, either a quanity or a value to be passed to Quantity with unit.
    
    unit: Unit
        Unit to initilize x, y and z with, if specified x, y and z will all be 
        converted to this unit.
    
    representation: BaseRepresentation
        A pre-existing Representation object to convert to Cartesian
    
    copy: bool
        If True arrays will be copied rather than referenced.
    """
    
    def __init__(self, x=None, y=None, z=None, unit=None, representation=None,
                 copy=True):
        
        if representation is not None:
            return
        
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

        self._x = u.Quantity(x, unit=unit, copy=copy)
        self._y = u.Quantity(y, unit=unit, copy=copy)
        self._z = u.Quantity(z, unit=unit, copy=copy)
        
    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z
        
    @classmethod
    def from_cartesian(self, other):
        return other

    def to_cartesian(self):
        return self


class SphericalRepresentation(BaseRepresentation):
    """
    Spherical Representation of a point based on Longitude, Latitude and Radius
    
    Parameters
    ----------
    lon: Londitude
        A Londitude object or a parameter to be parsed by Londitude.

    lat: Latitude
        A Latitude object or a parameter to be parsed by Latitude.

    distance: Distance
        Distance (Radius), either a Distance object or a parameter to be passed by Distance.

    representation: BaseRepresentation
        A pre-existing Representation object to convert to Spherical

    copy: bool
        If True arrays will be copied rather than referenced.
    """

    def __init__(self, lon=None, lat=None, distance=None, representation=None, copy=True):
        # Let the Longitude and Latitude classes deal with e.g. parsing
        self._lon = Longitude(lon, copy=copy)
        self._lat = Latitude(lat, copy=copy)
        
        if distance is not None:
            self._distance = Distance(distance, copy=copy)
        else:
            self._distance = None
    
    @property
    def lon(self):
        """ Longitude Value """
        return self._lon

    @property
    def lat(self):
        """ Latitude Value """
        return self._lat
    
    @property
    def distance(self):
        """ Radius Value """
        return self._distance
    
    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        
        if self.distance is None:
            raise ValueError("can only convert to cartesian coordinates if distance is set")

        x = self.distance * np.cos(self.lat) * np.cos(self.lon)
        y = self.distance * np.cos(self.lat) * np.sin(self.lon)
        z = self.distance * np.sin(self.lat)
        
        return CartesianRepresentation(x=x, y=y, z=z)
    
    def from_cartesian(self, cartesian_representation):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """
        
        xsq = cartesian_representation.x ** 2
        ysq = cartesian_representation.y ** 2
        zsq = cartesian_representation.z ** 2
        
        r = (xsq + ysq + zsq) ** 0.5
        s = (xsq + ysq) ** 0.5
        
        lon = np.arctan2(y, x)
        lat = np.arctan2(z, s)
        
        return SphericalRepresentation(lon=lon, lat=lat, distance=r)


class CylindricalRepresentation(BaseRepresentation):
    """
    Representation of a point in a Cylindrical system as rho, phi and z
    
    Parameters
    ----------
    rho: Distance
        The distance from the axis to the point. A Distance instance or a parameter to 
        be parsed by Distance.

    phi: Angle
        Angle around the axis. A Angle instance or a parameter to be parsed by 
        angle.

    z: Quantity
        Coordinate along the axis, a Quantity object. 

    representation: BaseRepresentation
        A pre-existing Representation object to convert to Cylindrical

    copy: bool
        If True arrays will be copied rather than referenced.
    """

    def __init__(self, rho=None, phi=None, z=None, representation=None, copy=True):
        
        self._rho = Distance(rho)
        self._phi = Angle(phi)
        
        if isinstance(z, u.Quantity) and z.unit.physical_type != 'length':
            raise u.UnitsError("z should have units of length")
        self._z = u.Quantity(z, copy=copy)
    
    @property
    def rho(self):
        return self._rho
    
    @property
    def phi(self):
        return self._phi
    
    @property
    def z(self):
        return self._z
    
    def from_cartesian(self, cartesian_representation):
        rho = np.sqrt(cartesian_representation.x**2 + cartesian_representation.y**2)

        if cartesian_representation.x == 0 and cartesian_representation.y == 0:
            phi = 0
        elif cartesian_representation.x >= 0:
            phi = np.arcsin(cartesian_representation.y / rho)
        else: #x < 0
            phi = -1 * np.arcsin(cartesian_representation.y / rho) + np.pi

        z = cartesian_representation.z
        
        return CylindricalRepresentation(rho=rho, phi=phi, z=z)

    def to_cartesian(self):
        x = self.rho * np.cos(self.phi)
        y = self.rho * np.sin(self.phi)
        z = self.z
        
        return CartesianRepresentation(x=x, y=y, z=z)
