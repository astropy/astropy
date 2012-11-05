# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the base classes and frameworks for coordinate objects.
"""
from abc import ABCMeta, abstractproperty, abstractmethod

import numpy as np

from .angles import RA, Dec, Angle, AngularSeparation
from .errors import UnitsError
from .. import units as u
from .. import cosmology

__all__ = ['SphericalCoordinatesBase', 'Coordinates', 'Distance',
           'CartesianPoint', 'cartesian_to_spherical', 'spherical_to_cartesian'
          ]


class SphericalCoordinatesBase(object):
    """
    Abstract superclass for all coordinate classes representing points in three
    dimensions.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, **kwargs):
        """
        Subclasses must override this, but they should also call this to set up
        internal state.
        """
        self._distance = None
        self._cartpoint = None

    _init_docstring_templ = """


    Parameters
    ----------
    {longnm} : `~astropy.coordinates.angle.Angle`, float, int, str
    {latnm} : `~astropy.coordinates.angle.Angle`, float, int, str
    unit : tuple
        If the units cannot be determined from the angle values provided, they must
        be specified as a tuple in the 'unit' parameter. The first value in the tuple
        is paired with the first angle provided, and the second with the second angle.
        (If the unit is specified in the first angle but not the second, the first
        value in the tuple may be 'None'.)
    """

    def _initialize_latlong(self, longname, latname, initargs, initkwargs):
        """
        Subclasses should use this to initialize standard lat/long-style
        coordinates.

        Parameters
        ----------
        longname : str
            The name of the longitude-like coordinate attribute
        latname : str
            The name of the latitude-like coordinate attribute
        initargs : list
            The ``*args`` from the initializer
        initkwargs : dict
            The ``**kwargs`` from the initializer
        """


        # Initialize values.
        # _ra, _dec are what we parse as potential values that still need validation
        _ra = None
        _dec = None
        self.ra = None
        self.dec = None

        if "unit" in kwargs:
            units = kwargs["unit"]
            del kwargs["unit"]
        else:
            units = list()

        if isinstance(units, tuple) or isinstance(units, list):
            pass  # good
        elif isinstance(units, u.Unit) or isinstance(units, str):
            # Only a single unit given, which is fine (assigned to 'ra').
            # The value, even if given as a tuple, is unpacked. Just make it
            # a tuple for consistency
            units = [units]
        else:
            raise ValueError("The value for units must be given as a tuple, e.g. "
                             "unit=(u.hour, u.degree). An object of type '{0}' "
                             "was given.".format(type(units).__name__))

        if len(args) == 0 and len(kwargs) == 0:
            raise ValueError("A coordinate object cannot be created without ra,dec values.")
        elif len(args) > 0 and len(kwargs) > 0:
            raise ValueError("The angle values can only be specified as keyword arguments "
                             "(e.g. ra=x, dec=y) or as a single value (e.g. a string) "
                             "not a combination.")
        elif len(args) == 0 and len(kwargs) > 0:
            # only "ra" and "dec" accepted as keyword arguments
            try:
                _ra = kwargs["ra"]
                _dec = kwargs["dec"]
            except KeyError:
                raise ValueError("When values are supplied as keyword arguments, both "
                                 "'ra' and 'dec' must be specified.")
            if isinstance(_ra, RA):
                self.ra = _ra
            if isinstance(_dec, Dec):
                self.dec = _dec

        elif len(args) == 1 and len(kwargs) == 0:
            # need to try to parge the coordinate from a single argument
            x = args[0]
            if isinstance(args[0], str):
                parsed = False
                if "," in x:
                    _ra, _dec = x.split(",")
                    parsed = True
                elif "\t" in x:
                    _ra, _dec = x.split("\t")
                    parsed = True
                elif len(x.split()) == 6:
                    _ra = " ".join(x.split()[0:3])
                    _dec = " ".join(x.split()[3:])
                    parsed = True
                elif len(x.split()) == 2:
                    _ra, _dec = x.split()
                    parsed = True

                if not parsed:
                    values = x.split()
                    i = 1
                    while i < len(values) and not parsed:
                        try:
                            self.ra = RA(" ".join(values[0:i]))
                            parsed = True
                        except:
                            i += 1

                    if parsed == True:
                        self.dec = Dec(" ".join(values[i:]))

                if not parsed:
                    raise ValueError("Could not parse ra,dec values from the string provided: '{0}'.".format(x))
            else:
                raise ValueError("A coordinate cannot be created with a value of type "
                                 "'{0}'.".format(type(args[0]).__name___))

        elif len(args) == 2 and len(kwargs) == 0:
            _ra = args[0]
            _dec = args[1]

        elif len(args) > 2 and len(kwargs) == 0:
            raise ValueError("More than two values were found where only ra and dec "
                             "were expected.")
        else:
            raise ValueError("Unable to create a coordinate using the values provided.")


#             # First try to see if RA, Dec objects were provided in the args.
#             for arg in args:
#                 if isinstance(arg, RA):
#                     _ra = arg
#                 elif isinstance(arg, Dec):
#                     _dec = arg
#
#             if None not in [_ra, _dec]:
#                 self.ra = _ra
#                 self.dec = _dec
#                 return
#             elif (_ra and not _dec) or (not _ra and _dec):
#                 raise ValueError("When an RA or Dec value is provided, the other "
#                                  "coordinate must also be given.")
#
#             # see if the whole coordinate might be parseable from arg[0]
#
#         try:
#             if isinstance(args[0], RA) and isinstance(args[1], Dec):
#                 _ra = args[0]
#                 _dec = args[1]
#             elif isinstance(args[1], RA) and isinstance(args[0], Dec):
#                 _ra = args[1]
#                 _dec = args[0]
#         except IndexError:
#             raise ValueError("Not enough parameters were provided.")

        if self.ra is None:
            self.ra = RA(_ra, unit=units[0]) if len(units) > 0 else RA(_ra)
        if self.dec is None:
            self.dec = Dec(_dec, unit=units[1]) if len(units) > 1 else Dec(_dec)


    @abstractproperty
    def latangle(self):
        """
        The latitudinal/elevation angle for these coordinates as an
        `~astropy.coorinates.angles.Angle` object.

        .. note ::
            This should be overridden in subclasses as a read-only property that
            just returns an attribute a way to abstract the exact choice of
            names for the coordiantes. E.g., `ICRSCoordinates` implements this
            by doing ``return self.ra``.
        """
        pass

    @abstractproperty
    def longangle(self):
        """
        The longitudinal/azimuthal angle for these coordinates as an
        `~astropy.coorinates.angles.Angle` object.

        .. note ::
            This should be overridden in subclasses as a read-only property that
            just returns an attribute a way to abstract the exact choice of
            names for the coordinates. E.g., `ICRSCoordinates` implements this
            by doing ``return self.dec``.
        """

    @property
    def distance(self):
        """
        The radial distance for this coordinate object as an
        `~astropy.coordinates.coordsystems.Distance` object.

        If set as a tuple, the tuple will be passed into the
        `~astropy.coordinates.coordsystems.Distance` constructor.

        Alternatively, this may be `None`, indicating an unknown/not given
        distance. Where necessary, this object will be interpreted as angles on
        the unit sphere.
        """
        return self._distance

    @distance.setter
    def distance(self, val):
        if val is None:
            self._distance = None
        elif isinstance(val, tuple):
            self._distance = Distance(*val)
        elif isinstance(val, Distance):
            self._distance = val
        else:
            raise TypeError('Spherical coordinate distance must be a ')

    @property
    def x(self):
        self._make_cart()
        return self._cartpoint.x

    @property
    def y(self):
        self._make_cart()
        return self._cartpoint.y

    @property
    def z(self):
        self._make_cart()
        return self._cartpoint.z

    @property
    def cartesian(self):
        self._make_cart()
        return self._cartpoint

    def _make_cart(self, override=False):
        if override or self._cartpoint is None:
            if self._distanceq is None:
                r = 1
                runit = None
            else:
                r = self._distance.value
                runit = self._distance.unit
            x, y, z = spherical_to_cartesian(r, self.latangle, self.longangle)
            self._cartpoint = CartesianPoint(x, y, z, runit)

    def separation(self, other):
        """
        Computes on-sky separation between this coordinate and another.

        Parameters
        ----------
        other : `~astropy.coordinates.coordsystems.SphericalCoordinatesBase`
            The coordinate system to get the separation to.

        Returns
        -------
        sep : `~astropy.coordinates.angles.AngularSeparation`
            The on-sky separation between this and the `other` coordinate.
        """
        other_in_self_system = other.convert_to(self.__class__)

        lat1 = self.latangle.radians
        long1 = other_in_self_system.latangle.radians
        lat2 = self.longangle.radians
        long2 = other_in_self_system.longangle.radians
        return AngularSeparation(lat1, long1, lat2, long2, u.radian)

    def separation3d(self, other):
        """
        Computes three dimensional separation between this coordinate
        and another.

        Parameters
        ----------
        other : `~astropy.coordinates.coordsystems.SphericalCoordinatesBase`
            The coordinate system to get the distance to.

        Returns
        -------
        sep : `~astropy.coordinates.coordsystems.Distance`
            The real-space distance between these two coordinates.

        Raises
        ------
        ValueError
            If this or the other coordinate do not have distances.
        """
        if self._distance is None:
            raise ValueError('This object does not have a distance; cannot '
                             'compute 3d separation.')

        # do this first just in case the conversion somehow creates a distance
        other_in_self_system = other.convert_to(self.__class__)

        if other_in_self_system._distance is None:
            raise ValueError('The other object does not have a distance; '
                             'cannot compute 3d separation.')

        dscale = other_in_self_system._distance.unit.to(self._distance.unit, 1)

        dx = self.x - other_in_self_system.x * dscale
        dy = self.y - other_in_self_system.y * dscale
        dz = self.z - other_in_self_system.z * dscale

        return (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5

    #<------------transformation-related stuff here-------------------->
    def transform_to(self, tosys):
        """
        Transform this coordinate to a new system.

        Parameters
        ----------
        tosys : class
            The system to transform this coordinate into.

        Returns
        -------
        transcoord
            A new object with this coordinate represented in the `tosys` system.

        Raises
        ------
        ValueError
            If there is no possible transformation route.
        """
        from .transformations import master_transform_graph

        trans = master_transform_graph.get_transform(self.__class__, tosys)
        if trans is None:
            raise ValueError('Cannot transform from {0} to '
                             '{1}'.format(self.__class__, tosys))
        return trans(self)

    def is_transformable_to(self, tosys):
        """
        Determines if this coordinate can be transformed to a particular system.

        Parameters
        ----------
        tosys : class
            The system to transform this coordinate into.

        Returns
        -------
        transformable : bool
            True if this can be trasnformed to `tosys`, False if not.
        """
        from .transformations import master_transform_graph

        trans = master_transform_graph.get_transform(self.__class__, tosys)
        return trans is not None

    def __getattr__(self, name):
        """
        Overrides getattr to return coordinates that this can be transformed
        to, based on the alias name in the master transform graph.
        """
        from .transformations import master_transform_graph

        nmsys = master_transform_graph.lookup_name(name)
        if nmsys is not None and self.is_transformable_to(nmsys):
            return self.transform_to(nmsys)
        else:
            objname = self.__class__.__name__
            raise AttributeError("'{0}' object has no attribute '{1}'".format(objname, name))

#FIXME: make this subclass Quantity once Quantity is in master
class Distance(object):
    """
    A one-dimensional distance.

    This can be initialized in one of two ways, using either a distance
    and a unit, or a redshift and (optionally) a cosmology.  `value`
    and `unit` may be provided as positional arguments, but `z` and
    `cosmology` are only valid as keyword arguments (see examples).

    Parameters
    ----------
    value : scalar
        The value of this distance
    unit : `~astropy.units.core.UnitBase`
        The units for this distance.  Must have dimensions of distance.
    z : float
        A redshift for this distance.  It will be converted to a distance
        by computing the luminosity distance for this redshift given the
        cosmology specified by `cosmology`.
    cosmology : `~astropy.cosmology.Cosmology` or None
        A cosmology that will be used to compute the distance from `z`.
        If None, the current cosmology will be used (see
        `astropy.cosmology` for details).

    Raises
    ------
    UnitsError
        If the `unit` is not a distance.

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.cosmology import WMAP3
    >>> d1 = Distance(10, u.Mpc)
    >>> d2 = Distance(40, unit=u.au)
    >>> d3 = Distance(value=5, unit=u.kpc)
    >>> d4 = Distance(z=0.23)
    >>> d5 = Distance(z=0.23, cosmology=WMAP3)
    """

    def __init__(self, *args, **kwargs):
        if 'z' in kwargs:
            z = kwargs.pop('z')
            cosmo = kwargs.pop('cosmology', None)
            if cosmo is None:
                cosmo = cosmology.get_current()

            if len(args) > 0 or len(kwargs) > 0:
                raise TypeError('Cannot give both distance and redshift')

            self._value = cosmo.luminosity_distance(z)
            self._unit = u.Mpc
        else:
            if len(args) == 0:
                value = kwargs.pop('value', None)
                unit = kwargs.pop('unit', None)
            if len(args) == 1:
                value = args[0]
                unit = kwargs.pop('unit', None)
            elif len(args) == 2:
                value, unit = args
            else:
                raise TypeError('Distance constructor cannot take more than 2 arguments')

            if len(kwargs) > 0:
                raise TypeError('Invalid keywords provided to Distance: ' +
                                str(kwargs.keys()))

            if value is None:
                raise ValueError('A value for the distance must be provided')
            if unit is None:
                raise UnitsError('A unit must be provided for distance.')

            if not unit.is_equivalent(u.m):
                raise UnitsError('provided unit for Distance is not a length')
            self._value = value
            self._unit = unit

    @property
    def lightyear(self):
        """
        The value of this distance in light years
        """
        return self._unit.to(u.lyr, self._value)

    @property
    def pc(self):
        """
        The value of this distance in parsecs
        """
        return self._unit.to(u.parsec, self._value)

    @property
    def kpc(self):
        """
        The value of this distance in kiloparsecs
        """
        return self._unit.to(u.kpc, self._value)

    @property
    def Mpc(self):
        """
        The value of this distance in megaparsecs
        """
        return self._unit.to(u.Mpc, self._value)

    @property
    def au(self):
        """
        The value of this distance in astronomical units
        """
        return self._unit.to(u.au, self._value)

    @property
    def m(self):
        """
        The value of this distance in meters
        """
        return self._unit.to(u.m, self._value)

    @property
    def km(self):
        """
        The value of this distance in kilometers
        """
        return self._unit.to(u.km, self._value)

    @property
    def z(self):
        """
        The redshift for this distance assuming its physical distance is
        a luminosity distance.

        .. note::
            This uses the "current" cosmology to determine the appropriate
            distance to redshift conversions.  See `astropy.cosmology`
            for details on how to change this.

        """
        return self.compute_z()

    def compute_z(self, cosmology=None):
        """
        The redshift for this distance assuming its physical distance is
        a luminosity distance.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.cosmology` or None
            The cosmology to assume for this calculation, or None to use the
            current cosmology.

        """
        from ..cosmology import luminosity_distance
        from scipy import optimize

        #FIXME: array: need to make this calculation more vector-friendly

        f = lambda z, d, cos: (luminosity_distance(z, cos) - d) ** 2
        return optimize.brent(f, (self.Mpc, cosmology))


class CartesianPoint(object):
    """
    A cartesian representation of a point in three-dimensional space.

    Attributes
    ----------
    x : number or array
        The first cartesian coordinate.
    y : number or array
        The second cartesian coordinate.
    z : number or array
        The third cartesian coordinate.
    unit : `~astropy.units.UnitBase` object or None
        The physical unit of the coordinate values.
    """

    def __init__(self, x, y, z, unit=None):
        self.x = x
        self.y = y
        self.z = z
        self.unit = unit

    def to_spherical(self):
        """
        Converts to the spherical representation of this point.

        Returns
        -------
        r : float or array
            The radial coordinate (in the same units as the inputs).
        lat : float or array
            The latitude in radians
        lng : float or array
            The longitude in radians

        """
        return cartesian_to_spherical(self.x, self.y, self.z)


class Coordinates(object):
    """
    A convenience factory class to create coordinate objects.

    This class can be used to create coordinate objects. The coordinate system is chosen
    based on the keywords used. For example, using the 'l' and 'b' keywords will return
    a `~astropy.coordinates.GalacticCoordinates` object. A "Coordinates" object cannot be
    created on its own.

    Parameters
    ----------
    (ra, dec) : `~astropy.coordinates.Angle`, str, float, int
        Right ascension and declination values. Returns an ICRSCoordinate object.
    (l, b) : `~astropy.coordinates.Angle`, str, float, int
        Galactic latitude and longitude. Returns a GalacticCoordinates object.
    (az, el) : `~astropy.coordinates.Angle`, str, float, int
        Azimuth and elevation values. Returns a HorizontaolCoordinates object.

    unit : `~astropy.units.Unit`, str, tuple
        Units must be provided for each of the angles provided. If the unit value can be
        determined from the Angle objects directly, `unit` does not need to be specified.
        If `unit` is a single value, it is applied to both of the given angles. If one angle
        requires a unit and the other does not, use `None` as a placeholder.
    """
    __meta__ = ABCMeta

    def __new__(self, *args, **kwargs):
        # coordinates, units=None, ra=None, dec=None, az=None, el=None, l=None, b=None):
        from .builtin_systems import ICRSCoordinates, GalacticCoordinates, HorizontalCoordinates

        #units = kwargs["unit"] if "unit" in kwargs.keys() else list()
        try:
            units = kwargs["unit"]
            if isinstance(units, u.Unit) or isinstance(units, str):
                units = (units, units)
        except KeyError:
            units = list()

        # first see if the keywords suggest what kind of coordinate is being requested.
        if "ra" in kwargs.keys() or "dec" in kwargs.keys():
            try:
                ra = kwargs["ra"]
                dec = kwargs["dec"]
            except KeyError:
                raise ValueError("When an 'ra' or 'dec' value is provided, the "
                                 "other coordinate must also be given.")
            for kw in ["l", "b", "az", "el"]:  # no others should be provided
                if kw in kwargs.keys():
                    raise ValueError("Conflicting coordinates were given.")
            ra = RA(ra, unit=units[0]) if len(units) > 0 else RA(ra)
            dec = Dec(dec, unit=units[1]) if len(units) > 1 else Dec(dec)
            return ICRSCoordinates(ra=ra, dec=dec)

        if "az" in kwargs.keys() or "el" in kwargs.keys():
            try:
                az = kwargs["az"]
                el = kwargs["el"]
            except KeyError:
                raise ValueError("When an 'az' or 'el' horizontal coordinates value "
                                 "is provided, the other coordinate must also be given.")
            for kw in ["ra", "dec", "l", "b"]:  # no others should be provided
                if kw in kwargs.keys():
                    raise ValueError("Conflicting coordinates were given.")
            az = Angle(az, unit=units[0]) if len(units) > 0 else Angle(az)
            el = Angle(el, unit=units[1]) if len(units) > 1 else Angle(el)
            return HorizontalCoordinates(az=az, el=el)

        if "l" in kwargs.keys() or "b" in kwargs.keys():
            try:
                l = kwargs["l"]
                b = kwargs["b"]
            except KeyError:
                raise ValueError("When an 'l' or 'b' galactic coordinates value is "
                                 "provided, the other coordinate must also be given.")
            for kw in ["ra", "dec", "az", "el"]:  # no others should be provided
                if kw in kwargs.keys():
                    raise ValueError("Conflicting coordinates were given.")
            l = Angle(l, unit=units[0]) if len(units) > 0 else Angle(l)
            b = Angle(b, unit=units[1]) if len(units) > 1 else Angle(b)
            return GalacticCoordinates(l=l, b=b)

        if len(args) == 1:
            x = args[0]

            if isinstance(x, str):
                raise ValueError("The coordinate system could not be determines from the value "
                                 "provided. Specify the system via keywords or use the "
                                 "corresponding class (e.g. GalacticCoordinate).")
            elif isinstance(x, list):
                return ValueError("Lists of coordinates are not yet supported")
            else:
                return ValueError("Could not create a Coordinate object from an object "
                                  "of type '{0}'.".format(type(x).__name__))
        if len(args) == 2:
            #a1, a2 = args[0:2]
            if isinstance(args[0], RA) and isinstance(args[1], Dec):
                return ICRSCoordinates(ra=args[0], dec=args[1])
            raise ValueError("Two angles were provided ('{0[0]}', '{0[1]}'), but the "
                             "coordinate system "
                             "was not provided. Specify the system via keywords or use the "
                             "corresponding class (e.g. GalacticCoordinate).".format(args))

        else:
            raise ValueError("Could not construct coordinates.")

        if False:  # old code - still useful?
            # determine units
            if units is not None:
                if isinstance(units, u.Unit):
                    pass  # great!
                elif isinstance(units, tuple):
                    if len(units) == 0 or len(units) > 2:
                        raise ValueError("The units parameter only accepts "
                                         "tuples with one or two values.")
                    else:
                        # validate them
                        for a in units:
                            if not isinstance(a, u.Unit):
                                raise ValueError("Units must be specified as u.degree, u.hour, etc.")
                    # units are valid
            else:
                #units were None - try to determine units from coordinate object
                if isinstance(coordinates, tuple):
                    if len(coordinates) != 2:
                        raise ValueError("Two coordinate values must be provided - '{0}' found.".format(len(coordinates)))
                    else:
                        # we have two values - the goal is to end up with two Angle
                        # objects.
                        pass


#<------------transformation-related utility functions----------------->

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
    import math

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
    import math

    if np.isscalar(r) and np.isscalar(lat) and np.isscalar(lng):
        x = r * math.cos(lat) * math.cos(lng)
        y = r * math.cos(lat) * math.sin(lng)
        z = r * math.sin(lat)
    else:
        x = r * np.cos(lat) * np.cos(lng)
        y = r * np.cos(lat) * np.sin(lng)
        z = r * np.sin(lat)

    return x, y, z
