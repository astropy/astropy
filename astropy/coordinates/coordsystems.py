# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the base classes and frameworks for coordinate objects.
"""
from abc import ABCMeta, abstractproperty, abstractmethod

from .. import units as u
from .angles import RA, Dec, Angle, AngularSeparation
from .distance import *

__all__ = ['SphericalCoordinatesBase', 'Coordinates']


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

    .. note::
        A range of different possible parameters are available to
        initialize this coordinate, but not all parameters can be used
        together.  See the examples in the documentation for more.

    Parameters
    ----------
    coordstr : str
        A single string with the coordinates.  Cannot be used with
        `{latnm}` and `{longnm}` nor `x`/`y`/`z`.
    {longnm} : `~astropy.coordinates.angle.Angle`, float, int, str
        This must be given with `{latnm}`.
    {latnm} : `~astropy.coordinates.angle.Angle`, float, int, str
        This must be given with `{longnm}`.
    distance : `~astropy.coordinates.coordsystems.Distance`, optional
        This may be given with `{latnm}` and `{longnm}` or `coordstr`
        and not `x`, `y`, or `z`.  If not given, `None` (unit sphere)
        will be assumed.
    x : number
        The first cartesian coordinate. Must be given with `y` and `z`
        and not with `{longnm}` or `{latnm}` nor `coordstr`.
    y : number
        The second cartesian coordinate. Must be given with `x` and `z`
        and not with `{longnm}` or `{latnm}` nor `coordstr`.
    z : number
        The third cartesian coordinate. Must be given with `x` and `y`
        and not with `{longnm}` or `{latnm}` nor `coordstr`.
    cartpoint : `~astropy.coordinates.distance.CartesianPoint`
        A cartesian point with the coordinates.  Cannot be used with
        any other arguments.
    unit : `~astropy.units.UnitBase` or tuple

        * If `{longnm}` and `{latnm}` or `coordstr` are given:
            If the units cannot be determined from the angle values
            provided, they must be specified as a tuple. The first value
            in the tuple is paired with `{longnm}`, and the second with
            `{latnm}`. If `coordstr` is applied or `{latnm}` is a string
            and a single unit is  given, it is assumed to apply to
            `{longnm}`. Otherwise, a single unit is applied to both.

        * If `x`, `y`, and `z` are given:
            `unit` must be present have dimensions of length"""

    def _initialize_latlong(self, longname, latname, useradec, initargs, initkwargs):
        """
        Subclasses should use this to initialize standard lat/long-style
        coordinates.

        This recognizes both the lat/long style and the cartesian form.

        Parameters
        ----------
        longname : str
            The name of the longitude-like coordinate attribute
        latname : str
            The name of the latitude-like coordinate attribute
        useradec : bool
            If True, the `RA` and `Dec` classes will be used for the
            angles.  Otherwise, a basic `Angle` will be used.
        initargs : list
            The ``*args`` from the initializer
        initkwargs : dict
            The ``**kwargs`` from the initializer
        """
        initkwargs = dict(initkwargs)  # copy
        nargs = len(initargs)
        if nargs == 1:
            if isinstance(initargs[0], CartesianPoint):
                initkwargs['cartpoint'] = initargs[0]
            else:
                initkwargs['coordstr'] = initargs[0]
        if nargs > 1:
            if longname in initkwargs:
                raise TypeError("_initialize_latlong() got multiple values for"
                                " keyword argument '{0}'".format(longname))
            initkwargs[longname] = initargs[0]
        if nargs >= 2:
            if latname in initkwargs:
                raise TypeError("_initialize_latlong() got multiple values for"
                                " keyword argument '{0}'".format(latname))
            initkwargs[latname] = initargs[1]
        if nargs == 3:
            if 'distance' in initkwargs:
                raise TypeError("_initialize_latlong() got multiple values for"
                                " keyword argument 'distance'")
            initkwargs['distance'] = initargs[2]
        if nargs > 3:
            raise TypeError('_initialize_latlong() takes up to 3 positional '
                            ' arguments ({0} given)'.format(len(initargs)))

        unit = initkwargs.pop('unit', None)
        coordstr = initkwargs.pop('coordstr', None)
        longval = initkwargs.pop(longname, None)
        latval = initkwargs.pop(latname, None)
        distval = initkwargs.pop('distance', None)
        cartpoint = initkwargs.pop('cartpoint', None)
        x = initkwargs.pop('x', None)
        y = initkwargs.pop('y', None)
        z = initkwargs.pop('z', None)

        if len(initkwargs) > 0:
            raise TypeError('_initialize_latlong() got unexpected keyword '
                            'arguments {0}'.format(initkwargs.keys()))

        ll = longval is not None and latval is not None
        xyz = x is not None or y is not None or z is not None

        if (ll or coordstr is not None) and not xyz and cartpoint is None:
            # lat/long-style initialization

            units = [] if unit is None else unit

            if isinstance(units, tuple) or isinstance(units, list):
                if len(units) > 2:
                    raise ValueError('Cannot give more than 2 units while '
                                     'initializing a coordinate')
            elif isinstance(units, u.UnitBase) or isinstance(units, str):
                # Only a single unit given, which is fine.  If the arguments are
                # strings, assign it to just the long, otherwise both
                if coordstr is not None or isinstance(latval, basestring):
                    units = (units, )
                else:
                    units = (units, units)
            else:
                raise ValueError("The value for units must be given as a tuple, e.g. "
                                 "unit=(u.hour, u.degree). An object of type '{0}' "
                                 "was given.".format(type(units).__name__))

            if coordstr is not None:
                # need to try to parse the coordinate from a single argument
                # populates latval and longval variables, which then get made
                # into coordinates below
                x = coordstr
                if isinstance(coordstr, str):
                    parsed = False
                    if "," in x:
                        longval, latval = x.split(",")
                        parsed = True
                    elif "\t" in x:
                        longval, latval = x.split("\t")
                        parsed = True
                    elif len(x.split()) == 6:
                        longval = " ".join(x.split()[0:3])
                        latval = " ".join(x.split()[3:])
                        parsed = True
                    elif len(x.split()) == 2:
                        longval, latval = x.split()
                        parsed = True

                    if not parsed:
                        values = x.split()
                        i = 1
                        while i < len(values) and not parsed:
                            try:
                                longval = " ".join(values[0:i])
                                parsed = True
                            except:
                                i += 1

                        if parsed == True:
                            latval = " ".join(values[i:])

                    if not parsed:
                        msg = ("Could not parse {longname}/{latname} values "
                               "from the string provided: '{coordstr}'.")
                        raise ValueError(msg.format(longname=longname,
                                                    latname=latname,
                                                    coordstr=coordstr))
                else:
                    raise ValueError("A coordinate cannot be created with a value of type "
                                     "'{0}'.".format(type(coordstr).__name__))
            if useradec:
                longang = RA(longval, unit=units[0]) if len(units) > 0 else RA(longval)
                latang = Dec(latval, unit=units[1]) if len(units) > 1 else Dec(latval)
            else:
                if isinstance(longval, RA):
                    raise TypeError('Cannot provide an RA object to a non-RA/Dec system')
                if isinstance(latval, Dec):
                    raise TypeError('Cannot provide a Dec object to a non-RA/Dec system')
                longang = Angle(longval, unit=units[0]) if len(units) > 0 else Angle(longval)
                latang = Angle(latval, unit=units[1]) if len(units) > 1 else Angle(latval)
            dist = None if distval is None else Distance(distval)  # copy

        elif (xyz or cartpoint is not None) and not ll and distval is None and coordstr is None:
            #cartesian-style initialization
            if cartpoint is not None:
                if xyz or unit is not None:
                    raise ValueError('Cannot give both a CartesianPoint and x/y/z/units.')
                x = cartpoint.x
                y = cartpoint.y
                z = cartpoint.z
                unit = cartpoint.unit
            r, latval, longval = cartesian_to_spherical(x, y, z)

            if useradec:
                longang = RA(longval, unit=u.radian)
                latang = Dec(latval, unit=u.radian)
            else:
                longang = Angle(longval, unit=u.radian)
                latang = Angle(latval, unit=u.radian)

            dist = None if unit is None else Distance(r, unit)

        else:
            raise TypeError('Must initialize coordinates with '
                            '{latname}/{longname}/(distance) or x/y/z '
                            ''.format(latname=latname, longname=longname))
        setattr(self, longname, longang)
        setattr(self, latname, latang)
        self._distance = dist


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
            if self._distance is None:
                r = 1
                runit = None
            else:
                r = self._distance._value
                runit = self._distance._unit
            x, y, z = spherical_to_cartesian(r, self.latangle.radians,
                                                self.longangle.radians)
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
        other_in_self_system = other.transform_to(self.__class__)

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
        other_in_self_system = other.transform_to(self.__class__)

        if other_in_self_system._distance is None:
            raise ValueError('The other object does not have a distance; '
                             'cannot compute 3d separation.')

        dscale = other_in_self_system._distance._unit.to(self._distance._unit, 1)

        dx = self.x - other_in_self_system.x * dscale
        dy = self.y - other_in_self_system.y * dscale
        dz = self.z - other_in_self_system.z * dscale

        distval = (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5
        return Distance(distval, self._distance._unit)

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
        from copy import deepcopy
        from .transformations import master_transform_graph
        from .errors import ConvertError

        if tosys is self.__class__:
            return deepcopy(self)

        trans = master_transform_graph.get_transform(self.__class__, tosys)
        if trans is None:
            raise ConvertError('Cannot transform from {0} to '
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
            msg = "'{0}' object has no attribute '{1}', nor a transform."
            raise AttributeError(msg.format(self.__class__.__name__, name))


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
            if isinstance(units, u.UnitBase) or isinstance(units, str):
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
                if isinstance(units, u.UnitBase):
                    pass  # great!
                elif isinstance(units, tuple):
                    if len(units) == 0 or len(units) > 2:
                        raise ValueError("The units parameter only accepts "
                                         "tuples with one or two values.")
                    else:
                        # validate them
                        for a in units:
                            if not isinstance(a, u.UnitBase):
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
