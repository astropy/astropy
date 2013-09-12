# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the base classes and frameworks for coordinate objects.
"""

import copy

from abc import ABCMeta, abstractproperty, abstractmethod

from ..extern import six
from .. import units as u
from .angles import Longitude, Latitude, Angle
from .distances import Distance, CartesianPoints, cartesian_to_spherical, spherical_to_cartesian
from ..utils.compat.misc import override__dir__
from . import angle_utilities
import numpy as np

__all__ = ['SphericalCoordinatesBase']


@six.add_metaclass(ABCMeta)
class SphericalCoordinatesBase(object):
    """
    Abstract superclass for all coordinate classes representing points
    in three dimensions.

    Notes
    -----
    Subclasses must implement `__init__`, and define the `latangle` and
    `lonangle` properties.  They may also override the `equinox`
    property, or leave it unaltered to indicate the coordinates are
    equinoxless.

    `_initialize_latlon` is provided to implement typical
    initialization features, and should be called from a subclass'
    `__init__`.  See the classes in
    `astropy.coordinates.builtin_systems` for examples of this.
    """

    @abstractmethod
    def __init__(self, *args, **kwargs):
        """
        Subclasses must override this, but they should also call this to set up
        internal state.
        """
        self._distance = None
        self._cartpoint = None

    def __eq__(self, other):
        try:
            return (self.latangle == other.latangle and
                    self.lonangle == other.lonangle and
                    self.distance == other.distance and
                    self.equinox == other.equinox)
        except AttributeError:
            return False

    _init_docstring_param_templ = """coordstr : str
        A single string with the coordinates.  Cannot be used with
        `{latnm}` and `{lonnm}` nor `x`/`y`/`z`.
    {lonnm} : `~astropy.coordinates.angle.Angle`, float, int, str
        This must be given with `{latnm}`.
    {latnm} : `~astropy.coordinates.angle.Angle`, float, int, str
        This must be given with `{lonnm}`.
    distance : `~astropy.coordinates.coordsystems.Distance`, optional
        This may be given with `{latnm}` and `{lonnm}` or `coordstr`
        and not `x`, `y`, or `z`.  If not given, `None` (unit sphere)
        will be assumed.
    x : number
        The first cartesian coordinate. Must be given with `y` and `z`
        and not with `{lonnm}` or `{latnm}` nor `coordstr`.
    y : number
        The second cartesian coordinate. Must be given with `x` and `z`
        and not with `{lonnm}` or `{latnm}` nor `coordstr`.
    z : number
        The third cartesian coordinate. Must be given with `x` and `y`
        and not with `{lonnm}` or `{latnm}` nor `coordstr`.
    cartpoint : `~astropy.coordinates.distance.CartesianPoints`
        A cartesian point with the coordinates.  Cannot be used with
        any other arguments.
    unit
        The `unit` parameter's interpretation depends on what other
        parameters are given:

        * If `{lonnm}` and `{latnm}` or `coordstr` are given:
            `unit` must be a length-2 sequence specifying the units of
            `{lonnm}` and `{latnm}`, respectively. They can be either
            `~astropy.units.UnitBase` objects or strings that will be
            converted using `~astropy.units.Unit`.  They can also be
            None to attempt to automatically interpret the units (see
            `~astropy.coordinates.angles.Angle` for details.) If `unit`
            is just `None`, this will be interpreted the same as
            ``(None, None)``.

        * If `x`, `y`, and `z` are given:
            `unit` must be a single unit with dimensions of length"""

    def _initialize_latlon(self, lonname, latname, initargs, initkwargs):
        """
        Subclasses should use this to initialize standard lat/lon-style
        coordinates.

        This recognizes both the lat/lon style and the cartesian form.

        Parameters
        ----------
        lonname : str
            The name of the longitude-like coordinate attribute
        latname : str
            The name of the latitude-like coordinate attribute
        initargs : list
            The ``*args`` from the initializer
        initkwargs : dict
            The ``**kwargs`` from the initializer
        """
        initkwargs = dict(initkwargs)  # copy
        nargs = len(initargs)
        sclsnm = self.__class__.__name__

        if nargs == 1:
            if isinstance(initargs[0], CartesianPoints):
                initkwargs['cartpoint'] = initargs[0]
            else:
                initkwargs['coordstr'] = initargs[0]
        if nargs > 1:
            if lonname in initkwargs:
                raise TypeError("{0} got multiple values for keyword argument "
                                "'{1}'".format(sclsnm, lonname))
            initkwargs[lonname] = initargs[0]
        if nargs >= 2:
            if latname in initkwargs:
                raise TypeError("{0} got multiple values for keyword argument "
                                "'{1}'".format(sclsnm, latname))
            initkwargs[latname] = initargs[1]
        if nargs > 2:
            raise TypeError('{0} takes up to 2 positional arguments '
                            '({1} given)'.format(sclsnm, len(initargs)))

        unit = initkwargs.pop('unit', None)
        coordstr = initkwargs.pop('coordstr', None)
        lonval = initkwargs.pop(lonname, None)
        latval = initkwargs.pop(latname, None)
        distval = initkwargs.pop('distance', None)
        cartpoint = initkwargs.pop('cartpoint', None)
        x = initkwargs.pop('x', None)
        y = initkwargs.pop('y', None)
        z = initkwargs.pop('z', None)

        if len(initkwargs) > 0:
            raise TypeError('{0} got unexpected keyword argument'
                            ' {1}'.format(sclsnm, initkwargs.keys()))

        angleinit = ((lonval is not None and latval is not None) or
                     coordstr is not None)
        cartinit = ((x is not None and y is not None and z is not None) or
                    cartpoint is not None)

        if angleinit and not cartinit:
            # lat/lon-style initialization
            for v in [x, y, z, cartpoint]:
                if v is not None:
                    raise ValueError('Cannot give both angular and cartesian '
                                     'coordinates while initializing ' + sclsnm)

            try:
                # this raises a TypeError if `unit` is not None or iterable
                units = [None, None] if unit is None else list(unit)
            except TypeError:
                raise ValueError('Must give a sequence of 2 units or None '
                                 'while initializing {0}. Instead got a '
                                 'non-sequence {1}'.format(sclsnm, unit))

            if len(units) == 2:
                try:
                    if units[0] is not None:
                        units[0] = u.Unit(units[0])
                    if units[1] is not None:
                        units[1] = u.Unit(units[1])
                except ValueError:
                    raise ValueError('Could not convert units to unit objects '
                                     'while initializing ' + sclsnm)
            else:
                raise ValueError('Must give a sequence of 2 units or None '
                             'while initializing {0}. Instead got a sequence '
                             'of {1}.'.format(sclsnm, len(units)))

            if coordstr is not None:
                # need to try to parse the coordinate from a single argument
                # populates latval and lonval variables, which then get made
                # into coordinates below
                if isinstance(coordstr, six.string_types):
                    if "," in coordstr:
                        lonval, latval = coordstr.split(",")
                    else:
                        coosplit = coordstr.split()
                        if len(coosplit) == 6:
                            lonval = " ".join(coosplit[0:3])
                            latval = " ".join(coosplit[3:])
                        elif len(coosplit) == 2:
                            lonval, latval = coosplit
                        else:
                            msg = ("Could not parse {lonname}/{latname} values "
                                   "from the string provided: '{coordstr}'.")
                            raise ValueError(msg.format(lonname=lonname,
                                                        latname=latname,
                                                        coordstr=coordstr))
                else:
                    raise ValueError("A {0} cannot be created with a single value of type "
                                     "'{1}', must be a string.".format(sclsnm, type(coordstr).__name__))

            # now actually create the angle objects
            lonang = Longitude(lonval, unit=units[0])
            latang = Latitude(latval, unit=units[1])

            dist = None if distval is None else Distance(distval)  # copy

        elif cartinit and not angleinit:
            # cartesian-style initialization
            for v in [coordstr, lonval, latval, distval]:
                if v is not None:
                    raise ValueError('Cannot give both angular and cartesian '
                                     'coordinates while initializing ' + sclsnm)

            if cartpoint is not None:
                for v in [x, y, z, unit]:
                    if v is not None:
                        raise ValueError('Cannot give both a CartesianPoints '
                                         'and x/y/z/unit parameters while '
                                         'initializing ' + sclsnm)
                x = cartpoint.x
                y = cartpoint.y
                z = cartpoint.z
                unit = cartpoint.unit
            r, latval, lonval = cartesian_to_spherical(x, y, z)

            lonang = Longitude(lonval, unit=u.radian)
            latang = Latitude(latval, unit=u.radian)

            dist = None if unit is None else Distance(r, unit)

        else:
            raise TypeError('Must initialize {coordnm} with '
                            '{latname}/{lonname}/(distance) or x/y/z '
                            ''.format(coordnm=sclsnm, latname=latname,
                                      lonname=lonname))

        # now actually set the values
        self._lonangle = lonang
        self._latangle = latang
        self._distance = dist


        #sanity-check that they are all consistent shapes
        if self.lonangle.shape != self.latangle.shape:
            raise ValueError('lonangle and latangle do not have matching shapes')

        if self._distance is not None and self._distance.shape != self.lonangle.shape:
            raise ValueError('distance and angles do not have matching shapes')

    def __repr__(self):
        if self.distance is not None:
            if self.isscalar:
                diststr = ', Distance={0:.2g} {1!s}'.format(self.distance.value, self.distance.unit)
            else:
                diststr = ', Distance={0} {1!s}'.format(self.distance.value, self.distance.unit)
        else:
            diststr = ''

        if self.isscalar:
            msg = "<{clsnm} {lonnm}={lonval:.5f} deg, {latnm}={latval:.5f} deg{diststr}>"
        else:
            msg = "<{clsnm} {lonnm}={lonval} deg, {latnm}={latval} deg{diststr}>"
        return msg.format(clsnm=self.__class__.__name__, lonval=self.lonangle.degree,
                          latval=self.latangle.degree, lonnm=self._repr_lon_name,
                          latnm=self._repr_lat_name, diststr=diststr)

    def __getitem__(self, key):
        from copy import deepcopy

        oldlat = self._latangle 
        oldlon = self._lonangle
        olddist = self._distance

        newlat = oldlat[key]
        newlon = oldlon[key]
        if olddist is not None:
            newdist = olddist[key]
        else:
            newdist = None

        try:
            #don't want to copy the old values, because we've already
            #copied them above as new*
            self._latangle = None
            self._lonangle = None
            self._distance = None

            newcoo =  deepcopy(self)

            newcoo._latangle = newlat
            newcoo._lonangle = newlon
            newcoo._distance = newdist

            return newcoo
        finally:
            self._latangle = oldlat
            self._lonangle = oldlon
            self._distance = olddist


    @property
    def latangle(self):
        """
        The latitudinal/elevation angle for these coordinates as an
        `~astropy.coordinates.angles.Angle` object.

        Subclasses will often provide properties returning the same object, but
        with a name more appropriate for the particular subclass.
        """
        return self._latangle

    @property
    def lonangle(self):
        """
        The longitudinal/azimuthal angle for these coordinates as an
        `~astropy.coordinates.angles.Angle` object.

        Subclasses will often provide properties returning the same object, but
        with a name more appropriate for the particular subclass.
        """
        return self._lonangle

    @property
    def equinox(self):
        """
        The equinox of this system, or None to indicate no equinox specified.
        """
        return None

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
        elif isinstance(val, u.Quantity):
            self._distance = Distance(val)
        else:
            raise TypeError(
                'Spherical coordinate distance must be a Distance object, a '
                'tuple that can be used to instantiate a Distance object, or '
                'None.')

        # must clear the old cached cartesian point, or it won't get updated
        # for the new distance
        self._cartpoint = None

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
                runit = u.dimensionless_unscaled
            else:
                r = self._distance.value
                runit = self._distance.unit
            x, y, z = spherical_to_cartesian(r, self.latangle.radian,
                                                self.lonangle.radian)
            self._cartpoint = CartesianPoints(x, y, z, runit)

    def separation(self, other):
        """
        Computes on-sky separation between this coordinate and another.

        Parameters
        ----------
        other : `~astropy.coordinates.coordsystems.SphericalCoordinatesBase`
            The coordinate to get the separation to.

        Returns
        -------
        sep : `~astropy.coordinates.angles.Angle`
            The on-sky separation between this and the `other` coordinate.

        Notes
        -----
        The separation is calculated using the Vincenty formula, which
        is stable at all locations, including poles and antipodes [1]_.

        .. [1] http://en.wikipedia.org/wiki/Great-circle_distance

        """
        other_in_self_system = other.transform_to(self.__class__)

        lon1 = self.lonangle
        lat1 = self.latangle
        lon2 = other_in_self_system.lonangle
        lat2 = other_in_self_system.latangle

        # Get the separation as a Quantity, convert to Angle in degrees
        sep = angle_utilities.angular_separation(lon1, lat1, lon2, lat2)
        return Angle(sep, unit=u.degree)

    def separation_3d(self, other):
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

        dx = self.x - other_in_self_system.x
        dy = self.y - other_in_self_system.y
        dz = self.z - other_in_self_system.z

        distval = (dx.value ** 2 + dy.value ** 2 + dz.value ** 2) ** 0.5
        return Distance(distval, dx.unit)

    def match_to_catalog_3d(self, catalogcoord, nthneighbor=1):
        """
        Finds the nearest 3-dimensional matches of this coordinate to a set
        of catalog coordinates.

        This finds the 3-dimensional closest neighbor, which is only different
        from the on-sky distance if `distance` is set in this object or the 
        `catalogcoord` object.
        
        Parameters
        ----------
        catalogcoord : `~astropy.coordinates.SphericalCoordinatesBase`
            The base catalog in which to search for matches. Typically this 
            will be a coordinate object that is an array (i.e., 
            ``catalogcoord.isscalar == False``)
        nthneighbor : int, optional
            Which closest neighbor to search for.  Typically ``1`` is desired here,
            as that is correct for matching one set of coordinates to another.
            The next likely use case is ``2``, for matching a coordinate catalog
            against *itself* (``1`` is inappropriate because each point will find 
            itself as the closest match).

        Returns
        -------
        idx : integer array
            Indecies into `catalogcoord` to get the matched points for each 
            `matchcoord`. Shape matches this coordinate.
        sep2d : `~astropy.units.quantity.Angle` 
            The on-sky separation between the closest match for each `matchcoord` and 
            the `matchcoord`. Shape matches `matchcoord`.
        dist3d : `~astropy.units.quantity.Quantity` 
            The 3D distance between the closest match for each `matchcoord` and 
            the `matchcoord`. Shape matches this coordinate.

        Notes
        -----
        This method requires `scipy` to be installed or it will fail.

        See Also
        --------
        astropy.coordinates.matching.match_coordinates_3d
        """
        from .matching import match_coordinates_3d

        return match_coordinates_3d(self, catalogcoord, nthneighbor=nthneighbor, storekdtree=True)

    def match_to_catalog_sky(self, catalogcoord, nthneighbor=1):
        """
        Finds the nearest on-sky matches of this coordinate in a set of 
        catalog coordinates.
        
        Parameters
        ----------
        catalogcoord : `~astropy.coordinates.SphericalCoordinatesBase`
            The base catalog in which to search for matches. Typically this 
            will be a coordinate object that is an array (i.e., 
            ``catalogcoord.isscalar == False``)
        nthneighbor : int, optional
            Which closest neighbor to search for.  Typically ``1`` is desired here,
            as that is correct for matching one set of coordinates to another.
            The next likely use case is ``2``, for matching a coordinate catalog
            against *itself* (``1`` is inappropriate because each point will find 
            itself as the closest match).

        Returns
        -------
        idx : integer array
            Indecies into `catalogcoord` to get the matched points for each 
            `matchcoord`. Shape matches this coordinate.
        sep2d : `~astropy.units.quantity.Angle` 
            The on-sky separation between the closest match for each `matchcoord` and 
            the `matchcoord`. Shape matches `matchcoord`.
        dist3d : `~astropy.units.quantity.Quantity` 
            The 3D distance between the closest match for each `matchcoord` and 
            the `matchcoord`. Shape matches this coordinate.

        Notes
        -----
        This method requires `scipy` to be installed or it will fail.
        
        See Also
        --------
        astropy.coordinates.matching.match_coordinates_sky
        """
        from .matching import match_coordinates_sky

        return match_coordinates_sky(self, catalogcoord, nthneighbor=nthneighbor, storekdtree=True)

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
        from .errors import ConvertError

        if tosys is self.__class__:
            return copy.deepcopy(self)

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
        transformable : bool or str
            True if this can be trasnformed to `tosys`, False if not. The
            string 'same' if `tosys` is the same system as this object
            (i.e. no transformation is needed).
        """
        from .transformations import master_transform_graph

        if self.__class__ is tosys:
            return 'same'
        else:
            trans = master_transform_graph.get_transform(self.__class__, tosys)
            return trans is not None

    def __getattr__(self, name):
        """
        Overrides getattr to return coordinates that this can be transformed
        to, based on the alias name in the master transform graph.
        """
        from .transformations import master_transform_graph

        nmsys = master_transform_graph.lookup_name(name)
        if self.__class__ is nmsys:
            return self
        if nmsys is not None and self.is_transformable_to(nmsys):
            return self.transform_to(nmsys)
        else:
            msg = "'{0}' object has no attribute '{1}', nor a transform."
            raise AttributeError(msg.format(self.__class__.__name__, name))

    @override__dir__
    def __dir__(self):
        """
        Overriding the builtin `dir` behavior allows us to add the
        transforms available by aliases.  This also allows ipython
        tab-completion to know about the transforms.
        """
        from .transformations import master_transform_graph

        dir_items = set()

        # determine the aliases that this can be transformed to.
        for alias in master_transform_graph.get_aliases():
            tosys = master_transform_graph.lookup_name(alias)
            if self.is_transformable_to(tosys):
                dir_items.add(alias)

        return dir_items

    @property
    def isscalar(self):
        """
        True if this coordinate contains scalar angles/distances, or False if
        they are array-like
        """
        #assumes input-validation occurs and thus lat/lon/dist consistent
        return self.lonangle.isscalar


    # Name resolve
    @classmethod
    def from_name(cls, name):
        """
        Given a name, query the CDS name resolver to attempt to retrieve
        coordinate information for that object. The search database, sesame
        url, and  query timeout can be set through configuration items in
        `astropy.coordinates.name_resolve` -- see docstring for
        `astropy.coordinates.name_resolve.get_icrs_coordinates` for more
        information.

        Parameters
        ----------
        name : str
            The name of the object to get coordinates for, e.g. m42.

        Returns
        -------
        coord : SphericalCoordinatesBase
            Instance of a Coordinates class, specified by the class this is
            called on, e.g. if `Galactic.from_name('m42')`, will
            get an instance of `Galactic` representing the
            position of M42.
        """

        from .name_resolve import get_icrs_coordinates

        icrs = get_icrs_coordinates(name)
        if cls == icrs.__class__:
            return icrs
        else:
            return icrs.transform_to(cls)

    _default_string_style = 'dms'

    def to_string(self, style=None, **kwargs):
        """
        A string representation of the coordinates.

        See :meth:`astropy.coordinates.Angle.to_string` for details and keyword
        arguments (the two angles forming the coordinates are are both
        :class:`astropy.coordinates.Angle` instances). Keyword arguments are passed to
        :meth:`astropy.coordinates.Angle.to_string`.

        Parameters
        ----------
        style : {'hmsdms', 'dms', 'decimal', None}
            The formatting specification to use. These encode the three most
            common ways to represent coordinates. If `None` is passed, the
            defaults for the current coordinate class is used.
        kwargs
            Keyword arguments are passed to :meth:`astropy.coordinates.Angle.to_string`.
        """

        if style is None:
            style = self._default_string_style

        styles = {
                  'hmsdms': {'lonargs': {'unit':u.hour},
                             'latargs': {'unit':u.degree}},
                  'dms':    {'lonargs': {'unit':u.degree},
                             'latargs': {'unit':u.degree}},
                  'decimal':{'lonargs': {'unit':u.degree,'decimal':True},
                             'latargs': {'unit':u.degree,'decimal':True}}
                 }

        lonargs = kwargs.copy()
        latargs = kwargs.copy()

        if style in styles:
            lonargs.update(styles[style]['lonargs'])
            latargs.update(styles[style]['latargs'])
        else:
            raise ValueError('Invalid style.  Valid options are: '+",".join(styles))

        if np.isscalar(self.lonangle.value):
            coord_string = (self.lonangle.to_string(**lonargs)
                            + " " +
                            self.latangle.to_string(**latargs))
        else:
            coord_string = []
            for lonangle, latangle in zip(self.lonangle, self.latangle):
                coord_string += [(lonangle.to_string(**lonargs)
                                 + " " +
                                 latangle.to_string(**latargs))]

        if hasattr(coord_string,'decode'):
            return coord_string.decode()

        return coord_string
