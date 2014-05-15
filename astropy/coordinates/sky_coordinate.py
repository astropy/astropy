from __future__ import (absolute_import, division, print_function, unicode_literals)

import collections

import numpy as np

from ..utils.compat.misc import override__dir__
from ..extern import six
from ..extern.six.moves import zip
from ..units import Unit
from .. import units as u

from .angles import Latitude, Longitude
from .baseframe import BaseCoordinateFrame, frame_transform_graph, GenericFrame
from .builtin_frames import NoFrame
from .representation import (BaseRepresentation, SphericalRepresentation,
                             UnitSphericalRepresentation)

__all__ = ['SkyCoord']

# Define a convenience mapping.  This is used like a module constants
# but is actually dynamically evaluated.
def FRAME_ATTR_NAMES_SET():
    """Set of all possible frame-specific attributes"""
    out = set()
    for frame_cls in frame_transform_graph.frame_set:
        for attr in frame_cls.frame_attr_names.keys():
            out.add(attr)
    return out


class SkyCoord(object):
    """
    A high-level celestial coordinate object.  This class represents celestial
    coordinates of various different types of coordinate systems using a single
    class.


    The actual frames/coordinate-systems available are the systems that are
    registered in the  `astropy.coordinates.frame_transform_graph` - E.g.,
    `~astropy.coordinates.ICRS`. see the documentation for specific frame classes
    to see the exact attributes available and how they are interpreted.

    Because it supports different kinds of frames, `SkyCoord` can be initialized
    in a variety of ways:

    * A single positional argument that is an
      `~astropy.coordinates.BaseCoordinateFrame` object.
    * A positional argument that is a string or list of strings of a form like
      "1h12m43.2s +1d12m43s", and a `frame` keyword (see Parameters section
      below)
    * A positional argument that is a string or list of strings of a form like
      "1:12:43.2s +1:12:43", a `unit` argument that is a 2-tuple (e.g.,
      ``(u.deg, u.deg)`` and a `frame` keyword (see Parameters section below)
    * A `frame` keyword (see Parameters section below), and keywords appropriate
      for that frame (matching the frame class' initializer).
    * A `~astropy.coordinates.BaseRepresentation` instance and a `frame` keyword
      (see Parameters section below).
    * Additional keywords will be interpreted as attributes of this `SkyCoord`,
      but will only be used by the containing coordinate if it is transformed.

    Parameters
    ----------
    frame : `~astropy.coordinates.BaseCoordinateFrame` class or string, optional
        The type of coordinate frame this `SkyCoord` should represent.  If a
        string, it should be one of the aliases available in
        `astropy.coordinates.frame_transform_graph`, typically an all-lower-case
        version of the system name.
    others
        Other parameters depend on the type of frame - see above.
    """

    def __init__(self, *args, **kwargs):

        # Parse the args and kwargs to assemble a sanitized and validated
        # kwargs dict for initializing attributes for this object and for
        # creating the internal self._sky_coord_frame object
        args = list(args)  # Make it mutable
        kwargs = self._parse_inputs(args, kwargs)

        # Set internal versions of object state attributes
        for attr in FRAME_ATTR_NAMES_SET():
            setattr(self, '_' + attr, kwargs[attr])

        # Set up the keyword args for creating the internal coordinate object.
        frame_cls = frame_transform_graph.lookup_name(kwargs['frame'])
        if frame_cls is None:
            frame_cls = NoFrame

        coord_kwargs = {}
        for attr, value in kwargs.items():
            if value is not None and (attr in frame_cls.preferred_attr_names
                                      or attr in frame_cls.frame_attr_names):
                coord_kwargs[attr] = value

        # Finally make the internal coordinate object.
        self._sky_coord_frame = frame_cls(**coord_kwargs)

    @property
    def frame_name(self):
        return self.frame.name  # This can be None so cannot use frame_transform_graph

    @property
    def frame(self):
        return self._sky_coord_frame

    def __len__(self):
        return len(self.frame)

    def __nonzero__(self):  # Py 2.x
        return self.frame.__nonzero__()

    def __bool__(self):  # Py 3.x
        return self.frame.__bool__()

    def __getitem__(self, item):
        self_frame = self._sky_coord_frame
        try:
            #first turn `self` into a mockup of the thing we want - we can copy
            # this to get all the right attributes
            self._sky_coord_frame = self_frame[item]
            return SkyCoord(self)
        finally:
            # now put back the right frame in self
            self._sky_coord_frame = self_frame

    def _parse_inputs(self, args, kwargs):
        """
        Assemble a validated and sanitized keyword args dict for instantiating a
        SkyCoord and coordinate object from the provided `args`, and `kwargs`.
        """
        valid_kwargs = {}

        # Put the SkyCoord attributes like frame, equinox, obstime, location
        # into valid_kwargs dict.  `Frame` could come from args or kwargs, so
        # set valid_kwargs['frame'] accordingly.  The others must be specified
        # by keyword args or else get a None default.  Pop them off of kwargs
        # in the process.
        frame = valid_kwargs['frame'] = _get_frame_name(args, kwargs)
        if frame is None:
            frame = NoFrame
        for attr in FRAME_ATTR_NAMES_SET():
            valid_kwargs[attr] = kwargs.pop(attr, None)

        # Get latitude and longitude units
        lon_unit, lat_unit = _get_units(args, kwargs)

        # Grab any frame-specific attr names like `ra` or `l` or `distance` from kwargs
        # and migrate to valid_kwargs.
        valid_kwargs.update(_get_preferred_attrs(frame, lon_unit, lat_unit, kwargs))

        # Error if anything is still left in kwargs
        if kwargs:
            raise ValueError('Unrecognized keyword argument(s) {0}'
                             .format(', '.join("'{0}'".format(key) for key in kwargs)))

        # Finally deal with the unnamed args.  This figures out what the arg[0] is
        # and returns a dict with appropriate key/values for initializing frame class.
        if args:
            if len(args) == 1:
                # One arg which must be a coordinate.  In this case
                # coord_kwargs will contain keys like 'ra', 'dec', 'distance'
                # along with any frame attributes like equinox or obstime which
                # were explicitly specified in the coordinate object (i.e. non-default).
                coord_kwargs = _parse_coordinate_arg(args[0], frame, lon_unit, lat_unit)

            elif len(args) == 2:
                # Must be longitude, latitude.
                frame_cls = frame_transform_graph.lookup_name(frame)
                if frame_cls is None:
                    frame_cls = NoFrame
                attr_name_for_type = dict((attr_type, name) for name, attr_type in
                                          frame_cls.preferred_attr_names.items())
                coord_kwargs = {}
                coord_kwargs[attr_name_for_type['lon']] = Longitude(args[0], unit=lon_unit)
                coord_kwargs[attr_name_for_type['lat']] = Latitude(args[1], unit=lat_unit)
            else:
                raise ValueError('Must supply no more than two positional arguments, got {}'
                                 .format(len(args)))

            # Copy the coord_kwargs into the final valid_kwargs dict.  For each
            # of the coord_kwargs ensure that there is no conflict with a value
            # specified by the user in the original kwargs.
            for attr, coord_value in coord_kwargs.items():
                if (attr in valid_kwargs
                        and valid_kwargs[attr] is not None
                        and valid_kwargs[attr] != coord_value):
                    raise ValueError("Coordinate attribute '{0}'={1!r} conflicts with "
                                     "keyword argument '{0}'={2!r}"
                                     .format(attr, coord_value, valid_kwargs[attr]))
                valid_kwargs[attr] = coord_value

        return valid_kwargs

    def transform_to(self, frame):
        """
        Transform this coordinate to a new frame.

        Parameters
        ----------
        frame : str or `BaseCoordinateFrame` class / instance or `SkyCoord` instance
            The frame to transform this coordinate into.

        Returns
        -------
        coord : `SkyCoord`
            A new object with this coordinate represented in the `frame` frame.

        Raises
        ------
        ValueError
            If there is no possible transformation route.
        """
        from astropy.coordinates.errors import ConvertError

        if frame is None or isinstance(self.frame, NoFrame):
            raise ValueError('Cannot transform to/from this SkyCoord because '
                             'the frame was not specified at creation.')

        frame_kwargs = {}

        # Frame name (string) or frame class?  Coerce into an instance.
        try:
            frame = _get_frame_class(frame)()
        except:
            pass

        if isinstance(frame, SkyCoord):
            frame = frame.frame  # Change to underlying coord frame instance

        if isinstance(frame, BaseCoordinateFrame):
            new_frame_cls = frame.__class__

            # Set the keyword args for making a new frame instance for the
            # transform.  If the supplied frame instance has a non-default
            # value set then use that, otherwise use the self attribute value
            # if it is not None.
            for attr in FRAME_ATTR_NAMES_SET():
                self_val = getattr(self, attr, None)
                frame_val = getattr(frame, attr, None)
                if (frame_val is not None and
                        attr not in frame._attr_names_with_defaults):
                    frame_kwargs[attr] = frame_val
                elif self_val is not None:
                    frame_kwargs[attr] = self_val
        else:
            raise ValueError('Transform `frame` must be a frame name, class, or instance')

        # Get the composite transform to the new frame
        trans = frame_transform_graph.get_transform(self.frame.__class__, new_frame_cls)
        if trans is None:
            raise ConvertError('Cannot transform from {0} to {1}'
                               .format(self.frame.__class__, new_frame_cls))

        # Make a generic frame which will accept all the frame kwargs that
        # are provided and allow for transforming through intermediate frames
        # which may require one or more of those kwargs.
        generic_frame = GenericFrame(frame_kwargs)

        # Do the transformation, returning a coordinate frame of the desired
        # final type (not generic).
        new_coord = trans(self.frame, generic_frame)

        # Finally make the new SkyCoord object from the `new_coord` and
        # remaining frame_kwargs that are not frame_attributes in `new_coord`.
        # We could remove overlaps here, but the init code is set up to accept
        # overlaps as long as the values are identical (which they must be).
        return self.__class__(new_coord, **frame_kwargs)

    def __getattr__(self, attr):
        """
        Overrides getattr to return coordinates that this can be transformed
        to, based on the alias attr in the master transform graph.
        """

        if self.frame_name == attr:
            return self  # Should this be a deepcopy of self?

        # Anything in the set of all possible frame_attr_names is handled
        # here. If the attr is relevant for the current frame then delegate
        # to self.frame otherwise get it from self._<attr>.
        if attr in FRAME_ATTR_NAMES_SET():
            if attr in self.frame.frame_attr_names:
                return getattr(self.frame, attr)
            else:
                return getattr(self, '_' + attr)

        # Some attributes might not fall in the above category but still
        # are available through self._sky_coord_frame.
        if not attr.startswith('_') and hasattr(self._sky_coord_frame, attr):
            return getattr(self._sky_coord_frame, attr)

        # Try to interpret as a new frame for transforming.
        frame_cls = frame_transform_graph.lookup_name(attr)
        if frame_cls is not None and self.frame.is_transformable_to(frame_cls):
            return self.transform_to(attr)

        # Fail
        raise AttributeError("'{0}' object has no attribute '{1}', nor a transform."
                             .format(self.__class__, attr))

    @override__dir__
    def __dir__(self):
        """
        Override the builtin `dir` behavior to include:
        - Transforms available by aliases
        - Attribute / methods of the underlying self.frame object
        """

        # determine the aliases that this can be transformed to.
        dir_values = set()
        for name in frame_transform_graph.get_names():
            frame_cls = frame_transform_graph.lookup_name(name)
            if self.frame.is_transformable_to(frame_cls):
                dir_values.add(name)

        # Add public attributes of self.frame
        dir_values.update(set(attr for attr in dir(self.frame) if not attr.startswith('_')))

        # Add all possible frame_attr_names
        dir_values.update(FRAME_ATTR_NAMES_SET())

        return dir_values

    def __repr__(self):
        clsnm = self.__class__.__name__
        coonm = self.frame.__class__.__name__

        s = '<{clsnm} ({coonm})'.format(**locals())
        crepr = repr(self.frame)
        return s + crepr[crepr.index(':'):]

    def to_string(self, style='decimal', **kwargs):
        """
        A string representation of the coordinates.

        The default styles definitions are::

          'decimal': 'lat': {'decimal': True, 'unit': "deg"}
                     'lon': {'decimal': True, 'unit': "deg"}
          'dms': 'lat': {'unit': "deg"}
                 'lon': {'unit': "deg"}
          'hmsdms': 'lat': {'alwayssign': True, 'pad': True, 'unit': "deg"}
                    'lon': {'pad': True, 'unit': "hour"}

        See :meth:`~astropy.coordinates.Angle.to_string` for details and
        keyword arguments (the two angles forming the coordinates are are
        both :class:`~astropy.coordinates.Angle` instances). Keyword
        arguments have precedence over the style defaults and are passed
        to :meth:`~astropy.coordinates.Angle.to_string`.

        Parameters
        ----------
        style : {'hmsdms', 'dms', 'decimal'}
            The formatting specification to use. These encode the three most
            common ways to represent coordinates. The default is `decimal`.
        kwargs
            Keyword args passed to :meth:`~astropy.coordinates.Angle.to_string`.
        """

        sph_coord = self.frame.represent_as(SphericalRepresentation)

        styles = {'hmsdms': {'lonargs': {'unit': u.hour, 'pad': True},
                             'latargs': {'unit': u.degree, 'pad': True, 'alwayssign': True}},
                  'dms': {'lonargs': {'unit': u.degree},
                          'latargs': {'unit': u.degree}},
                  'decimal': {'lonargs': {'unit': u.degree, 'decimal': True},
                              'latargs': {'unit': u.degree, 'decimal': True}}
                  }

        lonargs = {}
        latargs = {}

        if style in styles:
            lonargs.update(styles[style]['lonargs'])
            latargs.update(styles[style]['latargs'])
        else:
            raise ValueError('Invalid style.  Valid options are: {0}'.format(",".join(styles)))

        lonargs.update(kwargs)
        latargs.update(kwargs)

        if np.isscalar(sph_coord.lon.value):
            coord_string = (sph_coord.lon.to_string(**lonargs)
                            + " " +
                            sph_coord.lat.to_string(**latargs))
        else:
            coord_string = []
            for lonangle, latangle in zip(sph_coord.lon, sph_coord.lat):
                coord_string += [(lonangle.to_string(**lonargs)
                                 + " " +
                                 latangle.to_string(**latargs))]

        return coord_string

    # High-level convinience methods
    def separation(self, other):
        """
        Computes on-sky separation between this coordinate and another.

        Parameters
        ----------
        other : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate to get the separation to.

        Returns
        -------
        sep : `~astropy.coordinates.Angle`
            The on-sky separation between this and the ``other`` coordinate.

        Notes
        -----
        The separation is calculated using the Vincenty formula, which
        is stable at all locations, including poles and antipodes [1]_.

        .. [1] http://en.wikipedia.org/wiki/Great-circle_distance

        """
        from . import Angle
        from .angle_utilities import angular_separation

        if isinstance(other, SkyCoord):
            self_in_other_system = self.transform_to(other.frame)
        elif isinstance(other, BaseCoordinateFrame) and other.has_data:
            # it's a frame
            self_in_other_system = self.transform_to(other)
        else:
            raise TypeError('Can only get separation to another SkyCoord or a '
                            'coordinate frame with data')

        lon1 = self_in_other_system.spherical.lon
        lat1 = self_in_other_system.spherical.lat
        lon2 = other.spherical.lon
        lat2 = other.spherical.lat

        # Get the separation as a Quantity, convert to Angle in degrees
        sep = angular_separation(lon1, lat1, lon2, lat2)
        return Angle(sep, unit=u.degree)

    def separation_3d(self, other):
        """
        Computes three dimensional separation between this coordinate
        and another.

        Parameters
        ----------
        other : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate to get the separation to.

        Returns
        -------
        sep : `~astropy.coordinates.Distance`
            The real-space distance between these two coordinates.

        Raises
        ------
        ValueError
            If this or the other coordinate do not have distances.
        """
        from . import Distance

        if isinstance(other, SkyCoord):
            self_in_other_system = self.transform_to(other.frame)
        elif isinstance(other, BaseCoordinateFrame) and other.has_data:
            # it's a frame
            self_in_other_system = self.transform_to(other)
        else:
            raise TypeError('Can only get separation to another SkyCoord or a '
                            'coordinate frame with data')

        if self.data.__class__ == UnitSphericalRepresentation:
            raise ValueError('This object does not have a distance; cannot '
                             'compute 3d separation.')
        if other.data.__class__ == UnitSphericalRepresentation:
            raise ValueError('The other object does not have a distance; '
                             'cannot compute 3d separation.')

        dx = self_in_other_system.cartesian.x - other.cartesian.x
        dy = self_in_other_system.cartesian.y - other.cartesian.y
        dz = self_in_other_system.cartesian.z - other.cartesian.z

        distval = (dx.value ** 2 + dy.value ** 2 + dz.value ** 2) ** 0.5
        return Distance(distval, dx.unit)

    def match_to_catalog_sky(self, catalogcoord, nthneighbor=1):
        """
        Finds the nearest on-sky matches of this coordinate in a set of
        catalog coordinates.

        Parameters
        ----------
        catalogcoord : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The base catalog in which to search for matches. Typically this
            will be a coordinate object that is an array (i.e.,
            ``catalogcoord.isscalar == False``)
        nthneighbor : int, optional
            Which closest neighbor to search for.  Typically ``1`` is
            desired here, as that is correct for matching one set of
            coordinates to another. The next likely use case is ``2``,
            for matching a coordinate catalog against *itself* (``1``
            is inappropriate because each point will find itself as the
            closest match).

        Returns
        -------
        idx : integer array
            Indices into ``catalogcoord`` to get the matched points for
            each of this object's coordinates. Shape matches this
            object.
        sep2d : `~astropy.coordinates.Angle`
            The on-sky separation between the closest match for each
            element in this object in ``catalogcoord``. Shape matches
            this object.
        dist3d : `~astropy.units.Quantity`
            The 3D distance between the closest match for each element
            in this object in ``catalogcoord``. Shape matches this
            object.

        Notes
        -----
        This method requires `SciPy <http://www.scipy.org>`_ to be
        installed or it will fail.

        See Also
        --------
        astropy.coordinates.match_coordinates_sky
        """
        from .matching import match_coordinates_sky

        if (isinstance(catalogcoord, (SkyCoord, BaseCoordinateFrame))
                and catalogcoord.has_data):
            self_in_catalog_frame = self.transform_to(catalogcoord)
        else:
            raise TypeError('Can only get separation to another SkyCoord or a '
                            'coordinate frame with data')

        res = match_coordinates_sky(self_in_catalog_frame, catalogcoord,
                                    nthneighbor=nthneighbor,
                                    storekdtree='_kdtree_sky')
        return res

    def match_to_catalog_3d(self, catalogcoord, nthneighbor=1):
        """
        Finds the nearest 3-dimensional matches of this coordinate to a set
        of catalog coordinates.

        This finds the 3-dimensional closest neighbor, which is only different
        from the on-sky distance if ``distance`` is set in this object or the
        ``catalogcoord`` object.

        Parameters
        ----------
        catalogcoord : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The base catalog in which to search for matches. Typically this
            will be a coordinate object that is an array (i.e.,
            ``catalogcoord.isscalar == False``)
        nthneighbor : int, optional
            Which closest neighbor to search for.  Typically ``1`` is
            desired here, as that is correct for matching one set of
            coordinates to another.  The next likely use case is
            ``2``, for matching a coordinate catalog against *itself*
            (``1`` is inappropriate because each point will find
            itself as the closest match).

        Returns
        -------
        idx : integer array
            Indices into ``catalogcoord`` to get the matched points for
            each of this object's coordinates. Shape matches this
            object.
        sep2d : `~astropy.coordinates.Angle`
            The on-sky separation between the closest match for each
            element in this object in ``catalogcoord``. Shape matches
            this object.
        dist3d : `~astropy.units.Quantity`
            The 3D distance between the closest match for each element
            in this object in ``catalogcoord``. Shape matches this
            object.
        Notes
        -----
        This method requires `SciPy <http://www.scipy.org>`_ to be
        installed or it will fail.

        See Also
        --------
        astropy.coordinates.match_coordinates_3d
        """
        from .matching import match_coordinates_3d

        if (isinstance(catalogcoord, (SkyCoord, BaseCoordinateFrame))
                and catalogcoord.has_data):
            self_in_catalog_frame = self.transform_to(catalogcoord)
        else:
            raise TypeError('Can only get separation to another SkyCoord or a '
                            'coordinate frame with data')

        res = match_coordinates_3d(self_in_catalog_frame, catalogcoord,
                                   nthneighbor=nthneighbor,
                                   storekdtree='_kdtree_3d')

        return res

    def position_angle(self, other):
        """
        Computes the on-sky position angle (East of North) between this
        `SkyCoord` and another.

        Parameters
        ----------
        other : `SkyCoord`
            The other coordinate to compute the position angle to.  It is
            treated as the "head" of the vector of the position angle.

        Returns
        -------
        pa : `~astropy.coordinates.Angle`
            The (positive) position angle of the vector pointing from ``self``
            to ``other``.  If either ``self`` or ``other`` contain arrays, this
            will be an array following the appropriate `numpy` broadcasting
            rules.

        Examples
        --------

        >>> c1 = SkyCoord(0*u.deg, 0*u.deg)
        >>> c2 = SkyCoord(1*u.deg, 0*u.deg)
        >>> c1.position_angle(c2).degree
        90.0
        >>> c3 = SkyCoord(1*u.deg, 1*u.deg)
        >>> c1.position_angle(c3).degree
        44.9956...
        """
        from . import angle_utilities

        if self.frame_name == other.frame_name:
            other_in_self_frame = other
        else:
            other_in_self_frame = other.frame.transform_to(self.frame)

        slat = self.represent_as(UnitSphericalRepresentation).lat
        slon = self.represent_as(UnitSphericalRepresentation).lon
        olat = other_in_self_frame.represent_as(UnitSphericalRepresentation).lat
        olon = other_in_self_frame.represent_as(UnitSphericalRepresentation).lon

        return angle_utilities.position_angle(slon, slat, olon, olat)

    # Name resolve
    @classmethod
    def from_name(cls, name, frame='icrs'):
        """
        Given a name, query the CDS name resolver to attempt to retrieve
        coordinate information for that object. The search database, sesame
        url, and  query timeout can be set through configuration items in
        ``astropy.coordinates.name_resolve`` -- see docstring for
        `~astropy.coordinates.get_icrs_coordinates` for more
        information.

        Parameters
        ----------
        name : str
            The name of the object to get coordinates for, e.g. ``'M42'``.
        frame : str or `BaseCoordinateFrame` class or instance
            The frame to transform the object to.

        Returns
        -------
        coord : SkyCoord
            Instance of the SkyCoord class.
        """

        from .name_resolve import get_icrs_coordinates

        icrs_coord = get_icrs_coordinates(name)
        icrs_sky_coord = cls(icrs_coord)
        if frame in ('icrs', icrs_coord.__class__):
            return icrs_sky_coord
        else:
            return icrs_sky_coord.transform_to(frame)


#<----------------Private utility functions below here------------------------->


def _get_frame_class(frame):
    """
    Get a frame class from the input `frame`, which could be a frame name
    string, or frame class.
    """
    import inspect

    if isinstance(frame, six.string_types):
        frame_names = frame_transform_graph.get_names()
        if frame not in frame_names:
            raise ValueError('Coordinate frame {0} not in allowed values {1}'
                             .format(frame, sorted(frame_names)))
        frame_cls = frame_transform_graph.lookup_name(frame)

    elif inspect.isclass(frame) and issubclass(frame, BaseCoordinateFrame):
        frame_cls = frame

    else:
        raise ValueError('Coordinate frame {0} must be a frame name or frame class'
                         .format(frame))

    return frame_cls


def _get_frame_name(args, kwargs):
    """
    Determine the coordinate frame from input SkyCoord args and kwargs.  This
    modifies args and/or kwargs in-place to remove the item that provided
    `frame`.  It also infers the frame if an input coordinate was provided and
    checks for conflicts.

    This allows for frame to be specified as a string like 'icrs' or a frame
    class like ICRS, but not an instance ICRS() since the latter could have
    non-default preferred attributes which would require a three-way merge.
    """
    frame = kwargs.pop('frame', None)

    if frame is not None:
        # Frame was provided as kwarg so validate and coerce into corresponding frame.
        frame_cls = _get_frame_class(frame)
        frame_name = frame_cls.name
    else:
        # Look for the frame in args
        for arg in args:
            try:
                frame_cls = _get_frame_class(arg)
            except ValueError:
                pass
            else:
                frame_name = frame_cls.name
                args.remove(arg)
                break
        else:
            # Not in args nor kwargs
            frame_name = None

    # Check that the new frame doesn't conflict with existing coordinate frame
    # if a coordinate is supplied in the args list.  If the frame still had not
    # been set by this point and a coordinate was supplied, then use that frame.
    for arg in args:
        coord_frame_name = None
        if isinstance(arg, BaseCoordinateFrame):
            coord_frame_name = arg.__class__.name
        elif isinstance(arg, SkyCoord):
            coord_frame_name = arg.frame_name

        if coord_frame_name is not None:
            if frame_name is None:
                frame_name = coord_frame_name
            elif frame_name != coord_frame_name:
                raise ValueError("Cannot override frame='{0}' of input coordinate with "
                                 "new frame='{1}'.  Instead transform the coordinate."
                                 .format(coord_frame_name, frame_name))

    return frame_name


def _get_units(args, kwargs):
    """
    Get the longitude unit and latitude unit from kwargs.  Possible enhancement
    is to allow input from args as well.
    """
    if 'unit' not in kwargs:
        lon_unit, lat_unit = None, None

    else:
        units = kwargs.pop('unit')

        if isinstance(units, six.string_types):
            units = [x.strip() for x in units.split(',')]
            # Allow for input like `unit='deg'`
            if len(units) == 1:
                units = (units[0], units[0])
        elif isinstance(units, Unit):
            units = (units, units)

        try:
            lon_unit, lat_unit = [Unit(x) for x in units]
        except:
            raise ValueError('Unit keyword must have one unit value or two unit values as '
                             'tuple or comma-separated string')

    return lon_unit, lat_unit


def _parse_coordinate_arg(coords, frame, lon_unit, lat_unit):
    """
    Single unnamed arg supplied.  This must be:
    - Coordinate frame with data
    - Representation
    - List or tuple of:
      - String which splits into two values
      - Iterable with two values
    """
    is_scalar = False  # Differentiate between scalar and list input
    valid_kwargs = {}  # Returned dict of lon, lat, and distance (optional)

    # Get the mapping of attribute type (e.g. 'lat', 'lon', 'distance')
    # to corresponding class attribute name ('ra', 'dec', 'distance') for frame.
    frame_cls = frame_transform_graph.lookup_name(frame)
    if frame_cls is None:
        frame_cls = NoFrame
    attr_name_for_type = dict((attr_type, name) for name, attr_type in
                              frame_cls.preferred_attr_names.items())

    # Turn a single string into a list of strings for convenience
    if isinstance(coords, six.string_types):
        is_scalar = True
        coords = [coords]

    if isinstance(coords, (SkyCoord, BaseCoordinateFrame)):
        # Note that during parsing of `frame` it is checked that any coordinate
        # args have the same frame as explicitly supplied, so don't worry here.

        if not coords.has_data:
            raise ValueError('Cannot initialize from a frame without coordinate data')

        for attr, attr_type in coords.preferred_attr_names.items():
            value = getattr(coords, attr)
            if attr_type == 'lon':
                valid_kwargs[attr] = Longitude(value, unit=lon_unit)
            elif attr_type == 'lat':
                valid_kwargs[attr] = Latitude(value, unit=lat_unit)
            elif attr_type == 'distance':
                # don't pass on a distance if no distance was initially given
                if not isinstance(coords.data, UnitSphericalRepresentation):
                    valid_kwargs[attr] = value
            else:
                raise ValueError("Unexpected attribute type '{0}'".format(attr_type))

        for attr in FRAME_ATTR_NAMES_SET():
            value = getattr(coords, attr, None)
            use_value = (isinstance(coords, SkyCoord)
                         or attr not in coords._attr_names_with_defaults)
            if use_value and value is not None:
                valid_kwargs[attr] = value

    elif isinstance(coords, BaseRepresentation):
        sph_coords = coords.represent_as(SphericalRepresentation)
        valid_kwargs[attr_name_for_type['lon']] = sph_coords.lon
        valid_kwargs[attr_name_for_type['lat']] = sph_coords.lat
        valid_kwargs[attr_name_for_type['distance']] = sph_coords.distance

    elif isinstance(coords, collections.Sequence):
        # NOTE: we do not support SkyCoord((ra, dec)).  It has to be
        # SkyCoord(ra, dec) or SkyCoord([(ra, dec)])
        lons = []
        lats = []

        for coord in coords:
            # Each row must be either a 2-element iterable or a string that splits
            # into two elements.
            if isinstance(coord, six.string_types):
                coord = coord.split()
            try:
                lon, lat = coord
            except:
                raise ValueError('Cannot parse longitude and latitude from first argument')

            lons.append(lon)
            lats.append(lat)

        if is_scalar:
            lons, lats = lons[0], lats[0]

        try:
            valid_kwargs[attr_name_for_type['lon']] = Longitude(lons, unit=lon_unit)
            valid_kwargs[attr_name_for_type['lat']] = Latitude(lats, unit=lat_unit)
        except Exception as err:
            raise ValueError('Cannot parse longitude and latitude from first argument: {0}'
                             .format(err))
    else:
        raise ValueError('Cannot parse longitude and latitude from first argument')

    return valid_kwargs


def _get_preferred_attrs(frame, lon_unit, lat_unit, kwargs):
    """
    Find instances of the "preferred attributes" for specifying longitude,
    latitude, and distance for this frame.  Pop them off of kwargs, run
    through the appropriate class constructor (to validate and apply unit), and
    put into the output valid_kwargs.
    """
    valid_kwargs = {}
    frame_cls = frame_transform_graph.lookup_name(frame)
    if frame_cls is None:
        frame_cls = NoFrame
    for attr_cls, attr_type, unit in ((Longitude, 'lon', lon_unit),
                                      (Latitude, 'lat', lat_unit),
                                      (None, 'distance', None)):
        attr_name_for_type = dict((attr_type, name) for name, attr_type in
                                  frame_cls.preferred_attr_names.items())
        name = attr_name_for_type[attr_type]
        if name in kwargs:
            value = kwargs.pop(name)
            valid_kwargs[name] = attr_cls(value, unit=unit) if attr_cls else value

    return valid_kwargs
