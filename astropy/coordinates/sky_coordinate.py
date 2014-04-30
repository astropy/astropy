from __future__ import (absolute_import, division, print_function, unicode_literals)

import re
from copy import deepcopy
import collections

from ..utils.compat.misc import override__dir__
from ..extern import six
from ..units import Unit

from .angles import Latitude, Longitude
from .baseframe import BaseCoordinateFrame, frame_transform_graph

__all__ = ['SkyCoord']

# ALL THESE FIXED MAPPINGS NEED TO BE DYNAMIC

# Create FRAME_CLASSES dictionary mapping of frame name: frame class.
FRAME_NAMES = frame_transform_graph.get_aliases()
FRAME_CLASSES = dict((name, frame_transform_graph.lookup_name(name))
                     for name in FRAME_NAMES)
FRAME_CLASSES[None] = FRAME_CLASSES['icrs']

# Inverse mapping is useful
CLASS_TO_NAME_MAP = dict((cls, name) for name, cls in FRAME_CLASSES.items())


# Set of all possible frame-specific attributes
FRAME_ATTR_NAMES_SET = set()
for frame_cls in FRAME_CLASSES.values():
    for attr in frame_cls.frame_attr_names.keys():
        FRAME_ATTR_NAMES_SET.add(attr)


class SkyCoord(object):

    def __init__(self, *args, **kwargs):
        # Parse the args and kwargs to assemble a sanitized and validated
        # kwargs dict for initializing attributes for this object and for
        # creating the internal self._coords object
        args = list(args)  # Make it mutable
        kwargs = self._parse_inputs(args, kwargs)

        # Set internal versions of object state attributes
        self._frame = kwargs['frame']
        for attr in FRAME_ATTR_NAMES_SET:
            setattr(self, '_' + attr, kwargs[attr])

        # Set up the keyword args for creating the internal coordinate object.
        frame_cls = FRAME_CLASSES[self.frame]
        coord_kwargs = {}
        for attr, value in kwargs.items():
            if value is not None and (attr in frame_cls.preferred_attr_names
                                      or attr in frame_cls.frame_attr_names):
                coord_kwargs[attr] = value

        # Finally make the internal coordinate object.
        self._coord = frame_cls(**coord_kwargs)

    @property
    def frame(self):
        return self._frame

    def __len__(self):
        return len(self._coord)

    def _parse_inputs(self, args, kwargs):
        """
        Assemble a validated and sanitized keyword args dict for instantiating a
        SkyCoord and coordinate object from the provided `args`, and `kwargs`.
        """
        valid_kwargs = {}

        # Put the SkyCoord attributes frame, equinox, obstime, location into valid_kwargs
        # dict.  `Frame` could come from args or kwargs, so set valid_kwargs['frame']
        # accordingly.  The others must be specified by keyword args.  Pop them off
        # of kwargs in the process.
        frame = valid_kwargs['frame'] = _get_frame_name(args, kwargs)
        for attr in FRAME_ATTR_NAMES_SET:
            valid_kwargs[attr] = kwargs.pop(attr, None)

        # Get latitude and longitude units
        lon_unit, lat_unit = _get_units(args, kwargs)

        # Grab any frame-specific attr names like `ra` or `l` or `distance`
        valid_kwargs.update(_get_preferred_attrs(frame, lon_unit, lat_unit, kwargs))

        # Error if anything is still left in kwargs
        if kwargs:
            raise ValueError('Unrecognized keyword argument(s) {0}'
                             .format(', '.join("'{0}'".format(key) for key in kwargs)))

        # Finally deal with the unnamed args.  This figures out what the arg[0] is
        # and returns a dict with appropriate key/values for initializing frame class.
        if args:
            if len(args) == 1:
                # One arg which must be a coordinate
                coord_kwargs = _parse_coordinate_arg(args[0], frame, lon_unit, lat_unit)

            elif len(args) == 2:
                attr_name_for_type = dict((attr_type, name) for name, attr_type in
                                          FRAME_CLASSES[frame].preferred_attr_names.items())
                coord_kwargs = {}
                coord_kwargs[attr_name_for_type['lon']] = Longitude(args[0], unit=lon_unit)
                coord_kwargs[attr_name_for_type['lat']] = Latitude(args[1], unit=lat_unit)
            else:
                raise ValueError('Must supply no more than two positional arguments, got {}'
                                 .format(len(args)))

            for attr in coord_kwargs:
                if attr in valid_kwargs:
                    raise ValueError("Cannot supply a '{0}' keyword along with a coordinate input"
                                     .format(attr))
                valid_kwargs[attr] = coord_kwargs[attr]

        return valid_kwargs

    def transform_to(self, frame):
        """
        Transform this coordinate to a new frame.

        Parameters
        ----------
        frame : str
            The frame to transform this coordinate into.

        Returns
        -------
        coord
            A new object with this coordinate represented in the `frame` frame.

        Raises
        ------
        ValueError
            If there is no possible transformation route.
        """
        from astropy.coordinates.errors import ConvertError

        if frame is None or self.frame is None:
            raise ValueError('Cannot transform coordinates to/from `frame=None`')

        if frame not in FRAME_CLASSES:
            raise ValueError("Output coordinate frame '{0}' not in allowed values {1}"
                             .format(frame, sorted(FRAME_CLASSES)))

        # Make a frame instance for the new frame including the correct frame attributes.
        # Only define attributes that are not None.
        frame_cls = FRAME_CLASSES[frame]
        frame_kwargs = dict((attr, getattr(self, attr)) for attr in frame_cls.frame_attr_names
                            if getattr(self, attr) is not None)
        frame = frame_cls(**frame_kwargs)

        out = deepcopy(self)
        out._frame = frame
        out._coord = self._coord.transform_to(frame)
        if out._coord is None:
            raise ConvertError('Cannot transform from {0} to '
                               '{1}'.format(self.frame, frame))

        return out

    def __getattr__(self, attr):
        """
        Overrides getattr to return coordinates that this can be transformed
        to, based on the alias attr in the master transform graph.
        """

        if self.frame == attr:
            return self  # Should this be a deepcopy of self?

        # Anything in the set of all possible frame_attr_names is handled
        # here. If the attr is relevant for the current frame then delegate
        # to self._coord otherwise get it from self._<attr>.
        if attr in FRAME_ATTR_NAMES_SET:
            if attr in self._coord.frame_attr_names:
                return getattr(self._coord, attr)
            else:
                return getattr(self, '_' + attr)

        # Some attributes might not fall in the above category but still
        # are available through self._coord.
        if not attr.startswith('_') and hasattr(self._coord, attr):
            return getattr(self._coord, attr)

        # Try to interpret as a new frame for transforming.
        frame_cls = frame_transform_graph.lookup_name(attr)
        if frame_cls is not None and self._coord.is_transformable_to(frame_cls):
            return self.transform_to(attr)

        # Fail
        raise AttributeError("'{0}' object has no attribute '{1}', nor a transform."
                             .format(self.__class__, attr))

    @override__dir__
    def __dir__(self):
        """
        Override the builtin `dir` behavior to include:
        - Transforms available by aliases
        - Attribute / methods of the underlying self._coord object
        """

        # determine the aliases that this can be transformed to.
        dir_values = set()
        for name, frame_cls in FRAME_CLASSES.items():
            if self._coord.is_transformable_to(frame_cls):
                dir_values.add(name)

        # Add public attributes of self._coord
        dir_values.update(set(attr for attr in dir(self._coord) if not attr.startswith('_')))

        # Add all possible frame_attr_names
        dir_values.update(FRAME_ATTR_NAMES_SET)

        return dir_values

    def __repr__(self):
        out = re.sub('Coordinate:', self.__class__.__name__ + ':', repr(self._coord))
        return out


def _get_frame_class(frame):
    """
    Get a frame instance from the input `frame`, which could be a frame name
    string, frame class, or frame instance.
    """
    import inspect

    if isinstance(frame, six.string_types):
        if frame not in FRAME_NAMES:
            raise ValueError('Coordinate frame {0} not in allowed values {1}'
                             .format(frame, sorted(FRAME_NAMES)))
        frame_cls = FRAME_CLASSES[frame]

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
        frame_name = CLASS_TO_NAME_MAP[frame_cls]
    else:
        # Look for the frame in args
        for arg in args:
            try:
                frame_cls = _get_frame_class(arg)
            except ValueError:
                pass
            else:
                frame_name = CLASS_TO_NAME_MAP[frame_cls]
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
            coord_frame_name = CLASS_TO_NAME_MAP[arg.__class__]
        elif isinstance(arg, SkyCoord):
            coord_frame_name = arg.frame

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

        try:
            lon_unit, lat_unit = [Unit(x) for x in units]
        except:
            raise ValueError('Unit keyword must have two unit values as tuple or '
                             'comma-separated string')

    return lon_unit, lat_unit


def _parse_coordinate_arg(coords, frame, lon_unit, lat_unit):
    """
    Single unnamed arg supplied.  This must be:
    - Coordinate
    - List or tuple of:
      - String which splits into two values
      - Iterable with two values
    """
    is_scalar = False  # Differentiate between scalar and list input
    valid_kwargs = {}  # Returned dict of lon, lat, and distance (optional)

    # Get the mapping of attribute type (e.g. 'lat', 'lon', 'distance')
    # to corresponding class attribute name ('ra', 'dec', 'distance') for frame.
    frame_cls = FRAME_CLASSES[frame]
    attr_name_for_type = dict((attr_type, name) for name, attr_type in
                              frame_cls.preferred_attr_names.items())

    # Turn a single string into a list of strings for convenience
    if isinstance(coords, six.string_types):
        is_scalar = True
        coords = [coords]

    if isinstance(coords, (SkyCoord, BaseCoordinateFrame)):
        # Note that during parsing of `frame` it is checked that any coordinate
        # args have the same frame as explicitly supplied, so don't worry here.

        for attr, attr_type in coords.preferred_attr_names.items():
            value = getattr(coords, attr)
            if attr_type == 'lon':
                valid_kwargs[attr] = Longitude(value, unit=lon_unit)
            elif attr_type == 'lat':
                valid_kwargs[attr] = Latitude(value, unit=lat_unit)
            elif attr_type == 'distance':
                valid_kwargs[attr] = value
            else:
                raise ValueError("Unexpected attribute type '{0}'".format(attr_type))

        for attr in FRAME_ATTR_NAMES_SET:
            value = getattr(coords, attr, None)
            use_value = (isinstance(coords, SkyCoord)
                         or attr not in coords._attr_names_with_defaults)
            if use_value and value is not None:
                valid_kwargs[attr] = value

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
    frame_cls = FRAME_CLASSES[frame]
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
