from __future__ import (absolute_import, division, print_function, unicode_literals)

from copy import deepcopy
import collections

from ..utils.compat.misc import override__dir__
from ..extern import six
from ..units import Unit

from .angles import Latitude, Longitude
from .distances import Distance
from .baseframe import BaseCoordinateFrame, frame_transform_graph

__all__ = ['SkyCoord']

# Create FRAME_CLASSES dictionary mapping of system name: frame class.
FRAME_NAMES = frame_transform_graph.get_aliases()
FRAME_CLASSES = dict((name, frame_transform_graph.lookup_name(name))
                     for name in FRAME_NAMES)
FRAME_CLASSES[None] = FRAME_CLASSES['icrs']

# Inverse mapping is useful
CLASS_TO_NAME_MAP = dict((cls, name) for name, cls in FRAME_CLASSES.items())


# Map coordinate system names to allowed frame-specific attributes such as
# equinox or obstime.
FRAME_ATTR_NAMES = dict((name, tuple(FRAME_CLASSES[name].frame_attr_names.keys()))
                        for name in FRAME_NAMES)


class SkyCoord(object):

    def __init__(self, *args, **kwargs):
        # Parse the args and kwargs to assemble a sanitized and validated
        # kwargs dict for initializing attributes for this object and for
        # creating the internal self._coords object
        args = list(args)  # Make it mutable
        kwargs = self._parse_inputs(args, kwargs)

        # Set internal versions of object state attributes
        for attr in ('system', 'equinox', 'obstime', 'location'):
            setattr(self, '_' + attr, kwargs[attr])

        # Set up the keyword args for creating the internal coordinate object.
        frame_cls = FRAME_CLASSES[self.system]
        coord_kwargs = {}
        for attr, value in kwargs.items():
            if (attr in frame_cls.preferred_attr_names
                    or attr in frame_cls.frame_attr_names):
                coord_kwargs[attr] = value

        # Finally make the internal coordinate object.
        self._coord = frame_cls(**coord_kwargs)

    @property
    def system(self):
        return self._system

    @property
    def equinox(self):
        return self._equinox

    @property
    def location(self):
        return self._location

    @property
    def obstime(self):
        return self._obstime

    def __len__(self):
        return len(self._coord)

    def _parse_inputs(self, args, kwargs):
        """
        Assemble a validated and sanitized keyword args dict for instantiating a
        SkyCoord and coordinate object from the provided `args`, and `kwargs`.
        """
        valid_kwargs = {}

        # Put the SkyCoord attributes system, equinox, obstime, location into valid_kwargs
        # dict.  `System` could come from args or kwargs, so set valid_kwargs['system']
        # accordingly.  The others must be specified by keyword args.  Pop them off
        # of kwargs in the process.
        system = valid_kwargs['system'] = _get_system(args, kwargs)
        for attr in ('equinox', 'obstime', 'location'):
            valid_kwargs[attr] = kwargs.pop(attr, None)

        # Get latitude and longitude units
        lon_unit, lat_unit = _get_units(args, kwargs)

        # Grab any frame-specific attr names like `ra` or `l` or `distance`
        valid_kwargs.update(_get_preferred_attrs(system, lon_unit, lat_unit, kwargs))

        # Error if anything is still left in kwargs
        if kwargs:
            raise ValueError('Unrecognized keyword argument(s) {0}'
                             .format(', '.join("'{0}'".format(key) for key in kwargs)))

        # Finally deal with the unnamed args.  This figures out what the arg[0] is
        # and returns a dict with appropriate key/values for initializing frame class.
        if args:
            if len(args) == 1:
                # One arg which must be a coordinate
                coord_kwargs = _parse_coordinate_arg(args[0], system, lon_unit, lat_unit)
            elif len(args) == 2:
                raise NotImplementedError()
            else:
                raise ValueError()

            for attr in coord_kwargs:
                if attr in valid_kwargs:
                    raise ValueError("Cannot supply a '{0}' keyword along with a coordinate input"
                                     .format(attr))
                valid_kwargs[attr] = coord_kwargs[attr]

        return valid_kwargs

    def transform_to(self, system):
        """
        Transform this coordinate to a new system.

        Parameters
        ----------
        system : str
            The system to transform this coordinate into.

        Returns
        -------
        coord
            A new object with this coordinate represented in the `system` system.

        Raises
        ------
        ValueError
            If there is no possible transformation route.
        """
        from astropy.coordinates.errors import ConvertError
        import inspect

        if self.system is None:
            raise ValueError('Cannot transform coordinates if `system` is None')

        if isinstance(system, six.string_types):
            if system not in FRAME_CLASSES:
                raise ValueError('Coordinate system {0} not in allowed values {1}'
                                 .format(system, sorted(FRAME_CLASSES)))
            new_frame = FRAME_CLASSES[system]
        else:
            # Allow for a frame class or a frame class instance
            new_frame = system
            coord_cls = new_frame if inspect.isclass(new_frame) else new_frame.__class__
            if (not issubclass(coord_cls, BaseCoordinateFrame)
                    or coord_cls not in CLASS_TO_NAME_MAP):
                raise ValueError('Coordinate system {0} must be a frame class or '
                                 'frame class instance'.format(new_frame))
            system = CLASS_TO_NAME_MAP[coord_cls]

        out = deepcopy(self)
        out._system = system
        out._coord = self._coord.transform_to(new_frame)
        if out._coord is None:
            raise ConvertError('Cannot transform from {0} to '
                               '{1}'.format(self.system, system))

        return out

    def __getattr__(self, attr):
        """
        Overrides getattr to return coordinates that this can be transformed
        to, based on the alias attr in the master transform graph.
        """

        if self.system == attr:
            return self

        nmsys = frame_transform_graph.lookup_name(attr)
        if nmsys is not None and self._coord.is_transformable_to(nmsys):
            return self.transform_to(attr)
        else:
            if not attr.startswith('_') and hasattr(self._coord, attr):
                return getattr(self._coord, attr)
            else:
                msg = "'{0}' object has no attribute '{1}', nor a transform."
                raise AttributeError(msg.format(self.__class__, attr))

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

        return dir_values

    def __repr__(self):
        return '<{0} {1}'.format(self.__class__.__name__, repr(self._coord)[1:])


def _get_system(args, kwargs):
    """
    Determine the coordinate system from input SkyCoord args and kwargs.  This
    modifies args in-place to remove the item that provided `system`.  It
    also infers the system if an input coordinate was provided and checks
    for conflicts.

    """
    system = kwargs.pop('system', None)

    # If no system is provided via keyword then check for a positional arg
    # that is a valid system.
    if system is None:
        for arg in args:
            if arg in FRAME_NAMES:
                system = arg
                args.remove(system)
                break

    elif system not in FRAME_NAMES:
        # System was provided so check that it is allowed.
        raise ValueError('Coordinate system {0} not in allowed values {1}'
                         .format(system, sorted(FRAME_NAMES)))

    # Check that the new system doesn't conflict with existing coordinate system
    # if a coordinate is supplied in the args list.  If the system still had not
    # been set by this point and a coordinate was supplied, then use that system.
    for arg in args:
        coord_system = None
        if isinstance(arg, BaseCoordinateFrame):
            coord_system = CLASS_TO_NAME_MAP[arg.__class__]
        elif isinstance(arg, SkyCoord):
            coord_system = arg.system

        if coord_system is not None:
            if system is None:
                system = coord_system
            elif system != coord_system:
                raise ValueError("Cannot override system='{0}' of input coordinate with "
                                 "new system='{1}'.  Instead transform the coordinate."
                                 .format(coord_system, system))
    return system


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


def _parse_coordinate_arg(coords, system, lon_unit, lat_unit):
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
    frame_cls = FRAME_CLASSES[system]
    attr_name_for_type = dict((attr_type, name) for name, attr_type in
                              frame_cls.preferred_attr_names.items())

    # Turn a single string into a list of strings for convenience
    if isinstance(coords, six.string_types):
        is_scalar = True
        coords = [coords]

    if isinstance(coords, (SkyCoord, BaseCoordinateFrame)):
        # Note that during parsing of `system` it is checked that any coordinate
        # args have the same system as explicitly supplied, so don't worry here.

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


def _get_preferred_attrs(system, lon_unit, lat_unit, kwargs):
    """
    Find instances of the "preferred attributes" for specifying longitude,
    latitude, and distance for this system.  Pop them off of kwargs, run
    through the appropriate class constructor (to validate and apply unit), and
    put into the output valid_kwargs.
    """
    valid_kwargs = {}
    frame_cls = FRAME_CLASSES[system]
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
