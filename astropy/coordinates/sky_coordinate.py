from __future__ import (absolute_import, division, print_function, unicode_literals)

from copy import deepcopy
import collections

from ..utils.compat.misc import override__dir__
from ..extern import six

from .angles import Latitude, Longitude
from .baseframe import BaseCoordinateFrame, frame_transform_graph

__all__ = ['SkyCoord']

# Create FRAME_CLASSES dictionary mapping of system name: frame class.
FRAME_NAMES = frame_transform_graph.get_aliases()
FRAME_CLASSES = dict((name, frame_transform_graph.lookup_name(name))
                     for name in FRAME_NAMES)
FRAME_CLASSES[None] = FRAME_CLASSES['icrs']

# Map coordinate system names to allowed frame-specific attributes such as
# equinox or obstime.
FRAME_ATTR_NAMES = dict((name, tuple(FRAME_CLASSES[name].frame_attr_names.keys()))
                        for name in FRAME_NAMES)

# Get the preferred representation of longitude and latitude for each coord system.
PREFERRED_REPR_LON_LAT = {}
for name, frame_cls in FRAME_CLASSES.items():
    _rev_map = dict((val, key) for key, val in frame_cls.preferred_attr_names.items())
    PREFERRED_REPR_LON_LAT[name] = (_rev_map['lon'], _rev_map['lat'])


class SkyCoord(object):
    _sky_coordinate_attrs = ('system', 'equinox', 'obstime', 'location')

    def __init__(self, *args, **kwargs):
        # *args, **kwargs needed for desired flexibility in inputs

        # Get the coordinate system name from inputs
        system = self._get_system(args, kwargs)

        # Set self attributes from kwargs.  If the attr is an input to the
        # coordinate class for the `system` then leave it in kwargs, otherwise
        # pop it off.
        for attr in self._sky_coordinate_attrs:
            if attr in FRAME_ATTR_NAMES[system]:
                val = kwargs.get(attr)
            else:
                val = kwargs.pop(attr, None)
            setattr(self, attr, val)

        # Pull out lon, lat, and (optionally) distance from args and kwargs
        coord_kwargs = _get_coordinate_inputs(args, kwargs)

        # If found via `_get_coordinate_inputs()`, make new *args for creating coordinate
        args = [coord_kwargs.pop(axis) for axis in ('lon', 'lat') if axis in coord_kwargs]

        kwargs.update(coord_kwargs)  # Set `distance` (possibly)

        # Finally make the internal coordinate object.  This delegates much of
        # the input validation to the low-level coordinate class, including
        # potential conflicts like supplying an `ra` keyword along with an
        # inital coordinate arg, or missing `unit`, etc.
        self._coord = FRAME_CLASSES[system](*args, **kwargs)

    def _get_system(self, args, kwargs):
        """
        Determine the coordinate system from input args and kwargs.  This modifies
        args or kwargs in-place to remove the item that provided `system`.
        """
        system = kwargs.pop('system', None)

        if system is None:
            for arg in args:
                if arg in FRAME_NAMES:
                    system = arg
                    args.remove(system)
                    break

        if system not in FRAME_NAMES:
            raise ValueError('Coordinate system {0} not in allowed values {1}'
                             .format(system, sorted(FRAME_NAMES)))
        return system

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

        if system is None:
            raise ValueError('Cannot transform coordinates if `system` is None')

        if system not in FRAME_CLASSES:
            raise ValueError('Coordinate system {0} not in allowed values {1}'
                             .format(system, sorted(FRAME_CLASSES)))

        out = deepcopy(self)
        if system == self.system:
            return out

        out.system = system
        out._coord = self._coord.transform_to(FRAME_CLASSES[system])
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


def _get_coordinate_inputs(args, kwargs):
    """
    Figure out lat, lon, and distance from input args and kwargs.
    This updates `kwargs` in-place.
    """
    out = {}

    # Parse specified `unit` into a tuple and supply default
    unit = kwargs.get('unit', (None, None))
    if isinstance(unit, six.string_types):
        unit = unit.split(',')
        # Allow for input like `unit='deg'`
        if len(unit) == 1:
            unit = (unit, unit)
    try:
        lon_unit, lat_unit = unit
    except:
        raise ValueError('Unit keyword must have two values as tuple or '
                         'comma-separated string')

    # Convert any recognized kwargs like `ra` or `l` into the right Angle subclass
    for angle_class, attrs_index, angle_unit in ((Longitude, 0, lon_unit),
                                                 (Latitude, 1, lat_unit)):
        axis_attrs = set(lon_lat[attrs_index] for lon_lat in PREFERRED_REPR_LON_LAT.values())
        for key in kwargs:
            if key in axis_attrs:
                kwargs[key] = angle_class(kwargs[key], unit=angle_unit)

    # Finally deal with the unnamed args
    if len(args) == 1:
        out = _parse_one_arg(*args, **kwargs)
    elif len(args) == 2:
        out['lon'] = args[0]
        out['lat'] = args[1]
    return out


def _parse_one_arg(*args, **kwargs):
    """
    Single unnamed arg supplied.  This must be:
    - Coordinate
    - List or tuple of:
      - String which splits into two values
      - Iterable with two values
    """
    coords = args[0]
    is_scalar = False  # Differentiate between scalar and list input
    out = {}  # Returned dict of lon, lat, and distance (optional)

    # Turn a single string into a list of strings for convenience
    if isinstance(coords, six.string_types):
        is_scalar = True
        coords = [coords]

    if isinstance(coords, (SkyCoord, BaseCoordinateFrame)):
        out['lon'] = coords.lonangle
        out['lat'] = coords.latangle
        if coords.distance is not None:
            if 'distance' in kwargs:
                raise ValueError('Cannot supply a `distance` keyword that overrides existing '
                                 'coordinate distance')
            out['distance'] = coords.distance
        if 'unit' in kwargs and kwargs['unit'] != (None, None):
            raise ValueError('Cannot supply a `unit` keyword along with a coordinate input')

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

        out['lon'] = Longitude(lons, unit=kwargs['unit'][0])
        out['lat'] = Latitude(lats, unit=kwargs['unit'][1])
    else:
        raise ValueError('Cannot parse lon, lat from first argument')

    return out
