from __future__ import (absolute_import, division, print_function, unicode_literals)

from copy import deepcopy
import collections

from .. import units as u
from ..units import Quantity
from ..utils.compat.misc import override__dir__
from ..extern import six

from ..coordinates import Angle, SphericalCoordinatesBase, Latitude, Longitude
from .transformations import master_transform_graph

COORD_CLASSES = deepcopy(master_transform_graph._clsaliases)
COORD_CLASSES[None] = COORD_CLASSES['icrs']

# Stubs until coordinate-refactor is available
FRAME_ATTR_NAMES = {None: ('obstime'),
                    'icrs': ('obstime'),
                    'fk5': ('obstime', 'equinox'),
                    'fk4': ('obstime', 'equinox'),
                    'galactic': ('obstime'),
                    'altaz': ('obstime', 'equinox')}

PREFERRED_REPR_LON_LAT = {None: ('ra', 'dec'),
                          'icrs': ('ra', 'dec'),
                          'fk5': ('ra', 'dec'),
                          'fk4': ('ra', 'dec'),
                          'galactic': ('l', 'b'),
                          'altaz': ('az', 'alt')}


__all__ = ['SkyCoordinate']

class SkyCoordinate(object):
    _sky_coordinate_attrs = ('system', 'equinox', 'obstime', 'location')

    def __init__(self, *args, **kwargs):
        # *args, **kwargs needed for desired flexibility in inputs

        system = kwargs.get('system')

        if system not in COORD_CLASSES:
            raise ValueError('Coordinate system {0} not in allowed values {1}'
                             .format(system, sorted(COORD_CLASSES)))

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
        coord_kwargs = _get_coordinate_inputs(*args, **kwargs)

        # If found via `_get_coordinate_inputs()`, make new *args for creating coordinate
        args = [coord_kwargs.pop(axis) for axis in ('lon', 'lat') if axis in coord_kwargs]

        kwargs.update(coord_kwargs)  # Set `distance` (possibly)

        # Finally make the internal coordinate object.  This delegates much of
        # the input validation to the low-level coordinate class, including
        # potential conflicts like supplying an `ra` keyword along with an
        # inital coordinate arg, or missing `unit`, etc.
        self._coord = COORD_CLASSES[system](*args, **kwargs)

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

        if system not in COORD_CLASSES:
            raise ValueError('Coordinate system {0} not in allowed values {1}'
                             .format(system, sorted(COORD_CLASSES)))

        out = deepcopy(self)
        if system == self.system:
            return out

        out.system = system
        out._coord = self._coord.transform_to(COORD_CLASSES[system])
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

        nmsys = master_transform_graph.lookup_name(attr)
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
        for alias in master_transform_graph.get_aliases():
            tosys = master_transform_graph.lookup_name(alias)
            if self._coord.is_transformable_to(tosys):
                dir_values.add(alias)

        # Add public attributes of self._coord
        dir_values.update(set(attr for attr in dir(self._coord) if not attr.startswith('_')))

        return dir_values

    def __repr__(self):
        return '<{0} {1}'.format(self.__class__.__name__, repr(self._coord)[1:])
        

def _get_coordinate_inputs(*args, **kwargs):
    """
    Figure out lat, lon, and distance from input args and kwargs.

    """
    out = {}

    # One unnamed arg
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
        
    if isinstance(coords, (SkyCoordinate, SphericalCoordinatesBase)):
        out['lon'] = coords.lonangle
        out['lat'] = coords.latangle
        if coords.distance is not None:
            if 'distance' in kwargs:
                raise ValueError('Cannot supply a `distance` keyword that overrides existing '
                                 'coordinate distance')
            out['distance'] = coords.distance
        if 'unit' in kwargs:
            raise ValueError('Cannot supply a `unit` keyword along with a coordinate input')

    elif isinstance(coords, collections.Sequence):
        # NOTE: we do not support SkyCoordinate((ra, dec)).  It has to be 
        # SkyCoordinate(ra, dec) or SkyCoordinate([(ra, dec)])
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

        unit = kwargs.get('unit', (None, None))
        if isinstance(unit, six.string_types):
            unit = unit.split(',')
        try:
            lon_unit, lat_unit = unit
        except:
            raise ValueError('Unit keyword must have two values as tuple or '
                             'comma-separated string')

        if is_scalar:
            lons, lats = lons[0], lats[0]

        out['lon'] = Longitude(lons, unit=lon_unit)
        out['lat'] = Latitude(lats, unit=lat_unit)
    else:
        raise ValueError('Cannot parse lon, lat from first argument')

    return out
