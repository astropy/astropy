from __future__ import absolute_import

import base64
import numpy as np

from ..time import Time, TimeDelta
from .. import units as u
from .. import coordinates as coords
from ..coordinates.sky_coordinate import FRAME_ATTR_NAMES_SET
from .data_info import _get_obj_attrs_map

try:
    import yaml
except ImportError:
    raise ImportError('`import yaml` failed, PyYAML package is required for YAML')


__all__ = ['AstropyLoader', 'AstropyDumper', 'load', 'load_all', 'dump']

def _unit_representer(dumper, obj):
    out = {'name': str(obj.to_string())}
    return dumper.represent_mapping(u'!astropy.units.Unit', out)


def _unit_constructor(loader, node):
    map = loader.construct_mapping(node)
    return u.Unit(map['name'])


def _time_representer(dumper, obj):
    attrs = ('jd1', 'jd2', 'format', 'scale', 'precision', 'in_subfmt',
             'out_subfmt', 'location', '_delta_ut1_utc', '_delta_tdb_tt')
    out = _get_obj_attrs_map(obj, attrs)
    return dumper.represent_mapping(u'!astropy.time.Time', out)


def _time_constructor(loader, node):
    map = loader.construct_mapping(node)
    format = map.pop('format')
    delta_ut1_utc = map.pop('_delta_ut1_utc', None)
    delta_tdb_tt = map.pop('_delta_tdb_tt', None)

    map['format'] = 'jd'
    map['val'] = map.pop('jd1')
    map['val2'] = map.pop('jd2')

    out = Time(**map)
    out.format = format

    if delta_ut1_utc is not None:
        out._delta_ut1_utc = delta_ut1_utc
    if delta_tdb_tt is not None:
        out._delta_tdb_tt = delta_tdb_tt

    return out


def _timedelta_representer(dumper, obj):
    attrs = ('jd1', 'jd2', 'format', 'scale')
    out = _get_obj_attrs_map(obj, attrs)
    return dumper.represent_mapping(u'!astropy.time.TimeDelta', out)


def _timedelta_constructor(loader, node):
    map = loader.construct_mapping(node)
    format = map.pop('format')

    map['format'] = 'jd'
    map['val'] = map.pop('jd1')
    map['val2'] = map.pop('jd2')

    out = TimeDelta(**map)
    out.format = format

    return out


def _ndarray_representer(dumper, obj):
    if obj.flags['C_CONTIGUOUS']:
        obj_data = obj.data
    else:
        cont_obj = np.ascontiguousarray(obj)
        assert(cont_obj.flags['C_CONTIGUOUS'])
        obj_data = cont_obj.data
    data_b64 = base64.b64encode(bytes(obj_data))
    out = dict(__ndarray__=data_b64,
               dtype=str(obj.dtype),
               shape=obj.shape)

    return dumper.represent_mapping(u'!numpy.ndarray', out)


def _ndarray_constructor(loader, node):
    map = loader.construct_mapping(node)
    data = base64.b64decode(map['__ndarray__'])
    return np.frombuffer(data, map['dtype']).reshape(map['shape'])


# Define supported Quantity subclasses
QUANTITY_CLASSES = {cls.__name__: cls for cls in
                    (u.Quantity, coords.Angle, coords.Longitude, coords.Latitude)}

def _quantity_representer(dumper, obj):
    out = {'class': obj.__class__.__name__,
           'value': obj.value,
           'unit': obj.unit}
    if out['class'] not in QUANTITY_CLASSES:
        raise TypeError('cannot represent quantity subclass {}'
                        .format(out['class']))
    return dumper.represent_mapping(u'!astropy.units.Quantity', out)


def _quantity_constructor(loader, node):
    map = loader.construct_mapping(node)
    cls = map.pop('class')
    value = map.pop('value')
    return QUANTITY_CLASSES[cls](value, **map)


def _earthlocation_representer(dumper, obj):
    out = _get_obj_attrs_map(obj, ('x', 'y', 'z', 'ellipsoid'))
    return dumper.represent_mapping(u'!astropy.coordinates.earth.EarthLocation', out)


def _earthlocation_constructor(loader, node):
    map = loader.construct_mapping(node)
    ellipsoid = map.pop('ellipsoid')
    out = coords.EarthLocation(**map)
    out.ellipsoid = ellipsoid
    return out


def _skycoord_representer(dumper, obj):
    attrs = list(obj.representation_component_names)
    attrs += list(FRAME_ATTR_NAMES_SET())
    out = _get_obj_attrs_map(obj, attrs)

    # Don't output distance if it is all unitless 1.0
    if 'distance' in out and np.all(out['distance'] == 1.0):
        del out['distance']

    out['representation'] = obj.representation.get_name()
    out['frame'] = obj.frame.name

    return dumper.represent_mapping(u'!astropy.coordinates.sky_coordinate.SkyCoord',
                                    out)

def _skycoord_constructor(loader, node):
    map = loader.construct_mapping(node)
    out = coords.SkyCoord(**map)
    return out


class AstropyLoader(yaml.SafeLoader):
    """
    Custom Loader that constructs OrderedDict from an !!omap object.
    This does nothing but provide a namespace for adding the
    custom odict constructor.
    """
    def construct_python_tuple(self, node):
        return tuple(self.construct_sequence(node))

    def construct_python_unicode(self, node):
        return self.construct_scalar(node)

class AstropyDumper(yaml.SafeDumper):
    def represent_tuple(self, data):
        return self.represent_sequence(u'tag:yaml.org,2002:python/tuple', data)

AstropyDumper.add_representer(u.IrreducibleUnit, _unit_representer)
AstropyDumper.add_representer(u.CompositeUnit, _unit_representer)
AstropyDumper.add_multi_representer(u.Unit, _unit_representer)
AstropyDumper.add_representer(tuple, AstropyDumper.represent_tuple)
AstropyDumper.add_representer(np.ndarray, _ndarray_representer)
AstropyDumper.add_multi_representer(u.Quantity, _quantity_representer)
AstropyDumper.add_representer(Time, _time_representer)
AstropyDumper.add_representer(TimeDelta, _timedelta_representer)
AstropyDumper.add_representer(coords.EarthLocation, _earthlocation_representer)
AstropyDumper.add_representer(coords.SkyCoord, _skycoord_representer)

AstropyLoader.add_constructor(u'tag:yaml.org,2002:python/tuple',
                              AstropyLoader.construct_python_tuple)
AstropyLoader.add_constructor(u'tag:yaml.org,2002:python/unicode',
                              AstropyLoader.construct_python_unicode)
AstropyLoader.add_constructor('!astropy.units.Unit', _unit_constructor)
AstropyLoader.add_constructor('!numpy.ndarray', _ndarray_constructor)
AstropyLoader.add_constructor('!astropy.units.Quantity', _quantity_constructor)
AstropyLoader.add_constructor('!astropy.time.Time', _time_constructor)
AstropyLoader.add_constructor('!astropy.time.TimeDelta', _timedelta_constructor)
AstropyLoader.add_constructor('!astropy.coordinates.earth.EarthLocation',
                              _earthlocation_constructor)
AstropyLoader.add_constructor('!astropy.coordinates.sky_coordinate.SkyCoord',
                              _skycoord_constructor)


def load(stream):
    """Parse the first YAML document in a stream using the AstropyLoader and
    produce the corresponding Python object.

    Parameters
    ----------
    stream : str or file-like object
        YAML input

    Returns
    -------
    obj : object
        Object corresponding to YAML document
    """
    return yaml.load(stream, Loader=AstropyLoader)


def load_all(stream):
    """Parse the all YAML documents in a stream using the AstropyLoader class and
    produce the corresponding Python object.

    Parameters
    ----------
    stream : str or file-like object
        YAML input

    Returns
    -------
    obj : object
        Object corresponding to YAML document

    """
    return yaml.load_all(stream, Loader=AstropyLoader)


def dump(data, stream=None, **kwargs):
    """Serialize a Python object into a YAML stream using the AstropyDumper class.
    If stream is None, return the produced string instead.

    Parameters
    ----------
    data: object
        Object to serialize to YAML
    stream : file-like object, optional
        YAML output (if not supplied a string is returned)
    **kwargs
        Other keyword arguments that get passed to yaml.dump()

    Returns
    -------
    out : str or None
        If no ``stream`` is supplied then YAML output is returned as str

    """
    kwargs['Dumper'] = AstropyDumper
    return yaml.dump(data, stream=stream, **kwargs)
