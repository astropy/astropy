# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains functions for serializing core astropy objects via the
YAML protocol.

It provides functions `~astropy.io.misc.yaml.dump`,
`~astropy.io.misc.yaml.load`, and `~astropy.io.misc.yaml.load_all` which
call the corresponding functions in `PyYaml <http://pyyaml.org>`_ but use the
`~astropy.io.misc.yaml.AstropyDumper` and `~astropy.io.misc.yaml.AstropyLoader`
classes to define custom YAML tags for the following astropy classes:

- `astropy.units.Unit`
- `astropy.units.Quantity`
- `astropy.time.Time`
- `astropy.time.TimeDelta`
- `astropy.coordinates.SkyCoord`
- `astropy.coordinates.Angle`
- `astropy.coordinates.Latitude`
- `astropy.coordinates.Longitude`

.. Note ::

   This module requires PyYaml version 3.12 or later, which in turn requires
   Python 2.7 or Python 3.4 or later.

Example
=======
::

  >>> from astropy.io.misc import yaml
  >>> import astropy.units as u
  >>> from astropy.time import Time
  >>> from astropy.coordinates import EarthLocation

  >>> t = Time(2457389.0, format='mjd',
  ...          location=EarthLocation(1000, 2000, 3000, unit=u.km))
  >>> td = yaml.dump(t)

  >>> print(td)
  !astropy.time.Time
  format: mjd
  in_subfmt: '*'
  jd1: 4857389.5
  jd2: 0.0
  location: !astropy.coordinates.earth.EarthLocation
    ellipsoid: WGS84
    x: !astropy.units.Quantity
      unit: &id001 !astropy.units.Unit {unit: km}
      value: 1000.0
    y: !astropy.units.Quantity
      unit: *id001
      value: 2000.0
    z: !astropy.units.Quantity
      unit: *id001
      value: 3000.0
  out_subfmt: '*'
  precision: 3
  scale: utc

  >>> ty = yaml.load(td)
  >>> ty
  <Time object: scale='utc' format='mjd' value=2457389.0>

  >>> ty.location
  <EarthLocation ( 1000.,  2000.,  3000.) km>
"""

from __future__ import absolute_import, division

import base64
import numpy as np

from ...time import Time, TimeDelta
from ... import units as u
from ... import coordinates as coords
from ...utils import minversion
from ...extern import six


try:
    import yaml
except ImportError:
    raise ImportError('`import yaml` failed, PyYAML package is required for YAML')


YAML_LT_3_12 = not minversion(yaml, '3.12')


__all__ = ['AstropyLoader', 'AstropyDumper', 'load', 'load_all', 'dump']

def _unit_representer(dumper, obj):
    out = {'unit': str(obj.to_string())}
    return dumper.represent_mapping('!astropy.units.Unit', out)


def _unit_constructor(loader, node):
    map = loader.construct_mapping(node)
    return u.Unit(map['unit'])


def _time_representer(dumper, obj):
    out = obj.info._represent_as_dict()
    return dumper.represent_mapping('!astropy.time.Time', out)


def _time_constructor(loader, node):
    map = loader.construct_mapping(node)
    out = Time.info._construct_from_dict(map)
    return out


def _timedelta_representer(dumper, obj):
    out = obj.info._represent_as_dict()
    return dumper.represent_mapping('!astropy.time.TimeDelta', out)


def _timedelta_constructor(loader, node):
    map = loader.construct_mapping(node)
    out = TimeDelta.info._construct_from_dict(map)
    return out


def _ndarray_representer(dumper, obj):
    if not (obj.flags['C_CONTIGUOUS'] or obj.flags['F_CONTIGUOUS']):
        obj = np.ascontiguousarray(obj)

    if np.isfortran(obj):
        obj = obj.T
        order = 'F'
    else:
        order = 'C'

    data_b64 = base64.b64encode(obj.tostring())

    out = dict(buffer=data_b64,
               dtype=str(obj.dtype),
               shape=obj.shape,
               order=order)

    return dumper.represent_mapping('!numpy.ndarray', out)


def _ndarray_constructor(loader, node):
    map = loader.construct_mapping(node)
    map['buffer'] = base64.b64decode(map['buffer'])
    return np.ndarray(**map)


def _quantity_representer(tag):
    def representer(dumper, obj):
        out = obj.info._represent_as_dict()
        return dumper.represent_mapping(tag, out)
    return representer


def _quantity_constructor(cls):
    def constructor(loader, node):
        map = loader.construct_mapping(node)
        return cls.info._construct_from_dict(map)
    return constructor


def _skycoord_representer(dumper, obj):
    map = obj.info._represent_as_dict()
    out = dumper.represent_mapping('!astropy.coordinates.sky_coordinate.SkyCoord',
                                   map)
    return out

def _skycoord_constructor(loader, node):
    map = loader.construct_mapping(node)
    out = coords.SkyCoord.info._construct_from_dict(map)
    return out


class AstropyLoader(yaml.SafeLoader):
    """
    Custom SafeLoader that constructs astropy core objects as well
    as Python tuple and unicode objects.

    This class is not directly instantiated by user code, but instead is
    used to maintain the available constructor functions that are
    called when parsing a YAML stream.  See the `PyYaml documentation
    <http://pyyaml.org/wiki/PyYAMLDocumentation>`_ for details of the
    class signature.
    """
    def _construct_python_tuple(self, node):
        return tuple(self.construct_sequence(node))

    def _construct_python_unicode(self, node):
        return self.construct_scalar(node)

class AstropyDumper(yaml.SafeDumper):
    """
    Custom SafeDumper that represents astropy core objects as well
    as Python tuple and unicode objects.

    This class is not directly instantiated by user code, but instead is
    used to maintain the available representer functions that are
    called when generating a YAML stream from an object.  See the
    `PyYaml documentation <http://pyyaml.org/wiki/PyYAMLDocumentation>`_
    for details of the class signature.
    """
    def _represent_tuple(self, data):
        return self.represent_sequence('tag:yaml.org,2002:python/tuple', data)

    if YAML_LT_3_12:
        # pre-3.12, ignore-aliases could not deal with ndarray, so we backport
        # the more recent ignore_alises definition.
        def ignore_aliases(self, data):
            if data is None:
                return True
            if isinstance(data, tuple) and data == ():
                return True
            if isinstance(data, six.string_types + (bool, int, float)):
                return True

AstropyDumper.add_representer(u.IrreducibleUnit, _unit_representer)
AstropyDumper.add_representer(u.CompositeUnit, _unit_representer)
AstropyDumper.add_multi_representer(u.Unit, _unit_representer)
AstropyDumper.add_representer(tuple, AstropyDumper._represent_tuple)
AstropyDumper.add_representer(np.ndarray, _ndarray_representer)
AstropyDumper.add_representer(Time, _time_representer)
AstropyDumper.add_representer(TimeDelta, _timedelta_representer)
AstropyDumper.add_representer(coords.SkyCoord, _skycoord_representer)

AstropyLoader.add_constructor('tag:yaml.org,2002:python/tuple',
                              AstropyLoader._construct_python_tuple)
AstropyLoader.add_constructor('tag:yaml.org,2002:python/unicode',
                              AstropyLoader._construct_python_unicode)
AstropyLoader.add_constructor('!astropy.units.Unit', _unit_constructor)
AstropyLoader.add_constructor('!numpy.ndarray', _ndarray_constructor)
AstropyLoader.add_constructor('!astropy.time.Time', _time_constructor)
AstropyLoader.add_constructor('!astropy.time.TimeDelta', _timedelta_constructor)
AstropyLoader.add_constructor('!astropy.coordinates.sky_coordinate.SkyCoord',
                              _skycoord_constructor)

for cls, tag in ((u.Quantity, '!astropy.units.Quantity'),
                 (coords.Angle, '!astropy.coordinates.Angle'),
                 (coords.Latitude, '!astropy.coordinates.Latitude'),
                 (coords.Longitude, '!astropy.coordinates.Longitude'),
                 (coords.EarthLocation, '!astropy.coordinates.earth.EarthLocation')):
    AstropyDumper.add_multi_representer(cls, _quantity_representer(tag))
    AstropyLoader.add_constructor(tag, _quantity_constructor(cls))


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
