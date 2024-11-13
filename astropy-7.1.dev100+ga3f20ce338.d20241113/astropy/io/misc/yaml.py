# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Functions for serializing astropy objects to YAML.

It provides functions `~astropy.io.misc.yaml.dump`,
`~astropy.io.misc.yaml.load`, and `~astropy.io.misc.yaml.load_all` which
call the corresponding functions in `PyYaml <https://pyyaml.org>`_ but use the
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
- `astropy.coordinates.EarthLocation`
- `astropy.table.SerializedColumn`

Examples
--------
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
    jd1: 4857390.0
    jd2: -0.5
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
    >>> ty.location  # doctest: +FLOAT_CMP
    <EarthLocation (1000., 2000., 3000.) km>
"""

import base64

import numpy as np
import yaml

from astropy import coordinates as coords
from astropy import units as u
from astropy.table import SerializedColumn
from astropy.time import Time, TimeDelta

__all__ = ["AstropyLoader", "AstropyDumper", "load", "load_all", "dump"]


def _unit_representer(dumper, obj):
    out = {"unit": str(obj.to_string())}
    return dumper.represent_mapping("!astropy.units.Unit", out)


def _unit_constructor(loader, node):
    map = loader.construct_mapping(node)
    return u.Unit(map["unit"], parse_strict="warn")


def _serialized_column_representer(dumper, obj):
    out = dumper.represent_mapping("!astropy.table.SerializedColumn", obj)
    return out


def _serialized_column_constructor(loader, node):
    map = loader.construct_mapping(node)
    return SerializedColumn(map)


def _time_representer(dumper, obj):
    out = obj.info._represent_as_dict()
    return dumper.represent_mapping("!astropy.time.Time", out)


def _time_constructor(loader, node):
    map = loader.construct_mapping(node)
    out = Time.info._construct_from_dict(map)
    return out


def _timedelta_representer(dumper, obj):
    out = obj.info._represent_as_dict()
    return dumper.represent_mapping("!astropy.time.TimeDelta", out)


def _timedelta_constructor(loader, node):
    map = loader.construct_mapping(node)
    out = TimeDelta.info._construct_from_dict(map)
    return out


def _ndarray_representer(dumper, obj):
    if obj.dtype.hasobject:
        raise TypeError(f"cannot serialize numpy object array: {obj}")

    if not (obj.flags["C_CONTIGUOUS"] or obj.flags["F_CONTIGUOUS"]):
        obj = np.ascontiguousarray(obj)

    if np.isfortran(obj):
        obj = obj.T
        order = "F"
    else:
        order = "C"

    data_b64 = base64.b64encode(obj.tobytes())

    out = {
        "buffer": data_b64,
        "dtype": str(obj.dtype) if not obj.dtype.fields else obj.dtype.descr,
        "shape": obj.shape,
        "order": order,
    }

    return dumper.represent_mapping("!numpy.ndarray", out)


def _ndarray_constructor(loader, node):
    # Convert mapping to a dict useful for initializing ndarray.
    # Need deep=True since for structured dtype, the contents
    # include lists and tuples, which need recursion via
    # construct_sequence.
    map = loader.construct_mapping(node, deep=True)
    map["buffer"] = base64.b64decode(map["buffer"])

    if map["dtype"] == "object":
        raise TypeError("cannot load numpy array with dtype object")
    return np.ndarray(**map)


def _void_representer(dumper, obj):
    data_b64 = base64.b64encode(obj.tobytes())
    out = {
        "buffer": data_b64,
        "dtype": str(obj.dtype) if not obj.dtype.fields else obj.dtype.descr,
    }
    return dumper.represent_mapping("!numpy.void", out)


def _void_constructor(loader, node):
    # Interpret as node as an array scalar and then index to change to void.
    map = loader.construct_mapping(node, deep=True)
    map["buffer"] = base64.b64decode(map["buffer"])
    return np.ndarray(shape=(), **map)[()]


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
    out = dumper.represent_mapping("!astropy.coordinates.sky_coordinate.SkyCoord", map)
    return out


def _skycoord_constructor(loader, node):
    map = loader.construct_mapping(node)
    out = coords.SkyCoord.info._construct_from_dict(map)
    return out


# Straight from yaml's Representer
def _complex_representer(self, data):
    if data.imag == 0.0:
        data = f"{data.real!s}"
    elif data.real == 0.0:
        data = f"{data.imag!s}j"
    elif data.imag > 0:
        data = f"{data.real!s}+{data.imag!s}j"
    else:
        data = f"{data.real!s}{data.imag!s}j"
    return self.represent_scalar("tag:yaml.org,2002:python/complex", data)


def _complex_constructor(loader, node):
    map = loader.construct_scalar(node)
    return complex(map)


class AstropyLoader(yaml.SafeLoader):
    """
    Custom SafeLoader that constructs astropy core objects as well
    as Python tuple and unicode objects.

    This class is not directly instantiated by user code, but instead is
    used to maintain the available constructor functions that are
    called when parsing a YAML stream.  See the `PyYaml documentation
    <https://pyyaml.org/wiki/PyYAMLDocumentation>`_ for details of the
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
    `PyYaml documentation <https://pyyaml.org/wiki/PyYAMLDocumentation>`_
    for details of the class signature.
    """

    def _represent_tuple(self, data):
        return self.represent_sequence("tag:yaml.org,2002:python/tuple", data)

    def represent_float(self, data):
        # Override to change repr(data) to str(data) since otherwise all the
        # numpy scalars fail in not NUMPY_LT_2_0.
        # otherwise, this function is identical to yaml.SafeDumper.represent_float
        # (as of pyyaml 6.0.1)
        if data != data or (data == 0.0 and data == 1.0):
            value = ".nan"
        elif data == self.inf_value:
            value = ".inf"
        elif data == -self.inf_value:
            value = "-.inf"
        else:
            value = str(data).lower()
            # Note that in some cases `repr(data)` represents a float number
            # without the decimal parts.  For instance:
            #   >>> repr(1e17)
            #   '1e17'
            # Unfortunately, this is not a valid float representation according
            # to the definition of the `!!float` tag.  We fix this by adding
            # '.0' before the 'e' symbol.
            if "." not in value and "e" in value:
                value = value.replace("e", ".0e", 1)
        return self.represent_scalar("tag:yaml.org,2002:float", value)


AstropyDumper.add_multi_representer(u.UnitBase, _unit_representer)
AstropyDumper.add_multi_representer(u.FunctionUnitBase, _unit_representer)
AstropyDumper.add_multi_representer(u.StructuredUnit, _unit_representer)
AstropyDumper.add_representer(tuple, AstropyDumper._represent_tuple)
AstropyDumper.add_representer(np.ndarray, _ndarray_representer)
AstropyDumper.add_representer(np.void, _void_representer)
AstropyDumper.add_representer(Time, _time_representer)
AstropyDumper.add_representer(TimeDelta, _timedelta_representer)
AstropyDumper.add_representer(coords.SkyCoord, _skycoord_representer)
AstropyDumper.add_representer(SerializedColumn, _serialized_column_representer)

# Numpy dtypes
AstropyDumper.add_representer(np.bool_, yaml.representer.SafeRepresenter.represent_bool)
for np_type in [
    np.intc,
    np.intp,
    np.int8,
    np.int16,
    np.int32,
    np.int64,
    np.uint8,
    np.uint16,
    np.uint32,
    np.uint64,
]:
    AstropyDumper.add_representer(
        np_type, yaml.representer.SafeRepresenter.represent_int
    )
for np_type in [np.float16, np.float32, np.float64, np.longdouble]:
    AstropyDumper.add_representer(np_type, AstropyDumper.represent_float)
for np_type in [complex, np.complex64, np.complex128]:
    AstropyDumper.add_representer(np_type, _complex_representer)

AstropyLoader.add_constructor("tag:yaml.org,2002:python/complex", _complex_constructor)
AstropyLoader.add_constructor(
    "tag:yaml.org,2002:python/tuple", AstropyLoader._construct_python_tuple
)
AstropyLoader.add_constructor(
    "tag:yaml.org,2002:python/unicode", AstropyLoader._construct_python_unicode
)
AstropyLoader.add_constructor("!astropy.units.Unit", _unit_constructor)
AstropyLoader.add_constructor("!numpy.ndarray", _ndarray_constructor)
AstropyLoader.add_constructor("!numpy.void", _void_constructor)
AstropyLoader.add_constructor("!astropy.time.Time", _time_constructor)
AstropyLoader.add_constructor("!astropy.time.TimeDelta", _timedelta_constructor)
AstropyLoader.add_constructor(
    "!astropy.coordinates.sky_coordinate.SkyCoord", _skycoord_constructor
)
AstropyLoader.add_constructor(
    "!astropy.table.SerializedColumn", _serialized_column_constructor
)

for cls, tag in (
    (u.Quantity, "!astropy.units.Quantity"),
    (u.Magnitude, "!astropy.units.Magnitude"),
    (u.Dex, "!astropy.units.Dex"),
    (u.Decibel, "!astropy.units.Decibel"),
    (coords.Angle, "!astropy.coordinates.Angle"),
    (coords.Latitude, "!astropy.coordinates.Latitude"),
    (coords.Longitude, "!astropy.coordinates.Longitude"),
    (coords.EarthLocation, "!astropy.coordinates.earth.EarthLocation"),
):
    AstropyDumper.add_multi_representer(cls, _quantity_representer(tag))
    AstropyLoader.add_constructor(tag, _quantity_constructor(cls))

for cls in list(coords.representation.REPRESENTATION_CLASSES.values()) + list(
    coords.representation.DIFFERENTIAL_CLASSES.values()
):
    name = cls.__name__
    # Add representations/differentials defined in astropy.
    if name in coords.representation.__all__:
        tag = "!astropy.coordinates." + name
        AstropyDumper.add_multi_representer(cls, _quantity_representer(tag))
        AstropyLoader.add_constructor(tag, _quantity_constructor(cls))


def load(stream):
    """Parse the first YAML document in a stream using the AstropyLoader and
    produce the corresponding Python object.

    Parameters
    ----------
    stream : str or file-like
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
    stream : str or file-like
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
    data : object
        Object to serialize to YAML
    stream : file-like, optional
        YAML output (if not supplied a string is returned)
    **kwargs
        Other keyword arguments that get passed to yaml.dump()

    Returns
    -------
    out : str or None
        If no ``stream`` is supplied then YAML output is returned as str

    """
    kwargs["Dumper"] = AstropyDumper
    kwargs.setdefault("default_flow_style", None)
    return yaml.dump(data, stream=stream, **kwargs)
