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

NumPy arrays, masked arrays, and the `numpy.ma.masked` singleton are also
supported. Masked arrays retain their data, mask, fill value, and hard-mask
setting when serialized.

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
import itertools

import numpy as np
import yaml

from astropy import coordinates as coords
from astropy import units as u
from astropy.table import SerializedColumn
from astropy.time import Time, TimeDelta

__all__ = ["AstropyDumper", "AstropyLoader", "dump", "load", "load_all"]


def _unit_representer(dumper, obj):
    out = {"unit": str(obj.to_string())}
    return dumper.represent_mapping("!astropy.units.Unit", out)


def _unit_constructor(loader, node):
    map = loader.construct_mapping(node)
    return u.Unit(map["unit"], parse_strict="warn")


def _serialized_column_representer(dumper, obj):
    return dumper.represent_mapping("!astropy.table.SerializedColumn", obj)


def _serialized_column_constructor(loader, node):
    map = loader.construct_mapping(node)
    return SerializedColumn(map)


def _time_representer(dumper, obj):
    out = obj.info._represent_as_dict()
    return dumper.represent_mapping("!astropy.time.Time", out)


def _time_constructor(loader, node):
    map = loader.construct_mapping(node)
    return Time.info._construct_from_dict(map)


def _timedelta_representer(dumper, obj):
    out = obj.info._represent_as_dict()
    return dumper.represent_mapping("!astropy.time.TimeDelta", out)


def _timedelta_constructor(loader, node):
    map = loader.construct_mapping(node)
    return TimeDelta.info._construct_from_dict(map)


def _dtype_to_yaml(dtype):
    """Describe a dtype without turning padding into extra fields."""
    if dtype.subdtype is not None:
        base, shape = dtype.subdtype
        return _dtype_to_yaml(base), shape
    if dtype.names is None:
        return str(dtype)

    fields = [dtype.fields[name] for name in dtype.names]
    return {
        "names": list(dtype.names),
        "formats": [_dtype_to_yaml(field[0]) for field in fields],
        "offsets": [field[1] for field in fields],
        "titles": [field[2] if len(field) == 3 else None for field in fields],
        "itemsize": dtype.itemsize,
        "aligned": dtype.isalignedstruct,
    }


def _dtype_from_yaml(dtype):
    """Return numpy dtype from a YAML ``dtype`` value, refusing object dtypes.

    The dtype and buffer of a numpy array in YAML come straight from the file, which may
    be untrusted. For an object dtype the buffer would be interpreted as ``PyObject``
    pointers that get dereferenced on access, so such arrays must never be constructed.
    Using ``hasobject`` on the resolved dtype catches every spelling (``object``, ``O``)
    as well as object fields nested in a structured dtype.
    """
    if isinstance(dtype, dict) and "formats" in dtype:
        dtype = dict(dtype)
        # Resolve nested dtypes first: otherwise NumPy propagates the outer
        # alignment setting into fields that can themselves be packed.
        dtype["formats"] = [_dtype_from_yaml(item) for item in dtype["formats"]]
    elif isinstance(dtype, tuple):
        dtype = (_dtype_from_yaml(dtype[0]), dtype[1])

    dtype = np.dtype(dtype)
    if dtype.hasobject:
        raise TypeError(f"cannot load numpy array with object dtype {dtype}")
    return dtype


def _ndarray_representer(dumper, obj):
    if obj.dtype.hasobject:
        raise TypeError(f"cannot serialize numpy object array: {obj}")

    shape = obj.shape
    if obj.flags.c_contiguous:
        order = "C"
    elif obj.flags.f_contiguous:
        order = "F"
        obj = obj.T
    else:
        order = "C"
        obj = np.ascontiguousarray(obj)

    data_b64 = base64.b64encode(obj.tobytes())

    out = {
        "buffer": data_b64,
        "dtype": _dtype_to_yaml(obj.dtype),
        "shape": shape,
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
    map["dtype"] = _dtype_from_yaml(map["dtype"])
    return np.ndarray(**map)


def _void_representer(dumper, obj):
    if obj.dtype.hasobject:
        raise TypeError(f"cannot serialize numpy object array: {obj}")

    data_b64 = base64.b64encode(obj.tobytes())
    out = {
        "buffer": data_b64,
        "dtype": _dtype_to_yaml(obj.dtype),
    }
    return dumper.represent_mapping("!numpy.void", out)


def _void_constructor(loader, node):
    # Interpret as node as an array scalar and then index to change to void.
    map = loader.construct_mapping(node, deep=True)
    map["buffer"] = base64.b64decode(map["buffer"])
    map["dtype"] = _dtype_from_yaml(map["dtype"])
    return np.ndarray(shape=(), **map)[()]


def _masked_constant_representer(dumper, obj):
    return dumper.represent_scalar("!numpy.ma.masked", "")


def _masked_constant_constructor(loader, node):
    loader.construct_scalar(node)
    return np.ma.masked


def _masked_array_representer(dumper, obj):
    out = {
        "data": obj.data,
        "fill_value": np.asarray(obj.fill_value),
        "hard_mask": obj.hardmask,
    }
    if obj.mask is not np.ma.nomask:
        out["mask"] = obj.mask
    return dumper.represent_mapping("!numpy.ma.MaskedArray", out)


def _masked_array_constructor(loader, node):
    map = loader.construct_mapping(node, deep=True)
    fill_value = np.asarray(map.pop("fill_value"))
    default_fill = np.asarray(np.ma.array([], dtype=map["data"].dtype).fill_value)

    # NumPy's default fill value can exceed the data dtype's range. Passing
    # it explicitly to the constructor would cast it to that dtype.
    if fill_value.dtype != default_fill.dtype:
        map["fill_value"] = fill_value

    # Array buffers loaded from YAML are read-only. Copy them so that both
    # the data and the mask can be updated, including by masked assignment.
    result = np.ma.array(**map, copy=True)
    if fill_value.dtype == default_fill.dtype:
        # Serialization initialized the original fill value. Preserve that
        # state so a later dtype conversion carries the same fill value.
        result.get_fill_value()
        if not np.array_equal(fill_value, default_fill):
            result.fill_value = fill_value
    return result


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
    return dumper.represent_mapping("!astropy.coordinates.sky_coordinate.SkyCoord", map)


def _skycoord_constructor(loader, node):
    map = loader.construct_mapping(node)
    return coords.SkyCoord.info._construct_from_dict(map)


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
        # numpy scalars fail;
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
AstropyDumper.add_representer(type(np.ma.masked), _masked_constant_representer)
AstropyDumper.add_representer(np.ma.MaskedArray, _masked_array_representer)
AstropyDumper.add_representer(np.void, _void_representer)
AstropyDumper.add_representer(Time, _time_representer)
AstropyDumper.add_representer(TimeDelta, _timedelta_representer)
AstropyDumper.add_representer(coords.SkyCoord, _skycoord_representer)
AstropyDumper.add_representer(SerializedColumn, _serialized_column_representer)

# Numpy dtypes
AstropyDumper.add_representer(np.bool_, yaml.representer.SafeRepresenter.represent_bool)
AstropyDumper.add_representer(np.str_, yaml.representer.SafeRepresenter.represent_str)
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
AstropyLoader.add_constructor("!numpy.ma.masked", _masked_constant_constructor)
AstropyLoader.add_constructor("!numpy.ma.MaskedArray", _masked_array_constructor)
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

# Add representations, differentials, and built-in frames defined in astropy and in the
# ``astropy.coordinates`` public API.
cls_coords = [
    cls
    for cls in itertools.chain(
        coords.representation.REPRESENTATION_CLASSES.values(),
        coords.representation.DIFFERENTIAL_CLASSES.values(),
        coords.frame_transform_graph.frame_set,
    )
    if cls.__name__ in dir(coords)
]
for cls in cls_coords:
    name = cls.__name__
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
