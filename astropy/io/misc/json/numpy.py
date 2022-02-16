# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing core astropy objects via the
JSON protocol.
"""

import numpy as np

from .core import JSONExtendedDecoder, JSONExtendedEncoder, _json_base_encode


def register_json_extended():
    """:mod:`astropy.io.misc.json` entry points for :mod:`numpy`."""
    JSONExtendedEncoder.register_encoding(np.dtype)(encode_numpy_dtype)
    JSONExtendedDecoder.register_decoding(np.dtype)(decode_numpy_dtype)

    JSONExtendedEncoder.register_encoding(np.bool_)(encode_numpy_bool)
    JSONExtendedDecoder.register_decoding(np.bool_)(decode_numpy_bool)

    JSONExtendedEncoder.register_encoding(np.number)(encode_numpy_number)
    JSONExtendedDecoder.register_decoding(np.number)(decode_numpy_number)

    JSONExtendedEncoder.register_encoding(np.void)(encode_numpy_void)
    JSONExtendedDecoder.register_decoding(np.void)(decode_numpy_void)

    JSONExtendedEncoder.register_encoding(np.ndarray)(encode_ndarray)
    JSONExtendedDecoder.register_decoding(np.ndarray)(decode_ndarray)


# -------------------------------------------------------------------


def encode_numpy_dtype(obj):
    """Return a `numpy.dtype` as a JSON-able dictionary.

    `numpy.dtype` spans many different cases -- builtin types, aligned, shaped,
    structured, with metadata, etc. -- and the contents of the JSON will
    reflect this.
    The "metadata" field will only be included if there is metadata.
    If the type is a built-in (python or numpy), then "value" is just a string.
    If the type has a shape, then "value" is a list of the type and shape.
    Structured dtypes have field "align" and "value" is a nested dictionary
    of the contained dtypes.

    Examples
    --------
    >>> import json
    >>> from astropy.io.misc.json import JSONExtendedEncoder, JSONExtendedDecoder
    >>> def show(val):
    ...     print("value:", repr(val))
    ...     serialized = json.dumps(val, cls=JSONExtendedEncoder)
    ...     print("dump:", serialized)
    ...     out = json.loads(serialized, cls=JSONExtendedDecoder)
    ...     print("load:", repr(out))

    >>> show(np.dtype(float))
    value: dtype('float64')
    dump: {"!": "numpy.dtype", "value": "float64"}
    load: dtype('float64')

    >>> show(np.dtype(np.dtype("float", metadata={"a": 1})))
    value: dtype('float64')
    dump: {"!": "numpy.dtype", "value": "float64", "metadata": {"a": 1}}
    load: dtype('float64')

    >>> show(np.dtype("10float64", align=True))
    value: dtype(('<f8', (10,)))
    dump: {"!": "numpy.dtype", "value": "float64", "shape": [10]}
    load: dtype(('<f8', (10,)))
    """
    clsname = obj.__class__.__name__.split('[')[0]  # e.g. dtype[void] -> dtype
    code = {"!": f"{obj.__class__.__module__}.{clsname}"}

    # Built-in or only metadata making it not "isbuiltin"
    if obj.isbuiltin or (obj.fields is None and obj.subdtype is None):
        code["value"] = str(obj)
        # don't need `align` and copy is True
    elif obj.fields is None:  # Shaped, not structured
        code["value"] = str(obj.subdtype[0])
        code["shape"] = obj.subdtype[1]
    else:  # Structured
        code.update(value=dict(obj.fields),
                    align=False if obj.alignment == 1 else True)

    if obj.metadata is not None:
        code["metadata"] = dict(obj.metadata)

    # iterate through sub-dtypes (if structured), abbreviating.
    if isinstance(code["value"], dict):
        for k, v in code["value"].items():
            if not isinstance(v, tuple) or not isinstance(v[0], np.dtype):
                continue  # skip non dtypes

            dt = encode_numpy_dtype(v[0])
            if dt["!"] == "numpy.dtype":
                dt.pop("!")  # we know the dtype
            if dt.keys() == {"value"}:  # only value, no e.g. "metadata" or "align"
                dt = dt["value"]
            code["value"][k] = [dt, v[1]]  # (list not tuple b/c json)

    return code


def decode_numpy_dtype(constructor, value, code):
    """Return a `numpy.dtype` from an ``encode_numpy_dtype`` dictionary."""
    # Structured: iterate through sub-dtypes, decoding.
    if isinstance(value, dict):
        for k, v in value.items():  # need to convert lists to tuples for numpy
            if not isinstance(v, list):
                continue
            elif isinstance(v[0], np.dtype):
                value[k] = tuple(v)
            else:
                v_c = {"value": v[0]} if isinstance(v[0], str) else v[0]
                dt = decode_numpy_dtype(np.dtype, v_c.pop("value"), v_c)
                value[k] = (dt, v[1])
    elif "shape" in code:
        value = (value, tuple(code.pop("shape")))

    # Prune metadata, if None.
    if "metadata" in code and code["metadata"] is None:
        code.pop("metadata")

    return constructor(value, **code)


def _encode_abbreviate_dtype(dtype):
    """Abbreviate a dtype for nested structures, like `numpy.ndarray`"""
    dt = encode_numpy_dtype(dtype)
    dt.pop("!", None)  # we know the dtype

    if dt.keys() == {"value"}:  # only value, no e.g. "metadata" or "align"
        dt = dt["value"]
    return dt


def _decode_abbreviated_dtype(dtype):
    """Recreate a dtype from an abbreviated form, like in `numpy.ndarray`"""
    if isinstance(dtype, str):
        dtype = {"value": dtype}
    # iterate through subdtypes (if structured), similarly un-abbreviating.
    elif isinstance(dtype["value"], dict):
        for k, v in dtype["value"].items():
            if isinstance(v, list):
                dtype["value"][k] = (_decode_abbreviated_dtype(v[0]), v[1])

    # only if "!" wasn't a key, e.g. in non-structured ndarray
    if isinstance(dtype, dict):
        dtype = decode_numpy_dtype(np.dtype, dtype.pop("value"), dtype)

    return dtype


# ===================================================================
# Scalars

def encode_numpy_bool(obj):
    """Return `numpy.bool` as a JSON-able dictionary."""
    code = _json_base_encode(obj)
    code["value"] = bool(obj)
    return code


def decode_numpy_bool(constructor, value, code):
    """Return a `numpy.bool` from an ``encode_numpy_bool`` dictionary."""
    return constructor(value)


def encode_numpy_number(obj):
    """Return `numpy.number` as a JSON-able dictionary."""
    code = _json_base_encode(obj)
    code.update(value=obj.astype('U13'))  # TODO! dynamic length detection
    return code


def decode_numpy_number(constructor, value, code):
    """Return a `numpy.number` from an ``encode_numpy_number`` dictionary."""
    return constructor(value)


def encode_numpy_void(obj):
    """Encode `numpy.void` as a JSON-able dictionary.

    `~numpy.void` can be either a bytes or an item from a structured array.
    The former is a simple wrapper around the bytes encoding, the latter
    looks more like an encoded ndarray and has a "dtype" field.
    `numpy.void.flags` cannot be set, so are not included in this encoding.

    Returns
    -------
    dict
        With fields "value" (always) and "dtype" (if structured).
    """
    code = _json_base_encode(obj)
    value = tuple(map(str, obj))  # go through str to keep precision

    if value == ():  # it's actually bytes
        code["value"] = obj.tolist()
    else:
        code["value"] = value
        code["dtype"] = _encode_abbreviate_dtype(obj.dtype)
    return code


def decode_numpy_void(constructor, value, code):
    """Return a `numpy.void` from an ``encode_numpy_void`` dictionary."""
    if isinstance(value, bytes):
        # TODO? can't have field dtype, or needs to be "|V..."
        return constructor(value)

    dt = code.pop("dtype")
    dtype = decode_numpy_dtype(np.dtype, dt.pop("value"), dt)

    return np.array(tuple(value), dtype=dtype)[()]


# ===================================================================
# Arrays

def _check_flags(flags):
    """Checks that flags are not their default values."""
    out = {}
    if flags.writeable is False:
        out["write"] = False
    if flags.aligned is False:
        out["align"] = False
    if flags.writebackifcopy is True:
        out["uic"] = True
    return out


def encode_ndarray(obj):
    """Return a `numpy.ndarray` as a JSON-able dictionary.

    `numpy.ndarray` is a multipurpose object that can contain values of many
    different datatypes, even non-homogenous ones. Consequently the encoding
    can become quite verbose, but hopefully simple arrays are encoded simply.
    The "flags" field is included only if any flag is not the default value.
    The `numpy.dtype` is abbreviated and cannot be decoded separately.
    Structured arrays are encoded as a nested dictionary, with each field
    the encoded constituent array and the primary `numpy.dtype` used to
    coordinate the reconstruction. Below are examples.

    Examples
    --------
    >>> import json
    >>> from astropy.io.misc.json import JSONExtendedEncoder, JSONExtendedDecoder
    >>> def show(val):
    ...     print("value:", repr(val))
    ...     serialized = json.dumps(val, cls=JSONExtendedEncoder)
    ...     print("dump:", repr(serialized))
    ...     out = json.loads(serialized, cls=JSONExtendedDecoder)
    ...     print("load:", repr(out))

    We start with a scalar array.

    >>> show(np.array(np.float64(10)))
    value: array(10.)
    dump: '{"!": "numpy.ndarray", "value": "10.0", "dtype": "float64"}'
    load: array(10.)

    Actual arrays are econded very similarly. Note that the value is in a list.

    >>> show(np.array([3], dtype=float))
    value: array([3.])
    dump: '{"!": "numpy.ndarray", "value": ["3.0"], "dtype": "float64"}'
    load: array([3.])

    Non-trivial `numpy.dtype` are properly encoded. For further information see
    ``encode_numpy_dtype``.

    >>> show(np.array([3], dtype=np.dtype("int16", metadata={"a": 1})))
    value: array([3], dtype=int16)
    dump: '{"!": "numpy.ndarray", "value": ["3"],
            "dtype": {"value": "int16", "metadata": {"a": 1}}}'
    load: array([3], dtype=int16)

    Structured arrays are done in a nested structure so :mod:`json` properly
    decodes each component when it recurses through "value". The structured
    array is assembled using the "dtype".

    >>> show(np.array((0, 0.6), dtype=np.dtype([("f1", float), ("f2", np.float32)])))
    value: array((0., 0.6), dtype=[('f1', '<f8'), ('f2', '<f4')])
    dump: '{"!": "numpy.ndarray",
            "value": {"f1": {"!": "numpy.ndarray", "value": "0.0", "dtype": "float64"},             "f2": {"!": "numpy.ndarray", "value": "0.6", "dtype": "float32"}},
            "dtype": {"value": {"f1": ["float64", 0], "f2": ["float32", 8]},
                      "align": false}}'
    load: array([(0., 0.6)], dtype=[('f1', '<f8'), ('f2', '<f4')])

    A more complex example of the above:

    >>> show(np.array((0, (1, 0.6)),
    ...               dtype=np.dtype([("f1", float),
    ...                               ("f2", [("s1", np.float32), ("s2", np.float32)])])))
    value: array((0., (1., 0.6)),
                 dtype=[('f1', '<f8'), ('f2', [('s1', '<f4'), ('s2', '<f4')])])
    dump: '{"!": "numpy.ndarray",
            "value":
                {"f1": {"!": "numpy.ndarray", "value": "0.0", "dtype": "float64"},
                 "f2": {"!": "numpy.ndarray",
                        "value":
                            {"s1": {"!": "numpy.ndarray", "value": "1.0", "dtype": "float32"},
                             "s2": {"!": "numpy.ndarray", "value": "0.6", "dtype": "float32"}},
                        "dtype": {"value": {"s1": ["float32", 0], "s2": ["float32", 4]},
                                  "align": false}}},
            "dtype": {"value": {"f1": ["float64", 0],
                                "f2": [{"value": {"s1": ["float32", 0], "s2": ["float32", 4]}, "align": false}, 8]},
                      "align": false}}'
    load: array([(0., (1., 0.6))],
          dtype=[('f1', '<f8'), ('f2', [('s1', '<f4'), ('s2', '<f4')])])
    """
    # For convenience, check if actually a void
    if isinstance(obj, np.void):
        return encode_numpy_void(obj)

    code = _json_base_encode(obj)

    # Encode the actual value. Normal arrays are turned into a list, but first
    # routed through str, which keeps the precision and strips the numpy type.
    # Otherwise an array of e.g. int32 would become a list of dictionaries
    # each of which is the JSON serialization of one int32 number.
    # Structured arrays are serialized through the `fields`, which
    if obj.dtype.fields is not None:  # Let JSON detect and encode the values
        code["value"] = {k: obj[k] for k in obj.dtype.fields}
    else:
        code["value"] = obj.astype(str).tolist()

    # Encode the dtype. Normally JSON would do this automatically, using the
    # registerd dtype encoder. But we want to abbreviate the dtype, so must
    # encode it now, with a custom encoding function.
    code["dtype"] = _encode_abbreviate_dtype(obj.dtype)

    # flags
    if flags := _check_flags(obj.flags):
        code["flags"] = flags

    return code


def decode_ndarray(constructor, value, code):
    """Return a `numpy.ndarray` from an ``encode_ndarray`` dictionary."""
    if constructor is np.ndarray:
        constructor = np.array

    dtype = _decode_abbreviated_dtype(code.pop("dtype"))
    flags = code.pop("flags", {})

    # structured array
    if getattr(dtype, "fields", None) is not None:
        k1 = tuple(value.keys())[0]
        temp = np.empty(np.size(value[k1]), dtype=dtype)
        for k, val in value.items():
            temp[k] = val
        value = temp

    array = constructor(value, dtype=dtype, **code)
    array.setflags(**flags)
    return array
