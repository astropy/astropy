# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing core astropy objects via the
JSON protocol.
"""

import numpy as np

from .core import _json_base_encode, JSONExtendedEncoder, JSONExtendedDecoder


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
    >>> def show(val):
    ...     serialized = json.dumps(val, cls=JSONExtendedEncoder)
    ...     out = json.loads(serialized, cls=JSONExtendedDecoder)
    ...     return f"value: {val!r}\ndump: {serialized}\nload: {out!r}"

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
    dump: {"!": "numpy.dtype", "value": ["float64", [10]]}
    load: dtype(('<f8', (10,)))
    """
    code = {"!": "numpy.dtype"}
    if obj.isbuiltin:
        code["value"] = str(obj)
        # don't need `align` and copy is True
    elif obj.fields is None:  # not structured
        if obj.subdtype is None:  # only metadata
            code["value"] = str(obj)
        else:
            code["value"] = [str(obj.subdtype[0]), obj.subdtype[1]]
    else:  # structured
        code.update(value=dict(obj.fields),
                    align=False if obj.alignment == 1 else True)

    if obj.metadata is not None:
        code["metadata"] = dict(obj.metadata)

    return code


def decode_numpy_dtype(constructor, value, code):
    # not builtin, but scalar.
    if isinstance(value, list):
        value = (value[0], tuple(value[1]), *value[2:])
    elif isinstance(value, dict):  # structured
        for k, v in value.items():
            if isinstance(v, list):
                value[k] = tuple(v)

    # Prune metadata, if None.
    if "metadata" in code and code["metadata"] is None:
        code.pop("metadata")

    return constructor(value, **code)


def _abbreviate_dtype(dtype):
    dt = encode_numpy_dtype(dtype)
    dt.pop("!")  # we know the dtype
    if dt.keys() == {"value"}:  # only value, no e.g. metadata
        dt = dt["value"]
    return dt


def _unabbreviate_dtype(dtype):
    if isinstance(dtype, str):
        dtype = {"value": dtype}
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
        code["dtype"] = _abbreviate_dtype(obj.dtype)
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
    # For convenience, check if actually a void
    if isinstance(obj, np.void):
        return encode_numpy_void(obj)

    code = _json_base_encode(obj)

    # Encode the actual value. Normal arrays are turned into a list, but first
    # routed through str, which keeps the precision and strips the numpy type.
    # Otherwise an array of e.g. int32 would become a list of dictionaries
    # each of which is the JSON serialization of one int32 number.
    # Structured arrays are seriazed
    if obj.dtype.fields:
        code["value"] = dict()
        for k in obj.dtype.fields:
            code["value"][k] = encode_ndarray(obj[k])
    else:
        code["value"] = obj.astype(str).tolist()

    # dtype
    code["dtype"] = _abbreviate_dtype(obj.dtype)

    # flags
    if flags := _check_flags(obj.flags):
        code["flags"] = flags

    return code


def decode_ndarray(constructor, value, code):
    if constructor is np.ndarray:
        constructor = np.array

    dtype = _unabbreviate_dtype(code.pop("dtype"))
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
