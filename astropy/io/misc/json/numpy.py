# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing core astropy objects via the
JSON protocol.
"""

import numpy as np

from .core import _json_base_encode, JSONExtendedEncoder, JSONExtendedDecoder


def register_json_extended():
    """Register :mod:`numpy` into `JSON <http://json.org>`_."""
    JSONExtendedEncoder.register_encoding(np.number)(encode_numpy_number)
    JSONExtendedDecoder.register_decoding(np.number)(decode_numpy_number)

    JSONExtendedEncoder.register_encoding(np.dtype)(encode_numpy_dtype)
    JSONExtendedDecoder.register_decoding(np.dtype)(decode_numpy_dtype)

    JSONExtendedEncoder.register_encoding(np.void)(encode_numpy_void)
    JSONExtendedDecoder.register_decoding(np.void)(decode_numpy_void)

    JSONExtendedEncoder.register_encoding(np.ndarray)(encode_ndarray)
    JSONExtendedDecoder.register_decoding(np.ndarray)(decode_ndarray)


# -------------------------------------------------------------------


def encode_numpy_number(obj):
    code = _json_base_encode(obj)
    code.update(value=obj.astype('U13'))  # TODO! dynamic length detection
    return code


def decode_numpy_number(constructor, value, code):
    return constructor(value)


def encode_numpy_dtype(obj):
    code = {"__class__": "numpy.dtype"}
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


def encode_numpy_void(obj):
    code = _json_base_encode(obj)
    code["value"] = obj.tolist()
    
    dtype = encode_numpy_dtype(obj.dtype)
    code["dtype"] = dtype if not isinstance(dtype["value"], str) else dtype["value"]

    return code


def decode_numpy_void(constructor, value, code):
    return tuple(value)


def encode_ndarray(obj):
    if isinstance(obj, np.void):
        return encode_numpy_void(obj)

    code = _json_base_encode(obj)

    # TODO! structured ndarray need to be broken into each field and
    # serialized separately
    if obj.dtype.fields:
        code["value"] = dict()
        for k in obj.dtype.fields:
            code["value"][k] = encode_ndarray(obj[k])
    else:
        code["value"] = obj.tolist()

    dtype = encode_numpy_dtype(obj.dtype)
    code["dtype"] = dtype if not isinstance(dtype["value"], str) else dtype["value"]

    return code


def decode_ndarray(constructor, value, code):
    if constructor is np.ndarray:
        constructor = np.array

    dtype = code.pop("dtype")

    # structured array
    if getattr(dtype, "fields", None) is not None:
        k1 = tuple(value.keys())[0]
        temp = np.empty(len(value[k1]), dtype=dtype)        
        for k, val in value.items():
            temp[k] = val
        value = temp

    return constructor(value, dtype=dtype, **code)
