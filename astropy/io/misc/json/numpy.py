# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing core astropy objects via the
JSON protocol.
"""

import numpy as np

from .core import _json_base_encode, JSONExtendedEncoder, JSONExtendedDecoder


def register_json_extended():

    JSONExtendedEncoder.register_encoding(np.number)(encode_numpy_number)
    JSONExtendedDecoder.register_decoding(np.number)(decode_numpy_number)

    JSONExtendedEncoder.register_encoding(np.dtype)(encode_numpy_dtype)
    JSONExtendedDecoder.register_decoding(np.dtype)(decode_numpy_dtype)

    JSONExtendedEncoder.register_encoding(np.ndarray)(encode_ndarray)
    JSONExtendedDecoder.register_decoding(np.ndarray)(decode_ndarray)


# -------------------------------------------------------------------


def encode_numpy_number(obj):
    code = _json_base_encode(obj)
    code.update(value=obj.astype('U13'))
    return code


def decode_numpy_number(constructor, value, code):
    return constructor(value)


def encode_numpy_dtype(obj):
    code = {"__class__": "numpy.dtype", "value": str(obj)}
    # TODO! more sophisticated encoding for structured dtype
    return code


def decode_numpy_dtype(constructor, value, code):
    return constructor(value, **code)


def encode_ndarray(obj):
    code = _json_base_encode(obj)
    code.update(value=obj.tolist(), dtype=encode_numpy_dtype(obj.dtype))
    return code


def decode_ndarray(constructor, value, code):
    if constructor is np.ndarray:
        constructor = np.array
    return constructor(value, **code)
