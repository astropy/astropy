# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing core astropy objects via the
JSON protocol.
"""

import importlib
import json

import numpy as np

import astropy.units as u
import astropy.coordinates as coord
from astropy.io.utils import load_all_entry_points


__all__ = ['JSONExtendedEncoder', 'JSONExtendedDecoder']


class JSONExtendedEncoder(json.JSONEncoder):
    """Support for data types that JSON default encoder does not do."""

    _registry = []  # list of tuples

    def default(self, obj):
        for cls, func in self._registry:
            if isinstance(obj, cls):
                code = func(obj)
                break
        else:
            code = super().default(obj)

        return code

    @classmethod
    def register_encoding(cls, type):
        def register(func):
            # inserting subclasses before parent classes, so encountered first
            for i, (key, _) in enumerate(cls._registry):
                if issubclass(type, key):
                    cls._registry.insert(i, (type, func))
                    break
            else:  # put at the end
                cls._registry.append((type, func))
            return func

        return register


class JSONExtendedDecoder(json.JSONDecoder):

    _registry = []

    def __init__(self, *, parse_float=None, parse_int=None, parse_constant=None,
                 strict=True):
        super().__init__(object_hook=self.object_hook, parse_float=parse_float,
                         parse_int=parse_int, parse_constant=parse_constant,
                         strict=strict)

    @classmethod
    def object_hook(cls, code):
        try:
            qualname = code["__class__"].split(".")
            module = importlib.import_module(".".join(qualname[:-1]))
            constructor = getattr(module, qualname[-1])
        except ModuleNotFoundError as e:
            raise  # TODO!
        except KeyError:
            return code

        for key, func in cls._registry:
            if issubclass(constructor, key):
                code.pop("__class__", None)
                obj = func(constructor, code.pop("value"), code)
                break
        else:
            obj = code

        return obj

    @classmethod
    def register_decoding(cls, type):
        def register(func):
            # inserting subclasses before parent classes, so encountered first
            for i, (key, _) in enumerate(cls._registry):
                if issubclass(type, key):
                    cls._registry.insert(i, (type, func))
                    break
            else:  # put at the end
                cls._registry.append((type, func))
            return func

        return register


###############################################################################
# Register Encodings


def _json_base_encode(obj):
    qualname = obj.__class__.__module__ + "." + obj.__class__.__qualname__
    code = {"__class__": qualname}
    return code

# -------------------------------------------------------------------
# Builtin

@JSONExtendedEncoder.register_encoding(bytes)
def _encode_bytes(obj):
    code = _json_base_encode(obj)
    code.update(value=obj.decode("utf-8"))
    return code


@JSONExtendedEncoder.register_encoding(complex)
def _encode_complex(obj):
    code = _json_base_encode(obj)
    code.update(value=[obj.real, obj.imag])
    return code


@JSONExtendedEncoder.register_encoding(set)
def _encode_set(obj):
    code = _json_base_encode(obj)
    code.update(value=list(obj))
    return code


# -------------------------------------------------------------------
# NumPy

@JSONExtendedEncoder.register_encoding(np.number)
def _encode_numpy_number(obj):
    code = _json_base_encode(obj)
    code.update(value=obj.tolist())
    return code


@JSONExtendedEncoder.register_encoding(np.dtype)
def _encode_numpy_dtype(obj):
    code = {"__class__": "numpy.dtype", "value": str(obj)}
    # TODO! more sophisticated encoding for structured dtype
    return code


@JSONExtendedDecoder.register_decoding(np.dtype)
def _decode_numpy_dtype(constructor, value, code):
    return constructor(value, **code)


@JSONExtendedEncoder.register_encoding(np.ndarray)
def _encode_ndarray(obj):
    code = _json_base_encode(obj)
    code.update(value=obj.tolist(), dtype=_encode_numpy_dtype(obj.dtype))
    return code


@JSONExtendedDecoder.register_decoding(np.ndarray)
def _decode_ndarray(constructor, value, code):
    if constructor is np.ndarray:
        constructor = np.array
    return constructor(value, **code)


###############################################################################
# Load All Entry Points

load_all_entry_points('astropy_io_json_extensions')
