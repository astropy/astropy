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

    def __init__(self, *, parse_float=None, parse_int=None, parse_constant=None, strict=True):
        super().__init__(object_hook=self.object_hook, parse_float=parse_float,
                         parse_int=parse_int, parse_constant=parse_constant,
                         strict=strict)

    @classmethod
    def object_hook(cls, code):
        try:
            qualname = code.pop("__class__").split(".")
            module = importlib.import_module(".".join(qualname[:-1]))
            constructor = getattr(module, qualname[-1])
        except ModuleNotFoundError as e:
            raise  # TODO!

        for key, func in cls._registry:
            print(constructor, key)
            if issubclass(constructor, key):
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


def _base_encode(obj):
    qualname = obj.__class__.__module__ + "." + obj.__class__.__qualname__
    code = {"__class__": qualname}
    return code


@JSONExtendedEncoder.register_encoding(bytes)
def _encode_bytes(obj):
    code = _base_encode(obj)
    code.update(value=obj.decode("utf-8"))
    return code


@JSONExtendedEncoder.register_encoding(complex)
def _encode_complex(obj):
    code = _base_encode(obj)
    code.update(value=[obj.real, obj.imag])
    return code


@JSONExtendedEncoder.register_encoding(set)
def _encode_set(obj):
    code = _base_encode(obj)
    code.update(value=list(obj))
    return code


@JSONExtendedEncoder.register_encoding(np.number)
def _encode_numpy_number(obj):
    code = _base_encode(obj)
    code.update(value=obj.tolist())
    return code


@JSONExtendedEncoder.register_encoding(np.ndarray)
def _encode_ndarray(obj):
    code = _base_encode(obj)
    code.update(value=obj.tolist(), dtype=str(obj.dtype))  # TODO! encode dtype
    return code


@JSONExtendedEncoder.register_encoding(u.FunctionUnitBase)
@JSONExtendedEncoder.register_encoding(u.UnitBase)
def _encode_unit(obj):
    code = _base_encode(obj)
    if obj == u.dimensionless_unscaled:
        code.update(value="dimensionless_unit")
    else:
        code.update(value=obj.to_string())
    return code


@JSONExtendedEncoder.register_encoding(u.Quantity)
def _encode_quantity(obj):
    code = _base_encode(obj)
    code.update(value=obj.value.tolist(), dtype=str(obj.dtype))  # TODO! encode dtype
    code["unit"] = obj.unit.to_string()
    return code


@JSONExtendedEncoder.register_encoding(coord.Longitude)
def _encode_longitude(obj):
    code = _encode_quantity(obj)
    code["wrap_angle"] = _encode_quantity(obj.wrap_angle)
    return code


###############################################################################
# Register Decodings


@JSONExtendedDecoder.register_decoding(np.ndarray)
def _decode_ndarray(constructor, value, code):
    if constructor is np.ndarray:
        constructor = np.array
    return constructor(value, **code)


@JSONExtendedDecoder.register_decoding(u.Quantity)
def _decode_quantity(constructor, value, code):
    return constructor(value, **code)
