# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing core astropy objects via the
JSON protocol.
"""

import importlib
import inspect
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
            if isinstance(constructor, key) or (inspect.isclass(constructor) and issubclass(constructor, key)):
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


def _json_base_encode(obj):
    qualname = obj.__class__.__module__ + "." + obj.__class__.__qualname__
    code = {"__class__": qualname}
    return code
