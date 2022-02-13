# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing python builtins via JSON.
"""

from .core import JSONExtendedDecoder, JSONExtendedEncoder, _json_base_encode


def register_json_extended():
    """:mod:`astropy.io.misc.json` entry points for :mod:`builtins`."""
    JSONExtendedEncoder.register_encoding(bytes)(encode_bytes)
    JSONExtendedDecoder.register_decoding(bytes)(decode_bytes)

    JSONExtendedEncoder.register_encoding(complex)(encode_complex)
    JSONExtendedDecoder.register_decoding(complex)(decode_complex)

    JSONExtendedEncoder.register_encoding(set)(encode_set)
    JSONExtendedDecoder.register_decoding(set)(decode_set)

    JSONExtendedEncoder.register_encoding(type(NotImplemented))(encode_NotImplemented)
    JSONExtendedDecoder.register_decoding(type(NotImplemented))(decode_NotImplemented)

    JSONExtendedEncoder.register_encoding(type(Ellipsis))(encode_Ellipsis)
    JSONExtendedDecoder.register_decoding(type(Ellipsis))(decode_Ellipsis)


# -------------------------------------------------------------------


def encode_bytes(obj):
    """Return `bytes` as a JSON-able dictionary."""
    code = _json_base_encode(obj)
    code.update(value=obj.decode("utf-8"))
    return code


def decode_bytes(constructor, value, code):
    """Return a `bytes` from an ``encode_bytes`` dictionary."""
    return constructor(value, "utf-8")


def encode_complex(obj):
    """Return `complex` as a JSON-able dictionary."""
    code = _json_base_encode(obj)
    code.update(value=[obj.real, obj.imag])
    return code


def decode_complex(constructor, value, code):
    """Return a `complex` from an ``encode_complex`` dictionary."""
    return constructor(*value)


def encode_set(obj):
    """Return `set` as a JSON-able dictionary."""
    code = _json_base_encode(obj)
    code.update(value=list(obj))
    return code


def decode_set(constructor, value, code):
    """Return a `set` from an ``encode_set`` dictionary."""
    return constructor(value)


def encode_NotImplemented(obj):
    """Return `NotImplemented` as a JSON-able dictionary."""
    code = {"!": "builtins.NotImplemented", "value": str(obj)}
    return code


def decode_NotImplemented(constructor, value, code):
    """Return a `NotImplemented` from an ``encode_NotImplemented`` dictionary."""
    return NotImplemented


def encode_Ellipsis(obj):
    """Return `Ellipsis` as a JSON-able dictionary."""
    code = {"!": "builtins.Ellipsis", "value": str(obj)}
    return code


def decode_Ellipsis(constructor, value, code):
    """Return a `Ellipsis` from an ``encode_Ellipsis`` dictionary."""
    return Ellipsis
