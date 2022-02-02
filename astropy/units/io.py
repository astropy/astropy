# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""I/O for :mod:`astropy.units`"""

import astropy.units as u

__all__ = []  # Nothing is publicly scoped.


def register_json_extended():
    """Units entry points for JSONExtendedEncoder, JSONExtendedDecoder."""
    from astropy.io.misc.json import JSONExtendedDecoder, JSONExtendedEncoder

    # Unit
    JSONExtendedEncoder.register_encoding(u.UnitBase)(json_encode_unit)
    JSONExtendedDecoder.register_decoding(u.UnitBase)(json_decode_unit)
    
    JSONExtendedEncoder.register_encoding(u.FunctionUnitBase)(json_encode_unit)
    JSONExtendedDecoder.register_decoding(u.FunctionUnitBase)(json_decode_unit)

    JSONExtendedEncoder.register_encoding(u.StructuredUnit)(json_encode_unit)
    JSONExtendedDecoder.register_decoding(u.StructuredUnit)(json_decode_unit)

    JSONExtendedEncoder.register_encoding(u.CompositeUnit)(json_encode_composite_unit)
    JSONExtendedDecoder.register_decoding(u.CompositeUnit)(json_decode_composite_unit)

    # Quantity
    JSONExtendedEncoder.register_encoding(u.Quantity)(json_encode_quantity)
    JSONExtendedDecoder.register_decoding(u.Quantity)(json_decode_quantity)


def json_encode_unit(obj):  # FIXME so works with units defined outside units subpkg
    """Return `astropy.unit.Unit` as a JSON-able dictionary."""
    from astropy.io.misc.json import _json_base_encode

    if isinstance(obj, u.CompositeUnit):
        return json_encode_composite_unit(obj)

    code = _json_base_encode(obj)
    code.update(value=obj.to_string())
    return code


def json_decode_unit(constructor, value, code):
    """Return a |Unit| from an ``json_encode_unit`` dictionary."""
    if value == "dimensionless":
        value = ""
    return constructor(value, **code)


def json_encode_composite_unit(obj):
    """Return `astropy.unit.CompositeUnit` as a JSON-able dictionary."""
    from astropy.io.misc.json import _json_base_encode
    code = _json_base_encode(obj)

    if obj == u.dimensionless_unscaled:
        code.update(__class__= "astropy.units.core.Unit", value="dimensionless")
    else:
        code.update(value=None, scale=obj.scale, bases=obj.bases, powers=obj.powers)
    return code


def json_decode_composite_unit(constructor, value, code):
    """Return a `astropy.unit.CompositeUnit` from an ``json_encode_composite_unit`` dictionary."""
    return constructor(**code)


def json_encode_quantity(obj):
    """Return a |Quantity| as a JSON-able dictionary."""
    from astropy.io.misc.json.core import _json_base_encode
    from astropy.io.misc.json.numpy import encode_ndarray

    code = _json_base_encode(obj)
    value = encode_ndarray(obj.value)
    code["value"] = value

    # code["dtype"] = code["value"].pop("dtype") # move up a level

    unit = json_encode_unit(obj.unit)
    if (obj.unit == u.one or not isinstance(obj.unit, u.CompositeUnit)):
        code["unit"] = unit["value"]
    else:
        code["unit"] = unit

    return code


def json_decode_quantity(constructor, value, code):
    """Return a |Quantity| from an ``json_encode_quantity`` dictionary."""
    code["unit"] = u.Unit(code["unit"]) if code["unit"] != "dimensionless" else u.one
    return constructor(value, **code)
