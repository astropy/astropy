# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""I/O for :mod:`astropy.units`"""

import astropy.units as u

__all__ = []  # Nothing is publicly scoped.


def register_json_extended():
    """Units entry points for JSONExtendedEncoder, JSONExtendedDecoder."""
    from astropy.io.misc.json import JSONExtendedDecoder, JSONExtendedEncoder
    from astropy.io.misc.json.core import QUALNAME_SUBSTITUTIONS as QS

    # Unit
    JSONExtendedEncoder.register_encoding(u.UnitBase)(json_encode_unit)
    JSONExtendedDecoder.register_decoding(u.UnitBase)(json_decode_unit)

    JSONExtendedEncoder.register_encoding(u.FunctionUnitBase)(json_encode_unit)
    JSONExtendedDecoder.register_decoding(u.FunctionUnitBase)(json_decode_unit)

    JSONExtendedEncoder.register_encoding(u.StructuredUnit)(json_encode_unit)
    JSONExtendedDecoder.register_decoding(u.StructuredUnit)(json_decode_unit)

    JSONExtendedEncoder.register_encoding(u.CompositeUnit)(json_encode_unit)
    JSONExtendedDecoder.register_decoding(u.CompositeUnit)(json_decode_composite_unit)

    QS["astropy.units.core.Unit"] = "astropy.units.Unit"
    QS["astropy.units.core.CompositeUnit"] = "astropy.units.CompositeUnit"
    QS["astropy.units.core.PrefixUnit"] = "astropy.units.PrefixUnit"
    QS["astropy.units.structured.StructuredUnit"] = "astropy.units.StructuredUnit"

    # Quantity
    JSONExtendedEncoder.register_encoding(u.Quantity)(json_encode_quantity)
    JSONExtendedDecoder.register_decoding(u.Quantity)(json_decode_quantity)

    QS["astropy.units.quantity.Quantity"] = "astropy.units.Quantity"


def json_encode_unit(obj):  # FIXME so works with units defined outside units subpkg
    """Return `astropy.unit.Unit` as a JSON-able dictionary."""
    from astropy.io.misc.json import _json_base_encode

    code = _json_base_encode(obj)

    if obj == u.dimensionless_unscaled:
        code["!"] = "astropy.units.core.Unit"  # TODO? should it remain Composite?
        code["value"] = "dimensionless"
    else:
        code["value"] = obj.to_string()

    return code


def json_decode_unit(constructor, value, code):
    """Return a |Unit| from an ``json_encode_unit`` dictionary."""
    if value == "dimensionless":
        value = ""
    return constructor(value, **code)


def json_decode_composite_unit(constructor, value, code):
    """Return a `astropy.unit.CompositeUnit` from an ``json_encode_composite_unit`` dictionary."""
    return u.Unit(value, **code)


def json_encode_quantity(obj):
    """Return a |Quantity| as a JSON-able dictionary."""
    from astropy.io.misc.json.core import _json_base_encode
    from astropy.io.misc.json.numpy import encode_ndarray

    code = _json_base_encode(obj)
    code["value"] = obj.value

    # code["dtype"] = code["value"].pop("dtype") # move up a level

    # simplify unit representation
    unit = json_encode_unit(obj.unit)
    code["unit"] = unit["value"]

    return code


def json_decode_quantity(constructor, value, code):
    """Return a |Quantity| from an ``json_encode_quantity`` dictionary."""
    code["unit"] = u.Unit(code["unit"]) if code["unit"] != "dimensionless" else u.one
    return constructor(value, **code)
