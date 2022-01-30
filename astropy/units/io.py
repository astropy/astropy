# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""I/O for :mod:`astropy.units`"""

import astropy.units as u

__all__ = []  # Nothing is publicly scoped.


def json_encode_unit(obj):  # FIXME so works with units defined outside units subpkg
    """Return `astropy.unit.Unit` as a JSON-able dictionary."""
    from astropy.io.misc.json import _json_base_encode

    code = _json_base_encode(obj)
    if obj == u.dimensionless_unscaled:
        code.update(value="dimensionless_unit")
    else:
        code.update(value=obj.to_string())
    return code


def json_encode_quantity(obj):
    """Return a |Quantity| as a JSON-able dictionary."""
    from astropy.io.misc.json.core import _json_base_encode
    from astropy.io.misc.json.numpy import encode_numpy_dtype

    code = _json_base_encode(obj)
    code.update(value=obj.value.tolist(), dtype=encode_numpy_dtype(obj.dtype))
    code["unit"] = obj.unit.to_string()
    return code


def json_decode_quantity(constructor, value, code):
    """Return a |Quantity| from an ``json_encode_quantity`` dictionary."""
    return constructor(value, **code)


def register_json_extended():
    """Units entry points for JSONExtendedEncoder, JSONExtendedDecoder."""
    from astropy.io.misc.json import JSONExtendedDecoder, JSONExtendedEncoder

    # Unit
    JSONExtendedEncoder.register_encoding(u.UnitBase)(json_encode_unit)
    JSONExtendedEncoder.register_encoding(u.FunctionUnitBase)(json_encode_unit)

    # Quantity
    JSONExtendedEncoder.register_encoding(u.Quantity)(json_encode_quantity)
    JSONExtendedDecoder.register_decoding(u.Quantity)(json_decode_quantity)
