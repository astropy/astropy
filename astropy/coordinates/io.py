# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""I/O for :mod:`astropy.coordinates`."""

import astropy.coordinates as coord

__all__ = []  # Nothing is publicly scoped.


def json_encode_longitude(obj):
    from astropy.units.io import json_encode_quantity

    code = json_encode_quantity(obj)
    code["wrap_angle"] = json_encode_quantity(obj.wrap_angle)
    return code


def register_json_extended():
    from astropy.io.misc.json import JSONExtendedEncoder

    # Quantity subclasses. Decoding handled by Quantity decoder.
    JSONExtendedEncoder.register_encoding(coord.Longitude)(json_encode_longitude)
