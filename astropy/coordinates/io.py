# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""I/O for :mod:`astropy.coordinates`."""

import astropy.coordinates as coord

__all__ = []  # Nothing is publicly scoped.


def json_encode_longitude(obj):
    """Return a `~astropy.coordinates.Longitude` as a JSON-able dictionary."""
    from astropy.units.io import json_encode_quantity

    code = json_encode_quantity(obj)
    code["wrap_angle"] = obj.wrap_angle
    return code


def register_json_extended():
    """Coordinates entry points for JSONExtendedEncoder, JSONExtendedDecoder."""
    from astropy.io.misc.json import JSONExtendedEncoder
    from astropy.io.misc.json.core import QUALNAME_SUBSTITUTIONS as QS

    # Quantity subclasses. Decoding handled by Quantity decoder.
    JSONExtendedEncoder.register_encoding(coord.Longitude)(json_encode_longitude)

    QS["astropy.coordinates.angles.Angle"] = "astropy.coordinates.Angle"
    QS["astropy.coordinates.angles.Latitude"] = "astropy.coordinates.Latitude"
    QS["astropy.coordinates.angles.Longitude"] = "astropy.coordinates.Longitude"
    QS["astropy.coordinates.distances.Distance"] = "astropy.coordinates.Distance"
