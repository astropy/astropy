# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing Cosmology objects with JSON.
"""

from astropy.cosmology.core import Cosmology
from astropy.io.misc.json import _base_encode, JSONExtendedEncoder, JSONExtendedDecoder


@JSONExtendedEncoder.register_encoding(Cosmology)
def _encode_cosmology(cosmo):
    code = _base_encode(cosmo)

    map = cosmo.to_format("mapping")
    map.pop("cosmology")
    meta = map.pop("meta")

    code["value"] = map
    code["meta"] = meta

    return code


@JSONExtendedDecoder.register_decoding(Cosmology)
def _decode_cosmology(cosmo_cls, parameters, meta):
    map = {**parameters, "meta": meta}
    return Cosmology.from_format(map, format="mapping", cosmology=cosmo_cls)
