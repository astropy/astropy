# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing Cosmology objects with JSON.
"""

from astropy.cosmology.core import Cosmology


def encode_cosmology(cosmo):
    from astropy.io.misc.json.core import _json_base_encode

    code = _json_base_encode(cosmo)  # get type code for decoding

    map = cosmo.to_format("mapping")  # dictionary representation
    map.pop("cosmology")  # extraneous given `_json_base_encode`
    meta = map.pop("meta")  # going to move to sub-dict of `code`

    code["value"] = map  # parameters + name
    code["meta"] = meta  # everything else

    return code


def decode_cosmology(cosmo_cls, parameters, meta):
    map = {**parameters, "meta": meta}  # mix back to one dict
    return Cosmology.from_format(map, format="mapping", cosmology=cosmo_cls)


def register_json_extended():
    r"""Entry point for JSONExtendedEncoder, JSONExtendedDecoder."""
    from astropy.io.misc.json.core import JSONExtendedEncoder, JSONExtendedDecoder

    JSONExtendedEncoder.register_encoding(Cosmology)(encode_cosmology)
    JSONExtendedDecoder.register_decoding(Cosmology)(decode_cosmology)
