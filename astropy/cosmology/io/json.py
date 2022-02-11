# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for serializing Cosmology objects with JSON.
"""

from astropy.cosmology.core import Cosmology

__all__ = []  # Nothing is publicly scoped


def json_encode_cosmology(cosmo):
    """Return |Cosmology| as a JSON-able dictionary."""
    from astropy.io.misc.json.core import _json_base_encode

    code = _json_base_encode(cosmo)  # get type code for decoding

    map = cosmo.to_format("mapping")  # dictionary representation
    map.pop("cosmology")  # extraneous given `_json_base_encode`
    meta = map.pop("meta")  # going to move to sub-dict of `code`

    code["value"] = map  # parameters + name
    code["meta"] = meta  # everything else

    return code


def json_decode_cosmology(cosmo_cls, parameters, code):
    """Return a |Cosmology| from an ``json_encode_cosmology`` dictionary."""
    map = {**parameters, "meta": code.pop("meta", None)}  # mix back to one dict
    return Cosmology.from_format(map, format="mapping", cosmology=cosmo_cls)


def register_json_extended():
    """|Cosmology| entry point for JSONExtendedEncoder, JSONExtendedDecoder."""
    from astropy.io.misc.json.core import JSONExtendedEncoder, JSONExtendedDecoder
    from astropy.io.misc.json.core import QUALNAME_SUBSTITUTIONS as QS

    JSONExtendedEncoder.register_encoding(Cosmology)(json_encode_cosmology)
    JSONExtendedDecoder.register_decoding(Cosmology)(json_decode_cosmology)

    # TODO! automate discovery of Cosmology subclasses
    QS["astropy.cosmology.core.Cosmology"] = "astropy.cosmology.Cosmology"
    QS["astropy.cosmology.flrw.FLRW"] = "astropy.cosmology.FLRW"
    QS["astropy.cosmology.flrw.LambdaCDM"] = "astropy.cosmology.LambdaCDM"
    QS["astropy.cosmology.flrw.FlatLambdaCDM"] = "astropy.cosmology.FlatLambdaCDM"
    QS["astropy.cosmology.flrw.wCDM"] = "astropy.cosmology.wCDM"
    QS["astropy.cosmology.flrw.FlatwCDM"] = "astropy.cosmology.FlatwCDM"
    QS["astropy.cosmology.flrw.w0waCDM"] = "astropy.cosmology.w0waCDM"
    QS["astropy.cosmology.flrw.Flatw0waCDM"] = "astropy.cosmology.Flatw0waCDM"
    QS["astropy.cosmology.flrw.wpwaCDM"] = "astropy.cosmology.wpwaCDM"
    QS["astropy.cosmology.flrw.w0wzCDM"] = "astropy.cosmology.w0wzCDM"
