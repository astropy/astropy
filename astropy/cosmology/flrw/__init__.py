# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy FLRW classes."""

__all__ = []  # No public API -- see cosmology/__init__.py


def __getattr__(attr):
    """Lazy import deprecated private API."""
    base_attrs = {
        "H0units_to_invs": "_H0units_to_invs",
        "a_B_c2": "_a_B_c2",
        "critdens_const": "_critdens_const",
        "kB_evK": "_kB_evK",
        "radian_in_arcmin": "_radian_in_arcmin",
        "radian_in_arcsec": "_radian_in_arcsec",
        "sec_to_Gyr": "_sec_to_Gyr",
        "quad": "quad",
    }
    lambdacdm_attrs = {
        "ellipkinc": "ellipkinc",
        "hyp2f1": "hyp2f1",
    }

    if attr in base_attrs or attr in lambdacdm_attrs:
        import warnings

        from astropy.utils.exceptions import AstropyDeprecationWarning

        msg = (
            f"`astropy.cosmology.flrw.{attr}` is a private variable (since "
            "v5.1) and in future will raise an exception."
        )
        warnings.warn(msg, AstropyDeprecationWarning)

        if attr in base_attrs:
            from . import base

            return getattr(base, base_attrs[attr])
        elif attr in lambdacdm_attrs:
            from . import lambdacdm

            return getattr(lambdacdm, lambdacdm_attrs[attr])

    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}.")
