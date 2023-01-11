# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy FLRW classes."""

from . import base, lambdacdm, w0cdm, w0wacdm, w0wzcdm, wpwazpcdm
from .base import *
from .lambdacdm import *
from .w0cdm import *
from .w0wacdm import *
from .w0wzcdm import *
from .wpwazpcdm import *

__all__ = (
    base.__all__
    + lambdacdm.__all__
    + w0cdm.__all__
    + w0wacdm.__all__
    + wpwazpcdm.__all__
    + w0wzcdm.__all__
)


def __getattr__(attr):
    """Lazy import deprecated private API."""
    base_attrs = (
        "H0units_to_invs",
        "a_B_c2",
        "critdens_const",
        "kB_evK",
        "radian_in_arcmin",
        "radian_in_arcsec",
        "sec_to_Gyr",
    )

    if attr in base_attrs + ("quad",) + ("ellipkinc", "hyp2f1"):
        import warnings

        from astropy.utils.exceptions import AstropyDeprecationWarning

        from . import base, lambdacdm

        msg = (
            f"`astropy.cosmology.flrw.{attr}` is a private variable (since "
            "v5.1) and in future will raise an exception."
        )
        warnings.warn(msg, AstropyDeprecationWarning)

        if attr in base_attrs:
            return getattr(base, "_" + attr)
        elif attr == "quad":
            return getattr(base, attr)
        elif attr in ("ellipkinc", "hyp2f1"):
            return getattr(lambdacdm, attr)

    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}.")
