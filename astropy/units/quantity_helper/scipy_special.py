# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Quantity helpers for the scipy.special ufuncs.

Available ufuncs in this module are at
https://docs.scipy.org/doc/scipy/reference/special.html
"""

import os

from astropy.units.core import dimensionless_unscaled
from astropy.units.errors import UnitsError, UnitTypeError
from astropy.utils import minversion

from . import UFUNC_HELPERS
from .helpers import (
    get_converter,
    helper_cbrt,
    helper_dimensionless_to_dimensionless,
    helper_two_arg_dimensionless,
)

dimensionless_to_dimensionless_sps_ufuncs = (
    "erf", "erfc", "erfcx", "erfi", "erfinv", "erfcinv",
    "gamma", "gammaln", "loggamma", "gammasgn", "psi", "rgamma", "digamma",
    "wofz", "dawsn", "entr", "exprel", "expm1", "log1p", "exp2", "exp10",
    "j0", "j1", "y0", "y1", "i0", "i0e", "i1", "i1e",
    "k0", "k0e", "k1", "k1e", "itj0y0", "it2j0y0", "iti0k0", "it2i0k0",
    "ndtr", "ndtri",
)  # fmt: skip


scipy_special_ufuncs = dimensionless_to_dimensionless_sps_ufuncs
# ufuncs that require input in degrees and give dimensionless output.
degree_to_dimensionless_sps_ufuncs = ("cosdg", "sindg", "tandg", "cotdg")
scipy_special_ufuncs += degree_to_dimensionless_sps_ufuncs
two_arg_dimensionless_sps_ufuncs = (
    "jv", "jn", "jve", "yn", "yv", "yve", "kn", "kv", "kve", "iv", "ive",
    "hankel1", "hankel1e", "hankel2", "hankel2e",
)  # fmt: skip
scipy_special_ufuncs += two_arg_dimensionless_sps_ufuncs
# ufuncs handled as special cases
scipy_special_ufuncs += ("cbrt", "radian")


def helper_degree_to_dimensionless(f, unit):
    from astropy.units.si import degree

    try:
        return [get_converter(unit, degree)], dimensionless_unscaled
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to quantities with angle units"
        )


def helper_degree_minute_second_to_radian(f, unit1, unit2, unit3):
    from astropy.units.si import arcmin, arcsec, degree, radian

    try:
        return [
            get_converter(unit1, degree),
            get_converter(unit2, arcmin),
            get_converter(unit3, arcsec),
        ], radian
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to quantities with angle units"
        )


def _has_scipy_array_api() -> bool:
    if minversion("scipy", "2.0.0.dev0"):  # SCIPY_LT_2_0
        return True

    if not minversion("scipy", "1.12.0"):
        # at the time of writing, astropy requires scipy>=1.13 and Python>=3.12
        # but versions as early as 1.11.2 are installable against Python 3.12
        return False

    SCIPY_ARRAY_API: bool | str = os.getenv("SCIPY_ARRAY_API", "0")
    if isinstance(SCIPY_ARRAY_API, bool):
        return SCIPY_ARRAY_API

    return SCIPY_ARRAY_API.lower() in {"1", "true", "all"}


def get_scipy_special_helpers():
    # On scipy >= 2.0, as well as earlier versions with
    # SCIPY_ARRAY_API=1 set, sps functions are wrappers that call different
    # functions depending on the input array type.  Quantity are based on
    # numpy arrays, so Quantity.__array_ufunc__ will get passed an entry in
    # sps._ufuncs, which is what UFUNC_HELPERS will be keyed with.
    if _has_scipy_array_api():
        import scipy.special._ufuncs as sps
    else:
        import scipy.special as sps

    SCIPY_HELPERS = {}
    for name in dimensionless_to_dimensionless_sps_ufuncs:
        # Some names are aliases defined only in the scipy.special namespace
        # (e.g., digamma for psi), so are not present in sps._ufuncs.
        if (ufunc := getattr(sps, name, None)) is not None:
            SCIPY_HELPERS[ufunc] = helper_dimensionless_to_dimensionless

    for ufunc in degree_to_dimensionless_sps_ufuncs:
        SCIPY_HELPERS[getattr(sps, ufunc)] = helper_degree_to_dimensionless

    for ufunc in two_arg_dimensionless_sps_ufuncs:
        SCIPY_HELPERS[getattr(sps, ufunc)] = helper_two_arg_dimensionless

    # ufuncs handled as special cases
    SCIPY_HELPERS[sps.cbrt] = helper_cbrt
    SCIPY_HELPERS[sps.radian] = helper_degree_minute_second_to_radian
    return SCIPY_HELPERS


UFUNC_HELPERS.register_module(
    "scipy.special", scipy_special_ufuncs, get_scipy_special_helpers
)
