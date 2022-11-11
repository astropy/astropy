# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Quantity helpers for the scipy.special ufuncs.

Available ufuncs in this module are at
https://docs.scipy.org/doc/scipy/reference/special.html
"""
import numpy as np

from astropy.units.core import UnitsError, UnitTypeError, dimensionless_unscaled

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


def get_scipy_special_helpers():
    import scipy.special as sps

    SCIPY_HELPERS = {}
    for name in dimensionless_to_dimensionless_sps_ufuncs:
        # In SCIPY_LT_1_5, erfinv and erfcinv are not ufuncs.
        ufunc = getattr(sps, name, None)
        if isinstance(ufunc, np.ufunc):
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
