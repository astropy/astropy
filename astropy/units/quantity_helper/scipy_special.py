# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Quantity helpers for the scipy.special ufuncs.

Available ufuncs in this module are at
https://docs.scipy.org/doc/scipy/reference/special.html
"""

from astropy.units.core import dimensionless_unscaled
from astropy.units.errors import UnitsError, UnitTypeError

from . import UFUNC_HELPERS
from .helpers import (
    get_converter,
    helper_cbrt,
    helper_dimensionless_to_dimensionless,
    helper_invariant,
    helper_two_arg_dimensionless,
)

# ufuncs that take dimensionless input and give dimensionless output (nin=1, nout=1)
dimensionless_to_dimensionless_sps_ufuncs = (
    "erf", "erfc", "erfcx", "erfi", "erfinv", "erfcinv",
    "gamma", "gammaln", "loggamma", "gammasgn", "psi", "rgamma", "digamma",
    "wofz", "dawsn", "entr", "exprel", "expm1", "log1p", "exp2", "exp10",
    "j0", "j1", "y0", "y1", "i0", "i0e", "i1", "i1e",
    "k0", "k0e", "k1", "k1e",
    "ndtr", "ndtri",
    # New additions:
    "bei", "beip", "ber", "berp",
    "cosm1",
    "ellipe", "ellipk", "ellipkm1",
    "exp1", "expi", "expit",
    "it2struve0", "itmodstruve0", "itstruve0",
    "kei", "keip", "ker", "kerp",
    "kolmogi", "kolmogorov",
    "log_expit", "log_ndtr", "logit",
    "ndtri_exp",
    "spence",
    "wrightomega",
    "zetac",
)  # fmt: skip

# ufuncs that require input in degrees and give dimensionless output
degree_to_dimensionless_sps_ufuncs = ("cosdg", "sindg", "tandg", "cotdg")

# ufuncs with nin=2, nout=1 that require dimensionless inputs
two_arg_dimensionless_sps_ufuncs = (
    "jv", "jn", "jve", "yn", "yv", "yve", "kn", "kv", "kve", "iv", "ive",
    "hankel1", "hankel1e", "hankel2", "hankel2e",
    # New additions:
    "agm",
    "beta", "betaln", "binom",
    "boxcox", "boxcox1p",
    "chdtr", "chdtrc", "chdtri", "chdtriv",
    "ellipeinc", "ellipkinc", "elliprc",
    "eval_chebyc", "eval_chebys", "eval_chebyt", "eval_chebyu",
    "eval_hermite", "eval_hermitenorm",
    "eval_laguerre", "eval_legendre",
    "eval_sh_chebyt", "eval_sh_chebyu", "eval_sh_legendre",
    "expn",
    "gammainc", "gammaincc", "gammainccinv", "gammaincinv",
    "huber",
    "hyp0f1",
    "inv_boxcox", "inv_boxcox1p",
    "kl_div",
    "mathieu_a", "mathieu_b",
    "modstruve",
    "owens_t",
    "pdtr", "pdtrc", "pdtri", "pdtrik",
    "poch", "powm1",
    "pseudo_huber",
    "rel_entr",
    "smirnov", "smirnovi",
    "stdtr", "stdtridf", "stdtrit",
    "struve",
    "tklmbda",
    "xlog1py", "xlogy",
)  # fmt: skip

# ufuncs with nin=3, nout=1 that require dimensionless inputs
three_arg_dimensionless_sps_ufuncs = (
    "bdtr", "bdtrc", "bdtri", "bdtrik", "bdtrin",
    "besselpoly",
    "betainc", "betaincc", "betainccinv", "betaincinv",
    "btdtria", "btdtrib",
    "chndtr", "chndtridf", "chndtrinc", "chndtrix",
    "elliprd", "elliprf", "elliprg",
    "eval_gegenbauer", "eval_genlaguerre",
    "fdtr", "fdtrc", "fdtri", "fdtridfd",
    "gdtr", "gdtrc", "gdtria", "gdtrib", "gdtrix",
    "hyp1f1", "hyperu",
    "log_wright_bessel",
    "lpmv",
    "nbdtr", "nbdtrc", "nbdtri", "nbdtrik", "nbdtrin",
    "nctdtr", "nctdtridf", "nctdtrinc", "nctdtrit",
    "nrdtrimn", "nrdtrisd",
    "obl_cv", "pro_cv",
    "voigt_profile",
    "wright_bessel",
)  # fmt: skip

# ufuncs with nin=4, nout=1 that require dimensionless inputs
four_arg_dimensionless_sps_ufuncs = (
    "elliprj",
    "eval_jacobi", "eval_sh_jacobi",
    "hyp2f1",
    "ncfdtr", "ncfdtri", "ncfdtridfd", "ncfdtridfn", "ncfdtrinc",
)  # fmt: skip

# ufuncs with nin=1, nout=2 that require dimensionless input
one_arg_two_out_sps_ufuncs = (
    "fresnel",
    "it2i0k0", "it2j0y0", "iti0k0", "itj0y0",
    "modfresnelm", "modfresnelp",
    "shichi", "sici",
)  # fmt: skip

# ufuncs with nin=1, nout=4 that require dimensionless input
one_arg_four_out_sps_ufuncs = (
    "airy", "airye",
    "itairy",
    "kelvin",
)  # fmt: skip

# ufuncs with nin=2, nout=2 that require dimensionless inputs
two_arg_two_out_sps_ufuncs = (
    "pbdv", "pbvv", "pbwa",
)  # fmt: skip

# ufuncs with nin=2, nout=4 that require dimensionless inputs
two_arg_four_out_sps_ufuncs = (
    "ellipj",
)  # fmt: skip

# ufuncs with nin=3, nout=2 that require dimensionless inputs
three_arg_two_out_sps_ufuncs = (
    "mathieu_cem",
    "mathieu_modcem1", "mathieu_modcem2",
    "mathieu_modsem1", "mathieu_modsem2",
    "mathieu_sem",
)  # fmt: skip

# ufuncs with nin=4, nout=2 that require dimensionless inputs
four_arg_two_out_sps_ufuncs = (
    "obl_ang1", "obl_rad1", "obl_rad2",
    "pro_ang1", "pro_rad1", "pro_rad2",
)  # fmt: skip

# ufuncs with nin=5, nout=2 that require dimensionless inputs
five_arg_two_out_sps_ufuncs = (
    "obl_ang1_cv", "obl_rad1_cv", "obl_rad2_cv",
    "pro_ang1_cv", "pro_rad1_cv", "pro_rad2_cv",
)  # fmt: skip


# Collect all scipy_special ufunc names for registration
scipy_special_ufuncs = (
    dimensionless_to_dimensionless_sps_ufuncs
    + degree_to_dimensionless_sps_ufuncs
    + two_arg_dimensionless_sps_ufuncs
    + three_arg_dimensionless_sps_ufuncs
    + four_arg_dimensionless_sps_ufuncs
    + one_arg_two_out_sps_ufuncs
    + one_arg_four_out_sps_ufuncs
    + two_arg_two_out_sps_ufuncs
    + two_arg_four_out_sps_ufuncs
    + three_arg_two_out_sps_ufuncs
    + four_arg_two_out_sps_ufuncs
    + five_arg_two_out_sps_ufuncs
    + ("cbrt", "radian", "round")
)


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


# --- New multi-argument and multi-output helpers ---


def _get_dim_converter(unit):
    """Get a converter to dimensionless for a single unit, or None if already."""
    if unit is None:
        return None
    return get_converter(unit, dimensionless_unscaled)


def helper_three_arg_dimensionless(f, unit1, unit2, unit3):
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2, unit3)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, dimensionless_unscaled


def helper_four_arg_dimensionless(f, unit1, unit2, unit3, unit4):
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2, unit3, unit4)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, dimensionless_unscaled


def helper_dimensionless_to_dimensionless_2out(f, unit):
    if unit is None:
        return [None], (dimensionless_unscaled, dimensionless_unscaled)
    try:
        return (
            [get_converter(unit, dimensionless_unscaled)],
            (dimensionless_unscaled, dimensionless_unscaled),
        )
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )


def helper_dimensionless_to_dimensionless_4out(f, unit):
    out_units = (dimensionless_unscaled,) * 4
    if unit is None:
        return [None], out_units
    try:
        return [get_converter(unit, dimensionless_unscaled)], out_units
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )


def helper_two_arg_dimensionless_2out(f, unit1, unit2):
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled, dimensionless_unscaled)


def helper_two_arg_dimensionless_4out(f, unit1, unit2):
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled,) * 4


def helper_three_arg_dimensionless_2out(f, unit1, unit2, unit3):
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2, unit3)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled, dimensionless_unscaled)


def helper_four_arg_dimensionless_2out(f, unit1, unit2, unit3, unit4):
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2, unit3, unit4)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled, dimensionless_unscaled)


def helper_five_arg_dimensionless_2out(f, unit1, unit2, unit3, unit4, unit5):
    try:
        converters = [
            _get_dim_converter(u) for u in (unit1, unit2, unit3, unit4, unit5)
        ]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled, dimensionless_unscaled)


def get_scipy_special_helpers():
    import scipy.special as sps

    SCIPY_HELPERS = {}

    # nin=1, nout=1 dimensionless
    for name in dimensionless_to_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_dimensionless_to_dimensionless

    # nin=1, degree input
    for name in degree_to_dimensionless_sps_ufuncs:
        SCIPY_HELPERS[getattr(sps, name)] = helper_degree_to_dimensionless

    # nin=2, nout=1 dimensionless
    for name in two_arg_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_two_arg_dimensionless

    # nin=3, nout=1 dimensionless
    for name in three_arg_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_three_arg_dimensionless

    # nin=4, nout=1 dimensionless
    for name in four_arg_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_four_arg_dimensionless

    # nin=1, nout=2 dimensionless
    for name in one_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_dimensionless_to_dimensionless_2out

    # nin=1, nout=4 dimensionless
    for name in one_arg_four_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_dimensionless_to_dimensionless_4out

    # nin=2, nout=2 dimensionless
    for name in two_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_two_arg_dimensionless_2out

    # nin=2, nout=4 dimensionless
    for name in two_arg_four_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_two_arg_dimensionless_4out

    # nin=3, nout=2 dimensionless
    for name in three_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_three_arg_dimensionless_2out

    # nin=4, nout=2 dimensionless
    for name in four_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_four_arg_dimensionless_2out

    # nin=5, nout=2 dimensionless
    for name in five_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_five_arg_dimensionless_2out

    # ufuncs handled as special cases
    SCIPY_HELPERS[sps.cbrt] = helper_cbrt
    SCIPY_HELPERS[sps.radian] = helper_degree_minute_second_to_radian
    # round preserves the input unit, like np.rint
    if hasattr(sps, "round"):
        SCIPY_HELPERS[sps.round] = helper_invariant

    return SCIPY_HELPERS


UFUNC_HELPERS.register_module(
    "scipy.special", scipy_special_ufuncs, get_scipy_special_helpers
)
