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
    get_converters_and_unit,
    helper_cbrt,
    helper_dimensionless_to_dimensionless,
    helper_invariant,
    helper_radian_to_dimensionless,
    helper_two_arg_dimensionless,
)

# ============================================================================
# Category 1: Dimensionless -> dimensionless (nin=1, nout=1).
#
# Pure mathematical functions where all inputs and outputs are dimensionless.
# Includes error functions, gamma functions, Bessel functions (power series
# in x, so x must be dimensionless), statistical distributions, etc.
# ============================================================================
dimensionless_to_dimensionless_sps_ufuncs = (
    "erf", "erfc", "erfcx", "erfi", "erfinv", "erfcinv",
    "gamma", "gammaln", "loggamma", "gammasgn", "psi", "rgamma", "digamma",
    "wofz", "dawsn", "entr", "exprel", "expm1", "log1p", "exp2", "exp10",
    "j0", "j1", "y0", "y1", "i0", "i0e", "i1", "i1e",
    "k0", "k0e", "k1", "k1e",
    "ndtr", "ndtri",
    "bei", "beip", "ber", "berp",
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

# Dimensionless -> dimensionless (nin=2, nout=1).
two_arg_dimensionless_sps_ufuncs = (
    "jv", "jn", "jve", "yn", "yv", "yve", "kn", "kv", "kve", "iv", "ive",
    "hankel1", "hankel1e", "hankel2", "hankel2e",
    "beta", "betaln", "binom",
    "boxcox", "boxcox1p",
    "chdtr", "chdtrc", "chdtri", "chdtriv",
    "elliprc",
    "eval_chebyc", "eval_chebys", "eval_chebyt", "eval_chebyu",
    "eval_hermite", "eval_hermitenorm",
    "eval_laguerre", "eval_legendre",
    "eval_sh_chebyt", "eval_sh_chebyu", "eval_sh_legendre",
    "expn",
    "gammainc", "gammaincc", "gammainccinv", "gammaincinv",
    # NOTE: huber and pseudo_huber are loss functions. Mathematically,
    # if delta and r share unit U, the output has unit U^2. However,
    # implementing U^2 output for a ufunc helper is impractical, and
    # these functions are almost always used with normalized inputs.
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

# Dimensionless -> dimensionless (nin=3, nout=1).
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

# Dimensionless -> dimensionless (nin=4, nout=1).
four_arg_dimensionless_sps_ufuncs = (
    "elliprj",
    "eval_jacobi", "eval_sh_jacobi",
    "hyp2f1",
    "ncfdtr", "ncfdtri", "ncfdtridfd", "ncfdtridfn", "ncfdtrinc",
)  # fmt: skip

# Dimensionless -> dimensionless (nin=1, nout=2).
one_arg_two_out_sps_ufuncs = (
    "fresnel",
    "it2i0k0", "it2j0y0", "iti0k0", "itj0y0",
    "modfresnelm", "modfresnelp",
    "shichi", "sici",
)  # fmt: skip

# Dimensionless -> dimensionless (nin=1, nout=4).
one_arg_four_out_sps_ufuncs = (
    "airy", "airye",
    "itairy",
    "kelvin",
)  # fmt: skip

# Dimensionless -> dimensionless (nin=2, nout=2).
two_arg_two_out_sps_ufuncs = (
    "pbdv", "pbvv", "pbwa",
)  # fmt: skip

# Dimensionless -> dimensionless (nin=2, nout=4).
two_arg_four_out_sps_ufuncs = (
    "ellipj",
)  # fmt: skip

# Dimensionless -> dimensionless (nin=4, nout=2).
four_arg_two_out_sps_ufuncs = (
    "obl_ang1", "obl_rad1", "obl_rad2",
    "pro_ang1", "pro_rad1", "pro_rad2",
)  # fmt: skip

# Dimensionless -> dimensionless (nin=5, nout=2).
five_arg_two_out_sps_ufuncs = (
    "obl_ang1_cv", "obl_rad1_cv", "obl_rad2_cv",
    "pro_ang1_cv", "pro_rad1_cv", "pro_rad2_cv",
)  # fmt: skip

# ============================================================================
# Category 2: Angle -> dimensionless.
#
# Functions where one or more inputs are angles. Helpers accept any angle
# Quantity and convert to the appropriate unit (radians or degrees) internally.
# ============================================================================

# Ufuncs that require input in degrees and give dimensionless output.
degree_to_dimensionless_sps_ufuncs = ("cosdg", "sindg", "tandg", "cotdg")

# NOTE: cosm1(x) computes cos(x) - 1. The input x is an angle in radians,
# just like numpy's cos. We accept any angle Quantity and convert to radians.
radian_to_dimensionless_sps_ufuncs = ("cosm1",)

# NOTE: ellipeinc(phi, m) and ellipkinc(phi, m) are incomplete elliptic
# integrals where phi is the amplitude angle (in radians) and m is a
# dimensionless parameter. We accept any angle Quantity for phi.
angle_dimensionless_to_dimensionless_sps_ufuncs = ("ellipeinc", "ellipkinc")

# NOTE: mathieu_cem(m, q, x), mathieu_sem(m, q, x), and the modified
# variants mathieu_modcem1/2(m, q, x), mathieu_modsem1/2(m, q, x) all
# take x in degrees per the SciPy documentation. We accept any angle
# Quantity for x and convert to degrees.
dimless_dimless_angle_to_2out_sps_ufuncs = (
    "mathieu_cem", "mathieu_sem",
    "mathieu_modcem1", "mathieu_modcem2",
    "mathieu_modsem1", "mathieu_modsem2",
)  # fmt: skip

# ============================================================================
# Category 3: Arithmetic-like (preserve units).
#
# Functions where the output has the same units as the inputs, like addition
# or averaging operations.
# ============================================================================

# NOTE: agm(a, b) computes the arithmetic-geometric mean, which iterates
# a_{n+1} = (a_n + b_n)/2 and b_{n+1} = sqrt(a_n * b_n). Both operations
# require compatible units, and the result has the same unit as the inputs.
arithmetic_sps_ufuncs = ("agm",)

# ============================================================================
# Category 4: Special cases.
# ============================================================================

# cbrt: cube root, transforms units like sqrt.
# radian: converts degrees/arcmin/arcsec to radians (angle input).
# round: preserves the input unit, like np.rint.

# Collect all scipy_special ufunc names for registration.
scipy_special_ufuncs = (
    dimensionless_to_dimensionless_sps_ufuncs
    + degree_to_dimensionless_sps_ufuncs
    + radian_to_dimensionless_sps_ufuncs
    + two_arg_dimensionless_sps_ufuncs
    + angle_dimensionless_to_dimensionless_sps_ufuncs
    + arithmetic_sps_ufuncs
    + three_arg_dimensionless_sps_ufuncs
    + four_arg_dimensionless_sps_ufuncs
    + one_arg_two_out_sps_ufuncs
    + one_arg_four_out_sps_ufuncs
    + two_arg_two_out_sps_ufuncs
    + two_arg_four_out_sps_ufuncs
    + dimless_dimless_angle_to_2out_sps_ufuncs
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


# --- Multi-argument and multi-output helpers ---


def _get_dim_converter(unit):
    """Return a converter to dimensionless for a single unit, or None if already."""
    if unit is None:
        return None
    return get_converter(unit, dimensionless_unscaled)


def helper_three_arg_dimensionless(f, unit1, unit2, unit3):
    """All three inputs must be dimensionless, output is dimensionless."""
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2, unit3)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, dimensionless_unscaled


def helper_four_arg_dimensionless(f, unit1, unit2, unit3, unit4):
    """All four inputs must be dimensionless, output is dimensionless."""
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2, unit3, unit4)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, dimensionless_unscaled


def helper_dimensionless_to_dimensionless_2out(f, unit):
    """Single dimensionless input, two dimensionless outputs."""
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
    """Single dimensionless input, four dimensionless outputs."""
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
    """Two dimensionless inputs, two dimensionless outputs."""
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled, dimensionless_unscaled)


def helper_two_arg_dimensionless_4out(f, unit1, unit2):
    """Two dimensionless inputs, four dimensionless outputs."""
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled,) * 4


def helper_three_arg_dimensionless_2out(f, unit1, unit2, unit3):
    """Three dimensionless inputs, two dimensionless outputs."""
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2, unit3)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled, dimensionless_unscaled)


def helper_four_arg_dimensionless_2out(f, unit1, unit2, unit3, unit4):
    """Four dimensionless inputs, two dimensionless outputs."""
    try:
        converters = [_get_dim_converter(u) for u in (unit1, unit2, unit3, unit4)]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled, dimensionless_unscaled)


def helper_five_arg_dimensionless_2out(f, unit1, unit2, unit3, unit4, unit5):
    """Five dimensionless inputs, two dimensionless outputs."""
    try:
        converters = [
            _get_dim_converter(u) for u in (unit1, unit2, unit3, unit4, unit5)
        ]
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return converters, (dimensionless_unscaled, dimensionless_unscaled)


# --- Angle-aware helpers ---


def helper_angle_dimensionless_to_dimensionless(f, unit1, unit2):
    """First arg is angle (converted to radians), second is dimensionless."""
    from astropy.units.si import radian

    try:
        conv1 = get_converter(unit1, radian) if unit1 is not None else None
        conv2 = _get_dim_converter(unit2)
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to "
            f"(angle, dimensionless) quantities"
        )
    return [conv1, conv2], dimensionless_unscaled


def helper_dimless_dimless_angle_to_2out(f, unit1, unit2, unit3):
    """First two args dimensionless, third arg is angle (converted to degrees).

    Returns two dimensionless outputs. Used for mathieu_cem and mathieu_sem,
    which expect x in degrees internally.
    """
    from astropy.units.si import degree

    try:
        conv1 = _get_dim_converter(unit1)
        conv2 = _get_dim_converter(unit2)
        conv3 = get_converter(unit3, degree) if unit3 is not None else None
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to "
            f"(dimensionless, dimensionless, angle) quantities"
        )
    return [conv1, conv2, conv3], (dimensionless_unscaled, dimensionless_unscaled)


def get_scipy_special_helpers():
    import scipy.special as sps

    SCIPY_HELPERS = {}

    # Dimensionless -> dimensionless (nin=1, nout=1).
    for name in dimensionless_to_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_dimensionless_to_dimensionless

    # Angle (degrees) -> dimensionless (nin=1, nout=1).
    for name in degree_to_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_degree_to_dimensionless

    # Angle (radians) -> dimensionless (nin=1, nout=1).
    for name in radian_to_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_radian_to_dimensionless

    # Dimensionless -> dimensionless (nin=2, nout=1).
    for name in two_arg_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_two_arg_dimensionless

    # Angle + dimensionless -> dimensionless (nin=2, nout=1).
    for name in angle_dimensionless_to_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_angle_dimensionless_to_dimensionless

    # Arithmetic-like: preserve units (nin=2, nout=1).
    for name in arithmetic_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = get_converters_and_unit

    # Dimensionless -> dimensionless (nin=3, nout=1).
    for name in three_arg_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_three_arg_dimensionless

    # Dimensionless -> dimensionless (nin=4, nout=1).
    for name in four_arg_dimensionless_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_four_arg_dimensionless

    # Dimensionless -> dimensionless (nin=1, nout=2).
    for name in one_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_dimensionless_to_dimensionless_2out

    # Dimensionless -> dimensionless (nin=1, nout=4).
    for name in one_arg_four_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_dimensionless_to_dimensionless_4out

    # Dimensionless -> dimensionless (nin=2, nout=2).
    for name in two_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_two_arg_dimensionless_2out

    # Dimensionless -> dimensionless (nin=2, nout=4).
    for name in two_arg_four_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_two_arg_dimensionless_4out

    # Dimensionless, dimensionless, angle (degrees) -> 2 dimensionless outputs.
    for name in dimless_dimless_angle_to_2out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_dimless_dimless_angle_to_2out

    # Dimensionless -> dimensionless (nin=4, nout=2).
    for name in four_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_four_arg_dimensionless_2out

    # Dimensionless -> dimensionless (nin=5, nout=2).
    for name in five_arg_two_out_sps_ufuncs:
        ufunc = getattr(sps, name, None)
        if ufunc is not None:
            SCIPY_HELPERS[ufunc] = helper_five_arg_dimensionless_2out

    # Special cases.
    SCIPY_HELPERS[sps.cbrt] = helper_cbrt
    SCIPY_HELPERS[sps.radian] = helper_degree_minute_second_to_radian
    if hasattr(sps, "round"):
        SCIPY_HELPERS[sps.round] = helper_invariant

    return SCIPY_HELPERS


UFUNC_HELPERS.register_module(
    "scipy.special", scipy_special_ufuncs, get_scipy_special_helpers
)
