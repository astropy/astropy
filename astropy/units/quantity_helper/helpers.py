# Licensed under a 3-clause BSD style license - see LICENSE.rst
# The idea for this module (but no code) was borrowed from the
# quantities (http://pythonhosted.org/quantities/) package.
"""Helper functions for Quantity.

In particular, this implements the logic that determines scaling and result
units for a given ufunc, given input units.
"""

from fractions import Fraction

import numpy as np

from astropy.units.core import (
    UnitConversionError,
    UnitsError,
    UnitTypeError,
    dimensionless_unscaled,
    unit_scale_converter,
)
from astropy.utils.compat.numpycompat import NUMPY_LT_2_0

if NUMPY_LT_2_0:
    from numpy.core import umath as np_umath
else:
    from numpy._core import umath as np_umath

from . import UFUNC_HELPERS, UNSUPPORTED_UFUNCS


def _d(unit):
    if unit is None:
        return dimensionless_unscaled
    else:
        return unit


def get_converter(from_unit, to_unit):
    """Like Unit._get_converter, except returns None if no scaling is needed,
    i.e., if the inferred scale is unity.
    """
    converter = from_unit._get_converter(to_unit)
    return None if converter is unit_scale_converter else converter


def get_converters_and_unit(f, unit1, unit2):
    converters = [None, None]
    # By default, we try adjusting unit2 to unit1, so that the result will
    # be unit1 as well. But if there is no second unit, we have to try
    # adjusting unit1 (to dimensionless, see below).
    if unit2 is None:
        if unit1 is None:
            # No units for any input -- e.g., np.add(a1, a2, out=q)
            return converters, dimensionless_unscaled
        changeable = 0
        # swap units.
        unit2 = unit1
        unit1 = None
    elif unit2 is unit1:
        # ensure identical units is fast ("==" is slow, so avoid that).
        return converters, unit1
    else:
        changeable = 1

    # Try to get a converter from unit2 to unit1.
    if unit1 is None:
        try:
            converters[changeable] = get_converter(unit2, dimensionless_unscaled)
        except UnitsError:
            # special case: would be OK if unitless number is zero, inf, nan
            converters[1 - changeable] = False
            return converters, unit2
        else:
            return converters, dimensionless_unscaled
    else:
        try:
            converters[changeable] = get_converter(unit2, unit1)
        except UnitsError:
            raise UnitConversionError(
                f"Can only apply '{f.__name__}' function to quantities "
                "with compatible dimensions"
            )

        return converters, unit1


# SINGLE ARGUMENT UFUNC HELPERS
#
# The functions below take a single argument, which is the quantity upon which
# the ufunc is being used. The output of the helper function should be two
# values: a list with a single converter to be used to scale the input before
# it is being passed to the ufunc (or None if no conversion is needed), and
# the unit the output will be in.


def helper_onearg_test(f, unit):
    return ([None], None)


def helper_invariant(f, unit):
    return ([None], _d(unit))


def helper_square(f, unit):
    return ([None], unit**2 if unit is not None else dimensionless_unscaled)


def helper_reciprocal(f, unit):
    return ([None], unit**-1 if unit is not None else dimensionless_unscaled)


one_half = 0.5  # faster than Fraction(1, 2)
one_third = Fraction(1, 3)


def helper_sqrt(f, unit):
    return ([None], unit**one_half if unit is not None else dimensionless_unscaled)


def helper_cbrt(f, unit):
    return ([None], (unit**one_third if unit is not None else dimensionless_unscaled))


def helper_modf(f, unit):
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


def helper__ones_like(f, unit):
    return [None], dimensionless_unscaled


def helper_dimensionless_to_dimensionless(f, unit):
    if unit is None:
        return [None], dimensionless_unscaled

    try:
        return ([get_converter(unit, dimensionless_unscaled)], dimensionless_unscaled)
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )


def helper_dimensionless_to_radian(f, unit):
    from astropy.units.si import radian

    if unit is None:
        return [None], radian

    try:
        return [get_converter(unit, dimensionless_unscaled)], radian
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )


def helper_degree_to_radian(f, unit):
    from astropy.units.si import degree, radian

    try:
        return [get_converter(unit, degree)], radian
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to quantities with angle units"
        )


def helper_radian_to_degree(f, unit):
    from astropy.units.si import degree, radian

    try:
        return [get_converter(unit, radian)], degree
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to quantities with angle units"
        )


def helper_radian_to_dimensionless(f, unit):
    from astropy.units.si import radian

    try:
        return [get_converter(unit, radian)], dimensionless_unscaled
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to quantities with angle units"
        )


def helper_frexp(f, unit):
    if not unit.is_unity():
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to unscaled dimensionless"
            " quantities"
        )
    return [None], (None, None)


# TWO ARGUMENT UFUNC HELPERS
#
# The functions below take a two arguments. The output of the helper function
# should be two values: a tuple of two converters to be used to scale the
# inputs before being passed to the ufunc (None if no conversion is needed),
# and the unit the output will be in.


def helper_multiplication(f, unit1, unit2):
    return [None, None], _d(unit1) * _d(unit2)


def helper_division(f, unit1, unit2):
    return [None, None], _d(unit1) / _d(unit2)


def helper_power(f, unit1, unit2):
    # TODO: find a better way to do this, currently need to signal that one
    # still needs to raise power of unit1 in main code
    if unit2 is None:
        return [None, None], False

    try:
        return [None, get_converter(unit2, dimensionless_unscaled)], False
    except UnitsError:
        raise UnitTypeError("Can only raise something to a dimensionless quantity")


def helper_ldexp(f, unit1, unit2):
    if unit2 is not None:
        raise TypeError("Cannot use ldexp with a quantity as second argument.")
    else:
        return [None, None], _d(unit1)


def helper_copysign(f, unit1, unit2):
    # if first arg is not a quantity, just return plain array
    if unit1 is None:
        return [None, None], None
    else:
        return [None, None], unit1


def helper_heaviside(f, unit1, unit2):
    try:
        converter2 = (
            get_converter(unit2, dimensionless_unscaled) if unit2 is not None else None
        )
    except UnitsError:
        raise UnitTypeError(
            "Can only apply 'heaviside' function with a dimensionless second argument."
        )
    return ([None, converter2], dimensionless_unscaled)


def helper_two_arg_dimensionless(f, unit1, unit2):
    try:
        converter1 = (
            get_converter(unit1, dimensionless_unscaled) if unit1 is not None else None
        )
        converter2 = (
            get_converter(unit2, dimensionless_unscaled) if unit2 is not None else None
        )
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to dimensionless quantities"
        )
    return ([converter1, converter2], dimensionless_unscaled)


# This used to be a separate function that just called get_converters_and_unit.
# Using it directly saves a few us; keeping the clearer name.
helper_twoarg_invariant = get_converters_and_unit


def helper_twoarg_comparison(f, unit1, unit2):
    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, None


def helper_twoarg_invtrig(f, unit1, unit2):
    from astropy.units.si import radian

    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, radian


def helper_twoarg_floor_divide(f, unit1, unit2):
    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, dimensionless_unscaled


def helper_divmod(f, unit1, unit2):
    converters, result_unit = get_converters_and_unit(f, unit1, unit2)
    return converters, (dimensionless_unscaled, result_unit)


def helper_clip(f, unit1, unit2, unit3):
    # Treat the array being clipped as primary.
    converters = [None]
    if unit1 is None:
        result_unit = dimensionless_unscaled
        try:
            converters += [
                (None if unit is None else get_converter(unit, dimensionless_unscaled))
                for unit in (unit2, unit3)
            ]
        except UnitsError:
            raise UnitConversionError(
                f"Can only apply '{f.__name__}' function to quantities with "
                "compatible dimensions"
            )

    else:
        result_unit = unit1
        for unit in unit2, unit3:
            try:
                converter = get_converter(_d(unit), result_unit)
            except UnitsError:
                if unit is None:
                    # special case: OK if unitless number is zero, inf, nan
                    converters.append(False)
                else:
                    raise UnitConversionError(
                        f"Can only apply '{f.__name__}' function to quantities with "
                        "compatible dimensions"
                    )
            else:
                converters.append(converter)

    return converters, result_unit


# list of ufuncs:
# https://numpy.org/doc/stable/reference/ufuncs.html#available-ufuncs

UNSUPPORTED_UFUNCS |= {
    np.bitwise_and,
    np.bitwise_or,
    np.bitwise_xor,
    np.invert,
    np.left_shift,
    np.right_shift,
    np.logical_and,
    np.logical_or,
    np.logical_xor,
    np.logical_not,
    np.isnat,
    np.gcd,
    np.lcm,
}

if not NUMPY_LT_2_0:
    # string utilities - make no sense for Quantity.
    UNSUPPORTED_UFUNCS |= {
        np.bitwise_count,
        np._core.umath.count,
        np._core.umath.isalpha,
        np._core.umath.isdigit,
        np._core.umath.isspace,
        np._core.umath.isnumeric,
        np._core.umath.isdecimal,
        np._core.umath.isalnum,
        np._core.umath.istitle,
        np._core.umath.islower,
        np._core.umath.isupper,
        np._core.umath.index,
        np._core.umath.rindex,
        np._core.umath.startswith,
        np._core.umath.endswith,
        np._core.umath.find,
        np._core.umath.rfind,
        np._core.umath.str_len,
        np._core.umath._strip_chars,
        np._core.umath._lstrip_chars,
        np._core.umath._rstrip_chars,
        np._core.umath._strip_whitespace,
        np._core.umath._lstrip_whitespace,
        np._core.umath._rstrip_whitespace,
        np._core.umath._replace,
        np._core.umath._expandtabs,
        np._core.umath._expandtabs_length,
        np._core.umath._ljust,
        np._core.umath._rjust,
        np._core.umath._center,
        np._core.umath._zfill,
    }

# SINGLE ARGUMENT UFUNCS

# ufuncs that do not care about the unit and do not return a Quantity
# (but rather a boolean, or -1, 0, or +1 for np.sign).
onearg_test_ufuncs = (np.isfinite, np.isinf, np.isnan, np.sign, np.signbit)
for ufunc in onearg_test_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_onearg_test

# ufuncs that return a value with the same unit as the input
invariant_ufuncs = (
    np.absolute,
    np.fabs,
    np.conj,
    np.conjugate,
    np.negative,
    np.spacing,
    np.rint,
    np.floor,
    np.ceil,
    np.trunc,
    np.positive,
)
for ufunc in invariant_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_invariant

# ufuncs that require dimensionless input and and give dimensionless output
dimensionless_to_dimensionless_ufuncs = (
    np.exp,
    np.expm1,
    np.exp2,
    np.log,
    np.log10,
    np.log2,
    np.log1p,
)
# Default numpy does not ship an "erf" ufunc, but some versions hacked by
# intel do.  This is bad, since it means code written for that numpy will
# not run on non-hacked numpy.  But still, we might as well support it.
if isinstance(getattr(np_umath, "erf", None), np.ufunc):
    dimensionless_to_dimensionless_ufuncs += (np_umath.erf,)

for ufunc in dimensionless_to_dimensionless_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_dimensionless_to_dimensionless

# ufuncs that require dimensionless input and give output in radians
dimensionless_to_radian_ufuncs = (
    np.arccos,
    np.arcsin,
    np.arctan,
    np.arccosh,
    np.arcsinh,
    np.arctanh,
)
for ufunc in dimensionless_to_radian_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_dimensionless_to_radian

# ufuncs that require input in degrees and give output in radians
degree_to_radian_ufuncs = (np.radians, np.deg2rad)
for ufunc in degree_to_radian_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_degree_to_radian

# ufuncs that require input in radians and give output in degrees
radian_to_degree_ufuncs = (np.degrees, np.rad2deg)
for ufunc in radian_to_degree_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_radian_to_degree

# ufuncs that require input in radians and give dimensionless output
radian_to_dimensionless_ufuncs = (np.cos, np.sin, np.tan, np.cosh, np.sinh, np.tanh)
for ufunc in radian_to_dimensionless_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_radian_to_dimensionless

# ufuncs handled as special cases
UFUNC_HELPERS[np.sqrt] = helper_sqrt
UFUNC_HELPERS[np.square] = helper_square
UFUNC_HELPERS[np.reciprocal] = helper_reciprocal
UFUNC_HELPERS[np.cbrt] = helper_cbrt
UFUNC_HELPERS[np_umath._ones_like] = helper__ones_like
UFUNC_HELPERS[np.modf] = helper_modf
UFUNC_HELPERS[np.frexp] = helper_frexp


# TWO ARGUMENT UFUNCS

# two argument ufuncs that require dimensionless input and and give
# dimensionless output
two_arg_dimensionless_ufuncs = (np.logaddexp, np.logaddexp2)
for ufunc in two_arg_dimensionless_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_two_arg_dimensionless

# two argument ufuncs that return a value with the same unit as the input
twoarg_invariant_ufuncs = (
    np.add,
    np.subtract,
    np.hypot,
    np.maximum,
    np.minimum,
    np.fmin,
    np.fmax,
    np.nextafter,
    np.remainder,
    np.mod,
    np.fmod,
)
for ufunc in twoarg_invariant_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_twoarg_invariant

# two argument ufuncs that need compatible inputs and return a boolean
twoarg_comparison_ufuncs = (
    np.greater,
    np.greater_equal,
    np.less,
    np.less_equal,
    np.not_equal,
    np.equal,
)
for ufunc in twoarg_comparison_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_twoarg_comparison

# two argument ufuncs that do inverse trigonometry
twoarg_invtrig_ufuncs = (np.arctan2,)
# another private function in numpy; use getattr in case it disappears
if isinstance(getattr(np_umath, "_arg", None), np.ufunc):
    twoarg_invtrig_ufuncs += (np_umath._arg,)
for ufunc in twoarg_invtrig_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_twoarg_invtrig

# ufuncs handled as special cases
UFUNC_HELPERS[np.multiply] = helper_multiplication
UFUNC_HELPERS[np.matmul] = helper_multiplication
if isinstance(getattr(np, "vecdot", None), np.ufunc):
    UFUNC_HELPERS[np.vecdot] = helper_multiplication
UFUNC_HELPERS[np.divide] = helper_division
UFUNC_HELPERS[np.true_divide] = helper_division
UFUNC_HELPERS[np.power] = helper_power
UFUNC_HELPERS[np.ldexp] = helper_ldexp
UFUNC_HELPERS[np.copysign] = helper_copysign
UFUNC_HELPERS[np.floor_divide] = helper_twoarg_floor_divide
UFUNC_HELPERS[np.heaviside] = helper_heaviside
UFUNC_HELPERS[np.float_power] = helper_power
UFUNC_HELPERS[np.divmod] = helper_divmod
# Check for clip ufunc; note that np.clip is a wrapper function, not the ufunc.
if isinstance(getattr(np_umath, "clip", None), np.ufunc):
    UFUNC_HELPERS[np_umath.clip] = helper_clip

del ufunc
