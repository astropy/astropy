# The idea for this module (but no code) was borrowed from the
# quantities (http://pythonhosted.org/quantities/) package.

import numpy as np
from .core import UnitsError, dimensionless_unscaled


def _d(unit):
    if unit is None:
        return dimensionless_unscaled
    else:
        return unit

UFUNC_HELPERS = {}

# In this file, we implement the logic that determines for a given ufunc and
# input how the input should be scaled and what unit the output will have.

# list of ufuncs:
# http://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs

UNSUPPORTED_UFUNCS = set([np.bitwise_and, np.bitwise_or,
                          np.bitwise_xor, np.invert, np.left_shift,
                          np.right_shift, np.logical_and, np.logical_or,
                          np.logical_xor, np.logical_not])

# SINGLE ARGUMENT UFUNCS

# The functions below take a single argument, which is the quantity upon which
# the ufunc is being used. The output of the function should be two values: the
# scale by which the input needs to be multiplied before being passed to the
# ufunc, and the unit the output will be in.

# ufuncs that return a boolean and do not care about the unit
helper_onearg_test = lambda f, unit: ([1.], None)

UFUNC_HELPERS[np.isfinite] = helper_onearg_test
UFUNC_HELPERS[np.isinf] = helper_onearg_test
UFUNC_HELPERS[np.isnan] = helper_onearg_test
UFUNC_HELPERS[np.sign] = helper_onearg_test
UFUNC_HELPERS[np.signbit] = helper_onearg_test

# ufuncs that return a value with the same unit as the input

helper_invariant = lambda f, unit: ([1.], _d(unit))

UFUNC_HELPERS[np.absolute] = helper_invariant
UFUNC_HELPERS[np.fabs] = helper_invariant
UFUNC_HELPERS[np.conj] = helper_invariant
UFUNC_HELPERS[np.conjugate] = helper_invariant
UFUNC_HELPERS[np.negative] = helper_invariant
UFUNC_HELPERS[np.spacing] = helper_invariant
UFUNC_HELPERS[np.rint] = helper_invariant
UFUNC_HELPERS[np.floor] = helper_invariant
UFUNC_HELPERS[np.ceil] = helper_invariant
UFUNC_HELPERS[np.trunc] = helper_invariant

# ufuncs handled as special cases

UFUNC_HELPERS[np.sqrt] = lambda f, unit: ([1.], unit ** 0.5 if unit is not None
                                          else dimensionless_unscaled)
UFUNC_HELPERS[np.square] = lambda f, unit: ([1.], unit ** 2 if unit is not None
                                            else dimensionless_unscaled)
UFUNC_HELPERS[np.reciprocal] = lambda f, unit: ([1.], unit ** -1
                                                if unit is not None
                                                else dimensionless_unscaled)
# ones_like was not private in numpy <= 1.6
if isinstance(getattr(np.core.umath, 'ones_like', None), np.ufunc):
    UFUNC_HELPERS[np.core.umath.ones_like] = (lambda f, unit:
                                              ([1.], dimensionless_unscaled))
if isinstance(getattr(np.core.umath, '_ones_like', None), np.ufunc):
    UFUNC_HELPERS[np.core.umath._ones_like] = (lambda f, unit:
                                              ([1.], dimensionless_unscaled))


# ufuncs that require dimensionless input and and give dimensionless output
def helper_dimensionless_to_dimensionless(f, unit):
    try:
        scale = unit.to(dimensionless_unscaled) if unit is not None else 1.
    except UnitsError:
        raise TypeError("Can only apply '{0}' function to "
                        "dimensionless quantities"
                        .format(f.__name__))
    return [scale], dimensionless_unscaled

UFUNC_HELPERS[np.exp] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.expm1] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.exp2] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log10] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log2] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log1p] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.modf] = helper_dimensionless_to_dimensionless


# ufuncs that require dimensionless input and give output in radians
def helper_dimensionless_to_radian(f, unit):
    from .si import radian
    try:
        scale = unit.to(dimensionless_unscaled) if unit is not None else 1.
    except UnitsError:
        raise TypeError("Can only apply '{0}' function to "
                        "dimensionless quantities"
                        .format(f.__name__))
    return [scale], radian

UFUNC_HELPERS[np.arccos] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arcsin] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arctan] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arccosh] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arcsinh] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arctanh] = helper_dimensionless_to_radian


# ufuncs that require input in degrees and give output in radians
def helper_degree_to_radian(f, unit):
    from .si import degree, radian
    try:
        scale = unit.to(degree)
    except UnitsError:
        raise TypeError("Can only apply '{0}' function to "
                        "quantities with angle units"
                        .format(f.__name__))
    return [scale], radian

UFUNC_HELPERS[np.radians] = helper_degree_to_radian
UFUNC_HELPERS[np.deg2rad] = helper_degree_to_radian


# ufuncs that require input in radians and give output in degrees
def helper_radian_to_degree(f, unit):
    from .si import degree, radian
    try:
        scale = unit.to(radian)
    except UnitsError:
        raise TypeError("Can only apply '{0}' function to "
                        "quantities with angle units"
                        .format(f.__name__))
    return [scale], degree

UFUNC_HELPERS[np.degrees] = helper_radian_to_degree
UFUNC_HELPERS[np.rad2deg] = helper_radian_to_degree


# ufuncs that require input in radians and give dimensionless output
def helper_radian_to_dimensionless(f, unit):
    from .si import radian
    try:
        scale = unit.to(radian)
    except UnitsError:
        raise TypeError("Can only apply '{0}' function to "
                        "quantities with angle units"
                        .format(f.__name__))
    return [scale], dimensionless_unscaled

UFUNC_HELPERS[np.cos] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.sin] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.tan] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.cosh] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.sinh] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.tanh] = helper_radian_to_dimensionless


# ufuncs that require dimensionless_unscaled input and return non-quantities
def helper_dimensionless_to_none(f, unit):
    if not unit.is_unity():
        raise TypeError("Can only apply '{0}' function to "
                        "unscaled dimensionless quantities"
                        .format(f.__name__))
    return [1.], None

UFUNC_HELPERS[np.frexp] = helper_dimensionless_to_none

# TWO ARGUMENT UFUNCS

UFUNC_HELPERS[np.multiply] = lambda f, unit1, unit2: (
    [1., 1.], _d(unit1) * _d(unit2))

helper_division = lambda f, unit1, unit2: ([1., 1.], _d(unit1) / _d(unit2))

UFUNC_HELPERS[np.divide] = helper_division
UFUNC_HELPERS[np.true_divide] = helper_division
UFUNC_HELPERS[np.floor_divide] = helper_division


def helper_power(f, unit1, unit2):
    if unit2 is not None:
        try:
            scale2 = unit2.to(dimensionless_unscaled)
        except UnitsError:
            raise TypeError("Can only raise something to a "
                            "dimensionless quantity")
    else:
        scale2 = 1.

    # TODO: find a better way to do this, currently
    # need to raise power of unit1 in main code
    return [1., scale2], _d(unit1)

UFUNC_HELPERS[np.power] = helper_power


def helper_ldexp(f, unit1, unit2):
    if unit2 is not None:
        raise TypeError("Cannot use ldexp with a quantity "
                        "as second argument.")
    else:
        return [1., 1.], _d(unit1)

UFUNC_HELPERS[np.ldexp] = helper_ldexp


def helper_copysign(f, unit1, unit2):
    # if first arg is not a quantity, just return plain array
    if unit1 is None:
        return [1., 1.], None
    else:
        return [1., 1.], unit1

UFUNC_HELPERS[np.copysign] = helper_copysign


def helper_two_arg_dimensionless(f, unit1, unit2):
    try:
        scale1 = unit1.to(dimensionless_unscaled) if unit1 is not None else 1.
        scale2 = unit2.to(dimensionless_unscaled) if unit2 is not None else 1.
    except UnitsError:
        raise TypeError("Can only apply '{0}' function to "
                        "dimensionless quantities"
                        .format(f.__name__))
    return [scale1, scale2], dimensionless_unscaled

UFUNC_HELPERS[np.logaddexp] = helper_two_arg_dimensionless
UFUNC_HELPERS[np.logaddexp2] = helper_two_arg_dimensionless


def find_scales(f, *units):

    scales = [1., 1.]
    # no units for any input -- e.g., np.add(a1, a2, out=q)
    if all(unit is None for unit in units):
        return scales, dimensionless_unscaled

    fixed, changeable = (1, 0) if units[1] is None else (0, 1)
    if units[fixed] is None:
        try:
            scales[changeable] = units[changeable].to(dimensionless_unscaled)
        except UnitsError:
            # special case: would be OK if unitless number is zero, inf, nan
            scales[fixed] = 0.

        return scales, dimensionless_unscaled

    else:
        try:
            scales[changeable] = units[changeable].to(units[fixed])
        except UnitsError:
            raise UnitsError(
                "Can only apply '{0}' function to quantities "
                "with compatible dimensions"
                .format(f.__name__))

        return scales, units[fixed]


def helper_twoarg_invariant(f, unit1, unit2):
    return find_scales(f, unit1, unit2)

UFUNC_HELPERS[np.add] = helper_twoarg_invariant
UFUNC_HELPERS[np.subtract] = helper_twoarg_invariant
UFUNC_HELPERS[np.hypot] = helper_twoarg_invariant
UFUNC_HELPERS[np.maximum] = helper_twoarg_invariant
UFUNC_HELPERS[np.minimum] = helper_twoarg_invariant
UFUNC_HELPERS[np.fmin] = helper_twoarg_invariant
UFUNC_HELPERS[np.fmax] = helper_twoarg_invariant
UFUNC_HELPERS[np.nextafter] = helper_twoarg_invariant
UFUNC_HELPERS[np.remainder] = helper_twoarg_invariant
UFUNC_HELPERS[np.mod] = helper_twoarg_invariant
UFUNC_HELPERS[np.fmod] = helper_twoarg_invariant


def helper_twoarg_comparison(f, unit1, unit2):
    scales, _ = find_scales(f, unit1, unit2)
    return scales, None

UFUNC_HELPERS[np.greater] = helper_twoarg_comparison
UFUNC_HELPERS[np.greater_equal] = helper_twoarg_comparison
UFUNC_HELPERS[np.less] = helper_twoarg_comparison
UFUNC_HELPERS[np.less_equal] = helper_twoarg_comparison
UFUNC_HELPERS[np.not_equal] = helper_twoarg_comparison
UFUNC_HELPERS[np.equal] = helper_twoarg_comparison


def helper_twoarg_invtrig(f, unit1, unit2):
    from .si import radian
    scales, _ = find_scales(f, unit1, unit2)
    return scales, radian

UFUNC_HELPERS[np.arctan2] = helper_twoarg_invtrig
# another private function in numpy; use getattr in case it disappears
if isinstance(getattr(np.core.umath, '_arg', None), np.ufunc):
    UFUNC_HELPERS[np.core.umath._arg] = helper_twoarg_invtrig
