# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# The idea for this module (but no code) was borrowed from the
# quantities (http://pythonhosted.org/quantities/) package.
"""Helper functions for Quantity.

In particular, this implements the logic that determines scaling and result
units for a given ufunc, given input units.
"""

from fractions import Fraction

import numpy as np
from .core import (UnitsError, UnitConversionError, UnitTypeError,
                   dimensionless_unscaled, get_current_unit_registry)


def _d(unit):
    if unit is None:
        return dimensionless_unscaled
    else:
        return unit


def get_converter(from_unit, to_unit):
    """Like Unit._get_converter, except returns None if no scaling is needed,
    i.e., if the inferred scale is unity."""
    try:
        scale = from_unit._to(to_unit)
    except UnitsError:
        return from_unit._apply_equivalencies(
                from_unit, to_unit, get_current_unit_registry().equivalencies)
    except AttributeError:
        raise UnitTypeError("Unit '{0}' cannot be converted to '{1}'"
                            .format(from_unit, to_unit))
    if scale == 1.:
        return None
    else:
        return lambda val: scale * val


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
            converters[changeable] = get_converter(unit2,
                                                   dimensionless_unscaled)
        except UnitsError:
            # special case: would be OK if unitless number is zero, inf, nan
            converters[1-changeable] = False
            return converters, unit2
        else:
            return converters, dimensionless_unscaled
    else:
        try:
            converters[changeable] = get_converter(unit2, unit1)
        except UnitsError:
            raise UnitConversionError(
                "Can only apply '{0}' function to quantities "
                "with compatible dimensions"
                .format(f.__name__))

        return converters, unit1


def can_have_arbitrary_unit(value):
    """Test whether the items in value can have arbitrary units

    Numbers whose value does not change upon a unit change, i.e.,
    zero, infinity, or not-a-number

    Parameters
    ----------
    value : number or array

    Returns
    -------
    `True` if each member is either zero or not finite, `False` otherwise
    """
    return np.all(np.logical_or(np.equal(value, 0.), ~np.isfinite(value)))


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


def helper_sqrt(f, unit):
    return ([None], unit ** Fraction(1, 2) if unit is not None
            else dimensionless_unscaled)


def helper_square(f, unit):
    return ([None], unit ** 2 if unit is not None else dimensionless_unscaled)


def helper_reciprocal(f, unit):
    return ([None], unit ** -1 if unit is not None else dimensionless_unscaled)


def helper_cbrt(f, unit):
    return ([None], (unit ** Fraction(1, 3) if unit is not None
                     else dimensionless_unscaled))


def helper_modf(f, unit):
    if unit is None:
        return [None], (dimensionless_unscaled, dimensionless_unscaled)

    try:
        return ([get_converter(unit, dimensionless_unscaled)],
                (dimensionless_unscaled, dimensionless_unscaled))
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "dimensionless quantities"
                            .format(f.__name__))


def helper__ones_like(f, unit):
    return [None], dimensionless_unscaled


def helper_dimensionless_to_dimensionless(f, unit):
    if unit is None:
        return [None], dimensionless_unscaled

    try:
        return ([get_converter(unit, dimensionless_unscaled)],
                dimensionless_unscaled)
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "dimensionless quantities"
                            .format(f.__name__))


def helper_dimensionless_to_radian(f, unit):
    from .si import radian
    if unit is None:
        return [None], radian

    try:
        return [get_converter(unit, dimensionless_unscaled)], radian
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "dimensionless quantities"
                            .format(f.__name__))


def helper_degree_to_radian(f, unit):
    from .si import degree, radian
    try:
        return [get_converter(unit, degree)], radian
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


def helper_radian_to_degree(f, unit):
    from .si import degree, radian
    try:
        return [get_converter(unit, radian)], degree
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


def helper_radian_to_dimensionless(f, unit):
    from .si import radian
    try:
        return [get_converter(unit, radian)], dimensionless_unscaled
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


def helper_frexp(f, unit):
    if not unit.is_unity():
        raise UnitTypeError("Can only apply '{0}' function to "
                            "unscaled dimensionless quantities"
                            .format(f.__name__))
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
        raise UnitTypeError("Can only raise something to a "
                            "dimensionless quantity")


def helper_ldexp(f, unit1, unit2):
    if unit2 is not None:
        raise TypeError("Cannot use ldexp with a quantity "
                        "as second argument.")
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
        converter2 = (get_converter(unit2, dimensionless_unscaled)
                      if unit2 is not None else None)
    except UnitsError:
        raise UnitTypeError("Can only apply 'heaviside' function with a "
                            "dimensionless second argument.")
    return ([None, converter2], dimensionless_unscaled)


def helper_two_arg_dimensionless(f, unit1, unit2):
    try:
        converter1 = (get_converter(unit1, dimensionless_unscaled)
                      if unit1 is not None else None)
        converter2 = (get_converter(unit2, dimensionless_unscaled)
                      if unit2 is not None else None)
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "dimensionless quantities"
                            .format(f.__name__))
    return ([converter1, converter2], dimensionless_unscaled)


# This used to be a separate function that just called get_converters_and_unit.
# Using it directly saves a few us; keeping the clearer name.
helper_twoarg_invariant = get_converters_and_unit


def helper_twoarg_comparison(f, unit1, unit2):
    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, None


def helper_twoarg_invtrig(f, unit1, unit2):
    from .si import radian
    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, radian


def helper_twoarg_floor_divide(f, unit1, unit2):
    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, dimensionless_unscaled


def helper_divmod(f, unit1, unit2):
    converters, result_unit = get_converters_and_unit(f, unit1, unit2)
    return converters, (dimensionless_unscaled, result_unit)


def helper_degree_to_dimensionless(f, unit):
    from .si import degree
    try:
        return [get_converter(unit, degree)], dimensionless_unscaled
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


def helper_degree_minute_second_to_radian(f, unit1, unit2, unit3):
    from .si import degree, arcmin, arcsec, radian
    try:
        return [get_converter(unit1, degree),
                get_converter(unit2, arcmin),
                get_converter(unit3, arcsec)], radian
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


# list of ufuncs:
# http://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs
UFUNC_HELPERS = {}

UNSUPPORTED_UFUNCS = {
    np.bitwise_and, np.bitwise_or, np.bitwise_xor, np.invert, np.left_shift,
    np.right_shift, np.logical_and, np.logical_or, np.logical_xor,
    np.logical_not}
for name in 'isnat', 'gcd', 'lcm':
    # isnat was introduced in numpy 1.14, gcd+lcm in 1.15
    ufunc = getattr(np, name, None)
    if isinstance(ufunc, np.ufunc):
        UNSUPPORTED_UFUNCS |= {ufunc}

# SINGLE ARGUMENT UFUNCS

# ufuncs that return a boolean and do not care about the unit
onearg_test_ufuncs = (np.isfinite, np.isinf, np.isnan, np.sign, np.signbit)
for ufunc in onearg_test_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_onearg_test

# ufuncs that return a value with the same unit as the input
invariant_ufuncs = (np.absolute, np.fabs, np.conj, np.conjugate, np.negative,
                    np.spacing, np.rint, np.floor, np.ceil, np.trunc)
for ufunc in invariant_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_invariant
# positive was added in numpy 1.13
if isinstance(getattr(np, 'positive', None), np.ufunc):
    UFUNC_HELPERS[np.positive] = helper_invariant

# ufuncs that require dimensionless input and and give dimensionless output
dimensionless_to_dimensionless_ufuncs = (np.exp, np.expm1, np.exp2, np.log,
                                         np.log10, np.log2, np.log1p)
for ufunc in dimensionless_to_dimensionless_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_dimensionless_to_dimensionless

# ufuncs that require dimensionless input and give output in radians
dimensionless_to_radian_ufuncs = (np.arccos, np.arcsin, np.arctan, np.arccosh,
                                  np.arcsinh, np.arctanh)
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
radian_to_dimensionless_ufuncs = (np.cos, np.sin, np.tan, np.cosh, np.sinh,
                                  np.tanh)
for ufunc in radian_to_dimensionless_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_radian_to_dimensionless

# ufuncs handled as special cases
UFUNC_HELPERS[np.sqrt] = helper_sqrt
UFUNC_HELPERS[np.square] = helper_square
UFUNC_HELPERS[np.reciprocal] = helper_reciprocal
UFUNC_HELPERS[np.cbrt] = helper_cbrt
UFUNC_HELPERS[np.core.umath._ones_like] = helper__ones_like
UFUNC_HELPERS[np.modf] = helper_modf
UFUNC_HELPERS[np.frexp] = helper_frexp


# TWO ARGUMENT UFUNCS

# two argument ufuncs that require dimensionless input and and give
# dimensionless output
two_arg_dimensionless_ufuncs = (np.logaddexp, np.logaddexp2)
for ufunc in two_arg_dimensionless_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_two_arg_dimensionless

# two argument ufuncs that return a value with the same unit as the input
twoarg_invariant_ufuncs = (np.add, np.subtract, np.hypot, np.maximum,
                           np.minimum, np.fmin, np.fmax, np.nextafter,
                           np.remainder, np.mod, np.fmod)
for ufunc in twoarg_invariant_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_twoarg_invariant

# two argument ufuncs that need compatible inputs and return a boolean
twoarg_comparison_ufuncs = (np.greater, np.greater_equal, np.less,
                            np.less_equal, np.not_equal, np.equal)
for ufunc in twoarg_comparison_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_twoarg_comparison

# two argument ufuncs that do inverse trigonometry
twoarg_invtrig_ufuncs = (np.arctan2,)
# another private function in numpy; use getattr in case it disappears
if isinstance(getattr(np.core.umath, '_arg', None), np.ufunc):
    twoarg_invtrig_ufuncs += (np.core.umath._arg,)
for ufunc in twoarg_invtrig_ufuncs:
    UFUNC_HELPERS[ufunc] = helper_twoarg_invtrig

# ufuncs handled as special cases
UFUNC_HELPERS[np.multiply] = helper_multiplication
UFUNC_HELPERS[np.divide] = helper_division
UFUNC_HELPERS[np.true_divide] = helper_division
UFUNC_HELPERS[np.power] = helper_power
UFUNC_HELPERS[np.ldexp] = helper_ldexp
UFUNC_HELPERS[np.copysign] = helper_copysign
UFUNC_HELPERS[np.floor_divide] = helper_twoarg_floor_divide
# heaviside only was added in numpy 1.13
if isinstance(getattr(np, 'heaviside', None), np.ufunc):
    UFUNC_HELPERS[np.heaviside] = helper_heaviside
# float_power was added in numpy 1.12
if isinstance(getattr(np, 'float_power', None), np.ufunc):
    UFUNC_HELPERS[np.float_power] = helper_power
# divmod only was added in numpy 1.13
if isinstance(getattr(np, 'divmod', None), np.ufunc):
    UFUNC_HELPERS[np.divmod] = helper_divmod


# UFUNCS FROM SCIPY.SPECIAL
# available ufuncs in this module are at
# https://docs.scipy.org/doc/scipy/reference/special.html

try:
    import scipy
    import scipy.special as sps
except ImportError:
    pass
else:
    from ..utils import minversion

    # ufuncs that require dimensionless input and give dimensionless output
    dimensionless_to_dimensionless_sps_ufuncs = [
        sps.erf, sps.gamma, sps.gammasgn,
        sps.psi, sps.rgamma, sps.erfc, sps.erfcx, sps.erfi, sps.wofz,
        sps.dawsn, sps.entr, sps.exprel, sps.expm1, sps.log1p, sps.exp2,
        sps.exp10, sps.j0, sps.j1, sps.y0, sps.y1, sps.i0, sps.i0e, sps.i1,
        sps.i1e, sps.k0, sps.k0e, sps.k1, sps.k1e, sps.itj0y0,
        sps.it2j0y0, sps.iti0k0, sps.it2i0k0]

    # TODO: Revert https://github.com/astropy/astropy/pull/7219 when astropy
    #       requires scipy>=0.18.
    # See https://github.com/astropy/astropy/issues/7159
    if minversion(scipy, "0.18"):
        dimensionless_to_dimensionless_sps_ufuncs.append(sps.loggamma)

    for ufunc in dimensionless_to_dimensionless_sps_ufuncs:
        UFUNC_HELPERS[ufunc] = helper_dimensionless_to_dimensionless

    # ufuncs that require input in degrees and give dimensionless output
    degree_to_dimensionless_sps_ufuncs = (
        sps.cosdg, sps.sindg, sps.tandg, sps.cotdg)
    for ufunc in degree_to_dimensionless_sps_ufuncs:
        UFUNC_HELPERS[ufunc] = helper_degree_to_dimensionless

    # ufuncs that require 2 dimensionless inputs and give dimensionless output.
    # note: sps.jv and sps.jn are aliases in some scipy versions, which will
    # cause the same key to be written twice, but since both are handled by the
    # same helper there is no harm done.
    two_arg_dimensionless_sps_ufuncs = (
        sps.jv, sps.jn, sps.jve, sps.yn, sps.yv, sps.yve, sps.kn, sps.kv,
        sps.kve, sps.iv, sps.ive, sps.hankel1, sps.hankel1e, sps.hankel2,
        sps.hankel2e)
    for ufunc in two_arg_dimensionless_sps_ufuncs:
        UFUNC_HELPERS[ufunc] = helper_two_arg_dimensionless

    # ufuncs handled as special cases
    UFUNC_HELPERS[sps.cbrt] = helper_cbrt
    UFUNC_HELPERS[sps.radian] = helper_degree_minute_second_to_radian


def converters_and_unit(function, method, *args):
    """Determine the required converters and the unit of the ufunc result.

    Converters are functions required to convert to a ufunc's expected unit,
    e.g., radian for np.sin; or to ensure units of two inputs are consistent,
    e.g., for np.add.  In these examples, the unit of the result would be
    dimensionless_unscaled for np.sin, and the same consistent unit for np.add.

    Parameters
    ----------
    function : `~numpy.ufunc`
        Numpy universal function
    method : str
        Method with which the function is evaluated, e.g.,
        '__call__', 'reduce', etc.
    *args : Quantity or other ndarray subclass
        Input arguments to the function

    Raises
    ------
    TypeError : when the specified function cannot be used with Quantities
        (e.g., np.logical_or), or when the routine does not know how to handle
        the specified function (in which case an issue should be raised on
        https://github.com/astropy/astropy).
    UnitTypeError : when the conversion to the required (or consistent) units
        is not possible.
    """

    # Check whether we support this ufunc, by getting the helper function
    # (defined above) which returns a list of function(s) that convert the
    # input(s) to the unit required for the ufunc, as well as the unit the
    # result will have (a tuple of units if there are multiple outputs).
    try:
        ufunc_helper = UFUNC_HELPERS[function]
    except KeyError:
        if function in UNSUPPORTED_UFUNCS:
            raise TypeError("Cannot use function '{0}' with quantities"
                            .format(function.__name__))
        else:
            raise TypeError("Unknown ufunc {0}.  Please raise issue on "
                            "https://github.com/astropy/astropy"
                            .format(function.__name__))

    if method == '__call__' or (method == 'outer' and function.nin == 2):
        # Find out the units of the arguments passed to the ufunc; usually,
        # at least one is a quantity, but for two-argument ufuncs, the second
        # could also be a Numpy array, etc.  These are given unit=None.
        units = [getattr(arg, 'unit', None) for arg in args]

        # Determine possible conversion functions, and the result unit.
        converters, result_unit = ufunc_helper(function, *units)

        if any(converter is False for converter in converters):
            # for two-argument ufuncs with a quantity and a non-quantity,
            # the quantity normally needs to be dimensionless, *except*
            # if the non-quantity can have arbitrary unit, i.e., when it
            # is all zero, infinity or NaN.  In that case, the non-quantity
            # can just have the unit of the quantity
            # (this allows, e.g., `q > 0.` independent of unit)
            maybe_arbitrary_arg = args[converters.index(False)]
            try:
                if can_have_arbitrary_unit(maybe_arbitrary_arg):
                    converters = [None, None]
                else:
                    raise UnitsError("Can only apply '{0}' function to "
                                     "dimensionless quantities when other "
                                     "argument is not a quantity (unless the "
                                     "latter is all zero/infinity/nan)"
                                     .format(function.__name__))
            except TypeError:
                # _can_have_arbitrary_unit failed: arg could not be compared
                # with zero or checked to be finite. Then, ufunc will fail too.
                raise TypeError("Unsupported operand type(s) for ufunc {0}: "
                                "'{1}' and '{2}'"
                                .format(function.__name__,
                                        args[0].__class__.__name__,
                                        args[1].__class__.__name__))

        # In the case of np.power and np.float_power, the unit itself needs to
        # be modified by an amount that depends on one of the input values,
        # so we need to treat this as a special case.
        # TODO: find a better way to deal with this.
        if result_unit is False:
            if units[0] is None or units[0] == dimensionless_unscaled:
                result_unit = dimensionless_unscaled
            else:
                if units[1] is None:
                    p = args[1]
                else:
                    p = args[1].to(dimensionless_unscaled).value

                try:
                    result_unit = units[0] ** p
                except ValueError as exc:
                    # Changing the unit does not work for, e.g., array-shaped
                    # power, but this is OK if we're (scaled) dimensionless.
                    try:
                        converters[0] = units[0]._get_converter(
                            dimensionless_unscaled)
                    except UnitConversionError:
                        raise exc
                    else:
                        result_unit = dimensionless_unscaled

    else:  # methods for which the unit should stay the same
        nin = function.nin
        unit = getattr(args[0], 'unit', None)
        if method == 'at' and nin <= 2:
            if nin == 1:
                units = [unit]
            else:
                units = [unit, getattr(args[2], 'unit', None)]

            converters, result_unit = ufunc_helper(function, *units)

            # ensure there is no 'converter' for indices (2nd argument)
            converters.insert(1, None)

        elif method in {'reduce', 'accumulate', 'reduceat'} and nin == 2:
            converters, result_unit = ufunc_helper(function, unit, unit)
            converters = converters[:1]
            if method == 'reduceat':
                # add 'scale' for indices (2nd argument)
                converters += [None]

        else:
            if method in {'reduce', 'accumulate',
                          'reduceat', 'outer'} and nin != 2:
                raise ValueError("{0} only supported for binary functions"
                                 .format(method))

            raise TypeError("Unexpected ufunc method {0}.  If this should "
                            "work, please raise an issue on"
                            "https://github.com/astropy/astropy"
                            .format(method))

        # for all but __call__ method, scaling is not allowed
        if unit is not None and result_unit is None:
            raise TypeError("Cannot use '{1}' method on ufunc {0} with a "
                            "Quantity instance as the result is not a "
                            "Quantity.".format(function.__name__, method))

        if (converters[0] is not None or
            (unit is not None and unit is not result_unit and
             (not result_unit.is_equivalent(unit) or
              result_unit.to(unit) != 1.))):
            raise UnitsError("Cannot use '{1}' method on ufunc {0} with a "
                             "Quantity instance as it would change the unit."
                             .format(function.__name__, method))

    return converters, result_unit


def check_output(output, unit, inputs, function=None):
    """Check that function output can be stored in the output array given.

    Parameters
    ----------
    output : array or `~astropy.units.Quantity` or tuple
        Array that should hold the function output (or tuple of such arrays).
    unit : `~astropy.units.Unit` or None, or tuple
        Unit that the output will have, or `None` for pure numbers (should be
        tuple of same if output is a tuple of outputs).
    inputs : tuple
        Any input arguments.  These should be castable to the output.
    function : callable
        The function that will be producing the output.  If given, used to
        give a more informative error message.

    Returns
    -------
    arrays : `~numpy.ndarray` view of ``output`` (or tuple of such views).

    Raises
    ------
    UnitTypeError : If ``unit`` is inconsistent with the class of ``output``

    TypeError : If the ``inputs`` cannot be cast safely to ``output``.
    """
    if isinstance(output, tuple):
        return tuple(check_output(output_, unit_, inputs, function)
                     for output_, unit_ in zip(output, unit))

    # ``None`` indicates no actual array is needed.  This can happen, e.g.,
    # with np.modf(a, out=(None, b)).
    if output is None:
        return None

    if hasattr(output, '__quantity_subclass__'):
        # Check that we're not trying to store a plain Numpy array or a
        # Quantity with an inconsistent unit (e.g., not angular for Angle).
        if unit is None:
            raise TypeError("Cannot store non-quantity output{0} in {1} "
                            "instance".format(
                                (" from {0} function".format(function.__name__)
                                 if function is not None else ""),
                                type(output)))

        if output.__quantity_subclass__(unit)[0] is not type(output):
            raise UnitTypeError(
                "Cannot store output with unit '{0}'{1} "
                "in {2} instance.  Use {3} instance instead."
                .format(unit, (" from {0} function".format(function.__name__)
                               if function is not None else ""), type(output),
                        output.__quantity_subclass__(unit)[0]))

        # Turn into ndarray, so we do not loop into array_wrap/array_ufunc
        # if the output is used to store results of a function.
        output = output.view(np.ndarray)
    else:
        # output is not a Quantity, so cannot obtain a unit.
        if not (unit is None or unit is dimensionless_unscaled):
            raise UnitTypeError("Cannot store quantity with dimension "
                                "{0}in a non-Quantity instance."
                                .format("" if function is None else
                                        "resulting from {0} function "
                                        .format(function.__name__)))

    # check we can handle the dtype (e.g., that we are not int
    # when float is required).
    if not np.can_cast(np.result_type(*inputs), output.dtype,
                       casting='same_kind'):
        raise TypeError("Arguments cannot be cast safely to inplace "
                        "output with dtype={0}".format(output.dtype))
    return output
