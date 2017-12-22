# The idea for this module (but no code) was borrowed from the
# quantities (http://pythonhosted.org/quantities/) package.

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


UFUNC_HELPERS = {}

# In this file, we implement the logic that determines for a given ufunc and
# input how the input should be scaled and what unit the output will have.

# list of ufuncs:
# http://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs

UNSUPPORTED_UFUNCS = set([np.bitwise_and, np.bitwise_or,
                          np.bitwise_xor, np.invert, np.left_shift,
                          np.right_shift, np.logical_and, np.logical_or,
                          np.logical_xor, np.logical_not])
if isinstance(getattr(np, 'isnat', None), np.ufunc):
    UNSUPPORTED_UFUNCS |= {np.isnat}

# SINGLE ARGUMENT UFUNCS

# The functions below take a single argument, which is the quantity upon which
# the ufunc is being used. The output of the function should be two values: the
# scale by which the input needs to be multiplied before being passed to the
# ufunc, and the unit the output will be in.

# ufuncs that return a boolean and do not care about the unit
helper_onearg_test = lambda f, unit: ([None], None)

UFUNC_HELPERS[np.isfinite] = helper_onearg_test
UFUNC_HELPERS[np.isinf] = helper_onearg_test
UFUNC_HELPERS[np.isnan] = helper_onearg_test
UFUNC_HELPERS[np.sign] = helper_onearg_test
UFUNC_HELPERS[np.signbit] = helper_onearg_test

# ufuncs that return a value with the same unit as the input

helper_invariant = lambda f, unit: ([None], _d(unit))

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
# positive only was added in numpy 1.13
if isinstance(getattr(np, 'positive', None), np.ufunc):
    UFUNC_HELPERS[np.positive] = helper_invariant

# ufuncs handled as special cases

UFUNC_HELPERS[np.sqrt] = lambda f, unit: (
    [None], unit ** 0.5 if unit is not None else dimensionless_unscaled)
UFUNC_HELPERS[np.square] = lambda f, unit: (
    [None], unit ** 2 if unit is not None else dimensionless_unscaled)
UFUNC_HELPERS[np.reciprocal] = lambda f, unit: (
    [None], unit ** -1 if unit is not None else dimensionless_unscaled)

UFUNC_HELPERS[np.cbrt] = lambda f, unit: (
    [None], (unit ** Fraction(1, 3) if unit is not None
             else dimensionless_unscaled))
UFUNC_HELPERS[np.core.umath._ones_like] = (lambda f, unit:
                                           ([None], dimensionless_unscaled))

# ufuncs that require dimensionless input and and give dimensionless output


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


UFUNC_HELPERS[np.exp] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.expm1] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.exp2] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log10] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log2] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log1p] = helper_dimensionless_to_dimensionless


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


UFUNC_HELPERS[np.modf] = helper_modf


# ufuncs that require dimensionless input and give output in radians
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
        return [get_converter(unit, degree)], radian
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


UFUNC_HELPERS[np.radians] = helper_degree_to_radian
UFUNC_HELPERS[np.deg2rad] = helper_degree_to_radian


# ufuncs that require input in radians and give output in degrees
def helper_radian_to_degree(f, unit):
    from .si import degree, radian
    try:
        return [get_converter(unit, radian)], degree
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


UFUNC_HELPERS[np.degrees] = helper_radian_to_degree
UFUNC_HELPERS[np.rad2deg] = helper_radian_to_degree


# ufuncs that require input in radians and give dimensionless output
def helper_radian_to_dimensionless(f, unit):
    from .si import radian
    try:
        return [get_converter(unit, radian)], dimensionless_unscaled
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


UFUNC_HELPERS[np.cos] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.sin] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.tan] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.cosh] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.sinh] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.tanh] = helper_radian_to_dimensionless


# ufuncs that require dimensionless_unscaled input and return non-quantities
def helper_frexp(f, unit):
    if not unit.is_unity():
        raise UnitTypeError("Can only apply '{0}' function to "
                            "unscaled dimensionless quantities"
                            .format(f.__name__))
    return [None], (None, None)


UFUNC_HELPERS[np.frexp] = helper_frexp


# TWO ARGUMENT UFUNCS
def helper_multiplication(f, unit1, unit2):
    return [None, None], _d(unit1) * _d(unit2)


UFUNC_HELPERS[np.multiply] = helper_multiplication


def helper_division(f, unit1, unit2):
    return [None, None], _d(unit1) / _d(unit2)


UFUNC_HELPERS[np.divide] = helper_division
UFUNC_HELPERS[np.true_divide] = helper_division


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


UFUNC_HELPERS[np.power] = helper_power
# float_power was added in numpy 1.12
if isinstance(getattr(np, 'float_power', None), np.ufunc):
    UFUNC_HELPERS[np.float_power] = helper_power


def helper_ldexp(f, unit1, unit2):
    if unit2 is not None:
        raise TypeError("Cannot use ldexp with a quantity "
                        "as second argument.")
    else:
        return [None, None], _d(unit1)


UFUNC_HELPERS[np.ldexp] = helper_ldexp


def helper_copysign(f, unit1, unit2):
    # if first arg is not a quantity, just return plain array
    if unit1 is None:
        return [None, None], None
    else:
        return [None, None], unit1


UFUNC_HELPERS[np.copysign] = helper_copysign

# heaviside only was added in numpy 1.13
if isinstance(getattr(np, 'heaviside', None), np.ufunc):
    def helper_heaviside(f, unit1, unit2):
        try:
            converter2 = (get_converter(unit2, dimensionless_unscaled)
                          if unit2 is not None else None)
        except UnitsError:
            raise UnitTypeError("Can only apply 'heaviside' function with a "
                                "dimensionless second argument.")
        return ([None, converter2], dimensionless_unscaled)

    UFUNC_HELPERS[np.heaviside] = helper_heaviside


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


UFUNC_HELPERS[np.logaddexp] = helper_two_arg_dimensionless
UFUNC_HELPERS[np.logaddexp2] = helper_two_arg_dimensionless


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


# This used to be a separate function that just called get_converters_and_unit.
# Using it directly saves a few us; keeping the clearer name.
helper_twoarg_invariant = get_converters_and_unit


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
    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, None


UFUNC_HELPERS[np.greater] = helper_twoarg_comparison
UFUNC_HELPERS[np.greater_equal] = helper_twoarg_comparison
UFUNC_HELPERS[np.less] = helper_twoarg_comparison
UFUNC_HELPERS[np.less_equal] = helper_twoarg_comparison
UFUNC_HELPERS[np.not_equal] = helper_twoarg_comparison
UFUNC_HELPERS[np.equal] = helper_twoarg_comparison


def helper_twoarg_invtrig(f, unit1, unit2):
    from .si import radian
    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, radian


UFUNC_HELPERS[np.arctan2] = helper_twoarg_invtrig
# another private function in numpy; use getattr in case it disappears
if isinstance(getattr(np.core.umath, '_arg', None), np.ufunc):
    UFUNC_HELPERS[np.core.umath._arg] = helper_twoarg_invtrig


def helper_twoarg_floor_divide(f, unit1, unit2):
    converters, _ = get_converters_and_unit(f, unit1, unit2)
    return converters, dimensionless_unscaled


UFUNC_HELPERS[np.floor_divide] = helper_twoarg_floor_divide

# divmod only was added in numpy 1.13
if isinstance(getattr(np, 'divmod', None), np.ufunc):
    def helper_divmod(f, unit1, unit2):
        converters, result_unit = get_converters_and_unit(f, unit1, unit2)
        return converters, (dimensionless_unscaled, result_unit)

    UFUNC_HELPERS[np.divmod] = helper_divmod


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
                # with zero or checked to be finite.  Then, ufunc will fail too.
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
        # output is not a Quantity, so cannot attain a unit.
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
