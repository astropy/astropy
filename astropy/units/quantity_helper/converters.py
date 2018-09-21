# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Converters for Quantity."""

import numpy as np

from ..core import (UnitsError, UnitConversionError, UnitTypeError,
                    dimensionless_unscaled)

__all__ = ['can_have_arbitrary_unit', 'converters_and_unit',
           'check_output', 'UFUNC_HELPERS', 'UNSUPPORTED_UFUNCS']


class UfuncHelpers(dict):
    """Registry of unit conversion functions to help ufunc evaluation.

    Based on dict for quick access, but with a missing method to load
    helpers for additional modules such as scipy.special and erfa.

    Such modules should be registered using ``register_module``.
    """
    UNSUPPORTED = set()

    def register_module(self, module, names, importer):
        """Register (but do not import) a set of ufunc helpers.

        Parameters
        ----------
        module : str
            Name of the module with the ufuncs (e.g., 'scipy.special').
        names : iterable of str
            Names of the module ufuncs for which helpers are available.
        importer : callable
            Function that imports the ufuncs and returns a dict of helpers
            keyed by those ufuncs.  If the value is `None`, the ufunc is
            explicitly *not* supported.
        """
        self.modules[module] = {'names': names,
                                'importer': importer}

    @property
    def modules(self):
        """Modules for which helpers are available (but not yet loaded)."""
        if not hasattr(self, '_modules'):
            self._modules = {}
        return self._modules

    def import_module(self, module):
        """Import the helpers from the given module using its helper function.

        Parameters
        ----------
        module : str
            Name of the module. Has to have been registered beforehand.
        """
        module_info = self.modules.pop(module)
        self.update(module_info['importer']())

    def __missing__(self, ufunc):
        """Called if a ufunc is not found.

        Check if the ufunc is in any of the available modules, and, if so,
        import the helpers for that module.
        """
        if ufunc in self.UNSUPPORTED:
            raise TypeError("Cannot use ufunc '{0}' with quantities"
                            .format(ufunc.__name__))

        for module, module_info in list(self.modules.items()):
            if ufunc.__name__ in module_info['names']:
                # A ufunc with the same name is supported by this module.
                # Of course, this doesn't necessarily mean it is the
                # right module. So, we try let the importer do its work.
                # If it fails (e.g., for `scipy.special`), then that's
                # fine, just raise the TypeError.  If it succeeds, but
                # the ufunc is not found, that is also fine: we will
                # enter __missing__ again and either find another
                # module or get the TypeError there.
                try:
                    self.import_module(module)
                except ImportError:
                    pass
                else:
                    return self[ufunc]

        raise TypeError("unknown ufunc {0}.  If you believe this ufunc "
                        "should be supported, please raise an issue on "
                        "https://github.com/astropy/astropy"
                        .format(ufunc.__name__))

    def __setitem__(self, key, value):
        # Implementation note: in principle, we could just let `None`
        # mean that something is not implemented, but this means an
        # extra if clause for the output, slowing down the common
        # path where a ufunc is supported.
        if value is None:
            self.UNSUPPORTED |= {key}
            self.pop(key, None)
        else:
            super().__setitem__(key, value)
            self.UNSUPPORTED -= {key}


UFUNC_HELPERS = UfuncHelpers()
UNSUPPORTED_UFUNCS = UFUNC_HELPERS.UNSUPPORTED


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
    # (defined in helpers) which returns a list of function(s) that convert the
    # input(s) to the unit required for the ufunc, as well as the unit the
    # result will have (a tuple of units if there are multiple outputs).
    ufunc_helper = UFUNC_HELPERS[function]

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
