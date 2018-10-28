# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Function Units and Quantities."""

from abc import ABCMeta, abstractmethod

import numpy as np

from .. import (Unit, UnitBase, UnitsError, UnitTypeError,
                dimensionless_unscaled, Quantity)

__all__ = ['FunctionUnitBase', 'FunctionQuantity']

SUPPORTED_UFUNCS = set(getattr(np.core.umath, ufunc) for ufunc in (
    'isfinite', 'isinf', 'isnan', 'sign', 'signbit',
    'rint', 'floor', 'ceil', 'trunc', 'power',
    '_ones_like', 'ones_like', 'positive') if hasattr(np.core.umath, ufunc))

# TODO: the following could work if helper changed relative to Quantity:
# - spacing should return dimensionless, not same unit
# - negative should negate unit too,
# - add, subtract, comparisons can work if units added/subtracted

SUPPORTED_FUNCTIONS = set(getattr(np, function) for function in
                          ('clip', 'trace', 'mean', 'min', 'max', 'round'))


# subclassing UnitBase or CompositeUnit was found to be problematic, requiring
# a large number of overrides. Hence, define new class.
class FunctionUnitBase(metaclass=ABCMeta):
    """Abstract base class for function units.

    Function units are functions containing a physical unit, such as dB(mW).
    Most of the arithmetic operations on function units are defined in this
    base class.

    While instantiation is defined, this class should not be used directly.
    Rather, subclasses should be used that override the abstract properties
    `_default_function_unit` and `_quantity_class`, and the abstract methods
    `from_physical`, and `to_physical`.

    Parameters
    ----------
    physical_unit : `~astropy.units.Unit` or `string`
        Unit that is encapsulated within the function unit.
        If not given, dimensionless.

    function_unit :  `~astropy.units.Unit` or `string`
        By default, the same as the function unit set by the subclass.
    """
    # ↓↓↓ the following four need to be set by subclasses
    # Make this a property so we can ensure subclasses define it.
    @property
    @abstractmethod
    def _default_function_unit(self):
        """Default function unit corresponding to the function.

        This property should be overridden by subclasses, with, e.g.,
        `~astropy.unit.MagUnit` returning `~astropy.unit.mag`.
        """

    # This has to be a property because the function quantity will not be
    # known at unit definition time, as it gets defined after.
    @property
    @abstractmethod
    def _quantity_class(self):
        """Function quantity class corresponding to this function unit.

        This property should be overridden by subclasses, with, e.g.,
        `~astropy.unit.MagUnit` returning `~astropy.unit.Magnitude`.
        """

    @abstractmethod
    def from_physical(self, x):
        """Transformation from value in physical to value in function units.

        This method should be overridden by subclasses.  It is used to
        provide automatic transformations using an equivalency.
        """

    @abstractmethod
    def to_physical(self, x):
        """Transformation from value in function to value in physical units.

        This method should be overridden by subclasses.  It is used to
        provide automatic transformations using an equivalency.
        """
    # ↑↑↑ the above four need to be set by subclasses

    # have priority over arrays, regular units, and regular quantities
    __array_priority__ = 30000

    def __init__(self, physical_unit=None, function_unit=None):
        if physical_unit is None:
            self._physical_unit = dimensionless_unscaled
        else:
            self._physical_unit = Unit(physical_unit)
            if (not isinstance(self._physical_unit, UnitBase) or
                self._physical_unit.is_equivalent(
                    self._default_function_unit)):
                raise ValueError("Unit {0} is not a physical unit."
                                 .format(self._physical_unit))

        if function_unit is None:
            self._function_unit = self._default_function_unit
        else:
            # any function unit should be equivalent to subclass default
            function_unit = Unit(getattr(function_unit, 'function_unit',
                                         function_unit))
            if function_unit.is_equivalent(self._default_function_unit):
                self._function_unit = function_unit
            else:
                raise ValueError("Cannot initialize '{0}' instance with "
                                 "function unit '{1}', as it is not "
                                 "equivalent to default function unit '{2}'."
                                 .format(self.__class__.__name__,
                                         function_unit,
                                         self._default_function_unit))

    def _copy(self, physical_unit=None):
        """Copy oneself, possibly with a different physical unit."""
        if physical_unit is None:
            physical_unit = self.physical_unit
        return self.__class__(physical_unit, self.function_unit)

    @property
    def physical_unit(self):
        return self._physical_unit

    @property
    def function_unit(self):
        return self._function_unit

    @property
    def equivalencies(self):
        """List of equivalencies between function and physical units.

        Uses the `from_physical` and `to_physical` methods.
        """
        return [(self, self.physical_unit,
                 self.to_physical, self.from_physical)]

    # ↓↓↓ properties/methods required to behave like a unit
    def decompose(self, bases=set()):
        """Copy the current unit with the physical unit decomposed.

        For details, see `~astropy.units.UnitBase.decompose`.
        """
        return self._copy(self.physical_unit.decompose(bases))

    @property
    def si(self):
        """Copy the current function unit with the physical unit in SI."""
        return self._copy(self.physical_unit.si)

    @property
    def cgs(self):
        """Copy the current function unit with the physical unit in CGS."""
        return self._copy(self.physical_unit.cgs)

    def _get_physical_type_id(self):
        """Get physical type corresponding to physical unit."""
        return self.physical_unit._get_physical_type_id()

    @property
    def physical_type(self):
        """Return the physical type of the physical unit (e.g., 'length')."""
        return self.physical_unit.physical_type

    def is_equivalent(self, other, equivalencies=[]):
        """
        Returns `True` if this unit is equivalent to ``other``.

        Parameters
        ----------
        other : unit object or string or tuple
            The unit to convert to. If a tuple of units is specified, this
            method returns true if the unit matches any of those in the tuple.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            This list is in addition to the built-in equivalencies between the
            function unit and the physical one, as well as possible global
            defaults set by, e.g., `~astropy.units.set_enabled_equivalencies`.
            Use `None` to turn off any global equivalencies.

        Returns
        -------
        bool
        """
        if isinstance(other, tuple):
            return any(self.is_equivalent(u, equivalencies=equivalencies)
                       for u in other)

        other_physical_unit = getattr(other, 'physical_unit', (
            dimensionless_unscaled if self.function_unit.is_equivalent(other)
            else other))

        return self.physical_unit.is_equivalent(other_physical_unit,
                                                equivalencies)

    def to(self, other, value=1., equivalencies=[]):
        """
        Return the converted values in the specified unit.

        Parameters
        ----------
        other : `~astropy.units.Unit` object, `~astropy.units.function.FunctionUnitBase` object or string
            The unit to convert to.

        value : scalar int or float, or sequence convertible to array, optional
            Value(s) in the current unit to be converted to the specified unit.
            If not provided, defaults to 1.0.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            This list is in meant to treat only equivalencies between different
            physical units; the build-in equivalency between the function
            unit and the physical one is automatically taken into account.

        Returns
        -------
        values : scalar or array
            Converted value(s). Input value sequences are returned as
            numpy arrays.

        Raises
        ------
        UnitsError
            If units are inconsistent.
        """
        # conversion to one's own physical unit should be fastest
        if other is self.physical_unit:
            return self.to_physical(value)

        other_function_unit = getattr(other, 'function_unit', other)
        if self.function_unit.is_equivalent(other_function_unit):
            # when other is an equivalent function unit:
            # first convert physical units to other's physical units
            other_physical_unit = getattr(other, 'physical_unit',
                                          dimensionless_unscaled)
            if self.physical_unit != other_physical_unit:
                value_other_physical = self.physical_unit.to(
                    other_physical_unit, self.to_physical(value),
                    equivalencies)
                # make function unit again, in own system
                value = self.from_physical(value_other_physical)

            # convert possible difference in function unit (e.g., dex->dB)
            return self.function_unit.to(other_function_unit, value)

        else:
            # when other is not a function unit
            return self.physical_unit.to(other, self.to_physical(value),
                                         equivalencies)

    def is_unity(self):
        return False

    def __eq__(self, other):
        return (self.physical_unit == getattr(other, 'physical_unit',
                                              dimensionless_unscaled) and
                self.function_unit == getattr(other, 'function_unit', other))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __mul__(self, other):
        if isinstance(other, (str, UnitBase, FunctionUnitBase)):
            if self.physical_unit == dimensionless_unscaled:
                # If dimensionless, drop back to normal unit and retry.
                return self.function_unit * other
            else:
                raise UnitsError("Cannot multiply a function unit "
                                 "with a physical dimension with any unit.")
        else:
            # Anything not like a unit, try initialising as a function quantity.
            try:
                return self._quantity_class(other, unit=self)
            except Exception:
                return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __div__(self, other):
        if isinstance(other, (str, UnitBase, FunctionUnitBase)):
            if self.physical_unit == dimensionless_unscaled:
                # If dimensionless, drop back to normal unit and retry.
                return self.function_unit / other
            else:
                raise UnitsError("Cannot divide a function unit "
                                 "with a physical dimension by any unit.")
        else:
            # Anything not like a unit, try initialising as a function quantity.
            try:
                return self._quantity_class(1./other, unit=self)
            except Exception:
                return NotImplemented

    def __rdiv__(self, other):
        if isinstance(other, (str, UnitBase, FunctionUnitBase)):
            if self.physical_unit == dimensionless_unscaled:
                # If dimensionless, drop back to normal unit and retry.
                return other / self.function_unit
            else:
                raise UnitsError("Cannot divide a function unit "
                                 "with a physical dimension into any unit")
        else:
            # Don't know what to do with anything not like a unit.
            return NotImplemented

    __truediv__ = __div__

    __rtruediv__ = __rdiv__

    def __pow__(self, power):
        if power == 0:
            return dimensionless_unscaled
        elif power == 1:
            return self._copy()

        if self.physical_unit == dimensionless_unscaled:
            return self.function_unit ** power

        raise UnitsError("Cannot raise a function unit "
                         "with a physical dimension to any power but 0 or 1.")

    def __pos__(self):
        return self._copy()

    def to_string(self, format='generic'):
        """
        Output the unit in the given format as a string.

        The physical unit is appended, within parentheses, to the function
        unit, as in "dB(mW)", with both units set using the given format

        Parameters
        ----------
        format : `astropy.units.format.Base` instance or str
            The name of a format or a formatter object.  If not
            provided, defaults to the generic format.
        """
        if format not in ('generic', 'unscaled', 'latex'):
            raise ValueError("Function units cannot be written in {0} format. "
                             "Only 'generic', 'unscaled' and 'latex' are "
                             "supported.".format(format))
        self_str = self.function_unit.to_string(format)
        pu_str = self.physical_unit.to_string(format)
        if pu_str == '':
            pu_str = '1'
        if format == 'latex':
            self_str += r'$\mathrm{{\left( {0} \right)}}$'.format(
                pu_str[1:-1])   # need to strip leading and trailing "$"
        else:
            self_str += '({0})'.format(pu_str)
        return self_str

    def __str__(self):
        """Return string representation for unit."""
        self_str = str(self.function_unit)
        pu_str = str(self.physical_unit)
        if pu_str:
            self_str += '({0})'.format(pu_str)
        return self_str

    def __repr__(self):
        # By default, try to give a representation using `Unit(<string>)`,
        # with string such that parsing it would give the correct FunctionUnit.
        if callable(self.function_unit):
            return 'Unit("{0}")'.format(self.to_string())

        else:
            return '{0}("{1}"{2})'.format(
                self.__class__.__name__, self.physical_unit,
                "" if self.function_unit is self._default_function_unit
                else ', unit="{0}"'.format(self.function_unit))

    def _repr_latex_(self):
        """
        Generate latex representation of unit name.  This is used by
        the IPython notebook to print a unit with a nice layout.

        Returns
        -------
        Latex string
        """
        return self.to_string('latex')

    def __hash__(self):
        return hash((self.function_unit, self.physical_unit))


class FunctionQuantity(Quantity):
    """A representation of a (scaled) function of a number with a unit.

    Function quantities are quantities whose units are functions containing a
    physical unit, such as dB(mW).  Most of the arithmetic operations on
    function quantities are defined in this base class.

    While instantiation is also defined here, this class should not be
    instantiated directly.  Rather, subclasses should be made which have
    ``_unit_class`` pointing back to the corresponding function unit class.

    Parameters
    ----------
    value : number, sequence of convertible items, `~astropy.units.Quantity`, or `~astropy.units.function.FunctionQuantity`
        The numerical value of the function quantity. If a number or
        a `~astropy.units.Quantity` with a function unit, it will be converted
        to ``unit`` and the physical unit will be inferred from ``unit``.
        If a `~astropy.units.Quantity` with just a physical unit, it will
        converted to the function unit, after, if necessary, converting it to
        the physical unit inferred from ``unit``.

    unit : string, `~astropy.units.UnitBase` or `~astropy.units.function.FunctionUnitBase` instance, optional
        For an `~astropy.units.function.FunctionUnitBase` instance, the
        physical unit will be taken from it; for other input, it will be
        inferred from ``value``. By default, ``unit`` is set by the subclass.

    dtype : `~numpy.dtype`, optional
        The dtype of the resulting Numpy array or scalar that will
        hold the value.  If not provided, it is determined from the input,
        except that any input that cannot represent float (integer and bool)
        is converted to float.

    copy : bool, optional
        If `True` (default), then the value is copied.  Otherwise, a copy will
        only be made if ``__array__`` returns a copy, if value is a nested
        sequence, or if a copy is needed to satisfy an explicitly given
        ``dtype``.  (The `False` option is intended mostly for internal use,
        to speed up initialization where a copy is known to have been made.
        Use with care.)

    order : {'C', 'F', 'A'}, optional
        Specify the order of the array.  As in `~numpy.array`.  Ignored
        if the input does not need to be converted and ``copy=False``.

    subok : bool, optional
        If `False` (default), the returned array will be forced to be of the
        class used.  Otherwise, subclasses will be passed through.

    ndmin : int, optional
        Specifies the minimum number of dimensions that the resulting array
        should have.  Ones will be pre-pended to the shape as needed to meet
        this requirement.  This parameter is ignored if the input is a
        `~astropy.units.Quantity` and ``copy=False``.

    Raises
    ------
    TypeError
        If the value provided is not a Python numeric type.
    TypeError
        If the unit provided is not a `~astropy.units.function.FunctionUnitBase`
        or `~astropy.units.Unit` object, or a parseable string unit.
    """

    _unit_class = None
    """Default `~astropy.units.function.FunctionUnitBase` subclass.

    This should be overridden by subclasses.
    """

    # Ensure priority over ndarray, regular Unit & Quantity, and FunctionUnit.
    __array_priority__ = 40000

    # Define functions that work on FunctionQuantity.
    _supported_ufuncs = SUPPORTED_UFUNCS
    _supported_functions = SUPPORTED_FUNCTIONS

    def __new__(cls, value, unit=None, dtype=None, copy=True, order=None,
                subok=False, ndmin=0):

        if unit is not None:
            # Convert possible string input to a (function) unit.
            unit = Unit(unit)

        if not isinstance(unit, FunctionUnitBase):
            # By default, use value's physical unit.
            value_unit = getattr(value, 'unit', None)
            if value_unit is None:
                # if iterable, see if first item has a unit
                # (mixed lists fail in super call below).
                try:
                    value_unit = getattr(value[0], 'unit')
                except Exception:
                    pass
            physical_unit = getattr(value_unit, 'physical_unit', value_unit)
            unit = cls._unit_class(physical_unit, function_unit=unit)

        # initialise!
        return super().__new__(cls, value, unit, dtype=dtype, copy=copy,
                               order=order, subok=subok, ndmin=ndmin)

    # ↓↓↓ properties not found in Quantity
    @property
    def physical(self):
        """The physical quantity corresponding the function one."""
        return self.to(self.unit.physical_unit)

    @property
    def _function_view(self):
        """View as Quantity with function unit, dropping the physical unit.

        Use `~astropy.units.quantity.Quantity.value` for just the value.
        """
        return self._new_view(unit=self.unit.function_unit)

    # ↓↓↓ methods overridden to change the behavior
    @property
    def si(self):
        """Return a copy with the physical unit in SI units."""
        return self.__class__(self.physical.si)

    @property
    def cgs(self):
        """Return a copy with the physical unit in CGS units."""
        return self.__class__(self.physical.cgs)

    def decompose(self, bases=[]):
        """Generate a new `FunctionQuantity` with the physical unit decomposed.

        For details, see `~astropy.units.Quantity.decompose`.
        """
        return self.__class__(self.physical.decompose(bases))

    # ↓↓↓ methods overridden to add additional behavior
    def __quantity_subclass__(self, unit):
        if isinstance(unit, FunctionUnitBase):
            return self.__class__, True
        else:
            return super().__quantity_subclass__(unit)[0], False

    def _set_unit(self, unit):
        if not isinstance(unit, self._unit_class):
            # Have to take care of, e.g., (10*u.mag).view(u.Magnitude)
            try:
                # "or 'nonsense'" ensures `None` breaks, just in case.
                unit = self._unit_class(function_unit=unit or 'nonsense')
            except Exception:
                raise UnitTypeError(
                    "{0} instances require {1} function units"
                    .format(type(self).__name__, self._unit_class.__name__) +
                    ", so cannot set it to '{0}'.".format(unit))

        self._unit = unit

    # ↓↓↓ methods overridden to change behavior
    def __mul__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self._function_view * other

        raise UnitTypeError("Cannot multiply function quantities which "
                            "are not dimensionless with anything.")

    def __truediv__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self._function_view / other

        raise UnitTypeError("Cannot divide function quantities which "
                            "are not dimensionless by anything.")

    def __rtruediv__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self._function_view.__rdiv__(other)

        raise UnitTypeError("Cannot divide function quantities which "
                            "are not dimensionless into anything.")

    def _comparison(self, other, comparison_func):
        """Do a comparison between self and other, raising UnitsError when
        other cannot be converted to self because it has different physical
        unit, and returning NotImplemented when there are other errors."""
        try:
            # will raise a UnitsError if physical units not equivalent
            other_in_own_unit = self._to_own_unit(other, check_precision=False)
        except UnitsError as exc:
            if self.unit.physical_unit != dimensionless_unscaled:
                raise exc

            try:
                other_in_own_unit = self._function_view._to_own_unit(
                    other, check_precision=False)
            except Exception:
                raise exc

        except Exception:
            return NotImplemented

        return comparison_func(other_in_own_unit)

    def __eq__(self, other):
        try:
            return self._comparison(other, self.value.__eq__)
        except UnitsError:
            return False

    def __ne__(self, other):
        try:
            return self._comparison(other, self.value.__ne__)
        except UnitsError:
            return True

    def __gt__(self, other):
        return self._comparison(other, self.value.__gt__)

    def __ge__(self, other):
        return self._comparison(other, self.value.__ge__)

    def __lt__(self, other):
        return self._comparison(other, self.value.__lt__)

    def __le__(self, other):
        return self._comparison(other, self.value.__le__)

    # Ensure Quantity methods are used only if they make sense.
    def _wrap_function(self, function, *args, **kwargs):
        if function in self._supported_functions:
            return super()._wrap_function(function, *args, **kwargs)

        # For dimensionless, we can convert to regular quantities.
        if all(arg.unit.physical_unit == dimensionless_unscaled
               for arg in (self,) + args
               if (hasattr(arg, 'unit') and
                   hasattr(arg.unit, 'physical_unit'))):
            args = tuple(getattr(arg, '_function_view', arg) for arg in args)
            return self._function_view._wrap_function(function, *args, **kwargs)

        raise TypeError("Cannot use method that uses function '{0}' with "
                        "function quantities that are not dimensionless."
                        .format(function.__name__))
