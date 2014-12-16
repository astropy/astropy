# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Function Units and Quantities."""
from __future__ import (absolute_import, unicode_literals,
                        division, print_function)

from abc import ABCMeta, abstractmethod, abstractproperty

import numpy as np
# LOCAL
from ...extern import six
from .. import (Unit, UnitBase, UnitsError,
                dimensionless_unscaled, Quantity, quantity_helper as qh)

__all__ = ['FunctionUnitBase', 'FunctionQuantityBase']

_supported_ufuncs = set([ufunc for ufunc, helper in qh.UFUNC_HELPERS.items()
                         if helper in (qh.helper_onearg_test,
                                       qh.helper_invariant,
                                       qh.helper_division,
                                       qh.helper_copysign)])
_supported_ufuncs |= set([getattr(np.core.umath, ufunc)
                          for ufunc in ('sqrt', 'cbrt', 'square', 'reciprocal',
                                        'multiply', 'power',
                                        '_ones_like', 'ones_like')
                          if hasattr(np.core.umath, ufunc)])


# subclassing UnitBase or CompositeUnit was found to be problematic, requiring
# a large number of overrides. Hence, define new class.
@six.add_metaclass(ABCMeta)
class FunctionUnitBase(object):
    """Abstract base class for function units

    Function units are functions containing a physical unit, such as dB(mW)
    Most of the arithmetic operations on function units are defined in this
    base class.

    While instantiation is defined, this class should not be used directly.
    Rather, subclasses should be used that define the properties
    `_default_function_unit` and `_quantity_class`, and the methods
    `from_physical`, and `to_physical`.

    Parameters
    ----------
    physical_unit : `~astropy.units.Unit` or `string`
        Unit that is encapsulated within the function unit.
        If not given, dimensionless.

    function_unit :  `~astropy.units.Unit` or `string`
        By default, the same as the function unit set by the subclass.

    """
    # vvvvv the following four need to be set by subclasses
    # Make this a property so we can ensure subclasses define it.
    @abstractproperty
    def _default_function_unit(self):
        """Default function unit corresponding to the function.

        This property should be overridden by subclasses, with, e.g.,
        `~astropy.unit.MagUnit` returning `~astropy.unit.mag`.
        """

    # This has to be a property because the function quantity will not be
    # known at unit definition time, as it gets defined after.
    @abstractproperty
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
    # ^^^^^ the above four need to be set by subclasses

    # have priority over arrays, regular units, and regular quantities
    __array_priority__ = 3000

    def __init__(self, physical_unit=None, function_unit=None):
        if physical_unit is None:
            self._physical_unit = dimensionless_unscaled
        else:
            self._physical_unit = Unit(physical_unit)

        if function_unit is None:
            self._function_unit = self._default_function_unit
        else:
            # any function unit should be equivalent to subclass default
            function_unit = Unit(getattr(function_unit, 'function_unit',
                                         function_unit))
            if function_unit.is_equivalent(self._default_function_unit):
                self._function_unit = function_unit
            else:
                raise ValueError("Cannot initiliaze '{0}' instance with "
                                 "function unit '{1}', as it is not "
                                 "equivalent to default unit '{2}'"
                                 .format(self.__class__.__name__,
                                         function_unit,
                                         self._default_function_unit))


    def _copy(self, physical_unit=None):
        """Copy oneself, possibly changing physical unit"""
        if physical_unit is None:
            physical_unit = self.physical_unit
        return self.__class__(physical_unit, self.function_unit)

    @property
    def physical_unit(self):
        return self._physical_unit

    @property
    def function_unit(self):
        return self._function_unit

    def function_quantity(self, value):
        return self._quantity_class(value, unit=self)

    @property
    def equivalencies(self):
        """List of equivalencies between physical and function units.

        Uses the `from_physical` and `to_physical` methods.
        """
        return [(self.physical_unit, self,
                 self.from_physical, self.to_physical)]

    # vvvv properties/methods required to behave like a unit
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
        """Get physical type corresponding to physical unit"""
        return self.physical_unit._get_physical_type_id()

    @property
    def physical_type(self):
        """Return the physical type of the physical unit (e.g., 'length')"""
        return self.physical_unit.physical_type

    def _to(self, other):
        """
        Returns the scale to the specified function unit, raising
        UnitsError is the physical units are not equivalent.

        This is required to mimic behaviour expected for any units, e.g.,
        in `~astropy.units.core.UnitBase.apply_equivalencies`.
        """
        if not self.physical_unit._is_equivalent(
                getattr(other, 'physical_unit', other)):
            raise UnitsError("'{0!r}' is not a scaled version of '{1!r}'"
                             .format(self, other))

        return self.function_unit._to(getattr(other, 'function_unit', other))

    def is_equivalent(self, other, equivalencies=[]):
        if isinstance(other, tuple):
            return any(self.is_equivalent(u, equivalencies=equivalencies)
                       for u in other)

        other_physical = getattr(other, 'physical_unit',
                                 (dimensionless_unscaled
                                  if self.function_unit.is_equivalent(other)
                                  else other))

        return self.physical_unit.is_equivalent(other_physical, equivalencies)
    is_equivalent.__doc__ = UnitBase.is_equivalent.__doc__

    def to(self, other, value=1., equivalencies=[]):
        """
        Return the converted values in the specified unit.

        Parameters
        ----------
        other : unit object or string
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
            If units are inconsistent
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
                value_other_phys = self.physical_unit.to(
                    other_physical_unit, self.to_physical(value),
                    equivalencies)
                # make function unit again, in own system
                value = self.from_physical(value_other_phys)

            # convert possible difference in function unit (e.g., dex->dB)
            return self.function_unit.to(other_function_unit, value)

        else:
            # when other is not a function unit
            return self.physical_unit.to(other, self.to_physical(value),
                                         equivalencies)

    def __eq__(self, other):
        return (self.physical_unit == getattr(other, 'physical_unit',
                                              dimensionless_unscaled) and
                self.function_unit == getattr(other, 'function_unit', other))

    def __ne__(self, other):
        return (self.physical_unit != getattr(other, 'physical_unit',
                                              dimensionless_unscaled) or
                self.funtional_unit != getattr(other, 'function_unit', other))

    def __mul__(self, other):
        # anything not like a unit, try initialising as a function quantity
        if not isinstance(other, (six.string_types, UnitBase,
                                  FunctionUnitBase)):
            try:
                return self._quantity_class(other, unit=self)
            except UnitsError:
                return NotImplemented

        if self.physical_unit == dimensionless_unscaled:
            # if dimensionless, drop back to normal unit and retry
            return self.function_unit * other

        raise UnitsError("Cannot multiply a function unit "
                         "with a physical dimension with any unit")

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        # anything not like a unit, try initialising as a function quantity
        if not isinstance(other, (six.string_types, UnitBase,
                                  FunctionUnitBase)):
            try:
                return self.function_quantity(1./other)
            except UnitsError:
                return NotImplemented

        if self.physical_unit == dimensionless_unscaled:
            # if dimensionless, drop back to normal unit and retry
            return self.function_unit / other

        raise UnitsError("Cannot divide a function unit "
                         "with a physical dimension by any unit")

    def __rtruediv__(self, other):
        if not isinstance(other, (six.string_types, UnitBase,
                                  FunctionUnitBase)):
            return NotImplemented

        if self.physical_unit == dimensionless_unscaled:
            # if dimensionless, drop back to normal unit and retry
            return other / self.function_unit

        raise UnitsError("Cannot divide a function unit "
                         "with a physical dimension into any unit")

    def __pow__(self, power):
        if power == 0:
            return dimensionless_unscaled
        elif power == 1:
            return self._copy()

        if self.physical_unit == dimensionless_unscaled:
            return self.function_unit ** power

        raise UnitsError("Cannot raise a function unit "
                         "with a physical dimension to any power but 0 or 1")

    def __pos__(self):
        return self._copy()

    def to_string(self, format='generic'):
        """
        Output the unit in the given format as a string

        The physical unit is appended, within parentheses, to the function
        unit, as in "dB(mW)", with both units set using the given format

        Parameters
        ----------
        format : `astropy.format.Base` instance or str
            The name of a format or a formatter object.  If not
            provided, defaults to the generic format.
        """
        self_str = self.function_unit.to_string(format)
        if self.physical_unit != dimensionless_unscaled:
            self_str += '({0})'.format(self.physical_unit.to_string(format))
        return self_str

    def __str__(self):
        """Return string representation for unit"""
        self_str = str(self.function_unit)
        if self.physical_unit != dimensionless_unscaled:
            self_str += '({0})'.format(self.physical_unit)

        return self_str

    def __repr__(self):
        return "{0}('{1}'{2})".format(
            self.__class__.__name__, self.physical_unit,
            "" if self.function_unit is self._default_function_unit
            else ", unit='{0}'".format(self.function_unit))

    def __hash__(self):
        return hash((self.function_unit, self.physical_unit))


class FunctionQuantityBase(Quantity):
    """A representation of a (scaled) function of a number with a unit

    Function quanties are quantities whose units are functions containing a
    physical unit, such as dB(mW).  Most of the arithmetic operations on
    function quantities are defined in this base class.

    While instantiation is also defined here, this class should not be
    instantiated by users directly.  Rather, subclasses should be made which
    have `_FunctionUnit` pointing back to the corresponding function unit class.

    Parameters
    ----------
    value : number, `~astropy.units.Quantity`,
            `~astropy.units.FunctionQuantity`,
            or sequence of convertible items.
        The numerical value of the function quantity. If a number or
        a `Quantity` with a function unit, it will be converted to `unit`
        and the physical unit will be inferred from `unit`.
        If a `Quantity` with just a physical unit, it will converted to
        the function unit, after, if necessary, converting it to the
        physical unit inferred from `unit`.

    unit : function unit, or `~astropy.units.function.FunctionUnit`,
            optional
        E.g., `~astropy.units.mag`, `~astropy.units.dex`, or equivalent
        function units `~astropy.units.MagUnit`, `astropy.units.DexUnit`, etc.
        For a `FunctionUnit` instance, the physical unit will be taken from
        it; for non-`FunctionUnit` input, it will be inferred from `value`.
        By default, `unit` is that defined by the subclass.

    dtype : ~numpy.dtype, optional
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
        `Quantity` and ``copy=False``.

    Raises
    ------
    TypeError
        If the value provided is not a Python numeric type.
    TypeError
        If the unit provided is not a :class:`~astropy.units.FunctionUnit`
        or :class:`~astropy.units.Unit` object, or a parseable string unit.
    """

    # default FunctionUnit; to be overridden by subclasses
    _FunctionUnit = NotImplemented

    # have priority over ndarray, regular Unit & Quantity, and FunctionUnit
    __array_priority__ = 4000

    # store unit in different private property than is done by Quantity,
    # to allow a Quantity view with just the function unit (stored in _unit)
    _full_unit = None

    def __new__(cls, value, unit=None, dtype=None, copy=True, order=None,
                subok=False, ndmin=0):

        value_unit = getattr(value, 'unit', None)
        if value_unit is None:
            # if iterable, see if first item has a unit; mixed lists fail below
            try:
                value_unit = getattr(value[0], 'unit')
            except:
                pass

        if isinstance(value_unit, FunctionUnitBase):
            if unit is None:
                unit = value_unit
                value = value.function_value  # for initialising Quantity
            else:
                # convert value to its quantity (will convert back below)
                value_unit = value.physical_unit
                value = value.quantity

        if not isinstance(unit, FunctionUnitBase):
            unit = cls._FunctionUnit(value_unit, function_unit=unit)

        if value_unit is not None and value_unit is not unit:
            # convert to target unit
            value = unit.from_physical(value.to(unit.physical_unit).value)

        # initialise Quantity
        self = super(FunctionQuantityBase, cls).__new__(
            cls, value, unit.function_unit, dtype=dtype, copy=copy,
            order=order, subok=subok, ndmin=ndmin)
        # set the full unit returned by self.unit
        self._full_unit = unit
        return self

    # vvvv properties not found in Quantity
    @property
    def quantity(self):
        """The physical quantity corresponding the function one."""
        return self.to(self.unit.physical_unit)

    @property
    def function_value(self):
        """Quantity with plain function unit, dropping the physical unit.

        Use `~astropy.units.quantity.Quantity.value` for just the value.
        """
        return self._new_view(self, self.unit.function_unit)

    # vvvv properties overridden to point to different places
    @property
    def unit(self):
        """Function unit of the quantity, containing the physical unit.

        Represented by a `~astropy.units.function.FunctionUnitBase` object.
        """
        return self._full_unit

    @property
    def equivalencies(self):
        """Equivalencies applied by default during unit conversions.

        Contains the list to convert between function and physical unit,
        as set by the `~astropy.units.function.FunctionUnitBase` unit.
        """
        return self.unit.equivalencies

    # vvvv methods overridden to change the behaviour
    @property
    def si(self):
        return self.__class__(self.quantity.si)

    @property
    def cgs(self):
        return self.__class__(self.quantity.cgs)

    def decompose(self, bases=[]):
        return self.__class__(self.quantity.decompose(bases))

    # vvvv methods overridden to add additional behaviour
    def to(self, unit, equivalencies=[]):
        result = super(FunctionQuantityBase, self).to(unit, equivalencies)
        if isinstance(unit, FunctionUnitBase):
            result._full_unit = unit
        return result

    def __array_finalize__(self, obj):
        if isinstance(obj, FunctionQuantityBase):
            self._full_unit = obj._full_unit

        super(FunctionQuantityBase, self).__array_finalize__(obj)

    def __array_prepare__(self, obj, context=None):
        """Check that the ufunc can deal with a FunctionQuantity"""

        # If no context is set, just return the input
        if context is None:
            return obj

        # Find out whether ufunc is supported
        function = context[0]
        if not (function in _supported_ufuncs or
                all(arg.unit.physical_unit == dimensionless_unscaled
                    for arg in context[1][:function.nin]
                    if (hasattr(arg, 'unit') and
                        hasattr(arg.unit, 'physical_unit')))):
            raise TypeError("Cannot use function '{0}' with function "
                            "quantities that are not dimensionless"
                            .format(context[0].__name__))

        return super(FunctionQuantityBase,
                     self).__array_prepare__(obj, context)

    def __array_wrap__(self, obj, context=None):
        """Reorder units properly after ufunc calculation"""

        obj = super(FunctionQuantityBase, self).__array_wrap__(obj, context)

        if isinstance(obj, FunctionQuantityBase):
            obj._full_unit = obj._unit
            obj._unit = obj._full_unit.function_unit
            # has the function unit changed? (e.g., "mag" -> "2 mag")
            if obj._unit != obj._full_unit._default_function_unit:
                if obj._full_unit.physical_type is not dimensionless_unscaled:
                    raise UnitsError(
                        "Unexpected production in "
                        "FunctionQuantityBase.__array_wrap__ "
                        "of '{0}' instance with a physical dimension yet "
                        "function_unit '{1}' when '{2}' was expected. "
                        "Please alert the astropy developers"
                        .format(obj.__class__.__name__, obj._unit,
                                obj._full_unit._default_function_unit))

                if not obj._unit.is_equivalent(
                        obj._full_unit._default_function_unit):
                    # for dimensionless, fall back to regular quantity
                    obj = self._new_view(obj, obj._unit)

        return obj

    def __quantity_subclass__(self, unit):
        if isinstance(unit, FunctionUnitBase):
            return self.__class__, True
        else:
            return super(FunctionQuantityBase,
                         self).__quantity_subclass__(unit)[0], False

    # vvvv methods overridden to change behaviour
    def __mul__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self.function_value * other

        raise UnitsError("Cannot multiply function quantities which "
                         "are not dimensionless with anything")

    def __div__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self.function_value / other

        raise UnitsError("Cannot divide function quantities which "
                         "are not dimensionless by anything")

    def __rdiv__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self.function_value.__rdiv__(other)

        raise UnitsError("Cannot divide function quantities which "
                         "are not dimensionless into anything")

    def _comparison(self, other, comparison_func):
        """Do a comparison between self and other, raising UnitsError when
        other cannot be converted to self because it has different physical
        unit, and returning NotImplemented when there are other errors."""
        try:
            # will raise a UnitsError if physical units not equivalent
            return comparison_func(other.to(self.unit).value)
        except UnitsError:
            raise
        except:
            return NotImplemented

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

    # vvvv override Quantity methods that make no sense for function units
    def _not_implemented(self, *args, **kwargs):
        raise NotImplementedError

    dot = nansum = sum = cumsum = prod = cumprod = _not_implemented
