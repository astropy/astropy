# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functional Units and Quantities."""
from __future__ import (absolute_import, unicode_literals,
                        division, print_function)

import numpy as np
# LOCAL
from ...extern import six
from .. import (Unit, UnitBase, UnitsError,
                dimensionless_unscaled, Quantity, format as unit_format,
                quantity_helper as qh)

__all__ = ['FunctionalUnitBase', 'FunctionalQuantityBase']

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
class FunctionalUnitBase(object):
    """Abstract base class for functional units

    Functional units are functions containing a physical unit, such as dB(mW)
    Most of the arithmetic operations on functional units are defined in this
    base class.

    While instantiation is defined, this class should not be used directly.
    Rather, subclasses should be used that define `_functional_unit`,
    `from_physical`, `to_physical` and `functional_quantity`.

    Parameters
    ----------
    physical_unit : `~astropy.units.Unit` or `string`
        Unit that is encapsulated within the logarithmic functional unit.
        If not given, dimensionless.

    functional_unit :  `~astropy.units.Unit` or `string`
        By default, the same as the logarithmic unit set by the subclass.

    """
    # vvvvv the following four need to be set by subclasses
    _functional_unit = NotImplementedError

    def from_physical(self, x):
        """Transformation from value in physical to value in functional units.
        Used in equivalency.  Should be overridden by subclass."""
        return NotImplementedError

    def to_physical(self, x):
        """Transformation from value in functional to value in physical units.
        Used in equivalency.  Should be overridden by subclass."""
        return NotImplementedError

    # conversion to a quantity through a method rather than by setting, e.g.,
    # _FunctionalQuantity = Magnitude, since at FunctionalUnit definition time,
    # the corresponding FunctionalQuantity is not yet known.
    def functional_quantity(self, value):
        """Convert `value` to a functional quantity `FunctionalQuantity`
        with one's unit.  This method is typically overridden by subclasses,
        such that, e.g., for `MagUnit` it would contain
        return Magnitude(value, unit=self).
        """
        return NotImplementedError
    # ^^^^^ the above four need to be set by subclasses

    # have priority over arrays, regular units, and regular quantities
    __array_priority__ = 3000

    _physical_unit = dimensionless_unscaled

    def __init__(self, physical_unit=None, functional_unit=None):
        if physical_unit is not None:
            self._physical_unit = Unit(physical_unit)

        if functional_unit is not None:
            # any functional unit should be equivalent to subclass default
            if(hasattr(self._functional_unit, 'is_equivalent') and
               not self._functional_unit.is_equivalent(
                   getattr(functional_unit, 'functional_unit',
                           functional_unit))):
                raise UnitsError(
                    "Cannot initiliaze {0} instance with functional unit "
                    "'{1}', as it is not equivalent to class unit '{2}'"
                    .format(self.__class__.__name__, functional_unit,
                            self._functional_unit))

            self._functional_unit = Unit(functional_unit)

    def _copy(self, physical_unit=None, functional_unit=None):
        """Copy oneself, possibly changing physical or functional unit"""
        if physical_unit is None:
            physical_unit = self._physical_unit
        if functional_unit is None:
            functional_unit = self._functional_unit
        return self.__class__(physical_unit, functional_unit)

    @property
    def physical_unit(self):
        return self._physical_unit

    @property
    def functional_unit(self):
        return self._functional_unit

    @property
    def equivalencies(self):
        """List of equivalencies between physical and logarithmic units.

        Uses the `from_physical` and `to_physical` methods.
        """
        return [(self._physical_unit, self,
                 self.from_physical, self.to_physical)]

    # vvvv properties/methods required to behave like a unit
    def _get_physical_type_id(self):
        """Get physical type corresponding to physical unit"""
        return self._physical_unit._get_physical_type_id()

    @property
    def physical_type(self):
        return self.physical_unit.physical_type

    def _to(self, other):
        """
        Returns the scale to the specified functional unit, raising
        UnitsError is the physical units are not equivalent.

        This is required to mimic behaviour expected for any units, e.g.,
        in `~astropy.units.core.UnitBase.apply_equivalencies`.
        """
        if not self.physical_unit._is_equivalent(
                getattr(other, 'physical_unit', other)):
            raise UnitsError("'{0!r}' is not a scaled version of '{1!r}'"
                             .format(self, other))

        return self._functional_unit._to(getattr(other, 'functional_unit',
                                                 other))

    def is_equivalent(self, other, equivalencies=[]):
        if isinstance(other, tuple):
            return any(self.is_equivalent(u, equivalencies=equivalencies)
                       for u in other)

        other_physical = getattr(other, 'physical_unit',
                                 (dimensionless_unscaled
                                  if self.functional_unit.is_equivalent(other)
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
           physical units; the build-in equivalency between the logarithmic
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
        if other is self._physical_unit:
            return self.to_physical(value)

        other_functional_unit = getattr(other, 'functional_unit', other)
        if self._functional_unit.is_equivalent(other_functional_unit):
            # when other is logarithmic:
            # first convert physical units to other's physical units
            other_physical_unit = getattr(other, 'physical_unit',
                                          dimensionless_unscaled)
            if self._physical_unit is not other_physical_unit:
                value_other_phys = self._physical_unit.to(
                    other_physical_unit, self.to_physical(value),
                    equivalencies)
                # make functional unit again, in own system
                value = self.from_physical(value_other_phys)

            # convert possible difference in functional scale (e.g., dex->dB)
            return self._functional_unit.to(other_functional_unit, value)

        else:
            # when other is not a logarithmic unit
            return self._physical_unit.to(other, self.to_physical(value),
                                          equivalencies)

    def __eq__(self, other):
        return (self._physical_unit == getattr(other, 'physical_unit',
                                               dimensionless_unscaled) and
                self._functional_unit == getattr(other, 'functional_unit',
                                                 other))

    def __ne__(self, other):
        return (self._physical_unit != getattr(other, 'physical_unit',
                                               dimensionless_unscaled) or
                self._funtional_unit != getattr(other, 'functional_unit',
                                                other))

    def __mul__(self, other):
        # anything not like a unit, try initialising as a functional quantity
        if not isinstance(other, (six.string_types, UnitBase,
                                  FunctionalUnitBase)):
            try:
                return self.functional_quantity(other)
            except UnitsError:
                return NotImplemented

        if self._physical_unit == dimensionless_unscaled:
            # if dimensionless, drop back to normal unit and retry
            return self._functional_unit * other

        raise UnitsError("Cannot multiply a functional unit "
                         "with a physical dimension with any unit")

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        # anything not like a unit, try initialising as a functional quantity
        if not isinstance(other, (six.string_types, UnitBase,
                                  FunctionalUnitBase)):
            try:
                return self.functional_quantity(1./other)
            except UnitsError:
                return NotImplemented

        if self._physical_unit == dimensionless_unscaled:
            # if dimensionless, drop back to normal unit and retry
            return self._functional_unit / other

        raise UnitsError("Cannot divide a functional unit "
                         "with a physical dimension by any unit")

    def __rtruediv__(self, other):
        if not isinstance(other, (six.string_types, UnitBase,
                                  FunctionalUnitBase)):
            return NotImplemented

        if self._physical_unit == dimensionless_unscaled:
            # if dimensionless, drop back to normal unit and retry
            return other / self._functional_unit

        raise UnitsError("Cannot divide a functional unit "
                         "with a physical dimension into any unit")

    def __pow__(self, power):
        if power == 0:
            return dimensionless_unscaled
        elif power == 1:
            return self._copy()

        if self._physical_unit == dimensionless_unscaled:
            return self._functional_unit ** power

        raise UnitsError("Cannot raise a functional unit "
                         "with a physical dimension to any power but 0 or 1")

    def __pos__(self):
        return self._copy()

    def to_string(self, format='generic'):
        """
        Output the unit in the given format as a string

        The physical unit is appended, within parentheses, to the functional
        unit, as in "dB(mW)", with both units set using the given format

        Parameters
        ----------
        format : `astropy.format.Base` instance or str
            The name of a format or a formatter object.  If not
            provided, defaults to the generic format.
        """
        f = unit_format.get_format(format)
        self_str = f.to_string(self._functional_unit)
        if self._physical_unit != dimensionless_unscaled:
            self_str += '({0})'.format(f.to_string(self._physical_unit))
        return self_str

    def __str__(self):
        """Return string representation for unit"""
        self_str = self._functional_unit.__str__()
        if self._physical_unit != dimensionless_unscaled:
            self_str += '({0})'.format(self._physical_unit.__str__())

        return self_str

    def __repr__(self):
        return "{0}('{1}'{2})".format(
            self.__class__.__name__, self._physical_unit.__str__(),
            '' if self._functional_unit is self.__class__._functional_unit
            else ", unit='{0}'".format(self._functional_unit.__str__()))

    def __hash__(self):
        return hash((self._functional_unit, self._physical_unit))


class FunctionalQuantityBase(Quantity):
    """A representation of a (scaled) function of a number with a unit

    Functional quanties are quantities whose units are functions containing a
    physical unit, such as dB(mW).  Most of the arithmetic operations on
    functional quantities are defined in this base class.

    While instantiation is also defined here, this class should not be
    instantiated by users directly.  Rather, subclasses should be made which
    define `_FunctionalUnit` to point back to the functional unit class.

    Parameters
    ----------
    value : number, `~astropy.units.Quantity`,
            `~astropy.units.logarithmic.LogQuantity`,
            or sequence of convertible items.
        The numerical value of the functional quantity. If a number or
        a `Quantity` with a functional unit, it will be converted to `unit`
        and the physical unit will be inferred from `unit`.
        If a `Quantity` with just a physical unit, it will converted to
        the functional unit, after, if necessary, converting it to the
        physical unit inferred from `unit`.

    unit : functional unit, or `~astropy.units.functional.FunctionalUnit`,
            optional
        E.g., `~astropy.units.functional.mag`, `astropy.units.functional.dB`,
        `~astropy.units.functional.MagUnit`, etc.
        For a `FunctionalUnit` instance, the physical unit will be taken from
        it; for non-`FunctionalUnit` input, it will be inferred from `value`.
        By default, `unit` is set by the subclass.

    dtype : `~numpy.dtype`, optional
        The ``dtype`` of the resulting Numpy array or scalar that will
        hold the value.  If not provided, is is determined automatically
        from the input value.

    copy : bool, optional
        If `True` (default), then the value is copied.  Otherwise, a copy
        will only be made if :func:`__array__` returns a copy, if obj is a
        nested sequence, or if a copy is needed to satisfy ``dtype``.
        (The `False` option is intended mostly for internal use, to speed
        up initialization where it is known a copy has been made already.
        Use with care.)

    """

    # default FunctionalUnit; to be overridden by subclasses
    _FunctionalUnit = NotImplementedError

    # have priority over ndarray, regular Unit & Quantity, and FunctionalUnit
    __array_priority__ = 4000

    # store unit in different private property than is done by Quantity,
    # to allow a Quantity view with just the functional unit (stored in _unit)
    _full_unit = None

    def __new__(cls, value, unit=None, dtype=None, copy=True):

        value_unit = getattr(value, 'unit', None)
        if value_unit is None:
            # if iterable, see if first item has a unit; mixed lists fail below
            try:
                value_unit = getattr(value[0], 'unit')
            except:
                pass

        if isinstance(value_unit, FunctionalUnitBase):
            if unit is None:
                unit = value_unit
                value = value.functional_value  # for initialising Quantity
            else:
                # convert value to its quantity (will convert back below)
                value_unit = value.physical_unit
                value = value.quantity
        else:
            if unit is None:
                unit = cls._FunctionalUnit._functional_unit

        if not isinstance(unit, FunctionalUnitBase):
            unit = cls._FunctionalUnit(value_unit, functional_unit=unit)

        if value_unit is not None and value_unit is not unit:
            # convert to target unit
            value = unit.from_physical(value.to(unit.physical_unit).value)

        # initialise Quantity
        self = super(FunctionalQuantityBase, cls).__new__(
            cls, value, unit.functional_unit, dtype=dtype, copy=copy)
        # change view, and set the full unit returned by self.unit
        self = self.view(cls)
        self._full_unit = unit

        return self

    # vvvv properties not found in Quantity
    @property
    def quantity(self):
        """The physical quantity corresponding the functional one."""
        return self.to(self.unit.physical_unit)

    @property
    def functional_value(self):
        """Quantity with plain functional unit, dropping the physical unit.

        Use `~astropy.units.quantity.Quantity.value` for just the value.
        """
        return self._new_view(self, self.unit.functional_unit)

    # vvvv properties overridden to point to different places
    @property
    def unit(self):
        """Functional unit of the quantity, containing the physical unit.

        Represented by a `~astropy.units.functional.FunctionalUnitBase` object.
        """
        return self._full_unit

    @property
    def equivalencies(self):
        """Equivalencies applied by default during unit conversions.

        Contains the list to convert between functional and physical unit,
        as set by the `~astropy.units.functional.FunctionalUnitBase` unit.
        """
        return self.unit.equivalencies

    # vvvv methods overridden to add additional behaviour
    def to(self, unit, equivalencies=[]):
        result = super(FunctionalQuantityBase, self).to(unit, equivalencies)
        if isinstance(unit, FunctionalUnitBase):
            result._full_unit = unit
        return result

    def __array_finalize__(self, obj):
        if isinstance(obj, FunctionalQuantityBase):
            self._full_unit = obj._full_unit

        super(FunctionalQuantityBase, self).__array_finalize__(obj)

    def __array_prepare__(self, obj, context=None):
        """Check that the ufunc can deal with a LogQuantity"""

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
            raise TypeError("Cannot use function '{0}' with functional "
                            "quantities that are not dimensionless"
                            .format(context[0].__name__))

        return super(FunctionalQuantityBase,
                     self).__array_prepare__(obj, context)

    def __array_wrap__(self, obj, context=None):
        """Reorder units properly after ufunc calculation"""

        obj = super(FunctionalQuantityBase, self).__array_wrap__(obj, context)

        if isinstance(obj, FunctionalQuantityBase):
            obj._full_unit = obj._unit
            obj._unit = obj._full_unit.functional_unit
            # has the functional unit changed? (e.g., "mag" -> "2 mag")
            if obj._unit != obj.__class__._functional_unit:
                if obj._full_unit.physical_type is not dimensionless_unscaled:
                    raise UnitsError(
                        "Unexpected production in "
                        "FunctionalQuantityBase.__array_wrap__ "
                        "of '{0}' instance with a physical dimension yet "
                        "functional_unit '{1}' when '{2}' was expected. "
                        "Please alert the astropy developers"
                        .format(obj.__class__.__name__, obj._unit,
                                obj.__class__._functional_unit))
                if not obj._unit.is_equivalent(obj.__class__._functional_unit):
                    # for dimensionless, fall back to regular quantity
                    obj = self._new_view(obj, obj._unit)

        return obj

    def __quantity_subclass__(self, unit):
        if isinstance(unit, FunctionalUnitBase):
            return self.__class__, True
        else:
            return super(FunctionalQuantityBase,
                         self).__quantity_subclass__(unit)[0], False

    # vvvv methods overridden to change behaviour
    def __mul__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self.functional_value * other

        raise UnitsError("Cannot multiply functional quantities with "
                         "anything unless they are dimensionless")

    def __div__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self.functional_value / other

        raise UnitsError("Cannot divide functional quantities which "
                         "are not dimensionless by anything")

    def __rdiv__(self, other):
        if self.unit.physical_unit == dimensionless_unscaled:
            return self.functional_value.__rdiv__(other)

        raise UnitsError("Cannot divide functional quantities which "
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

    # vvvv override Quantity methods that make no sense for functional units
    def _not_implemented(self, *args, **kwargs):
        raise NotImplementedError

    dot = nansum = sum = cumsum = prod = cumprod = _not_implemented
