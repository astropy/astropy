# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines the `Quantity` object, which represents a number with some
associated units. `Quantity` objects support operations like ordinary numbers,
but will deal with unit conversions internally.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import copy
import numbers
import warnings

import numpy as np

# AstroPy
from .core import Unit, UnitBase, CompositeUnit, UnitsException
from ..config import ConfigurationItem
from ..utils.compat.misc import override__dir__


WARN_IMPLICIT_NUMERIC_CONVERSION = ConfigurationItem(
                                          "warn_implicit_numeric_conversion",
                                          True,
                                          "Whether to show a warning message "
                                          "in the log when converting a "
                                          "Quantity to a float/int")

__all__ = ["Quantity"]

# Numpy ufuncs that return unitless values
DIMENSIONLESS_UFUNCS = set([np.exp, np.log, np.log1p, np.log2, np.log10])

TRIG_UFUNCS = set([np.cos, np.sin, np.tan])

INVTRIG_UFUNCS = set([np.arccos, np.arcsin, np.arctan])


def _is_unity(value):
    x = value.decompose()
    return (len(x.bases) == 0 and x.scale == 1.0)


def _validate_value(value):
    """ Make sure that the input is a Python or Numpy numeric type.

    Parameters
    ----------
    value : number
        An object that will be checked whether it is a numeric type or not.

    Returns
    -------
    newval
        The new value either as an array or a scalar
    """

    from ..utils.misc import isiterable

    if isinstance(value, (numbers.Number, np.number)):
        value_obj = value
    elif isiterable(value):
        value_obj = np.array(value, copy=True)
    elif isinstance(value, np.ndarray):
        # A length-0 numpy array (i.e. numpy scalar) which we accept as-is
        value_obj = np.array(value, copy=True)
    else:
        raise TypeError("The value must be a valid Python or Numpy numeric "
                        "type.")

    return value_obj


class Quantity(object):
    """ A `Quantity` represents a number with some associated unit.

    Parameters
    ----------
    value : number, `Quantity` object, or sequence of `Quantity` objects.
        The numerical value of this quantity in the units given by
        unit.  If a `Quantity` or sequence of them, creates a new
        `Quantity` object, converting to `unit` units as needed.

    unit : `~astropy.units.UnitBase` instance, str
        An object that represents the unit associated with the input value.
        Must be an `~astropy.units.UnitBase` object or a string parseable by
        the `units` package.

    equivalencies : list of equivalence pairs, optional
        A list of equivalence pairs. See :ref:`unit_equivalencies`.

    Raises
    ------
    TypeError
        If the value provided is not a Python numeric type.
    TypeError
        If the unit provided is not either a `Unit` object or a parseable
        string unit.
    """
    # Need to set a class-level default for _equivalencies, or
    # Constants can not initialize properly
    _equivalencies = []

    # Make sure that __rmul__ of units gets called over the __mul__ of Numpy
    # arrays to avoid element-wise multiplication.
    __array_priority__ = 1000.

    def __array_prepare__(self, array, context=None):
        return array

    def __array_wrap__(self, out, context=None):
        if context:
            # Numpy will take the output of the ufuncs and try and convert it
            # to an array, which means that it will convert it to an array of
            # Quantity objects (because __len__ is defined, Numpy thinks it can
            # do that). Kind of a waste of computation to try and piece back
            # together the Quantity array, so we do the following, which is at
            # least correct.
            from . import dimensionless_unscaled
            from .si import radian
            if context[0] is np.sqrt:
                return Quantity(self.value ** 0.5, self.unit ** 0.5)
            elif context[0] in DIMENSIONLESS_UFUNCS:
                return Quantity(context[0](self.value), dimensionless_unscaled)
            elif context[0] in TRIG_UFUNCS:
                return Quantity(context[0](self.to(radian).value), dimensionless_unscaled)
            elif context[0] in INVTRIG_UFUNCS:
                return Quantity(context[0](self.value), radian)
        return self.__class__(out, unit=self._unit)

    def __init__(self, value, unit, equivalencies=[]):
        from ..utils.misc import isiterable

        self._unit = Unit(unit)
        self._equivalencies = Unit._normalize_equivalencies(equivalencies)

        if isinstance(value, Quantity):
            self._value = _validate_value(value.to(self._unit).value)
        elif isiterable(value) and all([isinstance(v, Quantity) for v in value]):
            self._value = _validate_value([q.to(self._unit).value for q in value])
        else:
            self._value = _validate_value(value)

    def to(self, unit, equivalencies=None):
        """ Returns a new `Quantity` object with the specified units.

        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            An object that represents the unit to convert to. Must be
            an `~astropy.units.UnitBase` object or a string parseable
            by the `units` package.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.  If
            not provided, the equivalencies that were provided in the
            constructor will be used.
        """
        if equivalencies is None:
            equivalencies = self._equivalencies
        new_val = self.unit.to(unit, self.value, equivalencies=equivalencies)
        new_unit = Unit(unit)
        return Quantity(new_val, new_unit)

    @property
    def value(self):
        """ The numerical value of this quantity. """

        return self._value

    @property
    def unit(self):
        """
        A `~astropy.units.UnitBase` object representing the unit of this
        quantity.
        """

        return self._unit

    @property
    def equivalencies(self):
        """
        A list of equivalencies that will be applied implicitly during
        unit conversions.
        """

        return self._equivalencies

    @property
    def si(self):
        """
        Returns a copy of the current `Quantity` instance with SI units. The
        value of the resulting object will be scaled.
        """

        from . import si
        si_unit = self.unit.to_system(si)[0]
        return Quantity(self.value * si_unit.scale, si_unit / si_unit.scale)

    @property
    def cgs(self):
        """
        Returns a copy of the current `Quantity` instance with CGS units. The
        value of the resulting object will be scaled.
        """

        from . import cgs
        cgs_unit = self.unit.to_system(cgs)[0]
        return Quantity(self.value * cgs_unit.scale, cgs_unit / cgs_unit.scale)

    @property
    def isscalar(self):
        """
        True if the `value` of this quantity is a scalar, or False if it
        is an array-like object.

        .. note::
            This is subtly different from `numpy.isscalar` in that
            `numpy.isscalar` returns False for a zero-dimensional array
            (e.g. ``np.array(1)``), while this is True in that case.
        """

        from ..utils.misc import isiterable

        return not isiterable(self.value)

    def copy(self):
        """ Return a copy of this `Quantity` instance """

        return self.__class__(self.value, unit=self.unit)

    @override__dir__
    def __dir__(self):
        """
        Quantities are able to directly convert to other units that
        have the same physical type.  This function is implemented in
        order to make autocompletion still work correctly in IPython.
        """
        extra_members = set()
        for equivalent in self.unit._get_units_with_same_physical_type(
                self._equivalencies):
            if len(equivalent.aliases):
                name = equivalent.aliases[0]
            else:
                name = equivalent.name
            extra_members.add(name)
        return extra_members

    def __getattr__(self, attr):
        """
        Quantities are able to directly convert to other units that
        have the same physical type.
        """
        def get_virtual_unit_attribute():
            try:
                to_unit = Unit(attr)
            except ValueError:
                return None

            if len(to_unit.aliases):
                if to_unit.aliases[0] != attr:
                    return None
            else:
                if to_unit.name != attr:
                    return None

            try:
                return self.unit.to(
                    to_unit, self.value, equivalencies=self.equivalencies)
            except UnitsException:
                return None

        value = get_virtual_unit_attribute()

        if value is None:
            raise AttributeError(
                "{0} instance has no attribute '{1}'".format(
                    self.__class__.__name__, attr))
        else:
            return value

    # Arithmetic operations
    def __add__(self, other):
        """ Addition between `Quantity` objects and other objects.  If
        they are both `Quantity` objects, results in the units of the
        **left** object if they are compatible, otherwise this fails.
        """

        if isinstance(other, Quantity):
            return Quantity(self.value + other.to(self.unit).value,
                            unit=self.unit)
        else:
            try:
                return self.__add__(other * Unit(1))
            except TypeError:
                raise TypeError(
                    "Object of type '{0}' cannot be added with a Quantity "
                    "object. Addition is only supported between Quantity "
                    "objects with compatible units.".format(other.__class__))

    __radd__ = __add__

    def __sub__(self, other):
        """ Subtraction between `Quantity` objects and other objects.
        If they are both `Quantity` objects, results in the units of the
        **left** object if they are compatible, otherwise this fails.
        """

        if isinstance(other, Quantity):
            return Quantity(self.value - other.to(self.unit).value,
                            unit=self.unit)
        else:
            try:
                return self.__sub__(other * Unit(1))
            except TypeError:
                raise TypeError(
                    "Object of type '{0}' cannot be subtracted from a Quantity "
                    "object. Subtraction is only supported between Quantity "
                    "objects with compatible units.".format(other.__class__))

    def __rsub__(self, other):
        return -__sub__(self, other)

    def __mul__(self, other):
        """ Multiplication between `Quantity` objects and other objects."""

        if isinstance(other, Quantity):
            return Quantity(self.value * other.value,
                            unit=self.unit * other.unit)
        elif isinstance(other, basestring):
            return Quantity(self.value, unit=Unit(other) * self.unit)
        elif isinstance(other, UnitBase):
            return Quantity(self.value, unit=other * self.unit)
        else:
            try:
                return Quantity(other * self.value, unit=self.unit)
            except TypeError:
                raise TypeError(
                    "Object of type '{0}' cannot be multiplied with a "
                    "Quantity object.".format(other.__class__))

    def __rmul__(self, other):
        """ Right Multiplication between `Quantity` objects and other
        objects.
        """

        return self.__mul__(other)

    def __div__(self, other):
        """ Division between `Quantity` objects and other objects."""

        if isinstance(other, Quantity):
            return Quantity(self.value / other.value,
                            unit=self.unit / other.unit)
        elif isinstance(other, basestring):
            return Quantity(self.value, unit=self.unit / Unit(other))
        elif isinstance(other, UnitBase):
            return Quantity(self.value, unit=self.unit / other)
        else:
            try:
                return Quantity(self.value / other, unit=self.unit)
            except TypeError:
                raise TypeError(
                    "Object of type '{0}' cannot be diveded with a Quantity "
                    "object.".format(other.__class__))

    def __rdiv__(self, other):
        """ Right Division between `Quantity` objects and other objects."""

        if isinstance(other, Quantity):
            return Quantity(other.value / self.value,
                            unit=other.unit / self.unit)
        elif isinstance(other, basestring):
            return Quantity(self.value, unit=Unit(other) / self.unit)
        elif isinstance(other, UnitBase):
            return Quantity(1. / self.value, unit=other / self.unit)
        else:
            try:
                return Quantity(other / self.value, unit=1. / self.unit)
            except TypeError:
                raise TypeError(
                    "Object of type '{0}' cannot be diveded with a Quantity "
                    "object.".format(other.__class__))

    def __truediv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__div__(other)

    def __rtruediv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__rdiv__(other)

    def __pow__(self, p):
        """ Raise `Quantity` object to a power. """

        if hasattr(p, 'unit'):
            raise TypeError(
                'Cannot raise a Quantity object to a power of something '
                'with a unit')
        return Quantity(self.value ** p, unit=self.unit ** p)

    def __neg__(self):
        """
        Minus the quantity. This is useful for doing -q where q is a quantity.
        """

        return Quantity(-self.value, unit=self.unit)

    def __pos__(self):
        """
        Plus the quantity. This is implemented in case users use +q where q is
        a quantity.
        """

        return Quantity(self.value, unit=self.unit)

    def __abs__(self):
        """
        Absolute value of the quantity.
        """

        return Quantity(abs(self.value), unit=self.unit)

    # Statistical methods

    def mean(self, axis=None, dtype=None, out=None):
        value = np.mean(self.value, axis=axis, dtype=dtype, out=out)
        return Quantity(value, unit=self.unit)

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        value = np.std(self.value, axis=axis, dtype=dtype, out=out, ddof=ddof)
        return Quantity(value, self._unit)

    # Trigonometric methods (these get called when np.cos/np.sin/np.tan are used)

    def cos(self):
        from . import dimensionless_unscaled
        from .si import radian
        try:
            return Quantity(np.cos(self.to(radian).value), unit=dimensionless_unscaled)
        except UnitsException:
            raise TypeError("Can only apply trigonometric functions to quantities with angle units")

    def sin(self):
        from . import dimensionless_unscaled
        from .si import radian
        try:
            return Quantity(np.sin(self.to(radian).value), unit=dimensionless_unscaled)
        except UnitsException:
            raise TypeError("Can only apply trigonometric functions to quantities with angle units")

    def tan(self):
        from . import dimensionless_unscaled
        from .si import radian
        try:
            print(self.to(radian).value)
            print(np.tan(self.to(radian).value))
            return Quantity(np.tan(self.to(radian).value), unit=dimensionless_unscaled)
        except UnitsException:
            raise TypeError("Can only apply trigonometric functions to quantities with angle units")

    def arccos(self):
        from .si import radian
        if _is_unity(self.unit):
            return Quantity(np.arccos(self.value), unit=radian)
        else:
            raise TypeError("Can only apply inverse trigonometric functions to dimensionless and unscaled quantities")

    def arcsin(self):
        from .si import radian
        if _is_unity(self.unit):
            return Quantity(np.arcsin(self.value), unit=radian)
        else:
            raise TypeError("Can only apply inverse trigonometric functions to dimensionless and unscaled quantities")

    def arctan(self):
        from .si import radian
        if _is_unity(self.unit):
            return Quantity(np.arctan(self.value), unit=radian)
        else:
            raise TypeError("Can only apply inverse trigonometric functions to dimensionless and unscaled quantities")

    # Mathematical methods (these get called when np.sqrt/np.exp/etc. are used)

    def sqrt(self):
        return Quantity(self.value ** 0.5, self.unit ** 0.5)

    def exp(self):
        if _is_unity(self._unit):
            return np.exp(self.value)
        else:
            raise TypeError("Can only apply exp function to dimensionless quantities")

    def log(self):
        if _is_unity(self._unit):
            return np.log(self.value)
        else:
            raise TypeError("Can only apply log function to dimensionless quantities")

    def log2(self):
        if _is_unity(self._unit):
            return np.log2(self.value)
        else:
            raise TypeError("Can only apply log2 function to dimensionless quantities")

    def log10(self):
        if _is_unity(self._unit):
            return np.log10(self.value)
        else:
            raise TypeError("Can only apply log10 function to dimensionless quantities")

    def log1p(self):
        if _is_unity(self._unit):
            return np.log1p(self.value)
        else:
            raise TypeError("Can only apply log1p function to dimensionless quantities")

    # Comparison operations
    def __eq__(self, other):
        if hasattr(other, 'value') and hasattr(other, 'to'):
            return self.value == other.to(self.unit).value
        else:
            return False

    def __ne__(self, other):
        if hasattr(other, 'value') and hasattr(other, 'to'):
            return self.value != other.to(self.unit).value
        else:
            return True

    def __lt__(self, other):
        if isinstance(other, Quantity):
            return self.value < other.to(self.unit).value
        else:
            raise TypeError("Quantity object cannot be compared to an object "
                            "of type {0}".format(other.__class__))

    def __le__(self, other):
        if isinstance(other, Quantity):
            return self.value <= other.to(self.unit).value
        else:
            raise TypeError("Quantity object cannot be compared to an object "
                            "of type {0}".format(other.__class__))

    def __gt__(self, other):
        if isinstance(other, Quantity):
            return self.value > other.to(self.unit).value
        else:
            raise TypeError("Quantity object cannot be compared to an object "
                            "of type {0}".format(other.__class__))

    def __ge__(self, other):
        if isinstance(other, Quantity):
            return self.value >= other.to(self.unit).value
        else:
            raise TypeError("Quantity object cannot be compared to an object "
                            "of type {0}".format(other.__class__))

    #other overrides of special functions
    def __hash__(self):
        return hash(self.value) ^ hash(self.unit)

    def __iter__(self):
        if self.isscalar:
            raise TypeError(
                "'{cls}' object with a scalar value is not iterable"
                .format(cls=self.__class__.__name__))

        # Otherwise return a generator
        def quantity_iter():
            for val in self.value:
                yield Quantity(val, unit=self.unit)

        return quantity_iter()

    def __getitem__(self, key):
        if self.isscalar:
            raise TypeError(
                "'{cls}' object with a scalar value does not support "
                "indexing".format(cls=self.__class__.__name__))
        else:
            return Quantity(self.value[key], unit=self.unit)

    def __nonzero__(self):
        """Quantities should always be treated as non-False; there is too much
        potential for ambiguity otherwise.
        """

        return True

    def __len__(self):
        if self.isscalar:
            raise TypeError("'{cls}' object with a scalar value has no "
                            "len()".format(cls=self.__class__.__name__))
        else:
            return len(self.value)

    # Numerical types
    def __float__(self):
        if not self.isscalar or not _is_unity(self.unit):
            raise TypeError('Only dimensionless scalar quantities can be '
                            'converted to Python scalars')
        else:
            return float(self.value)

    def __int__(self):
        if not self.isscalar or not _is_unity(self.unit):
            raise TypeError('Only dimensionless scalar quantities can be '
                            'converted to Python scalars')
        else:
            return int(self.value)

    def __long__(self):
        if not self.isscalar or not _is_unity(self.unit):
            raise TypeError('Only dimensionless scalar quantities can be '
                            'converted to Python scalars')
        else:
            return long(self.value)

    # Display
    # TODO: we may want to add a hook for dimensionless quantities?
    def __str__(self):
        return "{0} {1:s}".format(self.value, self.unit.to_string())

    def __repr__(self):
        return "<Quantity {0} {1:s}>".format(self.value, self.unit.to_string())

    def _repr_latex_(self):
        """
        Generate latex representation of unit name.  This is used by
        the IPython notebook to show it all latexified.

        Returns
        -------
        lstr
            LaTeX string
        """

        # Format value
        latex_value = "{0:g}".format(self.value)
        if "e" in latex_value:
            latex_value = latex_value.replace('e', '\\times 10^{') + '}'

        # Format unit
        # [1:-1] strips the '$' on either side needed for math mode
        latex_unit = self.unit._repr_latex_()[1:-1]  # note this is unicode

        return u'${0} \; {1}$'.format(latex_value, latex_unit)

    def decompose(self, bases=[]):
        """
        Generates a new `Quantity` with the units
        decomposed. Decomposed units have only irreducible units in
        them (see `astropy.units.UnitBase.decompose`).

        Parameters
        ----------
        bases : sequence of UnitBase, optional
            The bases to decompose into.  When not provided,
            decomposes down to any irreducible units.  When provided,
            the decomposed result will only contain the given units.
            This will raises a `UnitsException` if it's not possible
            to do so.

        Returns
        -------
        newq : `~astropy.units.quantity.Quantity`
            A new object equal to this quantity with units decomposed.
        """
        return self._decompose(False, bases=bases)

    def _decompose(self, allowscaledunits=False, bases=[]):
        """
        Generates a new `Quantity` with the units decomposed. Decomposed
        units have only irreducible units in them (see
        `astropy.units.UnitBase.decompose`).

        Parameters
        ----------
        allowscaledunits : bool
            If True, the resulting `Quantity` may have a scale factor
            associated with it.  If False, any scaling in the unit will
            be subsumed into the value of the resulting `Quantity`

        bases : sequence of UnitBase, optional
            The bases to decompose into.  When not provided,
            decomposes down to any irreducible units.  When provided,
            the decomposed result will only contain the given units.
            This will raises a `UnitsException` if it's not possible
            to do so.

        Returns
        -------
        newq : `~astropy.units.quantity.Quantity`
            A new object equal to this quantity with units decomposed.

        """
        newu = self.unit.decompose(bases=bases)
        newval = self.value
        if not allowscaledunits and hasattr(newu, 'scale'):
            newval *= newu.scale
            newu = newu / Unit(newu.scale)

        return Quantity(newval, newu)
