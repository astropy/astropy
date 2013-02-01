# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines the `Quantity` object, which represents a number with some associated
units. `Quantity` objects support operations like ordinary numbers, but will deal with
unit conversions internally.
"""

from __future__ import absolute_import, unicode_literals, division, print_function

# Standard library
import copy
import numbers

import numpy as np

# AstroPy
from .core import Unit, UnitBase, CompositeUnit
from .. import log
from ..config import ConfigurationItem

WARN_IMPLICIT_NUMERIC_CONVERSION = ConfigurationItem(
                                          "warn_implicit_numeric_conversion",
                                          True,
                                          "Whether to show a warning message "
                                          "in the log when converting a "
                                          "Quantity to a float/int")

__all__ = ["Quantity"]


def _is_unity(value):
    return not value.bases and value.scale == 1.


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
        raise TypeError("The value must be a valid Python or Numpy numeric type.")

    return value_obj


class Quantity(object):
    """ A `Quantity` represents a number with some associated unit.

    Parameters
    ----------
    value : number
        The numerical value of this quantity in the units given by unit.
    unit : `~astropy.units.UnitBase` instance, str
        An object that represents the unit associated with the input value. Must be an `~astropy.units.UnitBase`
        object or a string parseable by the `units` package.

    Raises
    ------
    TypeError
        If the value provided is not a Python numeric type.
    TypeError
        If the unit provided is not either a `Unit` object or a parseable string unit.
    """

    def __init__(self, value, unit):
        self._value = _validate_value(value)
        self._unit = Unit(unit)

    def to(self, unit, equivalencies=[]):
        """ Returns a new `Quantity` object with the specified units.

        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            An object that represents the unit to convert to. Must be
            an `~astropy.units.UnitBase` object or a string parseable
            by the `units` package.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
        """
        new_val = self.unit.to(unit, self.value, equivalencies=equivalencies)
        new_unit = Unit(unit)
        return Quantity(new_val, new_unit)

    @property
    def value(self):
        """ The numerical value of this quantity. """
        return self._value

    @value.setter
    def value(self, obj):
        """ Setter for the value attribute. We allow the user to change the value by setting this attribute,
        so this will validate the new object.

        Parameters
        ----------
        obj : number
            The numerical value of this quantity in the same units as stored internally.

        Raises
        ------
        TypeError
            If the value provided is not a Python numeric type.
        """
        self._value = _validate_value(obj)

    @property
    def unit(self):
        """ A `~astropy.units.UnitBase` object representing the unit of this quantity. """
        return self._unit

    @property
    def si(self):
        """ Returns a copy of the current `Quantity` instance with SI
            units. The value of the resulting object will be scaled.
        """
        from . import si
        si_unit = self.unit.to_system(si)[0]
        return Quantity(self.value * si_unit.scale, si_unit / si_unit.scale)

    @property
    def cgs(self):
        """ Returns a copy of the current `Quantity` instance with CGS
            units. The value of the resulting object will be scaled.
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
        return Quantity(self.value, unit=self.unit)

    # Arithmetic operations
    def __add__(self, other):
        """ Addition between `Quantity` objects and other objects.  If
        they are both `Quantity` objects, results in the units of the
        **left** object if they are compatible, otherwise this fails.
        """
        if isinstance(other, Quantity):
            return Quantity(self.value + other.to(self.unit).value, unit=self.unit)
        else:
            raise TypeError("Object of type '{0}' cannot be added with a "
                "Quantity object. Addition is only supported between Quantity "
                "objects with compatible units.".format(other.__class__))

    def __sub__(self, other):
        """ Subtraction between `Quantity` objects and other objects.
        If they are both `Quantity` objects, results in the units of the
        **left** object if they are compatible, otherwise this fails.
        """
        if isinstance(other, Quantity):
            return Quantity(self.value - other.to(self.unit).value, unit=self.unit)
        else:
            raise TypeError("Object of type '{0}' cannot be added with a "
                "Quantity object. Addition is only supported between Quantity "
                "objects with compatible units.".format(other.__class__))

    def __mul__(self, other):
        """ Multiplication between `Quantity` objects and other objects.
        """
        if isinstance(other, Quantity):
            return Quantity(self.value * other.value, unit=self.unit * other.unit)
        elif isinstance(other, UnitBase):
            return Quantity(self.value, unit=other * self.unit)
        else:
            try:
                return Quantity(other * self.value, unit=self.unit)
            except TypeError:
                raise TypeError("Object of type '{0}' cannot be multiplied with a Quantity object.".format(other.__class__))

    def __rmul__(self, other):
        """ Right Multiplication between `Quantity` objects and other
        objects.
        """
        return self.__mul__(other)

    def __div__(self, other):
        """ Division between `Quantity` objects and other objects.
        """
        if isinstance(other, Quantity):
            return Quantity(self.value / other.value, unit=self.unit / other.unit)
        elif isinstance(other, UnitBase):
            return Quantity(self.value, unit=self.unit / other)
        else:
            try:
                return Quantity(self.value / other, unit=self.unit)
            except TypeError:
                raise TypeError("Object of type '{0}' cannot be diveded with a Quantity object.".format(other.__class__))

    def __rdiv__(self, other):
        """ Right Division between `Quantity` objects and other objects.
        """
        if isinstance(other, Quantity):
            return Quantity(other.value / self.value, unit=other.unit / self.unit)
        elif isinstance(other, UnitBase):
            return Quantity(1. / self.value, unit=other / self.unit)
        else:
            try:
                return Quantity(other / self.value, unit=1. / self.unit)
            except TypeError:
                raise TypeError("Object of type '{0}' cannot be diveded with a Quantity object.".format(other.__class__))

    def __truediv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__div__(other)

    def __rtruediv__(self, other):
        """ Division between `Quantity` objects. """
        return self.__rdiv__(other)

    def __pow__(self, p):
        """ Raise `Quantity` object to a power. """
        if hasattr(p, 'unit'):
            raise TypeError('Cannot raise a Quantity object to a power of something with a unit')
        return Quantity(self.value ** p, unit=self.unit ** p)

    def __neg__(self):
        """
        Minus the quantity. This is useful for doing -q where q is a quantity.
        """
        return Quantity(-self.value, unit=self.unit)

    def __pos__(self):
        """
        Plus the quantity. This is implemented in case users use +q where q is a quantity.
        """
        return Quantity(self.value, unit=self.unit)

    def __abs__(self):
        """
        Absolute value of the quantity.
        """
        return Quantity(abs(self.value), unit=self.unit)

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
            raise TypeError("Quantity object cannot be compared to an object of type {0}".format(other.__class__))

    def __le__(self, other):
        if isinstance(other, Quantity):
            return self.value <= other.to(self.unit).value
        else:
            raise TypeError("Quantity object cannot be compared to an object of type {0}".format(other.__class__))

    def __gt__(self, other):
        if isinstance(other, Quantity):
            return self.value > other.to(self.unit).value
        else:
            raise TypeError("Quantity object cannot be compared to an object of type {0}".format(other.__class__))

    def __ge__(self, other):
        if isinstance(other, Quantity):
            return self.value >= other.to(self.unit).value
        else:
            raise TypeError("Quantity object cannot be compared to an object of type {0}".format(other.__class__))

    #other overrides of special functions
    def __hash__(self):
        return hash(self.value) ^ hash(self.unit)

    def __getitem__(self, key):
        if self.isscalar:
            raise TypeError("'{cls}' object with a scalar value does not support indexing".format(cls=self.__class__.__name__))
        else:
            return Quantity(self.value[key], unit=self.unit)

    def __len__(self):
        if self.isscalar:
            raise TypeError("'{cls}' object with a scalar value has no len()".format(cls=self.__class__.__name__))
        else:
            return len(self.value)

    # Numerical types
    def __float__(self):
        if not self.isscalar:
            raise TypeError('Only scalar quantities can be converted to Python scalars')
        # In the following, we check that unit is not the dimensionless unit
        # (note that this is different from ``is_dimensionless``, which checks
        # if the decomposed - not current - units are dimensionless.
        if not _is_unity(self.unit) and WARN_IMPLICIT_NUMERIC_CONVERSION():
            log.warn("Converting Quantity object in units '{0}' to a Python scalar".format(self.unit))
        return float(self.value)

    def __int__(self):
        if not self.isscalar:
            raise TypeError('Only scalar quantities can be converted to Python scalars')
        # In the following, we check that unit is not the dimensionless unit
        # (note that this is different from ``is_dimensionless``, which checks
        # if the decomposed - not current - units are dimensionless.
        if not _is_unity(self.unit) and WARN_IMPLICIT_NUMERIC_CONVERSION():
            log.warn("Converting Quantity object in units '{0}' to a Python scalar".format(self.unit))
        return int(self.value)

    def __long__(self):
        if not self.isscalar:
            raise TypeError('Only scalar quantities can be converted to Python scalars')
        # In the following, we check that unit is not the dimensionless unit
        # (note that this is different from ``is_dimensionless``, which checks
        # if the decomposed - not current - units are dimensionless.
        if not _is_unity(self.unit) and WARN_IMPLICIT_NUMERIC_CONVERSION():
            log.warn("Converting Quantity object in units '{0}' to a Python scalar".format(self.unit))
        return long(self.value)

    # Array types
    def __array__(self):
        # In the following, we check that unit is not the dimensionless unit
        # (note that this is different from ``is_dimensionless``, which checks
        # if the decomposed - not current - units are dimensionless.
        if not _is_unity(self.unit) and WARN_IMPLICIT_NUMERIC_CONVERSION():
            log.warn("Converting Quantity object in units '{0}' to a Numpy array".format(self.unit))
        return np.array(self.value)

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
