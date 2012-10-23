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

# AstroPy
from .core import UnitBase, Unit

__all__ = ["Quantity", "IncompatibleUnitsError"]

def _validate_units(unit):
    """ Make sure that the input is a either a string parseable by the `units` package, or a
        `Unit` object.
        
    Parameters
    ----------
    unit : `~astropy.units.UnitBase` instance, str
        Must be an `~astropy.units.UnitBase` object or a string parseable by the `units` package.
    """
    
    if isinstance(unit, (str, unicode)):
        unit_obj = Unit(unit)
    elif isinstance(unit, UnitBase):
        unit_obj = unit
    else:
        raise ValueError("The unit must be a Python string that is parseable by the Units package, or a Unit object.")
    
    return unit_obj

class Quantity(object):
    """ A `Quantity` represents a number wth some associated unit. 
        
    Parameters
    ----------
    value : number
        Any Python numeric type.
    unit : `~astropy.units.UnitBase` instance, str
        An object that represents the unit associated with the input value. Must be an `~astropy.units.UnitBase`
        object or a string parseable by the `units` package.
    
    Raises
    ------
    ValueError
        If the unit provided is not either a `Unit` object or a parseable string unit.
    """
    
    def __init__(self, value, unit):
        # The user must pass in a Python numeric type
        if isinstance(value, numbers.Number):
            self.value = value
        else:
            raise ValueError("The value must be a valid Python numeric type.")
        
        self._unit = _validate_units(unit)
        
    def to(self, unit):
        """ Returns a new `Quantity` object with the specified units.
        
        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            An object that represents the unit to convert to. Must be an `~astropy.units.UnitBase`
            object or a string parseable by the `units` package.
        """
        new_quantity = self.copy()
        new_quantity.unit = unit
        return new_quantity
    
    @property
    def unit(self):
        return self._unit
    
    @unit.setter
    def unit(self, unit):
        """ Setter for the unit attribute. We allow the user to change units by setting this attribute,
        so this will validate the unit and internally change the value to the new unit.
        
        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            An object that represents the unit to internally convert to. Must be an `~astropy.units.UnitBase`
            object or a string parseable by the `units` package.
        
        Raises
        ------
        IncompatibleUnitsError
            If the unit to convert to is not 'equivalent' (see `Unit` documentation) to the original unit.
        """
        new_unit = _validate_units(unit)
        
        if not self.unit.is_equivalent(new_unit):
            raise IncompatibleUnitsError("This object has units of '{0:s}' and can not be converted to '{1:s}'.".format(self.unit, new_unit))
        
        self.value *= self.unit.to(new_unit)
        self._unit = new_unit
    
    def copy(self):
        """ Return a copy of this `Quantity` instance """
        return Quantity(self.value, unit=self.unit)
    
    # Arithmetic operations
    def __add__(self, other):
        """ Addition between `Quantity` objects. All operations return a new `Quantity` object
        with the units of the **left** object.
        """
        return Quantity(self.value + other.to(self.unit).value, unit=self.unit)
    
    def __radd__(self, other):
        """ Addition between `Quantity` objects. All operations return a new `Quantity` object
        with the units of the **left** object.
        """
        return Quantity(self.to(other.unit).value + other.value, unit=other.unit)
    
    def __sub__(self, other):
        """ Subtraction between `Quantity` objects. All operations return a new `Quantity` object
        with the units of the **left** object.
        """
        return Quantity(self.value - other.to(self.unit).value, unit=self.unit)
    
    def __rsub__(self, other):
        """ Subtraction between `Quantity` objects. All operations return a new `Quantity` object
        with the units of the **left** object.
        """
        return Quantity(self.to(other.unit).value - other.value, unit=other.unit)
        
    def __mul__(self, other):
        """ Multiplication between `Quantity` objects. All operations return a new `Quantity` object
        with the units of the **left** object.
        """
        if self.unit.is_equivalent(other.unit):
            return Quantity(self.value * other.to(self.unit).value, unit=self.unit*self.unit)
        else:
            return Quantity(self.value * other.value, unit=self.unit*other.unit)
    
    def __rmul__(self, other):
        """ Multiplication between `Quantity` objects. All operations return a new `Quantity` object
        with the units of the **left** object.
        """
        if self.unit.is_equivalent(other.unit):
            return Quantity(other.value * self.to(other.unit).value, unit=other.unit*other.unit)
        else:
            return Quantity(self.value * other.value, unit=other.unit * self.unit)
    
    # TODO: We should make the distinction between the __future__ division and old-style division
    def __div__(self, other):
        """ Division between `Quantity` objects. This operation returns a dimensionless object. """
        if self.unit.is_equivalent(other.unit):
            return Quantity(self.value / other.to(self.unit).value, unit=Unit(""))
        else:
            return Quantity(self.value / other.value, unit=self.unit/other.unit)
    
    def __rdiv__(self, other):
        """ Division between `Quantity` objects. This operation returns a dimensionless object. """
        if self.unit.is_equivalent(other.unit):
            return Quantity(other.value / self.to(other.unit).value, unit=Unit(""))
        else:
            return Quantity(other.value / self.value, unit=other.unit/self.unit)
    
    def __truediv__(self, other):
        """ Division between `Quantity` objects. This operation returns a dimensionless object. """
        if self.unit.is_equivalent(other.unit):
            return Quantity(self.value / other.to(self.unit).value, unit=Unit(""))
        else:
            return Quantity(self.value / other.value, unit=self.unit/other.unit)
    
    def __rtruediv__(self, other):
        """ Division between `Quantity` objects. This operation returns a dimensionless object. """
        if self.unit.is_equivalent(other.unit):
            return Quantity(other.value / self.to(other.unit).value, unit=Unit(""))
        else:
            return Quantity(other.value / self.value, unit=other.unit/self.unit)
    
    def __pow__(self, p):
        """ Raise quantity object to a power. """
        return Quantity(self.value**2, unit=(self.unit*self.unit).simplify())
    
    # Comparison operations
    def __eq__(self, other):
        return self.value == other.to(self.unit).value

    def __ne__(self, other):
        return self.value != other.to(self.unit).value
        
    def __lt__(self, other):
        return self.value < other.to(self.unit).value
    
    def __le__(self, other):
        return self.value <= other.to(self.unit).value
    
    def __gt__(self, other):
        return self.value > other.to(self.unit).value
    
    def __ge__(self, other):
        return self.value >= other.to(self.unit).value
    
    # Display
    def __str__(self):
        return "{0:g} {1:s}".format(self.value, self.unit.to_string())
    
    def __repr__(self):
        return "<Quantity value: {0:g} unit: {1:s}>".format(self.value, self.unit.to_string())
        
class IncompatibleUnitsError(Exception):
    pass