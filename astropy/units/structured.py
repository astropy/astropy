# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines structured units and quantities.
"""


# Standard library
import numpy as np
from numpy.ma.core import _replace_dtype_fields  # May need our own version

from .core import Unit, UnitBase
from .quantity import Quantity


__all__ = ['StructuredUnit', 'StructuredQuantity']


class StructuredUnit(np.void):

    def __new__(cls, units, dtype=None):
        if not isinstance(dtype, cls):
            dtype = _replace_dtype_fields(dtype, 'O')
            dtype = np.dtype((cls, dtype))
        self = np.array(units, dtype)[()]
        self._recursively_check()
        return self

    def _recursively_check(self):
        for field in self.dtype.names:
            unit = self[field]
            if isinstance(unit, type(self)):
                unit._recursively_check()
            else:
                self[field] = Unit(unit)

    @property
    def _quantity_class(self):
        return StructuredQuantity

    def __getitem__(self, item):
        val = super().__getitem__(item)
        if isinstance(val, np.void) and val.dtype.fields:
            return val.view((self.__class__, val.dtype.fields))
        else:
            return val

    # define this to pass through Unit initializer
    def _get_physical_type_id(self):
        return NotImplementedError

    def _get_converter(self, unit, equivalencies=[]):
        """Helper method for to and to_value."""
        if not isinstance(unit, type(self)):
            unit = self.__class__(unit, self.dtype)

        def converter(value):
            result = np.empty(value.shape, value.dtype)
            for name, r_name, u_name in zip(self.dtype.names,
                                            result.dtype.names,
                                            unit.dtype.names):
                result[r_name] = self[name].to(unit[u_name], value[r_name],
                                               equivalencies=equivalencies)
            return result
        return converter

    def to(self, other, value, equivalencies=[]):
        return self._get_converter(other, equivalencies=equivalencies)(value)


class StructuredQuantity(Quantity):
    def __new__(cls, value, unit=None, dtype=None, *args, **kwargs):
        self = np.array(value, dtype, *args, **kwargs).view(cls)
        if not isinstance(unit, StructuredUnit):
            unit = StructuredUnit(unit, self.dtype)
        self._unit = unit
        return self

    def __getitem__(self, item):
        out = super().__getitem__(item)
        if out.dtype is not self.dtype:
            # item caused dtype change -> indexed with string-like.
            out_unit = self.unit[item]
            if not out.dtype.fields:
                out = out.view(getattr(out_unit, '_quantity_cls',
                                       Quantity))
            out._unit = out_unit

        return out

    def _set_unit(self, unit):
        if not isinstance(unit, StructuredUnit):
            unit = StructuredUnit(unit)
        self._unit = unit

    def to(self, unit, equivalencies=[]):
        if not isinstance(unit, StructuredUnit):
            unit = StructuredUnit(unit, self.dtype)
        result = self._to_value(unit, equivalencies=equivalencies)
        result = result.view(self.__class__)
        result._set_unit(unit)
        return result

    def to_value(self, unit=None, equivalencies=[]):
        if unit is None or unit is self.unit:
            result = self.view(np.ndarray)
        else:
            result = self._to_value(unit, equivalencies)
        # [()] to return a void rather than a tuple.
        return result if result.shape else result[()]

    value = property(to_value,
                     doc="""The numerical value of this instance.

    See also
    --------
    to_value : Get the numerical value in a given unit.
    """)
