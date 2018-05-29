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

    def _to_value(self, unit, equivalencies=[]):
        """Helper method for to and to_value."""
        if equivalencies == []:
            equivalencies = self._equivalencies
        result = np.empty(self.shape, self.dtype)
        for name, u_name in zip(self.dtype.names, unit.dtype.names):
            result[name] = self[name]._to_value(
                unit[u_name], equivalencies=equivalencies)
        return result

    def to_value(self, unit, equivalencies=[]):
        if not isinstance(unit, StructuredUnit):
            unit = StructuredUnit(unit, self.dtype)
        return self._to_value(self, unit, equivalencies)

    def to(self, unit, equivalencies=[]):
        if not isinstance(unit, StructuredUnit):
            unit = StructuredUnit(unit, self.dtype)
        result = self._to_value(unit, equivalencies=equivalencies)
        result = result.view(self.__class__)
        result._set_unit(unit)
        return result
