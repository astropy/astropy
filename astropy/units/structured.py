# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines structured units and quantities.
"""


# Standard library
import numpy as np

from .core import Unit, UnitBase
from .quantity import Quantity


__all__ = ['StructuredUnit', 'StructuredQuantity']


class StructuredUnit(np.void):

    def __new__(cls, units, dtype):
        if not isinstance(dtype, cls):
            dtype = np.dtype((cls, dtype))
        self = np.array(units, dtype)[()]
        return self

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
        if not isinstance(unit, StructuredUnit):
            unit = StructuredUnit(unit, dtype)
        self = np.array(value, dtype, *args, **kwargs).view(cls)
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
