# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines structured units and quantities.
"""


# Standard library
import numpy as np

from .core import Unit, UnitBase
from .quantity import Quantity


__all__ = ['StructuredUnit']


class StructuredUnit(np.void):

    def __new__(cls, units, dtype):
        if not isinstance(dtype, cls):
            dtype = np.dtype((cls, dtype))
        self = np.array(units, dtype)[()]
        return self

    def __getitem__(self, item):
        val = super().__getitem__(item)
        if isinstance(val, np.void) and val.dtype.fields:
            return val.view((self.__class__, val.dtype.fields))
        else:
            return val

    # define this to pass through Unit initializer
    def _get_physical_type_id(self):
        return NotImplementedError
