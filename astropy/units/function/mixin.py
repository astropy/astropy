# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ..core import IrreducibleUnit, Unit


class FunctionMixin:
    """Mixin class that makes UnitBase subclasses callable.

    Provides a __call__ method that passes on arguments to a FunctionUnit.
    Instances of this class should define ``_function_unit_class`` pointing
    to the relevant class.

    See units.py and logarithmic.py for usage.
    """
    def __call__(self, unit=None):
        return self._function_unit_class(physical_unit=unit,
                                         function_unit=self)


class IrreducibleFunctionUnit(FunctionMixin, IrreducibleUnit):
    pass


class RegularFunctionUnit(FunctionMixin, Unit):
    pass
