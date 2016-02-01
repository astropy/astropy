# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .. import core

class FunctionMixin(object):
    """Mixin class that makes UnitBase subclasses callable.

    Provides a __call__ method that passes on arguments to a FunctionUnit.

    See units.py and logarithmic.py for usage.
    """
    _function_unit_class = None

    def __call__(self, unit=None):
        if self._function_unit_class is None:
            raise TypeError("Cannot call unit since corresponding "
                            "function unit class is not defined.")
        return self._function_unit_class(physical_unit=unit,
                                         function_unit=self)


class IrreducibleFunctionUnit(FunctionMixin, core.IrreducibleUnit):
    pass


class RegularFunctionUnit(FunctionMixin, core.Unit):
    pass
