# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from ..units import Quantity
from .uncertainty import Variable

__all__ = ["UQuantity"]


class UQuantity(Variable, Quantity):

    def __new__(cls, value, uncertainty=None, unit=None, **kwargs):
        if uncertainty is None:
            uncertainty = getattr(value, 'uncertainty', None)
        value = getattr(value, 'nominal_value', value)
        subok = kwargs.pop('subok', True)
        if subok and isinstance(value, Quantity):
            q_cls = type(value)
            cls = _subclasses[q_cls]
        else:
            q_cls = Quantity
        value = q_cls(value, unit, **kwargs)
        if uncertainty is not None:
            uncertainty = q_cls(uncertainty, value.unit)
        return super().__new__(cls, value, uncertainty, **kwargs)

    @property
    def uncertainty(self):
        return (0 if self._uncertainty is None
                else self._uncertainty().to(self.unit).value)

    def __quantity_subclass__(self, unit):
        q_cls = self._nominal_value.__quantity_subclass__(unit)[0]
        return _subclasses[q_cls], True

    def __str__(self):
        return '{0}±{1}{2:s}'.format(self.value, self.uncertainty,
                                     self._unitstr)

    def __repr__(self):
        prefix1 = '<' + self.__class__.__name__ + ' '
        if self.isscalar:
            return '{0}{1}±{2}{3:s}>'.format(prefix1, self.value,
                                             self.uncertainty, self._unitstr)

        prefix2 = ' ' * (len(prefix1) - 1) + '±'
        valstr = np.array2string(self.value, separator=',', prefix=prefix1)
        errstr = np.array2string(self.uncertainty, separator=',',
                                 prefix=prefix2)
        return '{0}{1}\n{2}{3}{4:s}>'.format(prefix1, valstr,
                                             prefix2, errstr, self._unitstr)


class _SubclassDict(dict):
    def __getitem__(self, q_cls):
        if q_cls not in self:
            # Use str('U') to avoid unicode for class name in python2.
            self[q_cls] = type(str('U') + q_cls.__name__,
                               (UQuantity, Variable, q_cls), {})
        return super().__getitem__(q_cls)


_subclasses = _SubclassDict()
_subclasses[Quantity] = UQuantity
