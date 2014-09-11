"""
Classes that deal with stretching, i.e. mapping a range of [0:1] values onto
another set of [0:1] values with a transformation
"""

import abc

import numpy as np

from astropy.extern import six

__all__ = ["LinearStretch", "SqrtStretch", "PowerStretch"]


@six.add_metaclass(abc.ABCMeta)
class BaseStretch(object):
    pass


class LinearStretch(BaseStretch):
    """
    A linear stretch: y = x
    """

    def __call__(self, values):
        return values

    def inverted(self):
        return LinearStretch()


class SqrtStretch(BaseStretch):
    """
    A square root stretch: y = sqrt(x)
    """

    def __call__(self, values):
        return np.sqrt(values)

    def inverted(self):
        return PowerStretch(2)


class PowerStretch(BaseStretch):
    """
    A power-law stretch: y = x ** a
    """

    def __init__(self, a):
        self.power = a

    def __call__(self, values):
        return values ** self.power

    def inverted(self):
        return PowerStretch(1. / self.power)
