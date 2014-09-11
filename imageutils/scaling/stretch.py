"""
Classes that deal with stretching, i.e. mapping a range of [0:1] values onto
another set of [0:1] values with a transformation
"""

import abc

import numpy as np

from astropy.extern import six

__all__ = ["LinearStretch", "SqrtStretch", "PowerStretch", "HistEqStretch", "ContrastBiasStretch"]


@six.add_metaclass(abc.ABCMeta)
class BaseStretch(object):

    def __add__(self, other):
        return CompositeStretch(other, self)


class CompositeStretch(BaseStretch):
    """
    A combiantion of two stretches.

    Parameters
    ----------
    stretch_1:
        The first stretch to apply.
    stretch_2:
        The second stretch to apply.
    """

    def __init__(self, stretch_1, stretch_2):
        super(CompositeStretch, self).__init__()
        self.stretch_1 = stretch_1
        self.stretch_2 = stretch_2

    def __call__(self, values):
        return self.stretch_2(self.stretch_1(values))

    def inverted(self):
        return CompositeStretch(self.stretch_2.inverted(),
                                self.stretch_1.inverted())


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
        super(PowerStretch, self).__init__()
        self.power = a

    def __call__(self, values):
        return values ** self.power

    def inverted(self):
        return PowerStretch(1. / self.power)


class HistEqStretch(BaseStretch):
    """
    A histogram equalization stretch

    Parameters
    ----------
    data : float
        The data defining the equalization
    """

    def __init__(self, data, values=None):

        # Assume data is not necessarily normalized at this point
        self.data = np.sort(data.ravel())
        vmin = self.data.min()
        vmax = self.data.max()
        self.data = (self.data - vmin) / (vmax - vmin)

        # Compute relative position of each pixel
        if values is None:
            self.values = np.linspace(0., 1., len(self.data))
        else:
            self.values = values

    def __call__(self, values):
        return np.interp(values, self.data, self.values)

    def inverted(self):
        return InvertedHistEqStretch(self.data, values=self.values)


class InvertedHistEqStretch(BaseStretch):

    def __init__(self, data, values=None):
        self.data = data
        if values is None:
            self.values = np.linspace(0., 1., len(self.data))
        else:
            self.values = values

    def __call__(self, values):
        return np.interp(values, self.values, self.data)

    def inverted(self):
        return HistEqStretch(self.data, values=self.values)


class ContrastBiasStretch(BaseStretch):
    """
    A stretch that takes into account contrast and bias: y = clip((x - bias) * contrast + 0.5, 0, 1)
    """

    def __init__(self, contrast, bias):
        super(ContrastBiasStretch, self).__init__()
        self.contrast = contrast
        self.bias = bias

    def __call__(self, values):
        return np.clip((values - self.bias) * self.contrast + 0.5, 0., 1.)

    def inverted(self):
        return InvertedContrastBiasStretch(self.contrast, self.bias)


class InvertedContrastBiasStretch(BaseStretch):
    """
    A stretch that takes into account contrast and bias: y = clip((x - bias) * contrast + 0.5, 0, 1)
    """

    def __init__(self, contrast, bias):
        super(InvertedContrastBiasStretch, self).__init__()
        self.contrast = contrast
        self.bias = bias

    def __call__(self, values):
        return np.clip((values - 0.5) / self.contrast + self.bias, 0., 1.)

    def inverted(self):
        return ContrastBiasStretch(self.contrast, self.bias)
