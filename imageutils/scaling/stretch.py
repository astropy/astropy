"""
Classes that deal with stretching, i.e. mapping a range of [0:1] values onto
another set of [0:1] values with a transformation
"""

import abc

import numpy as np
from numpy.lib.scimath import logn

from astropy.extern import six

from .transform import BaseTransform

__all__ = ["LinearStretch", "SqrtStretch", "PowerStretch", "PowerDistStretch",
           "SquaredStretch", "LogStretch", "AsinhStretch", "SinhStretch",
           "HistEqStretch", "ContrastBiasStretch"]


@six.add_metaclass(abc.ABCMeta)
class BaseStretch(BaseTransform):
    pass


class LinearStretch(BaseStretch):
    """
    A linear stretch: y = x.
    """

    def __call__(self, values, out=None):
        return values

    def inverted(self):
        return LinearStretch()


class SqrtStretch(BaseStretch):
    """
    A square root stretch: y = sqrt(x).
    """

    def __call__(self, values, out=None):
        if out is None:
            return np.sqrt(values)
        else:
            return np.sqrt(values, out=out)

    def inverted(self):
        return PowerStretch(2)


class PowerStretch(BaseStretch):
    """
    A power-law stretch: y = x ** a.
    """

    def __init__(self, a):
        super(PowerStretch, self).__init__()
        self.power = a

    def __call__(self, values, out=None):
        if out is None:
            return np.power(values, self.power)
        else:
            return np.power(values, self.power, out=out)

    def inverted(self):
        return PowerStretch(1. / self.power)


class PowerDistStretch(BaseStretch):
    """
    An alternative power stretch: y = ((a ** x) - 1) / a.

    Notes
    -----

    This is the same as ds9's POW stretch, described at
    http://ds9.si.edu/doc/ref/how.html
    """

    def __init__(self, a=1000.0):
        super(PowerDistStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None):
        res = (self.exp ** values - 1.0) / self.exp
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return InvertedPowerDistStretch(a=self.exp)


class InvertedPowerDistStretch(BaseStretch):
    """
    Inverse transformation for PowerDistStretch.
    """

    def __init__(self, a=1000.0):
        super(InvertedPowerDistStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None):
        res = logn(self.exp, self.exp * values + 1.0)
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return PowerDistStretch(a=self.exp)


class SquaredStretch(PowerStretch):
    """
    A convenience class for a power stretch of 2.

    Notes
    -----

    This is the same as ds9's SQUARE stretch, described at
    http://ds9.si.edu/doc/ref/how.html
    """

    def __init__(self):
        super(SquaredStretch, self).__init__(2)

    def inverted(self):
        return SqrtStretch()


class LogStretch(BaseStretch):
    """
    A log stretch: y = log(a*x + 1) / log(a).

    Notes
    -----

    This is the same as ds9's LOG stretch, described here
        (http://ds9.si.edu/doc/ref/how.html)
    """

    def __init__(self, a=1000.0):
        super(LogStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None):
        res = np.log(self.exp * values + 1.0) / np.log(self.exp)
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return InvertedLogStretch(self.exp)


class InvertedLogStretch(BaseStretch):
    """
    Inverse transformation for LogStretch.
    """

    def __init__(self, a):
        super(InvertedLogStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None):
        res = (np.exp(values * np.log(self.exp)) - 1) / self.exp
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return LogStretch(self.exp)


class AsinhStretch(BaseStretch):
    """
    An asinh stretch: y = asinh(nonlinearity * x) / factor.

    Notes
    -----

    This is the same as ds9's ASINH stretch, described at
    http://ds9.si.edu/doc/ref/how.html
    """

    def __init__(self, factor=10, nonlinearity=3.0):
        super(AsinhStretch, self).__init__()
        self.factor = factor
        self.nonlinearity = nonlinearity

    def __call__(self, values, out=None):
        res = np.arcsinh(self.factor * values) / self.nonlinearity
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return InvertedAsinhStretch(factor=self.factor,
                                    nonlinearity=self.nonlinearity)


class InvertedAsinhStretch(BaseStretch):
    """
    Inverse transformation for AsinhStretch
    """

    def __init__(self, factor=10, nonlinearity=3.0):
        super(InvertedAsinhStretch, self).__init__()
        self.factor = factor
        self.nonlinearity = nonlinearity

    def __call__(self, values, out=None):
        res = np.sinh(self.nonlinearity * values) / self.factor
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return AsinhStretch(factor=self.factor,
                            nonlinearity=self.nonlinearity)


class SinhStretch(BaseStretch):
    """
    A sinh stretch: y = sinh(factor * x) / nonlinearity.
    
    Notes
    -----
    
    This is the same as ds9's SINH stretch, described at
    http://ds9.si.edu/doc/ref/how.html
    """

    def __init__(self, factor=3.0, nonlinearity=10.0):
        super(SinhStretch, self).__init__()
        self.factor = factor
        self.nonlinearity = nonlinearity

    def __call__(self, values, out=None):
        res = np.sinh(self.factor * values) / self.nonlinearity
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return InvertedSinhStretch(factor=self.factor,
                                   nonlinearity=self.nonlinearity)


class InvertedSinhStretch(BaseStretch):
    """
    Inverse transformation for SinhStretch.
    """

    def __init__(self, factor=3.0, nonlinearity=10.0):
        super(InvertedSinhStretch, self).__init__()
        self.factor = factor
        self.nonlinearity = nonlinearity

    def __call__(self, values, out=None):
        res = np.arcsinh(self.nonlinearity * values) / self.factor
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return SinhStretch(factor=self.factor,
                           nonlinearity=self.nonlinearity)


class HistEqStretch(BaseStretch):
    """
    A histogram equalization stretch.

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

    def __call__(self, values, out=None):
        if out is None:
            return np.interp(values, self.data, self.values)
        else:
            values[:] = np.interp(values, self.data, self.values)
            return values

    def inverted(self):
        return InvertedHistEqStretch(self.data, values=self.values)


class InvertedHistEqStretch(BaseStretch):
    """
    Inverse transformation for HistEqStretch
    """

    def __init__(self, data, values=None):
        self.data = data
        if values is None:
            self.values = np.linspace(0., 1., len(self.data))
        else:
            self.values = values

    def __call__(self, values, out=None):
        if out is None:
            return np.interp(values, self.values, self.data)
        else:
            values[:] = np.interp(values, self.values, self.data)
            return values

    def inverted(self):
        return HistEqStretch(self.data, values=self.values)


class ContrastBiasStretch(BaseStretch):
    """
    A stretch that takes into account contrast and bias.

    y = clip((x - bias) * contrast + 0.5, 0, 1).
    """

    def __init__(self, contrast, bias):
        super(ContrastBiasStretch, self).__init__()
        self.contrast = contrast
        self.bias = bias

    def __call__(self, values, out=None):

        if out is None:
            values = np.subtract(values, self.bias)
        else:
            np.subtract(values, self.bias, out=values)

        # Use in-place operations for the rest since they are faster
        np.multiply(values, self.contrast, out=values)
        np.add(values, 0.5, out=values)
        np.clip(values, 0, 1, out=values)

        return values

    def inverted(self):
        return InvertedContrastBiasStretch(self.contrast, self.bias)


class InvertedContrastBiasStretch(BaseStretch):
    """
    Inverse transformation for ContrastBiasStretch.
    """

    def __init__(self, contrast, bias):
        super(InvertedContrastBiasStretch, self).__init__()
        self.contrast = contrast
        self.bias = bias

    def __call__(self, values, out=None):

        if out is None:
            values = np.subtract(values, 0.5)
        else:
            np.subtract(values, 0.5, out=values)

        # Use in-place operations for the rest since they are faster
        np.divide(values, self.contrast, out=values)
        np.add(values, self.bias, out=values)
        np.clip(values, 0, 1, out=values)

        return values

    def inverted(self):
        return ContrastBiasStretch(self.contrast, self.bias)
