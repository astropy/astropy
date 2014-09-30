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
    A linear stretch.
    
    The stretch is given by:
    
    .. math::
        y = x
    """

    def __call__(self, values, out=None):
        return values

    def inverted(self):
        return LinearStretch()


class SqrtStretch(BaseStretch):
    r"""
    A square root stretch.
    
    The stretch is given by:
    
    .. math::
        y = \sqrt{x}
    """

    def __call__(self, values, out=None):
        if out is None:
            return np.sqrt(values)
        else:
            return np.sqrt(values, out=out)

    def inverted(self):
        return PowerStretch(2)


class PowerStretch(BaseStretch):
    r"""
    A power-law stretch.
    
    The stretch is given by:
    
    .. math::
        y = x^a
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
    r"""
    An alternative power stretch.
    
    The stretch is given by:
    
    .. math::
        y = \frac{a^x - 1}{a - 1}
    """

    def __init__(self, a=1000.0):
        if a == 1:  # singularity
            raise ValueError("a cannot be set to 1")
        super(PowerDistStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None):
        res = (self.exp ** values - 1.0) / (self.exp - 1.0)
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return InvertedPowerDistStretch(a=self.exp)


class InvertedPowerDistStretch(BaseStretch):
    """
    Inverse transformation for `~imageutils.scaling.PowerDistStretch`.
    """

    def __init__(self, a=1000.0):
        if a == 1:  # singularity
            raise ValueError("a cannot be set to 1")
        super(InvertedPowerDistStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None):
        res = logn(self.exp, (self.exp - 1.0) * values + 1.0)
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return PowerDistStretch(a=self.exp)


class SquaredStretch(PowerStretch):
    r"""
    A convenience class for a power stretch of 2.
    
    The stretch is given by:
    
    .. math::
        y = x^2
    """

    def __init__(self):
        super(SquaredStretch, self).__init__(2)

    def inverted(self):
        return SqrtStretch()


class LogStretch(BaseStretch):
    r"""
    A log stretch.
    
    The stretch is given by:
    
    .. math::
        y = \frac{\log{(a x + 1)}}{\log{(a + 1)}}.
    """

    def __init__(self, a=1000.0):
        super(LogStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None):
        res = np.log(self.exp * values + 1.0) / np.log(self.exp + 1)
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return InvertedLogStretch(self.exp)


class InvertedLogStretch(BaseStretch):
    """
    Inverse transformation for `~imageutils.scaling.LogStretch`.
    """

    def __init__(self, a):
        super(InvertedLogStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None):
        res = (np.exp(values * np.log(self.exp + 1.)) - 1) / self.exp
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return LogStretch(self.exp)


class AsinhStretch(BaseStretch):
    r"""
    An asinh stretch.
    
    The stretch is given by:

    .. math::
        y = \frac{{\rm asinh}(x / a)}{{\rm asinh}(1 / a)}.
    """

    def __init__(self, a=0.1):
        super(AsinhStretch, self).__init__()
        self.a = a

    def __call__(self, values, out=None):
        res = np.arcsinh(values / self.a) / np.arcsinh(1. / self.a)
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return InvertedAsinhStretch(a=self.a)


class InvertedAsinhStretch(BaseStretch):
    """
    Inverse transformation for `~imageutils.scaling.AsinhStretch`.
    """

    def __init__(self, a=0.1):
        super(InvertedAsinhStretch, self).__init__()
        self.a = a

    def __call__(self, values, out=None):
        res = np.sinh(np.arcsinh(1. / self.a) * values) * self.a
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return AsinhStretch(a=self.a)


class SinhStretch(BaseStretch):
    r"""
    A sinh stretch.
    
    The stretch is given by:

    .. math::
        y = \frac{{\rm sinh}(x / a)}{{\rm sinh}(1 / a)}
    """

    def __init__(self, a=1./3.):
        super(SinhStretch, self).__init__()
        self.a = a

    def __call__(self, values, out=None):
        res = np.sinh(values / self.a) / np.sinh(1. / self.a)
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return InvertedSinhStretch(a=self.a)


class InvertedSinhStretch(BaseStretch):
    """
    Inverse transformation for `~imageutils.scaling.SinhStretch`.
    """

    def __init__(self, a=1./3.):
        super(InvertedSinhStretch, self).__init__()
        self.a = a

    def __call__(self, values, out=None):
        res = np.arcsinh(np.sinh(1. / self.a) * values) * self.a
        return np.clip(res, 0.0, 1.0, out=out)

    def inverted(self):
        return SinhStretch(a=self.a)


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
    Inverse transformation for `~imageutils.scaling.HistEqStretch`.
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

    The stretch is given by:

    .. math::
        y = (x - {\\rm bias}) * {\\rm contrast} + 0.5
    
    and the output values are clipped to the [0:1] range.
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
