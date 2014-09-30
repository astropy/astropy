"""
Classes that deal with stretching, i.e. mapping a range of [0:1] values onto
another set of [0:1] values with a transformation
"""

import abc

import numpy as np

from astropy.extern import six

from .transform import BaseTransform

__all__ = ["LinearStretch", "SqrtStretch", "PowerStretch", "PowerDistStretch",
           "SquaredStretch", "LogStretch", "AsinhStretch", "SinhStretch",
           "HistEqStretch", "ContrastBiasStretch"]


def logn(n, x, out=None):
    # We define this because numpy.lib.scimath.logn doesn't support out=
    if out is None:
        return np.log(x) / np.log(n)
    else:
        np.log(x, out=out)
        np.true_divide(out, np.log(n), out=out)
        return out


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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = values.copy()
        else:
            out[:] = values
            values = out

        if clip:
            values = np.clip(values, 0., 1., out=values)

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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.sqrt(values)
        else:
            values = np.sqrt(values, out=out)

        if clip:
            values = np.clip(values, 0., 1., out=values)

        return values

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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.power(values, self.power)
        else:
            values = np.power(values, self.power, out=out)

        if clip:
            values = np.clip(values, 0., 1., out=values)

        return values

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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.power(self.exp, values)
        else:
            values = np.power(self.exp, values, out=out)

        np.subtract(values, 1, out=values)
        np.true_divide(values, self.exp - 1.0, out=values)

        if clip:
            np.clip(values, 0., 1., out=values)

        return values

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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.multiply(values, self.exp - 1.0)
        else:
            values = np.multiply(values, self.exp - 1.0, out=values)

        np.add(values, 1, out=values)
        logn(self.exp, values, out=values)

        if clip:
            np.clip(values, 0., 1., out=values)

        return values

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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.multiply(values, self.exp)
        else:
            values = np.multiply(values, self.exp, out=out)

        np.add(values, 1., out=values)
        np.log(values, out=values)
        np.true_divide(values, np.log(self.exp + 1.), out=values)

        if clip:
            np.clip(values, 0., 1., out=values)

        return values

    def inverted(self):
        return InvertedLogStretch(self.exp)


class InvertedLogStretch(BaseStretch):
    """
    Inverse transformation for `~imageutils.scaling.LogStretch`.
    """

    def __init__(self, a):
        super(InvertedLogStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.multiply(values, np.log(self.exp + 1.))
        else:
            values = np.multiply(values, np.log(self.exp + 1.), out=out)

        np.exp(values, out=values)
        np.subtract(values, 1., out=values)
        np.true_divide(values, self.exp, out=values)

        if clip:
            np.clip(values, 0., 1., out=values)

        return values

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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.true_divide(values, self.a)
        else:
            values = np.true_divide(values, self.a, out=out)

        np.arcsinh(values, out=values)
        np.true_divide(values, np.arcsinh(1. / self.a), out=values)

        if clip:
            np.clip(values, 0., 1., out=values)

        return values

    def inverted(self):
        return SinhStretch(a=1./np.arcsinh(1. / self.a))


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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.true_divide(values, self.a)
        else:
            values = np.true_divide(values, self.a, out=out)

        np.sinh(values, out=values)
        np.true_divide(values, np.sinh(1. / self.a), out=values)

        if clip:
            np.clip(values, 0., 1., out=values)

        return values

    def inverted(self):
        return AsinhStretch(a=1./np.sinh(1. / self.a))


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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.interp(values, self.data, self.values)
        else:
            out[:] = np.interp(values, self.data, self.values)
            values = out

        if clip:
            np.clip(values, 0., 1., out=values)

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

    def __call__(self, values, out=None, clip=False):
        if out is None:
            values = np.interp(values, self.values, self.data)
        else:
            values[:] = np.interp(values, self.values, self.data)

        if clip:
            np.clip(values, 0., 1., out=values)

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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.subtract(values, self.bias)
        else:
            values = np.subtract(values, self.bias, out=out)

        # Use in-place operations for the rest since they are faster
        np.multiply(values, self.contrast, out=values)
        np.add(values, 0.5, out=values)

        if clip:
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

    def __call__(self, values, out=None, clip=False):

        if out is None:
            values = np.subtract(values, 0.5)
        else:
            np.subtract(values, 0.5, out=values)

        # Use in-place operations for the rest since they are faster
        np.true_divide(values, self.contrast, out=values)
        np.add(values, self.bias, out=values)

        if clip:
            np.clip(values, 0, 1, out=values)

        return values

    def inverted(self):
        return ContrastBiasStretch(self.contrast, self.bias)
