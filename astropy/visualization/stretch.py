# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Classes that deal with stretching, i.e. mapping a range of [0:1] values onto
another set of [0:1] values with a transformation
"""

from __future__ import division, print_function

import numpy as np

from ..extern import six
from ..utils.misc import InheritDocstrings

from .transform import BaseTransform

__all__ = ["BaseStretch", "LinearStretch", "SqrtStretch", "PowerStretch",
           "PowerDistStretch", "SquaredStretch", "LogStretch", "AsinhStretch",
           "SinhStretch", "HistEqStretch", "ContrastBiasStretch"]


def logn(n, x, out=None):
    # We define this because numpy.lib.scimath.logn doesn't support out=
    if out is None:
        return np.log(x) / np.log(n)
    else:
        np.log(x, out=out)
        np.true_divide(out, np.log(n), out=out)
        return out


def _prepare(values, out=None, clip=True):
    """
    Prepare the data by optionally clipping and copying, and return the array
    that should be subsequently used for in-place calculations.
    """
    if clip:
        return np.clip(values, 0., 1., out=out)
    else:
        if out is None:
            return np.array(values, copy=True)
        else:
            out[:] = np.asarray(values)
            return out


@six.add_metaclass(InheritDocstrings)
class BaseStretch(BaseTransform):
    """
    Base class for the stretch classes, which, when called with an array of
    values in the range [0:1], return an transformed array of values, also in
    the range [0:1].
    """

    def __call__(self, values, out=None, clip=True):
        """
        Transform values using this stretch.

        Parameters
        ----------
        values : `~numpy.ndarray` or list
            The input values, which should already be normalized to the [0:1]
            range.
        out : `~numpy.ndarray`, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).
        clip : bool, optional
            If `True` (default), values outside the [0:1] range are clipped to
            the [0:1] range.

        Returns
        -------
        new_values : `~numpy.ndarray`
            The transformed values.
        """

    @property
    def inverse(self):
        """
        Return a stretch that performs the inverse operation.
        """


class LinearStretch(BaseStretch):
    """
    A linear stretch.

    The stretch is given by:

    .. math::
        y = x
    """

    def __call__(self, values, out=None, clip=True):

        return _prepare(values, out=out, clip=clip)

    @property
    def inverse(self):
        return LinearStretch()


class SqrtStretch(BaseStretch):
    r"""
    A square root stretch.

    The stretch is given by:

    .. math::
        y = \sqrt{x}
    """

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        np.sqrt(values, out=values)

        return values

    @property
    def inverse(self):
        return PowerStretch(2)


class PowerStretch(BaseStretch):
    r"""
    A power stretch.

    The stretch is given by:

    .. math::
        y = x^a
    """

    def __init__(self, a):
        super(PowerStretch, self).__init__()
        self.power = a

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        np.power(values, self.power, out=values)

        return values

    @property
    def inverse(self):
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

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        np.power(self.exp, values, out=values)
        np.subtract(values, 1, out=values)
        np.true_divide(values, self.exp - 1.0, out=values)

        return values

    @property
    def inverse(self):
        return InvertedPowerDistStretch(a=self.exp)


class InvertedPowerDistStretch(BaseStretch):
    """
    Inverse transformation for `~astropy.image.scaling.PowerDistStretch`.
    """

    def __init__(self, a=1000.0):
        if a == 1:  # singularity
            raise ValueError("a cannot be set to 1")
        super(InvertedPowerDistStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        np.multiply(values, self.exp - 1.0, out=values)
        np.add(values, 1, out=values)
        logn(self.exp, values, out=values)

        return values

    @property
    def inverse(self):
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

    @property
    def inverse(self):
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

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        np.multiply(values, self.exp, out=values)
        np.add(values, 1., out=values)
        np.log(values, out=values)
        np.true_divide(values, np.log(self.exp + 1.), out=values)

        return values

    @property
    def inverse(self):
        return InvertedLogStretch(self.exp)


class InvertedLogStretch(BaseStretch):
    """
    Inverse transformation for `~astropy.image.scaling.LogStretch`.
    """

    def __init__(self, a):
        super(InvertedLogStretch, self).__init__()
        self.exp = a

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        np.multiply(values, np.log(self.exp + 1.), out=values)
        np.exp(values, out=values)
        np.subtract(values, 1., out=values)
        np.true_divide(values, self.exp, out=values)

        return values

    @property
    def inverse(self):
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

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        np.true_divide(values, self.a, out=values)
        np.arcsinh(values, out=values)
        np.true_divide(values, np.arcsinh(1. / self.a), out=values)

        return values

    @property
    def inverse(self):
        return SinhStretch(a=1. / np.arcsinh(1. / self.a))


class SinhStretch(BaseStretch):
    r"""
    A sinh stretch.

    The stretch is given by:

    .. math::
        y = \frac{{\rm sinh}(x / a)}{{\rm sinh}(1 / a)}
    """

    def __init__(self, a=1. / 3.):
        super(SinhStretch, self).__init__()
        self.a = a

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        np.true_divide(values, self.a, out=values)
        np.sinh(values, out=values)
        np.true_divide(values, np.sinh(1. / self.a), out=values)

        return values

    @property
    def inverse(self):
        return AsinhStretch(a=1. / np.sinh(1. / self.a))


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

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        values[:] = np.interp(values, self.data, self.values)

        return values

    @property
    def inverse(self):
        return InvertedHistEqStretch(self.data, values=self.values)


class InvertedHistEqStretch(BaseStretch):
    """
    Inverse transformation for `~astropy.image.scaling.HistEqStretch`.
    """

    def __init__(self, data, values=None):
        self.data = data
        if values is None:
            self.values = np.linspace(0., 1., len(self.data))
        else:
            self.values = values

    def __call__(self, values, out=None, clip=True):

        values = _prepare(values, out=out, clip=clip)

        values[:] = np.interp(values, self.values, self.data)

        return values

    @property
    def inverse(self):
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

    def __call__(self, values, out=None, clip=True):

        # As a special case here, we only clip *after* the transformation since
        # it does not map [0:1] to [0:1]
        values = _prepare(values, out=out, clip=False)

        np.subtract(values, self.bias, out=values)
        np.multiply(values, self.contrast, out=values)
        np.add(values, 0.5, out=values)

        if clip:
            np.clip(values, 0, 1, out=values)

        return values

    @property
    def inverse(self):
        return InvertedContrastBiasStretch(self.contrast, self.bias)


class InvertedContrastBiasStretch(BaseStretch):
    """
    Inverse transformation for ContrastBiasStretch.
    """

    def __init__(self, contrast, bias):
        super(InvertedContrastBiasStretch, self).__init__()
        self.contrast = contrast
        self.bias = bias

    def __call__(self, values, out=None, clip=True):

        # As a special case here, we only clip *after* the transformation since
        # it does not map [0:1] to [0:1]

        values = _prepare(values, out=out, clip=False)

        np.subtract(values, 0.5, out=values)
        np.true_divide(values, self.contrast, out=values)
        np.add(values, self.bias, out=values)

        if clip:
            np.clip(values, 0, 1, out=values)

        return values

    @property
    def inverse(self):
        return ContrastBiasStretch(self.contrast, self.bias)
