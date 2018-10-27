# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Classes that deal with stretching, i.e. mapping a range of [0:1] values onto
another set of [0:1] values with a transformation
"""


import numpy as np

from ..utils.misc import InheritDocstrings
from .transform import BaseTransform


__all__ = ["BaseStretch", "LinearStretch", "SqrtStretch", "PowerStretch",
           "PowerDistStretch", "SquaredStretch", "LogStretch", "AsinhStretch",
           "SinhStretch", "HistEqStretch", "ContrastBiasStretch"]


def _logn(n, x, out=None):
    """Calculate the log base n of x."""
    # We define this because numpy.lib.scimath.logn doesn't support out=
    if out is None:
        return np.log(x) / np.log(n)
    else:
        np.log(x, out=out)
        np.true_divide(out, np.log(n), out=out)
        return out


def _prepare(values, clip=True, out=None):
    """
    Prepare the data by optionally clipping and copying, and return the
    array that should be subsequently used for in-place calculations.
    """

    if clip:
        return np.clip(values, 0., 1., out=out)
    else:
        if out is None:
            return np.array(values, copy=True)
        else:
            out[:] = np.asarray(values)
            return out


class BaseStretch(BaseTransform, metaclass=InheritDocstrings):
    """
    Base class for the stretch classes, which, when called with an array
    of values in the range [0:1], return an transformed array of values,
    also in the range [0:1].
    """

    def __call__(self, values, clip=True, out=None):
        """
        Transform values using this stretch.

        Parameters
        ----------
        values : array-like
            The input values, which should already be normalized to the
            [0:1] range.
        clip : bool, optional
            If `True` (default), values outside the [0:1] range are
            clipped to the [0:1] range.
        out : `~numpy.ndarray`, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).

        Returns
        -------
        result : `~numpy.ndarray`
            The transformed values.
        """

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""


class LinearStretch(BaseStretch):
    """
    A linear stretch.

    The stretch is given by:

    .. math::
        y = x
    """

    def __call__(self, values, clip=True, out=None):
        return _prepare(values, clip=clip, out=out)

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return LinearStretch()


class SqrtStretch(BaseStretch):
    r"""
    A square root stretch.

    The stretch is given by:

    .. math::
        y = \sqrt{x}
    """

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        with np.errstate(invalid='ignore'):
            np.sqrt(values, out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return PowerStretch(2)


class PowerStretch(BaseStretch):
    r"""
    A power stretch.

    The stretch is given by:

    .. math::
        y = x^a

    Parameters
    ----------
    a : float
        The power index (see the above formula).
    """

    def __init__(self, a):
        super().__init__()
        self.power = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.power(values, self.power, out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return PowerStretch(1. / self.power)


class PowerDistStretch(BaseStretch):
    r"""
    An alternative power stretch.

    The stretch is given by:

    .. math::
        y = \frac{a^x - 1}{a - 1}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  Default is 1000.
        ``a`` cannot be set to 1.
    """

    def __init__(self, a=1000.0):
        if a == 1:  # singularity
            raise ValueError("a cannot be set to 1")
        super().__init__()
        self.exp = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.power(self.exp, values, out=values)
        np.subtract(values, 1, out=values)
        np.true_divide(values, self.exp - 1.0, out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return InvertedPowerDistStretch(a=self.exp)


class InvertedPowerDistStretch(BaseStretch):
    r"""
    Inverse transformation for
    `~astropy.image.scaling.PowerDistStretch`.

    The stretch is given by:

    .. math::
        y = \frac{\log(y (a-1) + 1)}{\log a}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  Default is 1000.
        ``a`` cannot be set to 1.
    """

    def __init__(self, a=1000.0):
        if a == 1:  # singularity
            raise ValueError("a cannot be set to 1")
        super().__init__()
        self.exp = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.multiply(values, self.exp - 1.0, out=values)
        np.add(values, 1, out=values)
        _logn(self.exp, values, out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return PowerDistStretch(a=self.exp)


class SquaredStretch(PowerStretch):
    r"""
    A convenience class for a power stretch of 2.

    The stretch is given by:

    .. math::
        y = x^2
    """

    def __init__(self):
        super().__init__(2)

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return SqrtStretch()


class LogStretch(BaseStretch):
    r"""
    A log stretch.

    The stretch is given by:

    .. math::
        y = \frac{\log{(a x + 1)}}{\log{(a + 1)}}.

    Parameters
    ----------
    a : float
        The ``a`` parameter used in the above formula.  Default is 1000.
    """

    def __init__(self, a=1000.0):
        super().__init__()
        self.exp = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.multiply(values, self.exp, out=values)
        np.add(values, 1., out=values)
        np.log(values, out=values)
        np.true_divide(values, np.log(self.exp + 1.), out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return InvertedLogStretch(self.exp)


class InvertedLogStretch(BaseStretch):
    r"""
    Inverse transformation for `~astropy.image.scaling.LogStretch`.

    The stretch is given by:

    .. math::
        y = \frac{e^{y} (a + 1) -1}{a}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  Default is 1000.
    """

    def __init__(self, a):
        super().__init__()
        self.exp = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.multiply(values, np.log(self.exp + 1.), out=values)
        np.exp(values, out=values)
        np.subtract(values, 1., out=values)
        np.true_divide(values, self.exp, out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return LogStretch(self.exp)


class AsinhStretch(BaseStretch):
    r"""
    An asinh stretch.

    The stretch is given by:

    .. math::
        y = \frac{{\rm asinh}(x / a)}{{\rm asinh}(1 / a)}.

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  The value of
        this parameter is where the asinh curve transitions from linear
        to logarithmic behavior, expressed as a fraction of the
        normalized image.  Must be in the range between 0 and 1.
        Default is 0.1
    """

    def __init__(self, a=0.1):
        super().__init__()
        self.a = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.true_divide(values, self.a, out=values)
        np.arcsinh(values, out=values)
        np.true_divide(values, np.arcsinh(1. / self.a), out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return SinhStretch(a=1. / np.arcsinh(1. / self.a))


class SinhStretch(BaseStretch):
    r"""
    A sinh stretch.

    The stretch is given by:

    .. math::
        y = \frac{{\rm sinh}(x / a)}{{\rm sinh}(1 / a)}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  Default is 1/3.
    """

    def __init__(self, a=1./3.):
        super().__init__()
        self.a = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.true_divide(values, self.a, out=values)
        np.sinh(values, out=values)
        np.true_divide(values, np.sinh(1. / self.a), out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return AsinhStretch(a=1. / np.sinh(1. / self.a))


class HistEqStretch(BaseStretch):
    """
    A histogram equalization stretch.

    Parameters
    ----------
    data : array-like
        The data defining the equalization.
    values : array-like, optional
        The input image values, which should already be normalized to
        the [0:1] range.
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

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        values[:] = np.interp(values, self.data, self.values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return InvertedHistEqStretch(self.data, values=self.values)


class InvertedHistEqStretch(BaseStretch):
    """
    Inverse transformation for `~astropy.image.scaling.HistEqStretch`.

    Parameters
    ----------
    data : array-like
        The data defining the equalization.
    values : array-like, optional
        The input image values, which should already be normalized to
        the [0:1] range.
    """

    def __init__(self, data, values=None):
        self.data = data
        if values is None:
            self.values = np.linspace(0., 1., len(self.data))
        else:
            self.values = values

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        values[:] = np.interp(values, self.values, self.data)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return HistEqStretch(self.data, values=self.values)


class ContrastBiasStretch(BaseStretch):
    r"""
    A stretch that takes into account contrast and bias.

    The stretch is given by:

    .. math::
        y = (x - {\rm bias}) * {\rm contrast} + 0.5

    and the output values are clipped to the [0:1] range.

    Parameters
    ----------
    contrast : float
        The contrast parameter (see the above formula).

    bias : float
        The bias parameter (see the above formula).
    """

    def __init__(self, contrast, bias):
        super().__init__()
        self.contrast = contrast
        self.bias = bias

    def __call__(self, values, clip=True, out=None):
        # As a special case here, we only clip *after* the
        # transformation since it does not map [0:1] to [0:1]
        values = _prepare(values, clip=False, out=out)

        np.subtract(values, self.bias, out=values)
        np.multiply(values, self.contrast, out=values)
        np.add(values, 0.5, out=values)

        if clip:
            np.clip(values, 0, 1, out=values)

        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return InvertedContrastBiasStretch(self.contrast, self.bias)


class InvertedContrastBiasStretch(BaseStretch):
    """
    Inverse transformation for ContrastBiasStretch.

    Parameters
    ----------
    contrast : float
        The contrast parameter (see
        `~astropy.visualization.ConstrastBiasStretch).

    bias : float
        The bias parameter (see
        `~astropy.visualization.ConstrastBiasStretch).
    """

    def __init__(self, contrast, bias):
        super().__init__()
        self.contrast = contrast
        self.bias = bias

    def __call__(self, values, clip=True, out=None):
        # As a special case here, we only clip *after* the
        # transformation since it does not map [0:1] to [0:1]
        values = _prepare(values, clip=False, out=out)
        np.subtract(values, 0.5, out=values)
        np.true_divide(values, self.contrast, out=values)
        np.add(values, self.bias, out=values)

        if clip:
            np.clip(values, 0, 1, out=values)

        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return ContrastBiasStretch(self.contrast, self.bias)
