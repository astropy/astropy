# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Classes that deal with stretching, i.e. mapping a range of [0:1] values onto
another set of [0:1] values with a transformation.
"""

import numpy as np

from .transform import BaseTransform, CompositeTransform

__all__ = [
    "AsinhStretch",
    "BaseStretch",
    "CompositeStretch",
    "ContrastBiasStretch",
    "HistEqStretch",
    "LinearStretch",
    "LogStretch",
    "PowerDistStretch",
    "PowerStretch",
    "SinhStretch",
    "SqrtStretch",
    "SquaredStretch",
]


def _logn(n, x, out=None):
    """Calculate the log base n of x."""
    # We define this because numpy.emath.logn doesn't support the out
    # keyword.
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
        return np.clip(values, 0.0, 1.0, out=out)
    else:
        if out is None:
            return np.array(values, copy=True)
        else:
            out[:] = np.asarray(values)
            return out


class BaseStretch(BaseTransform):
    """
    Base class for the stretch classes, which when called with an array
    of values in the range [0:1], returns an transformed array of values
    also in the range [0:1].
    """

    @property
    def _supports_invalid_kw(self):
        return False

    def __add__(self, other):
        return CompositeStretch(other, self)

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
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).

        Returns
        -------
        result : ndarray
            The transformed values.
        """

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""


class LinearStretch(BaseStretch):
    """
    A linear stretch with a slope and offset.

    The stretch is given by:

    .. math::
        y = slope * x + intercept

    Parameters
    ----------
    slope : float, optional
        The ``slope`` parameter used in the above formula.  Default is 1.
    intercept : float, optional
        The ``intercept`` parameter used in the above formula.  Default is 0.

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import LinearStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        slopes = [1, 0.5, 1.3, 1.4, 2.0]
        intercepts = [0, 0.0, -0.4, 0., 0.2]
        for slope, intercept in zip(slopes, intercepts):
            stretch = LinearStretch(slope, intercept)
            label = f'{slope=}, {intercept=}'
            ax.plot(x, stretch(x, clip=True), label=label)

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='lower right', fontsize=8)
    """

    def __init__(self, slope=1, intercept=0):
        super().__init__()
        self.slope = slope
        self.intercept = intercept

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        if self.slope != 1:
            np.multiply(values, self.slope, out=values)
        if self.intercept != 0:
            np.add(values, self.intercept, out=values)

        if clip:
            np.clip(values, 0, 1, out=values)

        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return LinearStretch(1.0 / self.slope, -self.intercept / self.slope)


class SqrtStretch(BaseStretch):
    r"""
    A square root stretch.

    The stretch is given by:

    .. math::
        y = \sqrt{x}

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import SqrtStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        stretch = SqrtStretch()
        ax.plot(x, stretch(x, clip=True))

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='lower right', fontsize=8)
    """

    @property
    def _supports_invalid_kw(self):
        return True

    def __call__(self, values, clip=True, out=None, invalid=None):
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
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).
        invalid : None or float, optional
            Value to assign NaN values generated by this class.  NaNs in
            the input ``values`` array are not changed.  This option is
            generally used with matplotlib normalization classes, where
            the ``invalid`` value should map to the matplotlib colormap
            "under" value (i.e., any finite value < 0).  If `None`, then
            NaN values are not replaced.  This keyword has no effect if
            ``clip=True``.

        Returns
        -------
        result : ndarray
            The transformed values.
        """
        values = _prepare(values, clip=clip, out=out)
        replace_invalid = not clip and invalid is not None
        with np.errstate(invalid="ignore"):
            if replace_invalid:
                idx = values < 0
            np.sqrt(values, out=values)

        if replace_invalid:
            # Assign new NaN (i.e., NaN not in the original input
            # values, but generated by this class) to the invalid value.
            values[idx] = invalid

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
        The power index (see the above formula). ``a`` must be greater
        than 0.

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import PowerStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        a_vals = (0.3, 0.5, 0.7, 1, 1.5, 2, 3)
        for a in a_vals:
            stretch = PowerStretch(a)
            label = f'{a=}'
            ax.plot(x, stretch(x, clip=True), label=label)

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='lower right', fontsize=8)
    """

    @property
    def _supports_invalid_kw(self):
        return True

    def __init__(self, a):
        super().__init__()
        if a <= 0:
            raise ValueError("a must be > 0")
        self.a = a

    def __call__(self, values, clip=True, out=None, invalid=None):
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
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).
        invalid : None or float, optional
            Value to assign NaN values generated by this class.  NaNs in
            the input ``values`` array are not changed.  This option is
            generally used with matplotlib normalization classes, where
            the ``invalid`` value should map to the matplotlib colormap
            "under" value (i.e., any finite value < 0).  If `None`, then
            NaN values are not replaced.  This keyword has no effect if
            ``clip=True``.

        Returns
        -------
        result : ndarray
            The transformed values.
        """
        values = _prepare(values, clip=clip, out=out)
        replace_invalid = (
            not clip and invalid is not None and ((-1 < self.a < 0) or (0 < self.a < 1))
        )
        with np.errstate(invalid="ignore"):
            if replace_invalid:
                idx = values < 0
            np.power(values, self.a, out=values)

        if replace_invalid:
            # Assign new NaN (i.e., NaN not in the original input
            # values, but generated by this class) to the invalid value.
            values[idx] = invalid

        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return PowerStretch(1.0 / self.a)


class PowerDistStretch(BaseStretch):
    r"""
    An alternative power stretch.

    The stretch is given by:

    .. math::
        y = \frac{a^x - 1}{a - 1}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula. The stretch
        becomes more linear as ``a`` approaches 1, more exponential for
        ``a`` values greater than 1, and more logarithmic for ``a``
        values less than 1. ``a`` must be greater than 0, but cannot be
        set to 1. Default is 1000.

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import PowerDistStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        a_vals = (0.001, 0.05, 0.3, 0.8, 1.2, 3, 10, 30, 100, 1000)
        for a in a_vals:
            if a == 1000:
                lw = 3
            else:
                lw = 1
            stretch = PowerDistStretch(a)
            label = f'{a=}'
            ax.plot(x, stretch(x, clip=True), label=label, lw=lw)

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='upper left', fontsize=8)
    """

    def __init__(self, a=1000.0):
        if a <= 0 or a == 1:  # singularity
            raise ValueError("a must be > 0, but cannot be set to 1")
        super().__init__()
        self.a = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.power(self.a, values, out=values)
        np.subtract(values, 1, out=values)
        np.true_divide(values, self.a - 1.0, out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return InvertedPowerDistStretch(a=self.a)


class InvertedPowerDistStretch(BaseStretch):
    r"""
    Inverse transformation for
    `~astropy.visualization.PowerDistStretch`.

    The stretch is given by:

    .. math::
        y = \frac{\log(x (a - 1) + 1)}{\log a}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula. The stretch
        becomes more linear as ``a`` approaches 1, more logarithmic for
        ``a`` values greater than 1, and more exponential for ``a``
        values less than 1. ``a`` must be greater than 0, but cannot be
        set to 1. Default is 1000.
    """

    def __init__(self, a=1000.0):
        if a <= 0 or a == 1:  # singularity
            raise ValueError("a must be > 0, but cannot be set to 1")
        super().__init__()
        self.a = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.multiply(values, self.a - 1.0, out=values)
        np.add(values, 1, out=values)
        _logn(self.a, values, out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return PowerDistStretch(a=self.a)


class SquaredStretch(PowerStretch):
    r"""
    A convenience class for a power stretch of 2.

    The stretch is given by:

    .. math::
        y = x^2

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import SquaredStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        stretch = SquaredStretch()
        ax.plot(x, stretch(x, clip=True))

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='lower right', fontsize=8)
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
        y = \frac{\log{(a x + 1)}}{\log{(a + 1)}}

    Parameters
    ----------
    a : float
        The ``a`` parameter used in the above formula. The stretch
        becomes more linear for small ``a`` values. ``a`` must be
        greater than 0. Default is 1000.

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import LogStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        a_vals = (0.1, 1, 3, 10, 30, 100, 1000, 10000)
        for a in a_vals:
            if a == 1000:
                lw = 3
            else:
                lw = 1
            stretch = LogStretch(a)
            label = f'{a=}'
            ax.plot(x, stretch(x, clip=True), label=label, lw=lw)

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='lower right', fontsize=8)
    """

    @property
    def _supports_invalid_kw(self):
        return True

    def __init__(self, a=1000.0):
        super().__init__()
        if a <= 0:  # singularity
            raise ValueError("a must be > 0")
        self.a = a

    def __call__(self, values, clip=True, out=None, invalid=None):
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
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).
        invalid : None or float, optional
            Value to assign NaN values generated by this class.  NaNs in
            the input ``values`` array are not changed.  This option is
            generally used with matplotlib normalization classes, where
            the ``invalid`` value should map to the matplotlib colormap
            "under" value (i.e., any finite value < 0).  If `None`, then
            NaN values are not replaced.  This keyword has no effect if
            ``clip=True``.

        Returns
        -------
        result : ndarray
            The transformed values.
        """
        values = _prepare(values, clip=clip, out=out)
        replace_invalid = not clip and invalid is not None
        with np.errstate(invalid="ignore"):
            if replace_invalid:
                idx = values < 0
            np.multiply(values, self.a, out=values)
            np.add(values, 1.0, out=values)
            np.log(values, out=values)
            np.true_divide(values, np.log(self.a + 1.0), out=values)

        if replace_invalid:
            # Assign new NaN (i.e., NaN not in the original input
            # values, but generated by this class) to the invalid value.
            values[idx] = invalid

        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return InvertedLogStretch(self.a)


class InvertedLogStretch(BaseStretch):
    r"""
    Inverse transformation for `~astropy.visualization.LogStretch`.

    The stretch is given by:

    .. math::
        y = \frac{e^{x \log{a + 1}} - 1}{a} = \frac{(a + 1)^x - 1}{a}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula. The stretch
        becomes more linear for small ``a`` values and more exponential
        for large ``a`` values. ``a`` must be greater than 0. Default is
        1000.
    """

    def __init__(self, a):
        super().__init__()
        if a <= 0:  # singularity
            raise ValueError("a must be > 0")
        self.a = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.multiply(values, np.log(self.a + 1.0), out=values)
        np.exp(values, out=values)
        np.subtract(values, 1.0, out=values)
        np.true_divide(values, self.a, out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return LogStretch(self.a)


class AsinhStretch(BaseStretch):
    r"""
    An asinh stretch.

    The stretch is given by:

    .. math::
        y = \frac{{\rm asinh}(x / a)}{{\rm asinh}(1 / a)}.

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula. The value of this
        parameter is where the asinh curve transitions from linear to
        logarithmic behavior, expressed as a fraction of the normalized
        image. The stretch becomes more linear for larger ``a`` values
        and more logarithmic for smaller ``a`` values. ``a`` must be
        greater than 0. Default is 0.1.

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import AsinhStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        a_vals = (0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.9, 3.0)
        for a in a_vals:
            if a == 0.1:
                lw = 3
            else:
                lw = 1
            stretch = AsinhStretch(a)
            label = f'{a=}'
            ax.plot(x, stretch(x, clip=True), label=label, lw=lw)

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='lower right', fontsize=8)
    """

    def __init__(self, a=0.1):
        super().__init__()
        if a <= 0:
            raise ValueError("a must be > 0")
        self.a = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.true_divide(values, self.a, out=values)
        np.arcsinh(values, out=values)
        np.true_divide(values, np.arcsinh(1.0 / self.a), out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return SinhStretch(a=1.0 / np.arcsinh(1.0 / self.a))


class SinhStretch(BaseStretch):
    r"""
    A sinh stretch.

    The stretch is given by:

    .. math::
        y = \frac{{\rm sinh}(x / a)}{{\rm sinh}(1 / a)}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula. The stretch
        becomes more linear for larger ``a`` values and more exponential
        for smaller ``a`` values. ``a`` must be greater than 0. Default
        is 1/3.

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import SinhStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        a_vals = (0.1, 0.2, 0.3333, 0.5, 0.9, 3)
        for a in a_vals:
            if a == 0.3333:
                lw = 3
            else:
                lw = 1
            stretch = SinhStretch(a)
            label = f'{a=}'
            ax.plot(x, stretch(x, clip=True), label=label, lw=lw)

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='upper left', fontsize=8)
    """

    def __init__(self, a=1.0 / 3.0):
        super().__init__()
        if a <= 0:
            raise ValueError("a must be > 0")
        self.a = a

    def __call__(self, values, clip=True, out=None):
        values = _prepare(values, clip=clip, out=out)
        np.true_divide(values, self.a, out=values)
        np.sinh(values, out=values)
        np.true_divide(values, np.sinh(1.0 / self.a), out=values)
        return values

    @property
    def inverse(self):
        """A stretch object that performs the inverse operation."""
        return AsinhStretch(a=1.0 / np.sinh(1.0 / self.a))


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
        self.data = self.data[np.isfinite(self.data)]
        vmin = self.data.min()
        vmax = self.data.max()
        self.data = (self.data - vmin) / (vmax - vmin)

        # Compute relative position of each pixel
        if values is None:
            self.values = np.linspace(0.0, 1.0, len(self.data))
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
    Inverse transformation for `~astropy.visualization.HistEqStretch`.

    Parameters
    ----------
    data : array-like
        The data defining the equalization.
    values : array-like, optional
        The input image values, which should already be normalized to
        the [0:1] range.
    """

    def __init__(self, data, values=None):
        self.data = data[np.isfinite(data)]
        if values is None:
            self.values = np.linspace(0.0, 1.0, len(self.data))
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

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import ContrastBiasStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))

        x = np.linspace(0, 1, 100)
        contrasts = [1.0, 2.0, 0.7, 1.0, 1.0, 2.0]
        biases = [0.5, 0.5, 0.5, 0.3, 0.7, 0.3]
        for contrast, bias in zip(contrasts, biases):
            stretch = ContrastBiasStretch(contrast, bias)
            ax.plot(x, stretch(x, clip=True), label=f'{contrast=}, {bias=}')

        ax.axis('equal')
        ax.plot(x, x, ls='dotted', color='k', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Input Value')
        ax.set_ylabel('Output Value')
        ax.set_title(stretch.__class__.__name__)
        ax.legend(loc='lower right', fontsize=8)
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
    r"""
    Inverse transformation for
    `~astropy.visualization.ContrastBiasStretch`.

    The stretch is given by:

    .. math::
        y = \frac{x - 0.5}{{\rm contrast}} + {\rm bias}

    Parameters
    ----------
    contrast : float
        The contrast parameter (see
        `~astropy.visualization.ContrastBiasStretch`).

    bias : float
        The bias parameter (see
        `~astropy.visualization.ContrastBiasStretch`).
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


class CompositeStretch(CompositeTransform, BaseStretch):
    """
    A combination of two stretches.

    Parameters
    ----------
    stretch_1 : :class:`astropy.visualization.BaseStretch`
        The first stretch to apply.
    stretch_2 : :class:`astropy.visualization.BaseStretch`
        The second stretch to apply.
    """

    def __call__(self, values, clip=True, out=None):
        return self.transform_2(
            self.transform_1(values, clip=clip, out=out), clip=clip, out=out
        )
