"""
Classes that deal with computing intervals from arrays of values based on
various criteria.
"""

import abc

import numpy as np

from astropy.extern import six

from .transform import BaseTransform

__all__ = ['ManualInterval', 'MinMaxInterval', 'PercentileInterval', 'AsymmetricPercentileInterval']


@six.add_metaclass(abc.ABCMeta)
class BaseInterval(BaseTransform):

    def __call__(self, values, clip=False, out=None):

        vmin, vmax = self.get_limits(values)

        if out is None:
            values = np.subtract(values, vmin)
        else:
            values = np.subtract(values, vmin, out=out)

        np.true_divide(values, vmax - vmin, out=values)

        if clip:
            np.clip(values, 0., 1., out=values)

        return values


class ManualInterval(BaseInterval):
    """
    Interval based on user-specified values.

    Parameters
    ----------
    vmin : float
        The minimum value in the scaling
    vmax : float
        The maximum value in the scaling
    """

    def __init__(self, vmin, vmax):
        self.vmin = vmin
        self.vmax = vmax

    def get_limits(self, values):
        return self.vmin, self.vmax


class MinMaxInterval(BaseInterval):
    """
    Interval based on the minimum and maximum values in the data.
    """

    def get_limits(self, values):
        return np.min(values), np.max(values)


class AsymmetricPercentileInterval(BaseInterval):
    """
    Interval based on a keeping a specified fraction of pixels (can be asymmetric).

    Parameters
    ----------
    lower_percentile : float
        The lower percentile below which to ignore pixels.
    upper_percentile : float
        The upper percentile above which to ignore pixels.
    n_samples : int, optional
        Maximum number of values to use. If this is specified, and there are
        more values in the dataset as this, then values are randomly sampled
        from the array (with replacement)
    """

    def __init__(self, lower_percentile, upper_percentile, n_samples=None):
        self.lower_percentile = lower_percentile
        self.upper_percentile = upper_percentile
        self.n_samples = n_samples

    def get_limits(self, values):

        # If needed, limit the number of samples. We sample with replacement
        # since this is much faster.
        if self.n_samples is not None and values.size > self.n_samples:
            values = np.random.choice(values, self.n_samples)

        # Filter out invalid values (inf, nan)
        values = values[np.isfinite(values)]

        # Determine values at percentiles
        vmin = np.percentile(values, self.lower_percentile)
        vmax = np.percentile(values, self.upper_percentile)

        return vmin, vmax


class PercentileInterval(AsymmetricPercentileInterval):
    """
    Interval based on a keeping a specified fraction of pixels.

    Parameters
    ----------
    percentile : float
        The fraction of pixels to keep. The same fraction of pixels is
        eliminated from both ends.
    n_samples : int, optional
        Maximum number of values to use. If this is specified, and there are
        more values in the dataset as this, then values are randomly sampled
        from the array (with replacement)
    """

    def __init__(self, percentile, n_samples=None):
        lower_percentile = (100 - percentile) * 0.5
        upper_percentile = 100 - lower_percentile
        super(PercentileInterval, self).__init__(lower_percentile, upper_percentile, n_samples=n_samples)
