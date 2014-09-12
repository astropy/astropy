"""
Classes that deal with computing intervals from arrays of values based on
various criteria.
"""

import abc

import numpy as np

from astropy.extern import six

__all__ = ['ManualInterval', 'MinMaxInterval', 'PercentileInterval']


@six.add_metaclass(abc.ABCMeta)
class BaseInterval(object):
    pass


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


class PercentileInterval(BaseInterval):
    """
    Interval based on a keeping a specified fraction of pixels

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
        self.percentile = percentile
        self.n_samples = n_samples

    def get_limits(self, values):

        # If needed, limmit the number of samples. We sample with replacement
        # since this is much faster.
        if self.n_samples is not None and values.size > self.n_samples:
            values = np.random.choice(values, self.n_samples)

        pmin = (100. - self.percentile) * 0.5
        pmax = 100. - pmin

        # Filter out invalid values (inf, nan)
        values = values[np.isfinite(values)]

        # Determine values at percentiles
        vmin = np.percentile(values, pmin)
        vmax = np.percentile(values, pmax)

        return vmin, vmax
