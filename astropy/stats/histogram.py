# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Methods for selecting the bin width of histograms

Ported from the astroML project: http://astroML.org/
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..extern import six

import numpy as np
from . import bayesian_blocks

__all__ = ['histogram', 'scott_bin_width', 'freedman_bin_width',
           'knuth_bin_width']


def histogram(a, bins=10, range=None, weights=None, **kwargs):
    """Enhanced histogram function, providing adaptive binnings

    This is a histogram function that enables the use of more sophisticated
    algorithms for determining bins.  Aside from the ``bins`` argument allowing
    a string specified how bins are computed, the parameters are the same
    as ``numpy.histogram()``.

    Parameters
    ----------
    a : array_like
        array of data to be histogrammed

    bins : int or list or str (optional)
        If bins is a string, then it must be one of:

        - 'blocks' : use bayesian blocks for dynamic bin widths

        - 'knuth' : use Knuth's rule to determine bins

        - 'scott' : use Scott's rule to determine bins

        - 'freedman' : use the Freedman-Diaconis rule to determine bins

    range : tuple or None (optional)
        the minimum and maximum range for the histogram.  If not specified,
        it will be (x.min(), x.max())

    weights : array_like, optional
        Not Implemented

    other keyword arguments are described in numpy.histogram().

    Returns
    -------
    hist : array
        The values of the histogram. See ``normed`` and ``weights`` for a
        description of the possible semantics.
    bin_edges : array of dtype float
        Return the bin edges ``(length(hist)+1)``.

    See Also
    --------
    numpy.histogram
    """
    # if bins is a string, first compute bin edges with the desired heuristic
    if isinstance(bins, six.string_types):
        a = np.asarray(a).ravel()

        # TODO: if weights is specified, we need to modify things.
        #       e.g. we could use point measures fitness for Bayesian blocks
        if weights is not None:
            raise NotImplementedError("weights are not yet supported "
                                      "for the enhanced histogram")

        # if range is specified, we need to truncate the data for
        # the bin-finding routines
        if range is not None:
            a = a[(a >= range[0]) & (a <= range[1])]

        if bins == 'blocks':
            bins = bayesian_blocks(a)
        elif bins == 'knuth':
            da, bins = knuth_bin_width(a, True)
        elif bins == 'scott':
            da, bins = scott_bin_width(a, True)
        elif bins == 'freedman':
            da, bins = freedman_bin_width(a, True)
        else:
            raise ValueError("unrecognized bin code: '{}'".format(bins))

    # Now we call numpy's histogram with the resulting bin edges
    return np.histogram(a, bins=bins, range=range, weights=weights, **kwargs)


def scott_bin_width(data, return_bins=False):
    r"""Return the optimal histogram bin width using Scott's rule

    Scott's rule is a normal reference rule: it minimizes the integrated
    mean squared error in the bin approximation under the assumption that the
    data is approximately Gaussian.

    Parameters
    ----------
    data : array-like, ndim=1
        observed (one-dimensional) data
    return_bins : bool (optional)
        if True, then return the bin edges

    Returns
    -------
    width : float
        optimal bin width using Scott's rule
    bins : ndarray
        bin edges: returned if ``return_bins`` is True

    Notes
    -----
    The optimal bin width is

    .. math::
        \Delta_b = \frac{3.5\sigma}{n^{1/3}}

    where :math:`\sigma` is the standard deviation of the data, and
    :math:`n` is the number of data points [1]_.

    References
    ----------
    .. [1] Scott, David W. (1979). "On optimal and data-based histograms".
       Biometricka 66 (3): 605-610

    See Also
    --------
    knuth_bin_width
    freedman_bin_width
    bayesian_blocks
    histogram
    """
    data = np.asarray(data)
    if data.ndim != 1:
        raise ValueError("data should be one-dimensional")

    n = data.size
    sigma = np.std(data)

    dx = 3.5 * sigma / (n ** (1 / 3))

    if return_bins:
        Nbins = np.ceil((data.max() - data.min()) / dx)
        Nbins = max(1, Nbins)
        bins = data.min() + dx * np.arange(Nbins + 1)
        return dx, bins
    else:
        return dx


def freedman_bin_width(data, return_bins=False):
    r"""Return the optimal histogram bin width using the Freedman-Diaconis rule

    The Freedman-Diaconis rule is a normal reference rule like Scott's
    rule, but uses rank-based statistics for results which are more robust
    to deviations from a normal distribution.

    Parameters
    ----------
    data : array-like, ndim=1
        observed (one-dimensional) data
    return_bins : bool (optional)
        if True, then return the bin edges

    Returns
    -------
    width : float
        optimal bin width using the Freedman-Diaconis rule
    bins : ndarray
        bin edges: returned if ``return_bins`` is True

    Notes
    -----
    The optimal bin width is

    .. math::
        \Delta_b = \frac{2(q_{75} - q_{25})}{n^{1/3}}

    where :math:`q_{N}` is the :math:`N` percent quartile of the data, and
    :math:`n` is the number of data points [1]_.

    References
    ----------
    .. [1] D. Freedman & P. Diaconis (1981)
       "On the histogram as a density estimator: L2 theory".
       Probability Theory and Related Fields 57 (4): 453-476

    See Also
    --------
    knuth_bin_width
    scott_bin_width
    bayesian_blocks
    histogram
    """
    data = np.asarray(data)
    if data.ndim != 1:
        raise ValueError("data should be one-dimensional")

    n = data.size
    if n < 4:
        raise ValueError("data should have more than three entries")

    v25, v75 = np.percentile(data, [25, 75])
    dx = 2 * (v75 - v25) / (n ** (1 / 3))

    if return_bins:
        dmin, dmax = data.min(), data.max()
        Nbins = max(1, np.ceil((dmax - dmin) / dx))
        bins = dmin + dx * np.arange(Nbins + 1)
        return dx, bins
    else:
        return dx


def knuth_bin_width(data, return_bins=False, quiet=True):
    r"""Return the optimal histogram bin width using Knuth's rule.

    Knuth's rule is a fixed-width, Bayesian approach to determining
    the optimal bin width of a histogram.

    Parameters
    ----------
    data : array-like, ndim=1
        observed (one-dimensional) data
    return_bins : bool (optional)
        if True, then return the bin edges
    quiet : bool (optional)
        if True (default) then suppress stdout output from scipy.optimize

    Returns
    -------
    dx : float
        optimal bin width. Bins are measured starting at the first data point.
    bins : ndarray
        bin edges: returned if ``return_bins`` is True

    Notes
    -----
    The optimal number of bins is the value M which maximizes the function

    .. math::
        F(M|x,I) = n\log(M) + \log\Gamma(\frac{M}{2})
        - M\log\Gamma(\frac{1}{2})
        - \log\Gamma(\frac{2n+M}{2})
        + \sum_{k=1}^M \log\Gamma(n_k + \frac{1}{2})

    where :math:`\Gamma` is the Gamma function, :math:`n` is the number of
    data points, :math:`n_k` is the number of measurements in bin :math:`k`
    [1]_.

    References
    ----------
    .. [1] Knuth, K.H. "Optimal Data-Based Binning for Histograms".
       arXiv:0605197, 2006

    See Also
    --------
    freedman_bin_width
    scott_bin_width
    bayesian_blocks
    histogram
    """
    # import here because of optional scipy dependency
    from scipy import optimize

    knuthF = _KnuthF(data)
    dx0, bins0 = freedman_bin_width(data, True)
    M = optimize.fmin(knuthF, len(bins0), disp=not quiet)[0]
    bins = knuthF.bins(M)
    dx = bins[1] - bins[0]

    if return_bins:
        return dx, bins
    else:
        return dx


class _KnuthF(object):
    r"""Class which implements the function minimized by knuth_bin_width

    Parameters
    ----------
    data : array-like, one dimension
        data to be histogrammed

    Notes
    -----
    the function F is given by

    .. math::
        F(M|x,I) = n\log(M) + \log\Gamma(\frac{M}{2})
        - M\log\Gamma(\frac{1}{2})
        - \log\Gamma(\frac{2n+M}{2})
        + \sum_{k=1}^M \log\Gamma(n_k + \frac{1}{2})

    where :math:`\Gamma` is the Gamma function, :math:`n` is the number of
    data points, :math:`n_k` is the number of measurements in bin :math:`k`.

    See Also
    --------
    knuth_bin_width
    """
    def __init__(self, data):
        self.data = np.array(data, copy=True)
        if self.data.ndim != 1:
            raise ValueError("data should be 1-dimensional")
        self.data.sort()
        self.n = self.data.size

        # import here rather than globally: scipy is an optional dependency.
        # Note that scipy is imported in the function which calls this,
        # so there shouldn't be any issue importing here.
        from scipy import special

        # create a reference to gammaln to use in self.eval()
        self.gammaln = special.gammaln

    def bins(self, M):
        """Return the bin edges given a width dx"""
        return np.linspace(self.data[0], self.data[-1], int(M) + 1)

    def __call__(self, M):
        return self.eval(M)

    def eval(self, M):
        """Evaluate the Knuth function

        Parameters
        ----------
        dx : float
            Width of bins

        Returns
        -------
        F : float
            evaluation of the negative Knuth likelihood function:
            smaller values indicate a better fit.
        """
        M = int(M)

        if M <= 0:
            return np.inf

        bins = self.bins(M)
        nk, bins = np.histogram(self.data, bins)

        return -(self.n * np.log(M) +
                 self.gammaln(0.5 * M) -
                 M * self.gammaln(0.5) -
                 self.gammaln(self.n + 0.5 * M) +
                 np.sum(self.gammaln(nk + 0.5)))
