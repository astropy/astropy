# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains simple statistical algorithms that are straightforwardly
implemented as a single python function (or family of functions).

This module should generally not be used directly.  Everything in `__all__` is
imported into `astropy.stats`, and hence that package should be used for
access.
"""
from __future__ import division

import numpy as np

__all__ = ['sigma_clip', 'binom_conf_interval', 'effhist']


def sigma_clip(data, sig=3, iters=1, cenfunc=np.median, varfunc=np.var,
               maout=False):
    """ Perform sigma-clipping on the provided data.

    This performs the sigma clipping algorithm - i.e. the data will be iterated
    over, each time rejecting points that are more than a specified number of
    standard deviations discrepant.

    .. note::
        `scipy.stats.sigmaclip` provides a subset of the functionality in this
        function.

    Parameters
    ----------
    data : array-like
        The data to be sigma-clipped (any shape).
    sig : float
        The number of standard deviations (*not* variances) to use as the
        clipping limit.
    iters : int or None
        The number of iterations to perform clipping for, or None to clip until
        convergence is achieved (i.e. continue until the last iteration clips
        nothing).
    cenfunc : callable
        The technique to compute the center for the clipping. Must be a
        callable that takes in a 1D data array and outputs the central value.
        Defaults to the median.
    varfunc : callable
        The technique to compute the variance about the center. Must be a
        callable that takes in a 1D data array and outputs the width estimator
        that will be interpreted as a variance. Defaults to the variance.
    maout : bool or 'copy'
        If True, a masked array will be returned. If the special string
        'inplace', the masked array will contain the same array as `data`,
        otherwise the array data will be copied.

    Returns
    -------
    filtereddata : `numpy.ndarray` or `numpy.masked.MaskedArray`
        If `maout` is True, this is a masked array with a shape matching the
        input that is masked where the algorithm has rejected those values.
        Otherwise, a 1D array of values including only those that are not
        clipped.
    mask : boolean array
        Only present if `maout` is False. A boolean array with a shape matching
        the input `data` that is False for rejected values and True for all
        others.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    only the points that are within 2 *sample* standard deviation from the
    median::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> data,mask = sigma_clip(randvar, 2, 1)

    This will clipping on a similar distribution, but for 3 sigma relative to
    the sample *mean*, will clip until converged, and produces a
    `numpy.masked.MaskedArray`::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> maskedarr = sigma_clip(randvar, 3, None, mean, maout=True)

    """

    data = np.array(data, copy=False)
    oldshape = data.shape
    data = data.ravel()

    mask = np.ones(data.size, bool)
    if iters is None:
        i = -1
        lastrej = sum(mask) + 1
        while(sum(mask) != lastrej):
            i += 1
            lastrej = sum(mask)
            do = data - cenfunc(data[mask])
            mask = do * do <= varfunc(data[mask]) * sig ** 2
        iters = i + 1
        #TODO: ?print iters to the log if iters was None?
    else:
        for i in range(iters):
            do = data - cenfunc(data[mask])
            mask = do * do <= varfunc(data[mask]) * sig ** 2

    if maout:
        return np.ma.MaskedArray(data, ~mask, copy=maout != 'inplace')
    else:
        return data[mask], mask.reshape(oldshape)


def binom_conf_interval(k, n, conf=0.68269):
    r"""Confidence interval on binomial probability given k successes,
    n trials.

    Parameters
    ----------
    k : int or numpy.ndarray
        Number of successes (0 <= `k` <= `n`).
    n : int or numpy.ndarray
        Number of trials (`n` > 0).
    conf : float, optional
        Desired probability content of interval. Default is 0.68269.

    Returns
    -------
    interval : numpy.ndarray
        `interval[0]` and `interval[1]` correspond to the lower and upper
        limits, respectively, for each element in `k`, `n`.

    Notes
    -----

    In situations where a probability of success is not known, it can be
    estimated from an observed number of trials and successes. For example,
    this is done in Monte Carlo experiments designed to estimate a detection
    efficiency. Given a *known* probability of success :math:`\epsilon`,
    and number of trials ``N``, the binomial distribution gives the
    probability of observing ``k`` successes:

    .. math::

        P(k | \epsilon, N) = \frac{N!}{k!(N-k)!} \epsilon^k (1-\epsilon)^{N-k}

    Using Bayes' theorem and a flat prior on :math:`\epsilon` gives
    the desired posterior probability distribtion for
    :math:`\epsilon`, given ``k`` observed successes out of ``N``
    trials:

    .. math::

        P(\epsilon|k,N) = \frac{\Gamma(N+2)}{\Gamma(k+1)\Gamma(N-k+1)}
        \epsilon^k (1-\epsilon)^{N-k}

    This probability density function peaks at :math:`\epsilon` = ``k / N``,
    as expected for the best estimate of :math:`\epsilon`. The
    PDF fully describes one's knowledge about :math:`\epsilon`, but it
    is often desirable to calculate an interval that contains a given
    fraction :math:`\lambda` of the total probability content (for
    example, 68.3%).

    There is not a unique interval that satisfies this criterion.
    One choice is an interval that contains :math:`\lambda` of the
    total probability content below the peak of the distribution and
    :math:`\lambda` of the total probability content above the peak.
    This has the advantage of being easy to calculate with existing tools,
    and is what this function does. (Another well-motivated choice is
    the *shortest* interval that contains :math:`\lambda` of the total
    probability content. However, this is more difficult to implement
    and it should not differ significantly from the above interval.)

    Examples
    --------

    >>> binom_conf_interval(4, 5)
    array([ 0.5831211 ,  0.90251087])

    >>> binom_conf_interval([4, 2, 0, 5], 5)
    array([[ 0.5831211 ,  0.23375493,  0.        ,  0.82587432],
           [ 0.90251087,  0.60473566,  0.17412568,  1.        ]])

    """

    from scipy.special import betainc, betaincinv

    if conf < 0. or conf > 1.:
        raise ValueError('conf must be between 0. and 1.')
    tails = 1. - conf

    k = np.array(k).astype(np.int)
    n = np.array(n).astype(np.int)
    if (n <= 0).any():
        raise ValueError('n must be positive')
    if (k < 0).any() or (k > n).any():
        raise ValueError('k must be in {0, 1, .., k}')

    p_peak = k / n  # Peak of PDF.

    # Probability content below and above the peak.
    p_below = betainc(k + 1, n - k + 1, p_peak)
    p_above = 1. - p_below

    interval0 = betaincinv(k + 1, n - k + 1, tails * p_below)
    interval1 = betaincinv(k + 1, n - k + 1, 1. - tails * p_above)

    return np.array([interval0, interval1])


def effhist(x, success, bins=10, range=None, conf=0.68269):
    """Estimate success probability (e.g., efficiency) and uncertainty in
    bins in `x`.

    Parameters
    ----------
    x : list_like
        Values
    success : list_like (bool)
        Success (True) or failure (False) corresponding to items in 'x'
    bins : int or sequence of scalars, optional
        If bins is an int, it defines the number of equal-width bins
        in the given range (10, by default). If bins is a sequence, it
        defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths.
    range : (float, float), optional
        The lower and upper range of the bins. If not provided, range
        is simply (a.min(), a.max()). Values outside the range are
        ignored.
    conf : float, optional
        Desired probability content contained in confidence interval
        (p - perr, p + perr). Default is 0.68269.

    Returns
    -------
    bin_ctr : numpy.ndarray
        Central value of bins. Bins without any entries are not returned.
    bin_halfwidth : numpy.ndarray
        Half-width of each bin such that `bin_ctr - bin_halfwidth` and
        `bin_ctr + bins_halfwidth` give the left and right side of each bin,
        respectively.
    p : numpy.ndarray
        Efficiency in each bin.
    perr : numpy.ndarray
        2-d array of shape (2, len(p)) representing the upper and lower
        uncertainty on p in each bin.

    See Also
    --------
    binom_conf_interval : Function used to estimate confidence interval in
                          each bin.


    Examples
    --------

    The true probability is 0.9 for x < 5 and 0.1 at x >= 5. Estimate
    the probability as a function of x based on 100 random trials spread
    over the range 0 < x < 10.

    >>> x = 10. * np.random.rand(100)
    >>> success = np.zeros(len(x), dtype=np.bool)
    >>> success[x < 5.] = np.random.rand((x < 5.).sum()) > 0.1
    >>> success[x >= 5.] = np.random.rand((x >= 5.).sum()) > 0.9
    >>> bin_ctr, bin_hw, p, perr = effhist(x, success, bins=20)

    Plot it.

    >>> import matplotlib.pyplot as plt
    >>> plt.errorbar(bin_ctr, p, xerr=bin_hw, yerr=perr, ls='none', marker='o',
    ...              label='estimate')

    .. image:: ../stats/images/effhist.png

    """

    x = np.ravel(x)
    success = np.ravel(success).astype(np.bool)
    if x.shape != success.shape:
        raise ValueError('Length of x and success must match')

    # Put all values into a histogram.
    hist_x, bin_edges = np.histogram(x, bins=bins, range=range)

    # Put only values where success = True into the *same* histogram
    hist_success, bin_edges = np.histogram(x[success], bins=bin_edges)

    # Calculate bin centers and half-widths.
    bin_ctr = (bin_edges[:-1] + bin_edges[1:]) / 2.
    bin_halfwidth = bins - bin_edges[:-1]

    # Limit to bins with more than zero entries.
    valid = hist_x > 0
    bin_ctr = bin_ctr[valid]
    bin_halfwidth = bin_halfwidth[valid]
    hist_x = hist_x[valid]
    hist_success = hist_success[valid]

    p = hist_success.astype(np.float) / hist_x.astype(np.float)
    interval = binom_conf_interval(hist_success, hist_x, conf=conf)
    perr = np.abs(interval - p)

    return bin_ctr, bin_halfwidth, p, perr
