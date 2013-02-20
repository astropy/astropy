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

__all__ = ['sigma_clip', 'binom_conf_interval', 'binned_efficiency']


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


#TODO Note scipy dependency
def binom_conf_interval(k, n, conf=0.68269, interval='wilson'):
    r"""Binomial proportion confidence interval given k successes,
    n trials.

    Parameters
    ----------
    k : int or numpy.ndarray
        Number of successes (0 <= `k` <= `n`).
    n : int or numpy.ndarray
        Number of trials (`n` > 0).
    conf : float in [0, 1], optional
        Desired probability content of interval. Default is 0.68269.
    interval : {'wilson', 'jeffreys'}, optional
        Formula used for confidence interval (see notes for details).
        Default is 'wilson'.

    Returns
    -------
    conf_interval : numpy.ndarray
        `conf_interval[0]` and `conf_interval[1]` correspond to the lower
        and upper limits, respectively, for each element in `k`, `n`.

    Notes
    -----
    In situations where a probability of success is not known, it can
    be estimated from a number of trials (N) and number of
    observed successes (k). For example, this is done in Monte
    Carlo experiments designed to estimate a detection efficiency. It
    is simple to take the sample proportion of successes (k/N)
    as a reasonable best estimate of the true probability
    :math:`\epsilon`. However, deriving an accurate confidence
    interval on :math:`\epsilon` is non-trivial. There are several
    formulas for this interval (see [1]_). Two intervals are implemented
    here:

    **1. The Wilson Interval.** This interval, attributed to Wilson [2]_,
    is given by

    .. math::

        CI_W = \frac{k + \kappa^2/2}{N + \kappa^2}
        \pm \frac{\kappa n^{1/2}}{n + \kappa^2}
        ((\hat{\epsilon}(1 - \hat{\epsilon}) + \kappa^2/(4n))^{1/2}

    where :math:`\hat{\epsilon} = k / N` and :math:`\kappa` is the
    number of standard deviations corresponding to the desired
    confidence interval for a *normal* distribution (for example,
    1.0 for a confidence interval of 68.269%). For a
    confidence interval of 100(1 - :math:`\alpha`)%,

    .. math::

        \kappa = \Phi^{-1}(1-\alpha/2) = \sqrt{2}{\rm erf}^{-1}(1-\alpha).

    **2. The Jeffreys Interval.** This interval is derived by applying
    Bayes' theorem to the binomial distribution with the
    noninformative Jeffreys prior [3]_, [4]_. The noninformative Jeffreys
    prior is the Beta distribution, Beta(1/2, 1/2), which has the density
    function

    .. math::

        f(\epsilon) = \pi^{-1} \epsilon^{-1/2}(1-\epsilon)^{-1/2}.

    The posterior density function is also a Beta distribution: Beta(k
    + 1/2, N - k + 1/2). The interval is then chosen so that it is
    *equal-tailed*: Each tail (outside the interval) contains
    :math:`\alpha`/2 of the posterior probability, and the interval
    itself contains 1 - :math:`\alpha`. This interval must be
    calculated numerically. Additionally, when k = 0 the lower limit
    is set to 0 and when k = N the upper limit is set to 1, so that in
    these cases, there is only one tail containing :math:`\alpha`/2
    and the interval itself contains 1 - :math:`\alpha`/2 rather than
    the nominal 1 - :math:`\alpha`.

    References
    ----------
    .. [1] Brown, Lawrence D.; Cai, T. Tony; DasGupta, Anirban (2001).
       "Interval Estimation for a Binomial Proportion". Statistical
       Science 16 (2): 101-133. doi:10.1214/ss/1009213286

    .. [2] Wilson, E. B. (1927). "Probable inference, the law of
       succession, and statistical inference". Journal of the American
       Statistical Association 22: 209-212.

    .. [3] Jeffreys, Harold (1946). "An Invariant Form for the Prior
       Probability in Estimation Problems". Proc. R. Soc. Lond.. A 24 186
       (1007): 453-461. doi:10.1098/rspa.1946.0056

    .. [4] Jeffreys, Harold (1998). Theory of Probability. Oxford
       University Press, 3rd edition. ISBN 978-0198503682

    Examples
    --------
    Integer inputs return an array with shape (2,):

    >>> binom_conf_interval(4, 5, interval='wilson')
    array([ 0.57921724,  0.92078259])

    Arrays of arbitrary dimension are supported. The Wilson and Jeffreys
    intervals give similar results, even for small k, N:

    >>> binom_conf_interval([0, 1, 2, 5], 5, interval='wilson')
    array([[ 0.        ,  0.07921741,  0.21597328,  0.83333304],
           [ 0.16666696,  0.42078276,  0.61736012,  1.        ]])

    >>> binom_conf_interval([0, 1, 2, 5], 5, interval='jeffreys')
    array([[ 0.        ,  0.0842525 ,  0.21789949,  0.82788246],
           [ 0.17211754,  0.42218001,  0.61753691,  1.        ]])
    """

    # Check conf.
    if conf < 0. or conf > 1.:
        raise ValueError('conf must be between 0. and 1.')
    alpha = 1. - conf

    # Check k and n.
    k = np.asarray(k).astype(np.int)
    n = np.asarray(n).astype(np.int)
    if (n <= 0).any():
        raise ValueError('n must be positive')
    if (k < 0).any() or (k > n).any():
        raise ValueError('k must be in {0, 1, .., n}')

    if interval == 'wilson':
        from scipy.special import erfinv

        kappa = np.sqrt(2.) * erfinv(conf)
        if kappa > 1.e10: kappa = 1.e10  # Avoid floating point overflows.
        k = k.astype(np.float)
        n = n.astype(np.float)
        p = k / n
        midpoint = (k + kappa ** 2 / 2.) / (n + kappa ** 2)
        halflength = (kappa * np.sqrt(n)) / (n + kappa ** 2) * \
            np.sqrt(p * (1 - p) + kappa ** 2 / (4 * n))
        conf_interval = np.array([midpoint - halflength,
                                  midpoint + halflength])

        # Interval can be out of range [0, 1]
        # due to floating point errors. Correct this.
        conf_interval[conf_interval < 0.] = 0.
        conf_interval[conf_interval > 1.] = 1.

        return conf_interval

    elif interval == 'jeffreys':
        from scipy.special import betainc, betaincinv

        lowerbound = betaincinv(k + 0.5, n - k + 0.5, alpha / 2.)
        upperbound = betaincinv(k + 0.5, n - k + 0.5, 1. - alpha / 2.)
    
        # Set lower or upper bound to k/n when k/n = 0 or 1. 
        lowerbound[k == 0] = 0.
        upperbound[k == n] = 1.

        return np.array([lowerbound, upperbound])
    
    else:
        raise ValueError('Unrecognized interval: {0:s}'.format(interval))


#TODO Note scipy dependency (needed in binom_conf_interval)
def binned_efficiency(x, success, bins=10, range=None, conf=0.68269,
                      interval='wilson'):
    """Binomial proportion (efficiency) and confidence interval in bins
    in `x`.

    For a set of success/failure trials, each corresponding to some
    value of the parameter `x`, place the trials into bins and return
    an estimate of the efficiency (number of successes / number of
    failures) and its uncertainty in each bin.

    Parameters
    ----------
    x : list_like
        Values.
    success : list_like (bool)
        Success (True) or failure (False) corresponding to each value
        in `x`.  Must be same length as `x`.
    bins : int or sequence of scalars, optional
        If bins is an int, it defines the number of equal-width bins
        in the given range (10, by default). If bins is a sequence, it
        defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths (in this case, 'range' is ignored).
    range : (float, float), optional
        The lower and upper range of the bins. If `None` (default),
        the range is set to (x.min(), x.max()). Values outside the
        range are ignored.
    conf : float in [0, 1], optional
        Desired probability content in the confidence
        interval (p - perr[0], p + perr[1]) in each bin. Default is
        0.68269.
    interval : {'wilson', 'jeffreys'}, optional
        Formula used to calculate confidence interval on the
        efficiency in each bin. The two intervals give similar
        results.  See `binom_conf_interval` for definition of the
        intervals. Default is 'wilson'.

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
    Suppose we wish to estimate the detection efficiency of an astronomical
    source as a function of magnitude using a Monte Carlo experiment with
    100 trials between magnitude 20. and 30. The true 50% detection magnitude
    is 25.
    
    >>> from scipy.special import erf
    >>> from scipy.stats.distributions import binom
    >>> def true_efficiency(x):
    ...     return 0.5 - 0.5 * erf((x - 25.) / 2.)
    >>> mag = 20. + 10. * np.random.rand(100)
    >>> detected = binom.rvs(1, true_efficiency(mag))
    >>> bins, binshw, p, perr = binned_efficiency(mag, detected, bins=20)
    >>> plt.errorbar(bins, p, xerr=binshw, yerr=perr, ls='none', marker='o',
    ...              label='estimate')

    .. plot::

       import numpy as np
       from scipy.special import erf
       from scipy.stats.distributions import binom
       import matplotlib.pyplot as plt
       from astropy.stats import binned_efficiency
       def true_efficiency(x):
           return 0.5 - 0.5 * erf((x - 25.) / 2.)
       mag = 20. + 10. * np.random.rand(100)
       detected = binom.rvs(1, true_efficiency(mag))
       bins, binshw, p, perr = binned_efficiency(mag, detected, bins=20)
       plt.errorbar(bins, p, xerr=binshw, yerr=perr, ls='none', marker='o',
                    label='estimate')
       X = np.linspace(20., 30., 1000)
       plt.plot(X, true_efficiency(X), ls='-', color='r',
                label='true efficiency')
       plt.ylim(0., 1.)
       plt.xlabel('Magnitude')
       plt.ylabel('Detection efficiency')
       plt.legend()
       plt.show()

    """

    x = np.ravel(x)
    success = np.ravel(success).astype(np.bool)
    if x.shape != success.shape:
        raise ValueError('sizes of x and success must match')

    # Put values into a histogram (`n`). Put "successful" values
    # into a second histogram (`k`) with identical binning.
    n, bin_edges = np.histogram(x, bins=bins, range=range)
    k, bin_edges = np.histogram(x[success], bins=bin_edges)
    bin_ctr = (bin_edges[:-1] + bin_edges[1:]) / 2.
    bin_halfwidth = bin_ctr - bin_edges[:-1]

    # Remove bins with zero entries.
    valid = n > 0
    bin_ctr = bin_ctr[valid]
    bin_halfwidth = bin_halfwidth[valid]
    n = n[valid]
    k = k[valid]

    p = k / n
    bounds = binom_conf_interval(k, n, conf=conf, interval=interval)
    perr = np.abs(bounds - p)

    return bin_ctr, bin_halfwidth, p, perr
