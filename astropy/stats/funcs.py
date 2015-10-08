# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains simple statistical algorithms that are straightforwardly
implemented as a single python function (or family of functions).

This module should generally not be used directly.  Everything in `__all__` is
imported into `astropy.stats`, and hence that package should be used for
access.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ..extern.six.moves import xrange


__all__ = ['binom_conf_interval', 'binned_binom_proportion',
           'poisson_conf_interval',
           'median_absolute_deviation', 'biweight_location',
           'biweight_midvariance', 'signal_to_noise_oir_ccd', 'bootstrap',
           'mad_std', 'gaussian_fwhm_to_sigma', 'gaussian_sigma_to_fwhm']

__doctest_skip__ = ['binned_binom_proportion']
__doctest_requires__ = {'binom_conf_interval': ['scipy.special'],
                        'poisson_conf_interval': ['scipy.special']}


gaussian_sigma_to_fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))
"""
Factor with which to multiply Gaussian 1-sigma standard deviation(s) to
convert them to full width at half maximum(s).
"""

gaussian_fwhm_to_sigma = 1. / gaussian_sigma_to_fwhm
"""
Factor with which to multiply Gaussian full width at half maximum(s) to
convert them to 1-sigma standard deviation(s).
"""


# TODO Note scipy dependency
def binom_conf_interval(k, n, conf=0.68269, interval='wilson'):
    r"""Binomial proportion confidence interval given k successes,
    n trials.

    Parameters
    ----------
    k : int or numpy.ndarray
        Number of successes (0 <= ``k`` <= ``n``).
    n : int or numpy.ndarray
        Number of trials (``n`` > 0).  If both ``k`` and ``n`` are arrays,
        they must have the same shape.
    conf : float in [0, 1], optional
        Desired probability content of interval. Default is 0.68269,
        corresponding to 1 sigma in a 1-dimensional Gaussian distribution.
    interval : {'wilson', 'jeffreys', 'flat', 'wald'}, optional
        Formula used for confidence interval. See notes for details.  The
        ``'wilson'`` and ``'jeffreys'`` intervals generally give similar
        results, while 'flat' is somewhat different, especially for small
        values of ``n``.  ``'wilson'`` should be somewhat faster than
        ``'flat'`` or ``'jeffreys'``.  The 'wald' interval is generally not
        recommended.  It is provided for comparison purposes.  Default is
        ``'wilson'``.

    Returns
    -------
    conf_interval : numpy.ndarray
        ``conf_interval[0]`` and ``conf_interval[1]`` correspond to the lower
        and upper limits, respectively, for each element in ``k``, ``n``.

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
    formulas for this interval (see [1]_). Four intervals are implemented
    here:

    **1. The Wilson Interval.** This interval, attributed to Wilson [2]_,
    is given by

    .. math::

        CI_{\rm Wilson} = \frac{k + \kappa^2/2}{N + \kappa^2}
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

    The justification for this prior is that it is invariant under
    reparameterizations of the binomial proportion.
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

    **3. A Flat prior.** This is similar to the Jeffreys interval,
    but uses a flat (uniform) prior on the binomial proportion
    over the range 0 to 1 rather than the reparametrization-invariant
    Jeffreys prior.  The posterior density function is a Beta distribution:
    Beta(k + 1, N - k + 1).  The same comments about the nature of the
    interval (equal-tailed, etc.) also apply to this option.

    **4. The Wald Interval.** This interval is given by

    .. math::

       CI_{\rm Wald} = \hat{\epsilon} \pm
       \kappa \sqrt{\frac{\hat{\epsilon}(1-\hat{\epsilon})}{N}}

    The Wald interval gives acceptable results in some limiting
    cases. Particularly, when N is very large, and the true proportion
    :math:`\epsilon` is not "too close" to 0 or 1. However, as the
    later is not verifiable when trying to estimate :math:`\epsilon`,
    this is not very helpful. Its use is not recommended, but it is
    provided here for comparison purposes due to its prevalence in
    everyday practical statistics.

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

    >>> binom_conf_interval([0, 1, 2, 5], 5, interval='flat')
    array([[ 0.        ,  0.12139799,  0.24309021,  0.73577037],
           [ 0.26422963,  0.45401727,  0.61535699,  1.        ]])

    In contrast, the Wald interval gives poor results for small k, N.
    For k = 0 or k = N, the interval always has zero length.

    >>> binom_conf_interval([0, 1, 2, 5], 5, interval='wald')
    array([[ 0.        ,  0.02111437,  0.18091075,  1.        ],
           [ 0.        ,  0.37888563,  0.61908925,  1.        ]])

    For confidence intervals approaching 1, the Wald interval for
    0 < k < N can give intervals that extend outside [0, 1]:

    >>> binom_conf_interval([0, 1, 2, 5], 5, interval='wald', conf=0.99)
    array([[ 0.        , -0.26077835, -0.16433593,  1.        ],
           [ 0.        ,  0.66077835,  0.96433593,  1.        ]])

    """

    if conf < 0. or conf > 1.:
        raise ValueError('conf must be between 0. and 1.')
    alpha = 1. - conf

    k = np.asarray(k).astype(np.int)
    n = np.asarray(n).astype(np.int)

    if (n <= 0).any():
        raise ValueError('n must be positive')
    if (k < 0).any() or (k > n).any():
        raise ValueError('k must be in {0, 1, .., n}')

    if interval == 'wilson' or interval == 'wald':
        from scipy.special import erfinv
        kappa = np.sqrt(2.) * min(erfinv(conf), 1.e10)  # Avoid overflows.
        k = k.astype(np.float)
        n = n.astype(np.float)
        p = k / n

        if interval == 'wilson':
            midpoint = (k + kappa ** 2 / 2.) / (n + kappa ** 2)
            halflength = (kappa * np.sqrt(n)) / (n + kappa ** 2) * \
                np.sqrt(p * (1 - p) + kappa ** 2 / (4 * n))
            conf_interval = np.array([midpoint - halflength,
                                      midpoint + halflength])

            # Correct intervals out of range due to floating point errors.
            conf_interval[conf_interval < 0.] = 0.
            conf_interval[conf_interval > 1.] = 1.
        else:
            midpoint = p
            halflength = kappa * np.sqrt(p * (1. - p) / n)
            conf_interval = np.array([midpoint - halflength,
                                      midpoint + halflength])

    elif interval == 'jeffreys' or interval == 'flat':
        from scipy.special import betaincinv

        if interval == 'jeffreys':
            lowerbound = betaincinv(k + 0.5, n - k + 0.5, 0.5 * alpha)
            upperbound = betaincinv(k + 0.5, n - k + 0.5, 1. - 0.5 * alpha)
        else:
            lowerbound = betaincinv(k + 1, n - k + 1, 0.5 * alpha)
            upperbound = betaincinv(k + 1, n - k + 1, 1. - 0.5 * alpha)

        # Set lower or upper bound to k/n when k/n = 0 or 1
        #  We have to treat the special case of k/n being scalars,
        #  which is an ugly kludge
        if lowerbound.ndim == 0:
            if k == 0:
                lowerbound = 0.
            elif k == n:
                upperbound = 1.
        else:
            lowerbound[k == 0] = 0
            upperbound[k == n] = 1

        conf_interval = np.array([lowerbound, upperbound])
    else:
        raise ValueError('Unrecognized interval: {0:s}'.format(interval))

    return conf_interval


# TODO Note scipy dependency (needed in binom_conf_interval)
def binned_binom_proportion(x, success, bins=10, range=None, conf=0.68269,
                            interval='wilson'):
    """Binomial proportion and confidence interval in bins of a continuous
    variable ``x``.

    Given a set of datapoint pairs where the ``x`` values are
    continuously distributed and the ``success`` values are binomial
    ("success / failure" or "true / false"), place the pairs into
    bins according to ``x`` value and calculate the binomial proportion
    (fraction of successes) and confidence interval in each bin.

    Parameters
    ----------
    x : list_like
        Values.
    success : list_like (bool)
        Success (`True`) or failure (`False`) corresponding to each value
        in ``x``.  Must be same length as ``x``.
    bins : int or sequence of scalars, optional
        If bins is an int, it defines the number of equal-width bins
        in the given range (10, by default). If bins is a sequence, it
        defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths (in this case, 'range' is ignored).
    range : (float, float), optional
        The lower and upper range of the bins. If `None` (default),
        the range is set to ``(x.min(), x.max())``. Values outside the
        range are ignored.
    conf : float in [0, 1], optional
        Desired probability content in the confidence
        interval ``(p - perr[0], p + perr[1])`` in each bin. Default is
        0.68269.
    interval : {'wilson', 'jeffreys', 'flat', 'wald'}, optional
        Formula used to calculate confidence interval on the
        binomial proportion in each bin. See `binom_conf_interval` for
        definition of the intervals.  The 'wilson', 'jeffreys',
        and 'flat' intervals generally give similar results.  'wilson'
        should be somewhat faster, while 'jeffreys' and 'flat' are
        marginally superior, but differ in the assumed prior.
        The 'wald' interval is generally not recommended.
        It is provided for comparison purposes. Default is 'wilson'.

    Returns
    -------
    bin_ctr : numpy.ndarray
        Central value of bins. Bins without any entries are not returned.
    bin_halfwidth : numpy.ndarray
        Half-width of each bin such that ``bin_ctr - bin_halfwidth`` and
        ``bin_ctr + bins_halfwidth`` give the left and right side of each bin,
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
    Suppose we wish to estimate the efficiency of a survey in
    detecting astronomical sources as a function of magnitude (i.e.,
    the probability of detecting a source given its magnitude). In a
    realistic case, we might prepare a large number of sources with
    randomly selected magnitudes, inject them into simulated images,
    and then record which were detected at the end of the reduction
    pipeline. As a toy example, we generate 100 data points with
    randomly selected magnitudes between 20 and 30 and "observe" them
    with a known detection function (here, the error function, with
    50% detection probability at magnitude 25):

    >>> from scipy.special import erf
    >>> from scipy.stats.distributions import binom
    >>> def true_efficiency(x):
    ...     return 0.5 - 0.5 * erf((x - 25.) / 2.)
    >>> mag = 20. + 10. * np.random.rand(100)
    >>> detected = binom.rvs(1, true_efficiency(mag))
    >>> bins, binshw, p, perr = binned_binom_proportion(mag, detected, bins=20)
    >>> plt.errorbar(bins, p, xerr=binshw, yerr=perr, ls='none', marker='o',
    ...              label='estimate')

    .. plot::

       import numpy as np
       from scipy.special import erf
       from scipy.stats.distributions import binom
       import matplotlib.pyplot as plt
       from astropy.stats import binned_binom_proportion
       def true_efficiency(x):
           return 0.5 - 0.5 * erf((x - 25.) / 2.)
       np.random.seed(400)
       mag = 20. + 10. * np.random.rand(100)
       np.random.seed(600)
       detected = binom.rvs(1, true_efficiency(mag))
       bins, binshw, p, perr = binned_binom_proportion(mag, detected, bins=20)
       plt.errorbar(bins, p, xerr=binshw, yerr=perr, ls='none', marker='o',
                    label='estimate')
       X = np.linspace(20., 30., 1000)
       plt.plot(X, true_efficiency(X), label='true efficiency')
       plt.ylim(0., 1.)
       plt.title('Detection efficiency vs magnitude')
       plt.xlabel('Magnitude')
       plt.ylabel('Detection efficiency')
       plt.legend()
       plt.show()

    The above example uses the Wilson confidence interval to calculate
    the uncertainty ``perr`` in each bin (see the definition of various
    confidence intervals in `binom_conf_interval`). A commonly used
    alternative is the Wald interval. However, the Wald interval can
    give nonsensical uncertainties when the efficiency is near 0 or 1,
    and is therefore **not** recommended. As an illustration, the
    following example shows the same data as above but uses the Wald
    interval rather than the Wilson interval to calculate ``perr``:

    >>> bins, binshw, p, perr = binned_binom_proportion(mag, detected, bins=20,
    ...                                                 interval='wald')
    >>> plt.errorbar(bins, p, xerr=binshw, yerr=perr, ls='none', marker='o',
    ...              label='estimate')

    .. plot::

       import numpy as np
       from scipy.special import erf
       from scipy.stats.distributions import binom
       import matplotlib.pyplot as plt
       from astropy.stats import binned_binom_proportion
       def true_efficiency(x):
           return 0.5 - 0.5 * erf((x - 25.) / 2.)
       np.random.seed(400)
       mag = 20. + 10. * np.random.rand(100)
       np.random.seed(600)
       detected = binom.rvs(1, true_efficiency(mag))
       bins, binshw, p, perr = binned_binom_proportion(mag, detected, bins=20,
                                                       interval='wald')
       plt.errorbar(bins, p, xerr=binshw, yerr=perr, ls='none', marker='o',
                    label='estimate')
       X = np.linspace(20., 30., 1000)
       plt.plot(X, true_efficiency(X), label='true efficiency')
       plt.ylim(0., 1.)
       plt.title('The Wald interval can give nonsensical uncertainties')
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


def poisson_conf_interval(n, interval='root-n', sigma=1):
    r"""Poisson parameter confidence interval given observed counts

    Parameters
    ----------
    n : int or numpy.ndarray
        Number of counts (0 <= ``n``).
    interval : {'root-n','root-n-0','pearson','sherpagehrels','frequentist-confidence'}, optional
        Formula used for confidence interval. See notes for details.
        Default is ``'root-n'``.
    sigma : float
        Number of sigma for confidence interval; only supported for
        the 'frequentist-confidence' mode.


    Returns
    -------
    conf_interval : numpy.ndarray
        ``conf_interval[0]`` and ``conf_interval[1]`` correspond to the lower
        and upper limits, respectively, for each element in ``n``.

    Notes
    -----

    The "right" confidence interval to use for Poisson data is a
    matter of debate. The CDF working group [recommends][pois_eb]
    using root-n throughout, largely in the interest of
    comprehensibility, but discusses other possibilities. The ATLAS
    group also [discusses][ErrorBars] several possibilities but
    concludes that no single representation is suitable for all cases.
    The suggestion has also been [floated][ac12] that error bars should be
    attached to theoretical predictions instead of observed data,
    which this function will not help with (but it's easy; then you
    really should use the square root of the theoretical prediction).

    The intervals implemented here are:

    **1. 'root-n'** This is a very widely used standard rule derived
    from the maximum-likelihood estimator for the mean of the Poisson
    process. While it produces questionable results for small n and
    outright wrong results for n=0, it is standard enough that people are
    (supposedly) used to interpreting these wonky values. The interval is

    .. math::

        CI = (n-\sqrt{n}, n+\sqrt{n})

    **2. 'root-n-0'** This is identical to the above except that where
    n is zero the interval returned is (0,1).

    **3. 'pearson'** This is an only-slightly-more-complicated rule
    based on Pearson's chi-squared rule (as [explained][pois_eb] by
    the CDF working group). It also has the nice feature that
    if your theory curve touches an endpoint of the interval, then your
    data point is indeed one sigma away. The interval is

    .. math::

        CI = (n+0.5-\sqrt{n+0.25}, n+0.5+\sqrt{n+0.25})

    **4. 'sherpagehrels'** This rule is used by default in the fitting
    package 'sherpa'. The [documentation][sherpa_gehrels] claims it is
    based on a numerical approximation published in
    [Gehrels 1986][gehrels86] but it does not actually appear there.
    It is symmetrical, and while the upper limits
    are within about 1% of those given by 'frequentist-confidence', the
    lower limits can be badly wrong. The interval is

    .. math::

        CI = (n-1-\sqrt{n+0.75}, n+1+\sqrt{n+0.75})

    **5. 'frequentist-confidence'** These are frequentist central
    confidence intervals:

    .. math::

        CI = (0.5 F_{\chi^2}^{-1}(\alpha;2n),
              0.5 F_{\chi^2}^{-1}(1-\alpha;2(n+1)))

    where :math:`F_{\chi^2}^{-1}` is the quantile of the chi-square
    distribution with the indicated number of degrees of freedom and
    :math:\alpha is the one-tailed probability of the normal
    distribution (at the point given by the parameter 'sigma'). See
    [Maxwell 2011][maxw11] for further details.

    Examples
    --------

    >>> poisson_conf_interval(np.arange(10), interval='root-n').T
    array([[  0.        ,   0.        ],
           [  0.        ,   2.        ],
           [  0.58578644,   3.41421356],
           [  1.26794919,   4.73205081],
           [  2.        ,   6.        ],
           [  2.76393202,   7.23606798],
           [  3.55051026,   8.44948974],
           [  4.35424869,   9.64575131],
           [  5.17157288,  10.82842712],
           [  6.        ,  12.        ]])

    >>> poisson_conf_interval(np.arange(10), interval='root-n-0').T
    array([[  0.        ,   1.        ],
           [  0.        ,   2.        ],
           [  0.58578644,   3.41421356],
           [  1.26794919,   4.73205081],
           [  2.        ,   6.        ],
           [  2.76393202,   7.23606798],
           [  3.55051026,   8.44948974],
           [  4.35424869,   9.64575131],
           [  5.17157288,  10.82842712],
           [  6.        ,  12.        ]])

    >>> poisson_conf_interval(np.arange(10), interval='pearson').T
    array([[  0.        ,   1.        ],
           [  0.38196601,   2.61803399],
           [  1.        ,   4.        ],
           [  1.69722436,   5.30277564],
           [  2.43844719,   6.56155281],
           [  3.20871215,   7.79128785],
           [  4.        ,   9.        ],
           [  4.8074176 ,  10.1925824 ],
           [  5.62771868,  11.37228132],
           [  6.45861873,  12.54138127]])

    >>> poisson_conf_interval(np.arange(10),
    ...                       interval='frequentist-confidence').T
    array([[  0.        ,   1.84102165],
           [  0.17275378,   3.29952656],
           [  0.70818544,   4.63785962],
           [  1.36729531,   5.91818583],
           [  2.08566081,   7.16275317],
           [  2.84030886,   8.38247265],
           [  3.62006862,   9.58364155],
           [  4.41852954,  10.77028072],
           [  5.23161394,  11.94514152],
           [  6.05653896,  13.11020414]])

    >>> poisson_conf_interval(7,
    ...                       interval='frequentist-confidence').T
    array([  4.41852954,  10.77028072])

    [pois_eb]: http://www-cdf.fnal.gov/physics/statistics/notes/pois_eb.txt

    [ErrorBars]: http://www.pp.rhul.ac.uk/~cowan/atlas/ErrorBars.pdf

    [ac12]: http://adsabs.harvard.edu/abs/2012EPJP..127...24A

    [maxw11]: http://adsabs.harvard.edu/abs/2011arXiv1102.0822M

    [gehrels86]: http://adsabs.harvard.edu/abs/1986ApJ...303..336G

    [sherpa_gehrels]: http://cxc.harvard.edu/sherpa4.4/statistics/#chigehrels

    """

    if not np.isscalar(n):
        n = np.asanyarray(n)

    if interval == 'root-n':
        if sigma!=1:
            raise ValueError("Only sigma=1 supported for interval %s" % interval)
        conf_interval = np.array([n-np.sqrt(n),
                                  n+np.sqrt(n)])
    elif interval == 'root-n-0':
        if sigma!=1:
            raise ValueError("Only sigma=1 supported for interval %s" % interval)
        conf_interval = np.array([n-np.sqrt(n),
                                  n+np.sqrt(n)])
        if np.isscalar(n):
            if n==0:
                conf_interval[1] = 1
        else:
            conf_interval[1,n==0] = 1
    elif interval == 'pearson':
        if sigma!=1:
            raise ValueError("Only sigma=1 supported for interval %s" % interval)
        conf_interval = np.array([n+0.5-np.sqrt(n+0.25),
                                  n+0.5+np.sqrt(n+0.25)])
    elif interval == 'sherpagehrels':
        if sigma!=1:
            raise ValueError("Only sigma=1 supported for interval %s" % interval)
        conf_interval = np.array([n-1-np.sqrt(n+0.75),
                                  n+1+np.sqrt(n+0.75)])
    elif interval == 'frequentist-confidence':
        import scipy.stats
        alpha = scipy.stats.norm.sf(sigma)
        conf_interval = np.array([0.5*scipy.stats.chi2(2*n).ppf(alpha),
                                  0.5*scipy.stats.chi2(2*n+2).isf(alpha)])
        if np.isscalar(n):
            if n==0:
                conf_interval[0] = 0
        else:
            conf_interval[0,n==0] = 0
    else:
        raise ValueError("Invalid method for Poisson confidence intervals: %s" % interval)
    return conf_interval

def median_absolute_deviation(a, axis=None):
    """Compute the median absolute deviation.

    Returns the median absolute deviation (MAD) of the array elements.
    The MAD is defined as ``median(abs(a - median(a)))``.

    Parameters
    ----------
    a : array-like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (axis=None)
        is to compute the median along a flattened version of the array.

    Returns
    -------
    median_absolute_deviation : ndarray
        A new array holding the result. If the input contains
        integers, or floats of smaller precision than 64, then the output
        data-type is float64.  Otherwise, the output data-type is the same
        as that of the input.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the median absolute deviation for that distribution::

        >>> from astropy.stats import median_absolute_deviation
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> mad = median_absolute_deviation(randvar)

    See Also
    --------
    numpy.median

    """

    a = np.array(a, copy=False)
    a_median = np.median(a, axis=axis)

    # re-broadcast the output median array to subtract it
    if axis is not None:
        a_median = np.expand_dims(a_median, axis=axis)

    # calculated the median average deviation
    return np.median(np.abs(a - a_median), axis=axis)


def biweight_location(a, c=6.0, M=None):
    """Compute the biweight location for an array.

    Returns the biweight location for the array elements.
    The biweight is a robust statistic for determining the central
    location of a distribution.

    The biweight location is given by the following equation

    .. math::

        C_{bl}= M+\\frac{\Sigma_{\|u_i\|<1} (x_i-M)(1-u_i^2)^2}
        {\Sigma_{\|u_i\|<1} (1-u_i^2)^2}

    where M is the sample mean or if run iterative the initial guess,
    and u_i is given by

    .. math::

      u_{i} = \\frac{(x_i-M)}{cMAD}

    where MAD is the median absolute deviation.

    For more details, see Beers, Flynn, and Gebhardt, 1990, AJ, 100, 32B

    Parameters
    ----------
    a : array-like
        Input array or object that can be converted to an array.
    c : float
        Tuning constant for the biweight estimator.  Default value is 6.0.
    M : float, optional
        Initial guess for the biweight location.

    Returns
    -------
    biweight_location : float
        Returns the biweight location for the array elements.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the biweight location of the distribution::

    >>> from astropy.stats.funcs import biweight_location
    >>> from numpy.random import randn
    >>> randvar = randn(10000)
    >>> cbl = biweight_location(randvar)

    See Also
    --------
    median_absolute_deviation, biweight_midvariance

    """

    a = np.array(a, copy=False)

    if M is None:
        M = np.median(a)

    # set up the difference
    d = a - M

    # set up the weighting
    u = d / c / median_absolute_deviation(a)

    # now remove the outlier points
    mask = np.abs(u) < 1

    u = (1 - u ** 2) ** 2
    return M + (d[mask] * u[mask]).sum() / u[mask].sum()


def biweight_midvariance(a, c=9.0, M=None):
    """Compute the biweight midvariance for an array.

    Returns the biweight midvariance for the array elements.
    The biweight midvariance is a robust statistic for determining
    the midvariance (i.e. the standard deviation) of a distribution.

    The biweight location is given by the following equation

    .. math::

      C_{bl}= (n')^{1/2} \\frac{[\Sigma_{|u_i|<1} (x_i-M)^2(1-u_i^2)^4]^{0.5}}
      {|\Sigma_{|u_i|<1} (1-u_i^2)(1-5u_i^2)|}

    where :math:`u_i` is given by

    .. math::

      u_{i} = \\frac{(x_i-M)}{cMAD}

    where MAD is the median absolute deviation.

    :math:`n'` is the number of data for which :math:`|u_i| < 1` holds, while the
    summations are over all i up to n:

    .. math::

        n' = \Sigma_{|u_i|<1}^n 1

    This is slightly different than given in the reference below, but
    results in a value closer to the true midvariance.

    The midvariance parameter c is typically 9.0.

    For more details, see Beers, Flynn, and Gebhardt, 1990, AJ, 100, 32B

    Parameters
    ----------
    a : array-like
        Input array or object that can be converted to an array.
    c : float
        Tuning constant for the biweight estimator.  Default value is 9.0.
    M : float, optional
        Initial guess for the biweight location.

    Returns
    -------
    biweight_midvariance : float
        Returns the biweight midvariance for the array elements.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the biweight midvariance of the distribution::

    >>> from astropy.stats.funcs import biweight_midvariance
    >>> from numpy.random import randn
    >>> randvar = randn(10000)
    >>> scl = biweight_midvariance(randvar)

    See Also
    --------
    median_absolute_deviation, biweight_location
    """

    a = np.array(a, copy=False)

    if M is None:
        M = np.median(a)

    # set up the difference
    d = a - M

    # set up the weighting
    u = d / c / median_absolute_deviation(a)

    # now remove the outlier points
    mask = np.abs(u) < 1

    u = u ** 2
    n = mask.sum()
    return n ** 0.5 * (d[mask] * d[mask] * (1 - u[mask]) ** 4).sum() ** 0.5\
        / np.abs(((1 - u[mask]) * (1 - 5 * u[mask])).sum())


def signal_to_noise_oir_ccd(t, source_eps, sky_eps, dark_eps, rd, npix,
                            gain=1.0):
    """Computes the signal to noise ratio for source being observed in the
    optical/IR using a CCD.

    Parameters
    ----------
    t : float or numpy.ndarray
        CCD integration time in seconds
    source_eps : float
        Number of electrons (photons) or DN per second in the aperture from the
        source. Note that this should already have been scaled by the filter
        transmission and the quantum efficiency of the CCD. If the input is in
        DN, then be sure to set the gain to the proper value for the CCD.
        If the input is in electrons per second, then keep the gain as its
        default of 1.0.
    sky_eps : float
        Number of electrons (photons) or DN per second per pixel from the sky
        background. Should already be scaled by filter transmission and QE.
        This must be in the same units as source_eps for the calculation to
        make sense.
    dark_eps : float
        Number of thermal electrons per second per pixel. If this is given in
        DN or ADU, then multiply by the gain to get the value in electrons.
    rd : float
        Read noise of the CCD in electrons. If this is given in
        DN or ADU, then multiply by the gain to get the value in electrons.
    npix : float
        Size of the aperture in pixels
    gain : float
        Gain of the CCD. In units of electrons per DN.

    Returns
    ----------
    SNR : float or numpy.ndarray
        Signal to noise ratio calculated from the inputs
    """
    signal = t * source_eps * gain
    noise = np.sqrt(t * (source_eps * gain + npix *
                         (sky_eps * gain + dark_eps)) + npix * rd ** 2)
    return signal / noise


def bootstrap(data, bootnum=100, samples=None, bootfunc=None):
    """Performs bootstrap resampling on numpy arrays.

    Bootstrap resampling is used to understand confidence intervals of sample
    estimates. This function returns versions of the dataset resampled with
    replacement ("case bootstrapping"). These can all be run through a function
    or statistic to produce a distribution of values which can then be used to
    find the confidence intervals.

    Parameters
    ----------
    data : numpy.ndarray
        N-D array. The bootstrap resampling will be performed on the first
        index, so the first index should access the relevant information
        to be bootstrapped.
    bootnum : int
        Number of bootstrap resamples
    samples : int
        Number of samples in each resample. The default `None` sets samples to
        the number of datapoints
    bootfunc : function
        Function to reduce the resampled data. Each bootstrap resample will
        be put through this function and the results returned. If `None`, the
        bootstrapped data will be returned

    Returns
    -------
    boot : numpy.ndarray
        Bootstrapped data. Each row is a bootstrap resample of the data. The
        columns will correspond to the outputs of bootfunc.

    Examples
    --------
    Obtain a twice resampled array:

    >>> from astropy.stats import bootstrap
    >>> import numpy as np
    >>> from astropy.utils import NumpyRNGContext
    >>> bootarr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 0])
    >>> with NumpyRNGContext(1):
    ...     bootresult = bootstrap(bootarr, 2)
    ...
    >>> bootresult
    array([[ 6.,  9.,  0.,  6.,  1.,  1.,  2.,  8.,  7.,  0.],
           [ 3.,  5.,  6.,  3.,  5.,  3.,  5.,  8.,  8.,  0.]])
    >>> bootresult.shape
    (2, 10)

    Obtain a statistic on the array

    >>> with NumpyRNGContext(1):
    ...     bootresult = bootstrap(bootarr, 2, bootfunc=np.mean)
    ...
    >>> bootresult
    array([ 4. ,  4.6])

    Obtain a statistic with two outputs on the array

    >>> test_statistic = lambda x: (np.sum(x), np.mean(x))
    >>> with NumpyRNGContext(1):
    ...     bootresult = bootstrap(bootarr, 3, bootfunc=test_statistic)
    >>> bootresult
    array([[ 40. ,   4. ],
           [ 46. ,   4.6],
           [ 35. ,   3.5]])
    >>> bootresult.shape
    (3, 2)

    Obtain a statistic with two outputs on the array, keeping only the first
    output

    >>> bootfunc = lambda x:test_statistic(x)[0]
    >>> with NumpyRNGContext(1):
    ...     bootresult = bootstrap(bootarr, 3, bootfunc=bootfunc)
    ...
    >>> bootresult
    array([ 40.,  46.,  35.])
    >>> bootresult.shape
    (3,)

    """
    if samples is None:
        samples = data.shape[0]

    # make sure the input is sane
    assert samples > 0, "samples cannot be less than one"
    assert bootnum > 0, "bootnum cannot be less than one"

    if bootfunc is None:
        resultdims = (bootnum,) + (samples,) + data.shape[1:]
    else:
        # test number of outputs from bootfunc, avoid single outputs which are
        # array-like
        try:
            resultdims = (bootnum, len(bootfunc(data)))
        except TypeError:
            resultdims = (bootnum,)

    # create empty boot array
    boot = np.empty(resultdims)

    for i in xrange(bootnum):
        bootarr = np.random.randint(low=0, high=data.shape[0], size=samples)
        if bootfunc is None:
            boot[i] = data[bootarr]
        else:
            boot[i] = bootfunc(data[bootarr])

    return boot

def mad_std(data):
    """
    Calculate a robust standard deviation using the `median absolute
    deviation (MAD)
    <http://en.wikipedia.org/wiki/Median_absolute_deviation>`_.

    The standard deviation estimator is given by:

    .. math::

        \\sigma \\approx \\frac{\\textrm{MAD}}{\Phi^{-1}(3/4)} \\approx 1.4826 \ \\textrm{MAD}

    where :math:`\Phi^{-1}(P)` is the normal inverse cumulative
    distribution function evaluated at probability :math:`P = 3/4`.

    Parameters
    ----------
    data : array-like
        Data array or object that can be converted to an array.

    Returns
    -------
    result : float
        The robust standard deviation of the data.

    Examples
    --------
    >>> from astropy.stats import mad_std
    >>> from astropy.utils import NumpyRNGContext
    >>> from numpy.random import normal
    >>> with NumpyRNGContext(12345):
    ...     data = normal(5, 2, size=(100, 100))
    ...     mad_std(data)    # doctest: +FLOAT_CMP
    2.02327646594
    """

    # NOTE: 1. / scipy.stats.norm.ppf(0.75) = 1.482602218505602
    return median_absolute_deviation(data) * 1.482602218505602
