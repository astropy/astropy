# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains simple statistical algorithms that are
straightforwardly implemented as a single python function (or family of
functions).

This module should generally not be used directly.  Everything in
`__all__` is imported into `astropy.stats`, and hence that package
should be used for access.
"""


import math
import itertools

import numpy as np

from warnings import warn

from ..utils.decorators import deprecated_renamed_argument
from ..utils import isiterable
from . import _stats


__all__ = ['gaussian_fwhm_to_sigma', 'gaussian_sigma_to_fwhm',
           'binom_conf_interval', 'binned_binom_proportion',
           'poisson_conf_interval', 'median_absolute_deviation', 'mad_std',
           'signal_to_noise_oir_ccd', 'bootstrap', 'kuiper', 'kuiper_two',
           'kuiper_false_positive_probability', 'cdf_from_intervals',
           'interval_overlap_length', 'histogram_intervals', 'fold_intervals']

__doctest_skip__ = ['binned_binom_proportion']
__doctest_requires__ = {'binom_conf_interval': ['scipy.special'],
                        'poisson_conf_interval': ['scipy.special',
                                                  'scipy.optimize',
                                                  'scipy.integrate']}


gaussian_sigma_to_fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))
"""
Factor with which to multiply Gaussian 1-sigma standard deviation to
convert it to full width at half maximum (FWHM).
"""

gaussian_fwhm_to_sigma = 1. / gaussian_sigma_to_fwhm
"""
Factor with which to multiply Gaussian full width at half maximum (FWHM)
to convert it to 1-sigma standard deviation.
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

    k = np.asarray(k).astype(int)
    n = np.asarray(n).astype(int)

    if (n <= 0).any():
        raise ValueError('n must be positive')
    if (k < 0).any() or (k > n).any():
        raise ValueError('k must be in {0, 1, .., n}')

    if interval == 'wilson' or interval == 'wald':
        from scipy.special import erfinv
        kappa = np.sqrt(2.) * min(erfinv(conf), 1.e10)  # Avoid overflows.
        k = k.astype(float)
        n = n.astype(float)
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
    success = np.ravel(success).astype(bool)
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


def _check_poisson_conf_inputs(sigma, background, conflevel, name):
    if sigma != 1:
        raise ValueError("Only sigma=1 supported for interval {0}"
                         .format(name))
    if background != 0:
        raise ValueError("background not supported for interval {0}"
                         .format(name))
    if conflevel is not None:
        raise ValueError("conflevel not supported for interval {0}"
                         .format(name))


def poisson_conf_interval(n, interval='root-n', sigma=1, background=0,
                          conflevel=None):
    r"""Poisson parameter confidence interval given observed counts

    Parameters
    ----------
    n : int or numpy.ndarray
        Number of counts (0 <= ``n``).
    interval : {'root-n','root-n-0','pearson','sherpagehrels','frequentist-confidence', 'kraft-burrows-nousek'}, optional
        Formula used for confidence interval. See notes for details.
        Default is ``'root-n'``.
    sigma : float, optional
        Number of sigma for confidence interval; only supported for
        the 'frequentist-confidence' mode.
    background : float, optional
        Number of counts expected from the background; only supported for
        the 'kraft-burrows-nousek' mode. This number is assumed to be determined
        from a large region so that the uncertainty on its value is negligible.
    conflevel : float, optional
        Confidence level between 0 and 1; only supported for the
        'kraft-burrows-nousek' mode.


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
    :math:`\alpha` is the one-tailed probability of the normal
    distribution (at the point given by the parameter 'sigma'). See
    [Maxwell 2011][maxw11] for further details.

    **6. 'kraft-burrows-nousek'** This is a Bayesian approach which allows
    for the presence of a known background :math:`B` in the source signal
    :math:`N`.
    For a given confidence level :math:`CL` the confidence interval
    :math:`[S_\mathrm{min}, S_\mathrm{max}]` is given by:

    .. math::

       CL = \int^{S_\mathrm{max}}_{S_\mathrm{min}} f_{N,B}(S)dS

    where the function :math:`f_{N,B}` is:

    .. math::

       f_{N,B}(S) = C \frac{e^{-(S+B)}(S+B)^N}{N!}

    and the normalization constant :math:`C`:

    .. math::

       C = \left[ \int_0^\infty \frac{e^{-(S+B)}(S+B)^N}{N!} dS \right] ^{-1}
       = \left( \sum^N_{n=0} \frac{e^{-B}B^n}{n!}  \right)^{-1}

    See [KraftBurrowsNousek][kbn1991] for further details.

    These formulas implement a positive, uniform prior.
    [KraftBurrowsNousek][kbn1991] discuss this choice in more detail and show
    that the problem is relatively insensitive to the choice of prior.

    This functions has an optional dependency: Either scipy or
    `mpmath <http://mpmath.org/>`_  need to be available. (Scipy only works for
    N < 100).


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

    >>> poisson_conf_interval(10, background=1.5, conflevel=0.95,
    ...                       interval='kraft-burrows-nousek').T
    array([  3.47894005, 16.113329533])   # doctest: +FLOAT_CMP

    [pois_eb]: http://www-cdf.fnal.gov/physics/statistics/notes/pois_eb.txt

    [ErrorBars]: http://www.pp.rhul.ac.uk/~cowan/atlas/ErrorBars.pdf

    [ac12]: http://adsabs.harvard.edu/abs/2012EPJP..127...24A

    [maxw11]: http://adsabs.harvard.edu/abs/2011arXiv1102.0822M

    [gehrels86]: http://adsabs.harvard.edu/abs/1986ApJ...303..336G

    [sherpa_gehrels]: http://cxc.harvard.edu/sherpa4.4/statistics/#chigehrels

    [kbn1991]: http://adsabs.harvard.edu/abs/1991ApJ...374..344K

    """

    if not np.isscalar(n):
        n = np.asanyarray(n)

    if interval == 'root-n':
        _check_poisson_conf_inputs(sigma, background, conflevel, interval)
        conf_interval = np.array([n - np.sqrt(n),
                                  n + np.sqrt(n)])
    elif interval == 'root-n-0':
        _check_poisson_conf_inputs(sigma, background, conflevel, interval)
        conf_interval = np.array([n - np.sqrt(n),
                                  n + np.sqrt(n)])
        if np.isscalar(n):
            if n == 0:
                conf_interval[1] = 1
        else:
            conf_interval[1, n == 0] = 1
    elif interval == 'pearson':
        _check_poisson_conf_inputs(sigma, background, conflevel, interval)
        conf_interval = np.array([n + 0.5 - np.sqrt(n + 0.25),
                                  n + 0.5 + np.sqrt(n + 0.25)])
    elif interval == 'sherpagehrels':
        _check_poisson_conf_inputs(sigma, background, conflevel, interval)
        conf_interval = np.array([n - 1 - np.sqrt(n + 0.75),
                                  n + 1 + np.sqrt(n + 0.75)])
    elif interval == 'frequentist-confidence':
        _check_poisson_conf_inputs(1., background, conflevel, interval)
        import scipy.stats
        alpha = scipy.stats.norm.sf(sigma)
        conf_interval = np.array([0.5 * scipy.stats.chi2(2 * n).ppf(alpha),
                                  0.5 * scipy.stats.chi2(2 * n + 2).isf(alpha)])
        if np.isscalar(n):
            if n == 0:
                conf_interval[0] = 0
        else:
            conf_interval[0, n == 0] = 0
    elif interval == 'kraft-burrows-nousek':
        if conflevel is None:
            raise ValueError('Set conflevel for method {0}. (sigma is '
                             'ignored.)'.format(interval))
        conflevel = np.asanyarray(conflevel)
        if np.any(conflevel <= 0) or np.any(conflevel >= 1):
            raise ValueError('Conflevel must be a number between 0 and 1.')
        background = np.asanyarray(background)
        if np.any(background < 0):
            raise ValueError('Background must be >= 0.')
        conf_interval = np.vectorize(_kraft_burrows_nousek,
                                     cache=True)(n, background, conflevel)
        conf_interval = np.vstack(conf_interval)
    else:
        raise ValueError("Invalid method for Poisson confidence intervals: "
                         "{}".format(interval))
    return conf_interval


@deprecated_renamed_argument('a', 'data', '2.0')
def median_absolute_deviation(data, axis=None, func=None, ignore_nan=False):
    """
    Calculate the median absolute deviation (MAD).

    The MAD is defined as ``median(abs(a - median(a)))``.

    Parameters
    ----------
    data : array-like
        Input array or object that can be converted to an array.
    axis : {int, sequence of int, None}, optional
        Axis along which the MADs are computed.  The default (`None`) is
        to compute the MAD of the flattened array.
    func : callable, optional
        The function used to compute the median. Defaults to `numpy.ma.median`
        for masked arrays, otherwise to `numpy.median`.
    ignore_nan : bool
        Ignore NaN values (treat them as if they are not in the array) when
        computing the median.  This will use `numpy.ma.median` if ``axis`` is
        specified, or `numpy.nanmedian` if ``axis==None`` and numpy's version
        is >1.10 because nanmedian is slightly faster in this case.

    Returns
    -------
    mad : float or `~numpy.ndarray`
        The median absolute deviation of the input array.  If ``axis``
        is `None` then a scalar will be returned, otherwise a
        `~numpy.ndarray` will be returned.

    Examples
    --------
    Generate random variates from a Gaussian distribution and return the
    median absolute deviation for that distribution::

        >>> import numpy as np
        >>> from astropy.stats import median_absolute_deviation
        >>> rand = np.random.RandomState(12345)
        >>> from numpy.random import randn
        >>> mad = median_absolute_deviation(rand.randn(1000))
        >>> print(mad)    # doctest: +FLOAT_CMP
        0.65244241428454486

    See Also
    --------
    mad_std
    """

    if func is None:
        # Check if the array has a mask and if so use np.ma.median
        # See https://github.com/numpy/numpy/issues/7330 why using np.ma.median
        # for normal arrays should not be done (summary: np.ma.median always
        # returns an masked array even if the result should be scalar). (#4658)
        if isinstance(data, np.ma.MaskedArray):
            is_masked = True
            func = np.ma.median
            if ignore_nan:
                data = np.ma.masked_invalid(data)
        elif ignore_nan:
            is_masked = False
            func = np.nanmedian
        else:
            is_masked = False
            func = np.median
    else:
        is_masked = None

    data = np.asanyarray(data)
    # np.nanmedian has `keepdims`, which is a good option if we're not allowing
    # user-passed functions here
    data_median = func(data, axis=axis)

    # broadcast the median array before subtraction
    if axis is not None:
        if isiterable(axis):
            for ax in sorted(list(axis)):
                data_median = np.expand_dims(data_median, axis=ax)
        else:
            data_median = np.expand_dims(data_median, axis=axis)

    result = func(np.abs(data - data_median), axis=axis, overwrite_input=True)

    if axis is None and np.ma.isMaskedArray(result):
        # return scalar version
        result = result.item()
    elif np.ma.isMaskedArray(result) and not is_masked:
        # if the input array was not a masked array, we don't want to return a
        # masked array
        result = result.filled(fill_value=np.nan)

    return result


def mad_std(data, axis=None, func=None, ignore_nan=False):
    r"""
    Calculate a robust standard deviation using the `median absolute
    deviation (MAD)
    <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_.

    The standard deviation estimator is given by:

    .. math::

        \sigma \approx \frac{\textrm{MAD}}{\Phi^{-1}(3/4)}
            \approx 1.4826 \ \textrm{MAD}

    where :math:`\Phi^{-1}(P)` is the normal inverse cumulative
    distribution function evaluated at probability :math:`P = 3/4`.

    Parameters
    ----------
    data : array-like
        Data array or object that can be converted to an array.
    axis : {int, sequence of int, None}, optional
        Axis along which the robust standard deviations are computed.
        The default (`None`) is to compute the robust standard deviation
        of the flattened array.
    func : callable, optional
        The function used to compute the median. Defaults to `numpy.ma.median`
        for masked arrays, otherwise to `numpy.median`.
    ignore_nan : bool
        Ignore NaN values (treat them as if they are not in the array) when
        computing the median.  This will use `numpy.ma.median` if ``axis`` is
        specified, or `numpy.nanmedian` if ``axis=None`` and numpy's version is
        >1.10 because nanmedian is slightly faster in this case.

    Returns
    -------
    mad_std : float or `~numpy.ndarray`
        The robust standard deviation of the input data.  If ``axis`` is
        `None` then a scalar will be returned, otherwise a
        `~numpy.ndarray` will be returned.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import mad_std
    >>> rand = np.random.RandomState(12345)
    >>> madstd = mad_std(rand.normal(5, 2, (100, 100)))
    >>> print(madstd)    # doctest: +FLOAT_CMP
    2.0232764659422626

    See Also
    --------
    biweight_midvariance, biweight_midcovariance, median_absolute_deviation
    """

    # NOTE: 1. / scipy.stats.norm.ppf(0.75) = 1.482602218505602
    MAD = median_absolute_deviation(
        data, axis=axis, func=func, ignore_nan=ignore_nan)
    return MAD * 1.482602218505602


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
    gain : float, optional
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
    bootnum : int, optional
        Number of bootstrap resamples
    samples : int, optional
        Number of samples in each resample. The default `None` sets samples to
        the number of datapoints
    bootfunc : function, optional
        Function to reduce the resampled data. Each bootstrap resample will
        be put through this function and the results returned. If `None`, the
        bootstrapped data will be returned

    Returns
    -------
    boot : numpy.ndarray

        If bootfunc is None, then each row is a bootstrap resample of the data.
        If bootfunc is specified, then the columns will correspond to the
        outputs of bootfunc.

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
    >>> bootresult  # doctest: +FLOAT_CMP
    array([[6., 9., 0., 6., 1., 1., 2., 8., 7., 0.],
           [3., 5., 6., 3., 5., 3., 5., 8., 8., 0.]])
    >>> bootresult.shape
    (2, 10)

    Obtain a statistic on the array

    >>> with NumpyRNGContext(1):
    ...     bootresult = bootstrap(bootarr, 2, bootfunc=np.mean)
    ...
    >>> bootresult  # doctest: +FLOAT_CMP
    array([4. , 4.6])

    Obtain a statistic with two outputs on the array

    >>> test_statistic = lambda x: (np.sum(x), np.mean(x))
    >>> with NumpyRNGContext(1):
    ...     bootresult = bootstrap(bootarr, 3, bootfunc=test_statistic)
    >>> bootresult  # doctest: +FLOAT_CMP
    array([[40. ,  4. ],
           [46. ,  4.6],
           [35. ,  3.5]])
    >>> bootresult.shape
    (3, 2)

    Obtain a statistic with two outputs on the array, keeping only the first
    output

    >>> bootfunc = lambda x:test_statistic(x)[0]
    >>> with NumpyRNGContext(1):
    ...     bootresult = bootstrap(bootarr, 3, bootfunc=bootfunc)
    ...
    >>> bootresult  # doctest: +FLOAT_CMP
    array([40., 46., 35.])
    >>> bootresult.shape
    (3,)

    """
    if samples is None:
        samples = data.shape[0]

    # make sure the input is sane
    if samples < 1 or bootnum < 1:
        raise ValueError("neither 'samples' nor 'bootnum' can be less than 1.")

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

    for i in range(bootnum):
        bootarr = np.random.randint(low=0, high=data.shape[0], size=samples)
        if bootfunc is None:
            boot[i] = data[bootarr]
        else:
            boot[i] = bootfunc(data[bootarr])

    return boot


def _scipy_kraft_burrows_nousek(N, B, CL):
    '''Upper limit on a poisson count rate

    The implementation is based on Kraft, Burrows and Nousek
    `ApJ 374, 344 (1991) <http://adsabs.harvard.edu/abs/1991ApJ...374..344K>`_.
    The XMM-Newton upper limit server uses the same formalism.

    Parameters
    ----------
    N : int
        Total observed count number
    B : float
        Background count rate (assumed to be known with negligible error
        from a large background area).
    CL : float
       Confidence level (number between 0 and 1)

    Returns
    -------
    S : source count limit

    Notes
    -----
    Requires `scipy`. This implementation will cause Overflow Errors for about
    N > 100 (the exact limit depends on details of how scipy was compiled).
    See `~astropy.stats.mpmath_poisson_upper_limit` for an implementation that
    is slower, but can deal with arbitrarily high numbers since it is based on
    the `mpmath <http://mpmath.org/>`_ library.
    '''

    from scipy.optimize import brentq
    from scipy.integrate import quad

    from math import exp

    def eqn8(N, B):
        n = np.arange(N + 1, dtype=np.float64)
        # Create an array containing the factorials. scipy.special.factorial
        # requires SciPy 0.14 (#5064) therefore this is calculated by using
        # numpy.cumprod. This could be replaced by factorial again as soon as
        # older SciPy are not supported anymore but the cumprod alternative
        # might also be a bit faster.
        factorial_n = np.ones(n.shape, dtype=np.float64)
        np.cumprod(n[1:], out=factorial_n[1:])
        return 1. / (exp(-B) * np.sum(np.power(B, n) / factorial_n))

    # The parameters of eqn8 do not vary between calls so we can calculate the
    # result once and reuse it. The same is True for the factorial of N.
    # eqn7 is called hundred times so "caching" these values yields a
    # significant speedup (factor 10).
    eqn8_res = eqn8(N, B)
    factorial_N = float(math.factorial(N))

    def eqn7(S, N, B):
        SpB = S + B
        return eqn8_res * (exp(-SpB) * SpB**N / factorial_N)

    def eqn9_left(S_min, S_max, N, B):
        return quad(eqn7, S_min, S_max, args=(N, B), limit=500)

    def find_s_min(S_max, N, B):
        '''
        Kraft, Burrows and Nousek suggest to integrate from N-B in both
        directions at once, so that S_min and S_max move similarly (see
        the article for details). Here, this is implemented differently:
        Treat S_max as the optimization parameters in func and then
        calculate the matching s_min that has has eqn7(S_max) =
        eqn7(S_min) here.
        '''
        y_S_max = eqn7(S_max, N, B)
        if eqn7(0, N, B) >= y_S_max:
            return 0.
        else:
            return brentq(lambda x: eqn7(x, N, B) - y_S_max, 0, N - B)

    def func(s):
        s_min = find_s_min(s, N, B)
        out = eqn9_left(s_min, s, N, B)
        return out[0] - CL

    S_max = brentq(func, N - B, 100)
    S_min = find_s_min(S_max, N, B)
    return S_min, S_max


def _mpmath_kraft_burrows_nousek(N, B, CL):
    '''Upper limit on a poisson count rate

    The implementation is based on Kraft, Burrows and Nousek in
    `ApJ 374, 344 (1991) <http://adsabs.harvard.edu/abs/1991ApJ...374..344K>`_.
    The XMM-Newton upper limit server used the same formalism.

    Parameters
    ----------
    N : int
        Total observed count number
    B : float
        Background count rate (assumed to be known with negligible error
        from a large background area).
    CL : float
       Confidence level (number between 0 and 1)

    Returns
    -------
    S : source count limit

    Notes
    -----
    Requires the `mpmath <http://mpmath.org/>`_ library.  See
    `~astropy.stats.scipy_poisson_upper_limit` for an implementation
    that is based on scipy and evaluates faster, but runs only to about
    N = 100.
    '''
    from mpmath import mpf, factorial, findroot, fsum, power, exp, quad

    N = mpf(N)
    B = mpf(B)
    CL = mpf(CL)

    def eqn8(N, B):
        sumterms = [power(B, n) / factorial(n) for n in range(int(N) + 1)]
        return 1. / (exp(-B) * fsum(sumterms))

    eqn8_res = eqn8(N, B)
    factorial_N = factorial(N)

    def eqn7(S, N, B):
        SpB = S + B
        return eqn8_res * (exp(-SpB) * SpB**N / factorial_N)

    def eqn9_left(S_min, S_max, N, B):
        def eqn7NB(S):
            return eqn7(S, N, B)
        return quad(eqn7NB, [S_min, S_max])

    def find_s_min(S_max, N, B):
        '''
        Kraft, Burrows and Nousek suggest to integrate from N-B in both
        directions at once, so that S_min and S_max move similarly (see
        the article for details). Here, this is implemented differently:
        Treat S_max as the optimization parameters in func and then
        calculate the matching s_min that has has eqn7(S_max) =
        eqn7(S_min) here.
        '''
        y_S_max = eqn7(S_max, N, B)
        if eqn7(0, N, B) >= y_S_max:
            return 0.
        else:
            def eqn7ysmax(x):
                return eqn7(x, N, B) - y_S_max
            return findroot(eqn7ysmax, (N - B) / 2.)

    def func(s):
        s_min = find_s_min(s, N, B)
        out = eqn9_left(s_min, s, N, B)
        return out - CL

    S_max = findroot(func, N - B, tol=1e-4)
    S_min = find_s_min(S_max, N, B)
    return float(S_min), float(S_max)


def _kraft_burrows_nousek(N, B, CL):
    '''Upper limit on a poisson count rate

    The implementation is based on Kraft, Burrows and Nousek in
    `ApJ 374, 344 (1991) <http://adsabs.harvard.edu/abs/1991ApJ...374..344K>`_.
    The XMM-Newton upper limit server used the same formalism.

    Parameters
    ----------
    N : int
        Total observed count number
    B : float
        Background count rate (assumed to be known with negligible error
        from a large background area).
    CL : float
       Confidence level (number between 0 and 1)

    Returns
    -------
    S : source count limit

    Notes
    -----
    This functions has an optional dependency: Either `scipy` or `mpmath
    <http://mpmath.org/>`_  need to be available. (Scipy only works for
    N < 100).
    '''
    try:
        import scipy
        HAS_SCIPY = True
    except ImportError:
        HAS_SCIPY = False

    try:
        import mpmath
        HAS_MPMATH = True
    except ImportError:
        HAS_MPMATH = False

    if HAS_SCIPY and N <= 100:
        try:
            return _scipy_kraft_burrows_nousek(N, B, CL)
        except OverflowError:
            if not HAS_MPMATH:
                raise ValueError('Need mpmath package for input numbers this '
                                 'large.')
    if HAS_MPMATH:
        return _mpmath_kraft_burrows_nousek(N, B, CL)

    raise ImportError('Either scipy or mpmath are required.')


def kuiper_false_positive_probability(D, N):
    """Compute the false positive probability for the Kuiper statistic.

    Uses the set of four formulas described in Paltani 2004; they report
    the resulting function never underestimates the false positive
    probability but can be a bit high in the N=40..50 range.
    (They quote a factor 1.5 at the 1e-7 level.)

    Parameters
    ----------
    D : float
        The Kuiper test score.
    N : float
        The effective sample size.

    Returns
    -------
    fpp : float
        The probability of a score this large arising from the null hypothesis.

    Notes
    -----
    Eq 7 of Paltani 2004 appears to incorrectly quote the original formula
    (Stephens 1965). This function implements the original formula, as it
    produces a result closer to Monte Carlo simulations.

    References
    ----------

    .. [1] Paltani, S., "Searching for periods in X-ray observations using
           Kuiper's test. Application to the ROSAT PSPC archive",
           Astronomy and Astrophysics, v.240, p.789-790, 2004.

    .. [2] Stephens, M. A., "The goodness-of-fit statistic VN: distribution
           and significance points", Biometrika, v.52, p.309, 1965.

    """
    try:
        from scipy.special import factorial, comb
    except ImportError:
        # Retained for backwards compatibility with older versions of scipy
        # (factorial appears to have moved here in 0.14)
        from scipy.misc import factorial, comb

    if D < 0. or D > 2.:
        raise ValueError("Must have 0<=D<=2 by definition of the Kuiper test")

    if D < 2. / N:
        return 1. - factorial(N) * (D - 1. / N)**(N - 1)
    elif D < 3. / N:
        k = -(N * D - 1.) / 2.
        r = np.sqrt(k**2 - (N * D - 2.)**2 / 2.)
        a, b = -k + r, -k - r
        return 1 - (factorial(N - 1) * (b**(N - 1) * (1 - a) - a**(N - 1) * (1 - b))
                    / N**(N - 2) / (b - a))
    elif (D > 0.5 and N % 2 == 0) or (D > (N - 1.) / (2. * N) and N % 2 == 1):
        # NOTE: the upper limit of this sum is taken from Stephens 1965
        t = np.arange(np.floor(N * (1 - D)) + 1)
        y = D + t / N
        Tt = y**(t - 3) * (y**3 * N
                           - y**2 * t * (3 - 2 / N)
                           + y * t * (t - 1) * (3 - 2 / N) / N
                           - t * (t - 1) * (t - 2) / N**2)
        term = Tt * comb(N, t) * (1 - D - t / N)**(N - t - 1)
        return term.sum()
    else:
        z = D * np.sqrt(N)
        # When m*z>18.82 (sqrt(-log(finfo(double))/2)), exp(-2m**2z**2)
        # underflows.  Cutting off just before avoids triggering a (pointless)
        # underflow warning if `under="warn"`.
        ms = np.arange(1, 18.82 / z)
        S1 = (2 * (4 * ms**2 * z**2 - 1) * np.exp(-2 * ms**2 * z**2)).sum()
        S2 = (ms**2 * (4 * ms**2 * z**2 - 3) * np.exp(-2 * ms**2 * z**2)).sum()
        return S1 - 8 * D / 3 * S2


def kuiper(data, cdf=lambda x: x, args=()):
    """Compute the Kuiper statistic.

    Use the Kuiper statistic version of the Kolmogorov-Smirnov test to
    find the probability that a sample like ``data`` was drawn from the
    distribution whose CDF is given as ``cdf``.

    .. warning::
        This will not work correctly for distributions that are actually
        discrete (Poisson, for example).

    Parameters
    ----------
    data : array-like
        The data values.
    cdf : callable
        A callable to evaluate the CDF of the distribution being tested
        against. Will be called with a vector of all values at once.
        The default is a uniform distribution.
    args : list-like, optional
        Additional arguments to be supplied to cdf.

    Returns
    -------
    D : float
        The raw statistic.
    fpp : float
        The probability of a D this large arising with a sample drawn from
        the distribution whose CDF is cdf.

    Notes
    -----
    The Kuiper statistic resembles the Kolmogorov-Smirnov test in that
    it is nonparametric and invariant under reparameterizations of the data.
    The Kuiper statistic, in addition, is equally sensitive throughout
    the domain, and it is also invariant under cyclic permutations (making
    it particularly appropriate for analyzing circular data).

    Returns (D, fpp), where D is the Kuiper D number and fpp is the
    probability that a value as large as D would occur if data was
    drawn from cdf.

    .. warning::
        The fpp is calculated only approximately, and it can be
        as much as 1.5 times the true value.

    Stephens 1970 claims this is more effective than the KS at detecting
    changes in the variance of a distribution; the KS is (he claims) more
    sensitive at detecting changes in the mean.

    If cdf was obtained from data by fitting, then fpp is not correct and
    it will be necessary to do Monte Carlo simulations to interpret D.
    D should normally be independent of the shape of CDF.

    References
    ----------

    .. [1] Stephens, M. A., "Use of the Kolmogorov-Smirnov, Cramer-Von Mises
           and Related Statistics Without Extensive Tables", Journal of the
           Royal Statistical Society. Series B (Methodological), Vol. 32,
           No. 1. (1970), pp. 115-122.


    """

    data = np.sort(data)
    cdfv = cdf(data, *args)
    N = len(data)
    D = (np.amax(cdfv - np.arange(N) / float(N)) +
         np.amax((np.arange(N) + 1) / float(N) - cdfv))

    return D, kuiper_false_positive_probability(D, N)


def kuiper_two(data1, data2):
    """Compute the Kuiper statistic to compare two samples.

    Parameters
    ----------
    data1 : array-like
        The first set of data values.
    data2 : array-like
        The second set of data values.

    Returns
    -------
    D : float
        The raw test statistic.
    fpp : float
        The probability of obtaining two samples this different from
        the same distribution.

    .. warning::
        The fpp is quite approximate, especially for small samples.

    """
    data1 = np.sort(data1)
    data2 = np.sort(data2)
    n1, = data1.shape
    n2, = data2.shape
    common_type = np.find_common_type([], [data1.dtype, data2.dtype])
    if not (np.issubdtype(common_type, np.number)
            and not np.issubdtype(common_type, np.complexfloating)):
        raise ValueError('kuiper_two only accepts real inputs')
    # nans, if any, are at the end after sorting.
    if np.isnan(data1[-1]) or np.isnan(data2[-1]):
        raise ValueError('kuiper_two only accepts non-nan inputs')
    D = _stats.ks_2samp(np.asarray(data1, common_type),
                        np.asarray(data2, common_type))
    Ne = len(data1) * len(data2) / float(len(data1) + len(data2))
    return D, kuiper_false_positive_probability(D, Ne)


def fold_intervals(intervals):
    """Fold the weighted intervals to the interval (0,1).

    Convert a list of intervals (ai, bi, wi) to a list of non-overlapping
    intervals covering (0,1). Each output interval has a weight equal
    to the sum of the wis of all the intervals that include it. All intervals
    are interpreted modulo 1, and weights are accumulated counting
    multiplicity. This is appropriate, for example, if you have one or more
    blocks of observation and you want to determine how much observation
    time was spent on different parts of a system's orbit (the blocks
    should be converted to units of the orbital period first).

    Parameters
    ----------
    intervals : list of three-element tuples (ai,bi,wi)
        The intervals to fold; ai and bi are the limits of the interval, and
        wi is the weight to apply to the interval.

    Returns
    -------
    breaks : array of floats length N
        The endpoints of a set of intervals covering [0,1]; breaks[0]=0 and
        breaks[-1] = 1
    weights : array of floats of length N-1
        The ith element is the sum of number of times the interval
        breaks[i],breaks[i+1] is included in each interval times the weight
        associated with that interval.

    """
    r = []
    breaks = set()
    tot = 0
    for (a, b, wt) in intervals:
        tot += (np.ceil(b) - np.floor(a)) * wt
        fa = a % 1
        breaks.add(fa)
        r.append((0, fa, -wt))
        fb = b % 1
        breaks.add(fb)
        r.append((fb, 1, -wt))

    breaks.add(0.)
    breaks.add(1.)
    breaks = sorted(breaks)
    breaks_map = dict([(f, i) for (i, f) in enumerate(breaks)])
    totals = np.zeros(len(breaks) - 1)
    totals += tot
    for (a, b, wt) in r:
        totals[breaks_map[a]:breaks_map[b]] += wt
    return np.array(breaks), totals


def cdf_from_intervals(breaks, totals):
    """Construct a callable piecewise-linear CDF from a pair of arrays.

    Take a pair of arrays in the format returned by fold_intervals and
    make a callable cumulative distribution function on the interval
    (0,1).

    Parameters
    ----------
    breaks : array of floats of length N
        The boundaries of successive intervals.
    totals : array of floats of length N-1
        The weight for each interval.

    Returns
    -------
    f : callable
        A cumulative distribution function corresponding to the
        piecewise-constant probability distribution given by breaks, weights

    """
    if breaks[0] != 0 or breaks[-1] != 1:
        raise ValueError("Intervals must be restricted to [0,1]")
    if np.any(np.diff(breaks) <= 0):
        raise ValueError("Breaks must be strictly increasing")
    if np.any(totals < 0):
        raise ValueError(
            "Total weights in each subinterval must be nonnegative")
    if np.all(totals == 0):
        raise ValueError("At least one interval must have positive exposure")
    b = breaks.copy()
    c = np.concatenate(((0,), np.cumsum(totals * np.diff(b))))
    c /= c[-1]
    return lambda x: np.interp(x, b, c, 0, 1)


def interval_overlap_length(i1, i2):
    """Compute the length of overlap of two intervals.

    Parameters
    ----------
    i1, i2 : pairs of two floats
        The two intervals.

    Returns
    -------
    l : float
        The length of the overlap between the two intervals.

    """
    (a, b) = i1
    (c, d) = i2
    if a < c:
        if b < c:
            return 0.
        elif b < d:
            return b - c
        else:
            return d - c
    elif a < d:
        if b < d:
            return b - a
        else:
            return d - a
    else:
        return 0


def histogram_intervals(n, breaks, totals):
    """Histogram of a piecewise-constant weight function.

    This function takes a piecewise-constant weight function and
    computes the average weight in each histogram bin.

    Parameters
    ----------
    n : int
        The number of bins
    breaks : array of floats of length N
        Endpoints of the intervals in the PDF
    totals : array of floats of length N-1
        Probability densities in each bin

    Returns
    -------
    h : array of floats
        The average weight for each bin

    """
    h = np.zeros(n)
    start = breaks[0]
    for i in range(len(totals)):
        end = breaks[i + 1]
        for j in range(n):
            ol = interval_overlap_length((float(j) / n,
                                          float(j + 1) / n), (start, end))
            h[j] += ol / (1. / n) * totals[i]
        start = end

    return h
