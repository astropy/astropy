"""
Utilities for computing periodogram statistics.
"""
import numpy as np
from scipy.special import gammaln

from .core import LombScargle


def _weighted_sum(val, dy):
    return (val / dy ** 2).sum()


def _weighted_mean(val, dy):
    return _weighted_sum(val, dy) / _weighted_sum(np.ones_like(val), dy)


def _weighted_var(val, dy):
    return _weighted_mean(val ** 2, dy) - _weighted_mean(val, dy) ** 2


def _gamma(N):
    # Note: this is closely approximated by (1 - 0.75 / N) for large N
    return np.sqrt(2 / N) * np.exp(gammaln(N / 2) - gammaln((N - 1) / 2))


def _log_gamma(N):
    return 0.5 * np.log(2 / N) + gammaln(N / 2) - gammaln((N - 1) / 2)


def pdf_single(z, N, normalization, dH=1, dK=3):
    """Probability density function for Lomb-Scargle periodogram

    Compute the expected probability density function of the periodogram
    for the null hypothesis - i.e. data consisting of Gaussian noise.

    Parameters
    ----------
    z : array-like
        the periodogram value
    N : int
        the number of data points from which the periodogram was computed
    normalization : string
        The periodogram normalization. Must be one of
        ['standard', 'model', 'log', 'psd']
    dH, dK : integers (optional)
        The number of parameters in the null hypothesis and the model

    Returns
    -------
    pdf : np.ndarray
        The expected probability density function

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.
    All expressions used here are adapted from Table 1 of Baluev 2008 [1]_.

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    if dK - dH != 2:
        raise NotImplementedError("Degrees of freedom != 2")
    Nk = N - dK

    if normalization == 'psd':
        return np.exp(-z)
    elif normalization == 'standard':
        return 0.5 * Nk * (1 - z) ** (0.5 * Nk - 1)
    elif normalization == 'model':
        return 0.5 * Nk * (1 + z) ** (-0.5 * Nk - 1)
    elif normalization == 'log':
        return 0.5 * Nk * np.exp(-0.5 * Nk * z)
    else:
        raise ValueError("normalization='{0}' is not recognized"
                         "".format(normalization))


def fap_single(z, N, normalization, dH=1, dK=3):
    """Single-frequency false alarm probability for the Lomb-Scargle periodogram

    This is equal to 1 - cdf, where cdf is the cumulative distribution.
    The single-frequency false alarm probability should not be confused with
    the false alarm probability for the largest peak.

    Parameters
    ----------
    z : array-like
        the periodogram value
    N : int
        the number of data points from which the periodogram was computed
    normalization : string
        The periodogram normalization. Must be one of
        ['standard', 'model', 'log', 'psd']
    dH, dK : integers (optional)
        The number of parameters in the null hypothesis and the model

    Returns
    -------
    fap : np.ndarray
        The expected cumulative distribution function

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.
    All expressions used here are adapted from Table 1 of Baluev 2008 [1]_.

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    if dK - dH != 2:
        raise NotImplementedError("Degrees of freedom != 2")
    Nk = N - dK

    if normalization == 'psd':
        return np.exp(-z)
    elif normalization == 'standard':
        return (1 - z) ** (0.5 * Nk)
    elif normalization == 'model':
        return (1 + z) ** (-0.5 * Nk)
    elif normalization == 'log':
        return np.exp(-0.5 * Nk * z)
    else:
        raise ValueError("normalization='{0}' is not recognized"
                         "".format(normalization))


def cdf_single(z, N, normalization, dH=1, dK=3):
    """Cumulative distribution for the Lomb-Scargle periodogram

    Compute the expected cumulative distribution of the periodogram
    for the null hypothesis - i.e. data consisting of Gaussian noise.

    Parameters
    ----------
    z : array-like
        the periodogram value
    N : int
        the number of data points from which the periodogram was computed
    normalization : string
        The periodogram normalization. Must be one of
        ['standard', 'model', 'log', 'psd']
    dH, dK : integers (optional)
        The number of parameters in the null hypothesis and the model

    Returns
    -------
    cdf : np.ndarray
        The expected cumulative distribution function

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.
    All expressions used here are adapted from Table 1 of Baluev 2008 [1]_.

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    return 1 - fap_single(z, N, normalization=normalization, dH=dH, dK=dK)


def tau_davies(Z, fmax, t, y, dy, normalization='standard', dH=1, dK=3):
    """tau factor for estimating Davies bound (Baluev 2008, Table 1)"""
    N = len(t)
    NH = N - dH  # DOF for null hypothesis
    NK = N - dK  # DOF for periodic hypothesis
    Dt = _weighted_var(t, dy)
    Teff = np.sqrt(4 * np.pi * Dt)
    W = fmax * Teff
    if normalization == 'psd':
        # 'psd' normalization is same as Baluev's z
        return W * np.exp(-Z) * np.sqrt(Z)
    elif normalization == 'standard':
        # 'standard' normalization is Z = 2/NH * z_1
        return (_gamma(NH) * W * (1 - Z) ** (0.5 * (NK - 1))
                * np.sqrt(0.5 * NH * Z))
    elif normalization == 'model':
        # 'model' normalization is Z = 2/NK * z_2
        return (_gamma(NK) * W * (1 + Z) ** (-0.5 * NK)
                * np.sqrt(0.5 * NK * Z))
    elif normalization == 'log':
        # 'log' normalization is Z = 2/NK * z_3
        return (_gamma(NK) * W * np.exp(-0.5 * Z * (NK - 0.5))
                * np.sqrt(NK * np.sinh(0.5 * Z)))
    else:
        raise NotImplementedError("normalization={0}".format(normalization))


def fap_simple(Z, fmax, t, y, dy, normalization='standard'):
    """False Alarm Probability based on estimated number of indep frequencies"""
    N = len(t)
    T = max(t) - min(t)
    N_eff = fmax * T
    p_s = cdf_single(Z, N, normalization=normalization)
    print(Z)
    print(p_s)
    return 1 - p_s ** N_eff


def fap_davies(Z, fmax, t, y, dy, normalization='standard'):
    """Davies upper-bound to the false alarm probability

    (Eqn 5 of Baluev 2008)
    """
    N = len(t)
    fap_s = fap_single(Z, N, normalization=normalization)
    tau = tau_davies(Z, fmax, t, y, dy, normalization=normalization)
    return fap_s + tau


def fap_baluev(Z, fmax, t, y, dy, normalization='standard'):
    """Alias-free approximation to false alarm probability

    (Eqn 6 of Baluev 2008)
    """
    cdf = cdf_single(Z, len(t), normalization)
    tau = tau_davies(Z, fmax, t, y, dy, normalization=normalization)
    return 1 - cdf * np.exp(-tau)


def fap_bootstrap(Z, fmax, t, y, dy, normalization='standard',
                  n_bootstraps=1000, random_seed=None):
    rng = np.random.RandomState(random_seed)

    def bootstrapped_power():
        resample = rng.randint(0, len(y), len(y))  # sample with replacement
        ls_boot = LombScargle(t, y[resample], dy[resample])
        freq, power = ls_boot.autopower(normalization=normalization,
                                        maximum_frequency=fmax)
        return power.max()

    pmax = np.array([bootstrapped_power() for i in range(n_bootstraps)])
    pmax.sort()
    return 1 - np.searchsorted(pmax, Z) / len(pmax)


METHODS = {'simple': fap_simple,
           'davies': fap_davies,
           'baluev': fap_baluev,
           'bootstrap': fap_bootstrap}


def false_alarm_probability(Z, fmax, t, y, dy, normalization,
                            method='baluev', method_kwds=None):
    """Approximate the False Alarm Probability

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO
    """
    if method not in METHODS:
        raise ValueError("Unrecognized method: {0}".format(method))
    method = METHODS[method]
    method_kwds = method_kwds or {}

    return method(Z, fmax, t, y, dy, normalization, **method_kwds)
