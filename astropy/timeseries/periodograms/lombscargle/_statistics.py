"""
Utilities for computing periodogram statistics.

This is an internal module; users should access this functionality via the
``false_alarm_probability`` and ``false_alarm_level`` methods of the
``astropy.timeseries.LombScargle`` API.
"""

from functools import wraps

import numpy as np

from astropy import units as u


def _weighted_sum(val, dy):
    if dy is not None:
        return (val / dy ** 2).sum()
    else:
        return val.sum()


def _weighted_mean(val, dy):
    if dy is None:
        return val.mean()
    else:
        return _weighted_sum(val, dy) / _weighted_sum(np.ones(val.shape), dy)


def _weighted_var(val, dy):
    return _weighted_mean(val ** 2, dy) - _weighted_mean(val, dy) ** 2


def _gamma(N):
    from scipy.special import gammaln
    # Note: this is closely approximated by (1 - 0.75 / N) for large N
    return np.sqrt(2 / N) * np.exp(gammaln(N / 2) - gammaln((N - 1) / 2))


def vectorize_first_argument(func):
    @wraps(func)
    def new_func(x, *args, **kwargs):
        x = np.asarray(x)
        return np.array([func(xi, *args, **kwargs)
                         for xi in x.flat]).reshape(x.shape)
    return new_func


def pdf_single(z, N, normalization, dH=1, dK=3):
    """Probability density function for Lomb-Scargle periodogram

    Compute the expected probability density function of the periodogram
    for the null hypothesis - i.e. data consisting of Gaussian noise.

    Parameters
    ----------
    z : array-like
        The periodogram value.
    N : int
        The number of data points from which the periodogram was computed.
    normalization : {'standard', 'model', 'log', 'psd'}
        The periodogram normalization.
    dH, dK : int, optional
        The number of parameters in the null hypothesis and the model.

    Returns
    -------
    pdf : np.ndarray
        The expected probability density function.

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.
    All expressions used here are adapted from Table 1 of Baluev 2008 [1]_.

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    z = np.asarray(z)
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
        raise ValueError(f"normalization='{normalization}' is not recognized")


def fap_single(z, N, normalization, dH=1, dK=3):
    """Single-frequency false alarm probability for the Lomb-Scargle periodogram

    This is equal to 1 - cdf, where cdf is the cumulative distribution.
    The single-frequency false alarm probability should not be confused with
    the false alarm probability for the largest peak.

    Parameters
    ----------
    z : array-like
        The periodogram value.
    N : int
        The number of data points from which the periodogram was computed.
    normalization : {'standard', 'model', 'log', 'psd'}
        The periodogram normalization.
    dH, dK : int, optional
        The number of parameters in the null hypothesis and the model.

    Returns
    -------
    false_alarm_probability : np.ndarray
        The single-frequency false alarm probability.

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.
    All expressions used here are adapted from Table 1 of Baluev 2008 [1]_.

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    z = np.asarray(z)
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
        raise ValueError(f"normalization='{normalization}' is not recognized")


def inv_fap_single(fap, N, normalization, dH=1, dK=3):
    """Single-frequency inverse false alarm probability

    This function computes the periodogram value associated with the specified
    single-frequency false alarm probability. This should not be confused with
    the false alarm level of the largest peak.

    Parameters
    ----------
    fap : array-like
        The false alarm probability.
    N : int
        The number of data points from which the periodogram was computed.
    normalization : {'standard', 'model', 'log', 'psd'}
        The periodogram normalization.
    dH, dK : int, optional
        The number of parameters in the null hypothesis and the model.

    Returns
    -------
    z : np.ndarray
        The periodogram power corresponding to the single-peak false alarm
        probability.

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.
    All expressions used here are adapted from Table 1 of Baluev 2008 [1]_.

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    fap = np.asarray(fap)
    if dK - dH != 2:
        raise NotImplementedError("Degrees of freedom != 2")
    Nk = N - dK

    # No warnings for fap = 0; rather, just let it give the right infinity.
    with np.errstate(divide='ignore'):
        if normalization == 'psd':
            return -np.log(fap)
        elif normalization == 'standard':
            return 1 - fap ** (2 / Nk)
        elif normalization == 'model':
            return -1 + fap ** (-2 / Nk)
        elif normalization == 'log':
            return -2 / Nk * np.log(fap)
        else:
            raise ValueError(f"normalization='{normalization}' is not recognized")


def cdf_single(z, N, normalization, dH=1, dK=3):
    """Cumulative distribution for the Lomb-Scargle periodogram

    Compute the expected cumulative distribution of the periodogram
    for the null hypothesis - i.e. data consisting of Gaussian noise.

    Parameters
    ----------
    z : array-like
        The periodogram value.
    N : int
        The number of data points from which the periodogram was computed.
    normalization : {'standard', 'model', 'log', 'psd'}
        The periodogram normalization.
    dH, dK : int, optional
        The number of parameters in the null hypothesis and the model.

    Returns
    -------
    cdf : np.ndarray
        The expected cumulative distribution function.

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
    Teff = np.sqrt(4 * np.pi * Dt)  # Effective baseline
    W = fmax * Teff
    Z = np.asarray(Z)
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
        raise NotImplementedError(f"normalization={normalization}")


def fap_naive(Z, fmax, t, y, dy, normalization='standard'):
    """False Alarm Probability based on estimated number of indep frequencies"""
    N = len(t)
    T = max(t) - min(t)
    N_eff = fmax * T
    fap_s = fap_single(Z, N, normalization=normalization)
    # result is 1 - (1 - fap_s) ** N_eff
    # this is much more precise for small Z / large N
    # Ignore divide by zero no np.log1p - fine to let it return -inf.
    with np.errstate(divide='ignore'):
        return -np.expm1(N_eff * np.log1p(-fap_s))


def inv_fap_naive(fap, fmax, t, y, dy, normalization='standard'):
    """Inverse FAP based on estimated number of indep frequencies"""
    fap = np.asarray(fap)
    N = len(t)
    T = max(t) - min(t)
    N_eff = fmax * T
    # fap_s = 1 - (1 - fap) ** (1 / N_eff)
    # Ignore divide by zero no np.log - fine to let it return -inf.
    with np.errstate(divide='ignore'):
        fap_s = -np.expm1(np.log(1 - fap) / N_eff)
    return inv_fap_single(fap_s, N, normalization)


def fap_davies(Z, fmax, t, y, dy, normalization='standard'):
    """Davies upper-bound to the false alarm probability

    (Eqn 5 of Baluev 2008)
    """
    N = len(t)
    fap_s = fap_single(Z, N, normalization=normalization)
    tau = tau_davies(Z, fmax, t, y, dy, normalization=normalization)
    return fap_s + tau


@vectorize_first_argument
def inv_fap_davies(p, fmax, t, y, dy, normalization='standard'):
    """Inverse of the davies upper-bound"""
    from scipy import optimize
    args = (fmax, t, y, dy, normalization)
    z0 = inv_fap_naive(p, *args)
    func = lambda z, *args: fap_davies(z, *args) - p
    res = optimize.root(func, z0, args=args, method='lm')
    if not res.success:
        raise ValueError(f'inv_fap_baluev did not converge for p={p}')
    return res.x


def fap_baluev(Z, fmax, t, y, dy, normalization='standard'):
    """Alias-free approximation to false alarm probability

    (Eqn 6 of Baluev 2008)
    """
    fap_s = fap_single(Z, len(t), normalization)
    tau = tau_davies(Z, fmax, t, y, dy, normalization=normalization)
    # result is 1 - (1 - fap_s) * np.exp(-tau)
    # this is much more precise for small numbers
    return -np.expm1(-tau) + fap_s * np.exp(-tau)


@vectorize_first_argument
def inv_fap_baluev(p, fmax, t, y, dy, normalization='standard'):
    """Inverse of the Baluev alias-free approximation"""
    from scipy import optimize
    args = (fmax, t, y, dy, normalization)
    z0 = inv_fap_naive(p, *args)
    func = lambda z, *args: fap_baluev(z, *args) - p
    res = optimize.root(func, z0, args=args, method='lm')
    if not res.success:
        raise ValueError(f'inv_fap_baluev did not converge for p={p}')
    return res.x


def _bootstrap_max(t, y, dy, fmax, normalization, random_seed, n_bootstrap=1000):
    """Generate a sequence of bootstrap estimates of the max"""
    from .core import LombScargle
    rng = np.random.default_rng(random_seed)
    power_max = []
    for _ in range(n_bootstrap):
        s = rng.integers(0, len(y), len(y))  # sample with replacement
        ls_boot = LombScargle(t, y[s], dy if dy is None else dy[s],
                              normalization=normalization)
        freq, power = ls_boot.autopower(maximum_frequency=fmax)
        power_max.append(power.max())

    power_max = u.Quantity(power_max)
    power_max.sort()

    return power_max


def fap_bootstrap(Z, fmax, t, y, dy, normalization='standard',
                  n_bootstraps=1000, random_seed=None):
    """Bootstrap estimate of the false alarm probability"""
    pmax = _bootstrap_max(t, y, dy, fmax, normalization, random_seed,
                          n_bootstraps)

    return 1 - np.searchsorted(pmax, Z) / len(pmax)


def inv_fap_bootstrap(fap, fmax, t, y, dy, normalization='standard',
                      n_bootstraps=1000, random_seed=None):
    """Bootstrap estimate of the inverse false alarm probability"""
    fap = np.asarray(fap)
    pmax = _bootstrap_max(t, y, dy, fmax, normalization, random_seed,
                          n_bootstraps)

    return pmax[np.clip(np.floor((1 - fap) * len(pmax)).astype(int),
                        0, len(pmax) - 1)]


METHODS = {'single': fap_single,
           'naive': fap_naive,
           'davies': fap_davies,
           'baluev': fap_baluev,
           'bootstrap': fap_bootstrap}


def false_alarm_probability(Z, fmax, t, y, dy, normalization='standard',
                            method='baluev', method_kwds=None):
    """Compute the approximate false alarm probability for periodogram peaks Z

    This gives an estimate of the false alarm probability for the largest value
    in a periodogram, based on the null hypothesis of non-varying data with
    Gaussian noise. The true probability cannot be computed analytically, so
    each method available here is an approximation to the true value.

    Parameters
    ----------
    Z : array-like
        The periodogram value.
    fmax : float
        The maximum frequency of the periodogram.
    t, y, dy : array-like
        The data times, values, and errors.
    normalization : {'standard', 'model', 'log', 'psd'}, optional
        The periodogram normalization.
    method : {'baluev', 'davies', 'naive', 'bootstrap'}, optional
        The approximation method to use.
    method_kwds : dict, optional
        Additional method-specific keywords.

    Returns
    -------
    false_alarm_probability : np.ndarray
        The false alarm probability.

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.

    See Also
    --------
    false_alarm_level : compute the periodogram level for a particular fap

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    if method == 'single':
        return fap_single(Z, len(t), normalization)
    elif method not in METHODS:
        raise ValueError(f"Unrecognized method: {method}")
    method = METHODS[method]
    method_kwds = method_kwds or {}

    return method(Z, fmax, t, y, dy, normalization, **method_kwds)


INV_METHODS = {'single': inv_fap_single,
               'naive': inv_fap_naive,
               'davies': inv_fap_davies,
               'baluev': inv_fap_baluev,
               'bootstrap': inv_fap_bootstrap}


def false_alarm_level(p, fmax, t, y, dy, normalization,
                      method='baluev', method_kwds=None):
    """Compute the approximate periodogram level given a false alarm probability

    This gives an estimate of the periodogram level corresponding to a specified
    false alarm probability for the largest peak, assuming a null hypothesis
    of non-varying data with Gaussian noise. The true level cannot be computed
    analytically, so each method available here is an approximation to the true
    value.

    Parameters
    ----------
    p : array-like
        The false alarm probability (0 < p < 1).
    fmax : float
        The maximum frequency of the periodogram.
    t, y, dy : arrays
        The data times, values, and errors.
    normalization : {'standard', 'model', 'log', 'psd'}, optional
        The periodogram normalization.
    method : {'baluev', 'davies', 'naive', 'bootstrap'}, optional
        The approximation method to use.
    method_kwds : dict, optional
        Additional method-specific keywords.

    Returns
    -------
    z : np.ndarray
        The periodogram level.

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.

    See Also
    --------
    false_alarm_probability : compute the fap for a given periodogram level

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    if method == 'single':
        return inv_fap_single(p, len(t), normalization)
    elif method not in INV_METHODS:
        raise ValueError(f"Unrecognized method: {method}")
    method = INV_METHODS[method]
    method_kwds = method_kwds or {}

    return method(p, fmax, t, y, dy, normalization, **method_kwds)
