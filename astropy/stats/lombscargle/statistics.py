"""
Periodogram Statistics
======================

This module implements various periodogram statistics such as probability
distribution functions, cumulative distributions, and false alarm
probabilities.
"""
import numpy as np

from ... import units


def pdf(z, N, normalization, dH=1, dK=3):
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

    if isinstance(z, units.Quantity):
        if z.unit == units.dimensionless_unscaled:
            z = z.value
        else:
            raise ValueError('The distribution can be computed only for '
                             'normalized (unitless) periodograms.')

    if normalization == 'psd':
        return np.exp(-z)
    elif normalization == 'standard':
        return 0.5 * Nk * (1 + z) ** (-0.5 * Nk - 1)
    elif normalization == 'model':
        return 0.5 * Nk * (1 - z) ** (0.5 * Nk - 1)
    elif normalization == 'log':
        return 0.5 * Nk * np.exp(-0.5 * Nk * z)
    else:
        raise ValueError("normalization='{0}' is not recognized"
                         "".format(normalization))


def cdf(z, N, normalization, dH=1, dK=3):
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
    if dK - dH != 2:
        raise NotImplementedError("Degrees of freedom != 2")
    Nk = N - dK

    if isinstance(z, units.Quantity):
        if z.unit == units.dimensionless_unscaled:
            z = z.value
        else:
            raise ValueError('The distribution can be computed only for '
                             'normalized (unitless) periodograms.')

    if normalization == 'psd':
        return 1 - np.exp(-z)
    elif normalization == 'standard':
        return 1 - (1 + z) ** (-0.5 * Nk)
    elif normalization == 'model':
        return 1 - (1 - z) ** (0.5 * Nk)
    elif normalization == 'log':
        return 1 - np.exp(-0.5 * Nk * z)
    else:
        raise ValueError("normalization='{0}' is not recognized"
                         "".format(normalization))
