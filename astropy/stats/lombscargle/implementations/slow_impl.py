from __future__ import print_function, division

import numpy as np


def lombscargle_slow(t, y, dy, frequency, normalization='normalized',
                     fit_bias=True, center_data=True):
    """Lomb-Scargle Periodogram

    This is a pure-python implemnetation of the original Lomb-Scargle formalism
    (e.g. [1]_, [2]_), with the addition of the floating mean (e.g. [3]_)

    Parameters
    ----------
    t, y, dy : array_like
        times, values, and errors of the data points. These should be
        broadcastable to the same shape.
    frequency : array_like
        frequencies (not angular frequencies) at which to calculate periodogram
    normalization : string (optional, default='normalized')
        Normalization to use for the periodogram
    fit_bias : bool (optional, default=True)
        if True, include a constant offet as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool (optional, default=True)
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if ``fit_bias = False``

    Returns
    -------
    power : array_like
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    References
    ----------
    .. [1] W. Press et al, Numerical Recipies in C (2002)
    .. [2] Scargle, J.D. 1982, ApJ 263:835-853
    .. [3] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    """
    if dy is None:
        dy = 1

    t, y, dy = np.broadcast_arrays(t, y, dy)
    frequency = np.asarray(frequency)
    assert t.ndim == 1
    assert frequency.ndim == 1

    w = dy ** -2.0
    w /= w.sum()

    # if fit_bias is true, centering the data now simplifies the math below.
    if fit_bias or center_data:
        y = y - np.dot(w, y)

    omega = 2 * np.pi * frequency
    omega = omega.ravel()[np.newaxis, :]

    # make following arrays into column vectors
    t, y, dy, w = map(lambda x: x[:, np.newaxis], (t, y, dy, w))

    sin_omega_t = np.sin(omega * t)
    cos_omega_t = np.cos(omega * t)

    # compute time-shift tau
    # S2 = np.dot(w.T, np.sin(2 * omega * t)
    S2 = 2 * np.dot(w.T, sin_omega_t * cos_omega_t)
    # C2 = np.dot(w.T, np.cos(2 * omega * t)
    C2 = 2 * np.dot(w.T, 0.5 - sin_omega_t ** 2)

    if fit_bias:
        S = np.dot(w.T, sin_omega_t)
        C = np.dot(w.T, cos_omega_t)

        S2 -= (2 * S * C)
        C2 -= (C * C - S * S)

    # compute components needed for the fit
    tan_2omega_tau = S2 / C2
    omega_t_tau = omega * t - 0.5 * np.arctan(tan_2omega_tau)

    sin_omega_t_tau = np.sin(omega_t_tau)
    cos_omega_t_tau = np.cos(omega_t_tau)

    Y = np.dot(w.T, y)

    wy = w * y

    YCtau = np.dot(wy.T, cos_omega_t_tau)
    YStau = np.dot(wy.T, sin_omega_t_tau)
    CCtau = np.dot(w.T, cos_omega_t_tau * cos_omega_t_tau)
    SStau = np.dot(w.T, sin_omega_t_tau * sin_omega_t_tau)

    if fit_bias:
        Ctau = np.dot(w.T, cos_omega_t_tau)
        Stau = np.dot(w.T, sin_omega_t_tau)

        YCtau -= Y * Ctau
        YStau -= Y * Stau
        CCtau -= Ctau * Ctau
        SStau -= Stau * Stau

    p = (YCtau * YCtau / CCtau + YStau * YStau / SStau)

    if normalization == 'normalized':
        p /= np.dot(w.T, y * y)
    elif normalization == 'unnormalized':
        p *= 0.5 * t.size
    else:
        raise ValueError("normalization='{0}' "
                         "not recognized".format(normalization))
    return p.ravel()
