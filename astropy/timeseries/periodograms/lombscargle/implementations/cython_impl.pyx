#cython: language_level=3

import numpy as np

cimport cython
cimport numpy as np

np.import_array()

cdef extern from "math.h":
    double sin(double)
    double cos(double)
    double atan2(double, double)

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.intp
ctypedef np.intp_t ITYPE_t


def lombscargle_cython(t, y, dy, frequency, normalization='standard',
                       fit_mean=True, center_data=True):
    """Lomb-Scargle Periodogram

    This is a pure-python implementation of the original Lomb-Scargle formalism
    (e.g. [1]_, [2]_), with the addition of the floating mean (e.g. [3]_)

    Parameters
    ----------
    t, y, dy : array-like
        times, values, and errors of the data points. These should be
        broadcastable to the same shape. None should be `~astropy.units.Quantity`.
    frequency : array-like
        frequencies (not angular frequencies) at which to calculate periodogram
    normalization : str, optional
        Normalization to use for the periodogram.
        Options are 'standard', 'model', 'log', or 'psd'.
    fit_mean : bool, optional
        if True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if ``fit_mean = False``

    Returns
    -------
    power : array-like
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    References
    ----------
    .. [1] W. Press et al, Numerical Recipes in C (2002)
    .. [2] Scargle, J.D. 1982, ApJ 263:835-853
    .. [3] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    """
    if dy is None:
        dy = 1

    t, y, dy = np.broadcast_arrays(t, y, dy)
    t = np.asarray(t, dtype=DTYPE, order='C')
    y = np.asarray(y, dtype=DTYPE, order='C')
    dy = np.asarray(dy, dtype=DTYPE, order='C')
    frequency = np.asarray(frequency, dtype=DTYPE, order='C')

    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if frequency.ndim != 1:
        raise ValueError("frequency should be one-dimensional")

    PLS = np.zeros(frequency.shape, dtype=DTYPE, order='C')

    # pre-center the data: not technically required if fit_mean=True,
    # but it simplifies the math.
    if fit_mean or center_data:
        # compute MLE for mean in the presence of noise.
        w = dy ** -2
        y = y - np.dot(w, y) / np.sum(w)

    if fit_mean:
        _generalized_lomb_scargle(t, y, dy, 2 * np.pi * frequency, PLS)
    else:
        _standard_lomb_scargle(t, y, dy, 2 * np.pi * frequency, PLS)

    if normalization == 'standard':
        pass
    elif normalization == 'model':
        return PLS / (1 - PLS)
    elif normalization == 'log':
        return -np.log(1 - PLS)
    elif normalization == 'psd':
        w = dy ** -2
        PLS *= 0.5 * np.dot(w, y * y)
    else:
        raise ValueError("normalization='{0}' "
                         "not recognized".format(normalization))
    return PLS.ravel()


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef _standard_lomb_scargle(const DTYPE_t[::1] t, const DTYPE_t[::1] y, const DTYPE_t[::1] dy,
                            const DTYPE_t[::1] omega, DTYPE_t[::1] PLS):
    cdef ITYPE_t N_freq = omega.shape[0]
    cdef ITYPE_t N_obs = t.shape[0]

    cdef DTYPE_t w, omega_t, sin_omega_t, cos_omega_t
    cdef DTYPE_t S2, C2, tau, Y, wsum, YY, YCtau, YStau, CCtau, SStau

    for i in range(N_freq):
        # first pass: determine tau
        S2 = 0
        C2 = 0
        for j in range(N_obs):
            w = 1. / dy[j]
            w *= w

            omega_t = omega[i] * t[j]
            sin_omega_t = sin(omega_t)
            cos_omega_t = cos(omega_t)

            S2 += 2 * w * sin_omega_t * cos_omega_t
            C2 += w * (1 - 2 * sin_omega_t * sin_omega_t)

        tau = 0.5 * atan2(S2, C2) / omega[i]

        wsum = 0
        Y = 0
        YY = 0
        YCtau = 0
        YStau = 0
        CCtau = 0
        SStau = 0

        # second pass: compute the power
        for j in range(N_obs):
            w = 1. / dy[j]
            w *= w
            wsum += w

            omega_t = omega[i] * (t[j] - tau)
            sin_omega_t = sin(omega_t)
            cos_omega_t = cos(omega_t)

            Y += w * y[j]
            YY += w * y[j] * y[j]
            YCtau += w * y[j] * cos_omega_t
            YStau += w * y[j] * sin_omega_t
            CCtau += w * cos_omega_t * cos_omega_t
            SStau += w * sin_omega_t * sin_omega_t

        Y /= wsum
        YY /= wsum
        YCtau /= wsum
        YStau /= wsum
        CCtau /= wsum
        SStau /= wsum

        PLS[i] = (YCtau * YCtau / CCtau + YStau * YStau / SStau) / YY


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef _generalized_lomb_scargle(const DTYPE_t[::1] t, const DTYPE_t[::1] y, const DTYPE_t[::1] dy,
                               const DTYPE_t[::1] omega, DTYPE_t[::1] PLS):
    cdef ITYPE_t N_freq = omega.shape[0]
    cdef ITYPE_t N_obs = t.shape[0]

    cdef DTYPE_t w, omega_t, sin_omega_t, cos_omega_t
    cdef DTYPE_t S, C, S2, C2, tau, Y, wsum, YY
    cdef DTYPE_t Stau, Ctau, YCtau, YStau, CCtau, SStau

    for i in range(N_freq):
        # first pass: determine tau
        wsum = 0
        S = 0
        C = 0
        S2 = 0
        C2 = 0
        for j in range(N_obs):
            w = 1. / dy[j]
            w *= w
            wsum += w

            omega_t = omega[i] * t[j]
            sin_omega_t = sin(omega_t)
            cos_omega_t = cos(omega_t)

            S += w * sin_omega_t
            C += w * cos_omega_t

            S2 += 2 * w * sin_omega_t * cos_omega_t
            C2 += w - 2 * w * sin_omega_t * sin_omega_t

        S2 /= wsum
        C2 /= wsum
        S /= wsum
        C /= wsum

        S2 -= (2 * S * C)
        C2 -= (C * C - S * S)

        tau = 0.5 * atan2(S2, C2) / omega[i]

        Y = 0
        YY = 0
        Stau = 0
        Ctau = 0
        YCtau = 0
        YStau = 0
        CCtau = 0
        SStau = 0

        # second pass: compute the power
        for j in range(N_obs):
            w = 1. / dy[j]
            w *= w

            omega_t = omega[i] * (t[j] - tau)
            sin_omega_t = sin(omega_t)
            cos_omega_t = cos(omega_t)

            Y += w * y[j]
            YY += w * y[j] * y[j]
            Ctau += w * cos_omega_t
            Stau += w * sin_omega_t
            YCtau += w * y[j] * cos_omega_t
            YStau += w * y[j] * sin_omega_t
            CCtau += w * cos_omega_t * cos_omega_t
            SStau += w * sin_omega_t * sin_omega_t

        Y /= wsum
        YY /= wsum
        Ctau /= wsum
        Stau /= wsum
        YCtau /= wsum
        YStau /= wsum
        CCtau /= wsum
        SStau /= wsum

        YCtau -= Y * Ctau
        YStau -= Y * Stau
        CCtau -= Ctau * Ctau
        SStau -= Stau * Stau

        YY -= Y * Y

        PLS[i] = (YCtau * YCtau / CCtau + YStau * YStau / SStau) / YY
