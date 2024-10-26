#cython: language_level=3

import numpy as np

cimport cython
cimport numpy as np

np.import_array()

cdef extern from "math.h":
    double sin(double)
    double cos(double)
    double sqrt(double)

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.intp
ctypedef np.intp_t ITYPE_t


def lombscargle_cython(t, y, dy, frequency, normalization='standard',
                       fit_mean=True, center_data=True, assume_regular_frequency=False):
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
    assume_regular_frequency : bool, optional
        if True, allocate the memory buffers to accelerate the computation

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
    w = dy ** -2
    w /= np.sum(w)
    if fit_mean or center_data:
        # compute MLE for mean in the presence of noise.
        y = y - np.dot(w, y) / np.sum(w)

    _lomb_scargle(t, y, w, 2 * np.pi * frequency, PLS, fit_mean = fit_mean, assume_regular_frequency = assume_regular_frequency)

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

from libc.stdlib cimport malloc, free
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)

cdef _lomb_scargle(const DTYPE_t[::1] t, const DTYPE_t[::1] y, const DTYPE_t[::1] w,
                            const DTYPE_t[::1] omega, DTYPE_t[::1] PLS, const bint fit_mean, const bint assume_regular_frequency):

    cdef int i, j
    cdef ITYPE_t N_freq = omega.shape[0]
    cdef ITYPE_t N_obs = t.shape[0]

    cdef DTYPE_t S, C, Sh, Ch, Sw, Cw, S2, C2, S2w, C2w, YC, YS, YY, tan_2omega_tau

    # allocate buffers used to speed-up the computation for regular grids
    cdef DTYPE_t* sin_buffer = <DTYPE_t*>malloc(N_obs * sizeof(DTYPE_t))
    cdef DTYPE_t* cos_buffer = <DTYPE_t*>malloc(N_obs * sizeof(DTYPE_t))
    cdef DTYPE_t* dsin_buffer = NULL
    cdef DTYPE_t* dcos_buffer = NULL
    if assume_regular_frequency:
        dsin_buffer = <DTYPE_t*>malloc(N_obs * sizeof(DTYPE_t))
        dcos_buffer = <DTYPE_t*>malloc(N_obs * sizeof(DTYPE_t))

    if assume_regular_frequency:
        #init the buffers to enable recursion
        for j in range(N_obs):
            tmp = t[j] * (omega[N_freq-1] - omega[0]) / (N_freq - 1)
            dsin_buffer[j] = sin(tmp)
            dcos_buffer[j] = cos(tmp)

    YY = 0

    for j in range(N_obs):
        YY += y[j] * y[j] * w[j]

    for i in range(N_freq):
        S = C = Sh = Ch = Sw = Cw = S2 = C2 = S2w = C2w = YS = YC = 0


        if not assume_regular_frequency or i % 64 == 0 :
             for j in range(N_obs):
                sin_buffer[j] = sin(omega[i] * t[j])
                cos_buffer[j] = cos(omega[i] * t[j])
        elif assume_regular_frequency :
            # exp(x + dx) = exp(x) * exp(dx)
            for j in range(N_obs):
                tmp = cos_buffer[j] * dcos_buffer[j] - sin_buffer[j] * dsin_buffer[j]
                sin_buffer[j] = cos_buffer[j] * dsin_buffer[j] + sin_buffer[j] * dcos_buffer[j]
                cos_buffer[j] = tmp

        for j in range(N_obs):

            Sh += y[j] * w[j] * sin_buffer[j]
            Ch += y[j] * w[j] * cos_buffer[j]
            S2 += w[j] * 2 * sin_buffer[j] * cos_buffer[j]
            C2 += w[j] * (1 - (2 * sin_buffer[j] * sin_buffer[j]))

        if fit_mean:
            for j in range(N_obs):
                S += w[j] * sin_buffer[j]
                C += w[j] * cos_buffer[j]

        if fit_mean:
            tan_2omega_tau = (S2 - (2 * S * C)) / (C2 - ((C * C) - (S * S)))
        else:
            tan_2omega_tau = S2 / C2

        S2w = tan_2omega_tau / sqrt(1 + tan_2omega_tau * tan_2omega_tau)
        C2w = 1 / sqrt(1 + tan_2omega_tau * tan_2omega_tau)

        Cw = sqrt(0.5 * (1 + C2w))
        Sw = sqrt(0.5 * (1 - C2w))
        if S2w < 0:
            Sw = - Sw

        YC = Ch * Cw + Sh * Sw
        YS = Sh * Cw - Ch * Sw

        CC = 0.5 * (1 + C2 * C2w + S2 * S2w)
        SS = 0.5 * (1 - C2 * C2w - S2 * S2w)

        if fit_mean:
            CC -= (C * Cw + S * Sw) * (C * Cw + S * Sw)
            SS -= (S * Cw - C * Sw) * (S * Cw - C * Sw)

        PLS[i] = ((YC * YC / CC) + (YS * YS / SS)) / YY

    free(sin_buffer)
    free(cos_buffer)
    if assume_regular_frequency:
        free(dsin_buffer)
        free(dcos_buffer)
