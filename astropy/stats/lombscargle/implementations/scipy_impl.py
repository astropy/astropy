from __future__ import print_function, division

import numpy as np

try:
    from scipy import signal
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def lombscargle_scipy(t, y, frequency, normalization='normalized',
                      center_data=True):
    """Lomb-Scargle Periodogram

    This is a wrapper of ``scipy.signal.lombscargle`` for computation of the
    Lomb-Scargle periodogram. This is a relatively fast version of the naive
    O[N^2] algorithm, but cannot handle heteroskedastic errors.

    Parameters
    ----------
    t, y: array_like  (NOT astropy.Quantities)
        times, values, and errors of the data points. These should be
        broadcastable to the same shape.
    frequency : array_like
        frequencies (not angular frequencies) at which to calculate periodogram
    normalization : string (optional, default='normalized')
        Normalization to use for the periodogram
    center_data : bool (optional, default=True)
        if True, pre-center the data by subtracting the weighted mean
        of the input data.

    Returns
    -------
    power : array_like
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    References
    ----------
    .. [1] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    .. [2] W. Press et al, Numerical Recipies in C (2002)
    .. [3] Scargle, J.D. 1982, ApJ 263:835-853
    """
    if not HAS_SCIPY:
        raise ValueError("scipy must be installed to use lombscargle_scipy")

    t, y = np.broadcast_arrays(t, y)
    frequency = np.asarray(frequency)

    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if frequency.ndim != 1:
        raise ValueError("frequency should be one-dimensional")

    if center_data:
        y = y - y.mean()

    # Note: scipy input accepts angular frequencies
    p = signal.lombscargle(t, y, 2 * np.pi * frequency)

    if normalization == 'unnormalized':
        pass
    elif normalization == 'normalized':
        p *= 2 / (t.size * np.mean(y ** 2))
    else:
        raise ValueError("normalization='{0}' "
                         "not recognized".format(normalization))
    return p
