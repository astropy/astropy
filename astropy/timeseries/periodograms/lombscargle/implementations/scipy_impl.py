
import numpy as np


def lombscargle_scipy(t, y, frequency, normalization='standard',
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
    normalization : string (optional, default='standard')
        Normalization to use for the periodogram.
        Options are 'standard', 'model', 'log', or 'psd'.
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
    .. [2] W. Press et al, Numerical Recipes in C (2002)
    .. [3] Scargle, J.D. 1982, ApJ 263:835-853
    """
    try:
        from scipy import signal
    except ImportError:
        raise ImportError("scipy must be installed to use lombscargle_scipy")

    t, y = np.broadcast_arrays(t, y)

    # Scipy requires floating-point input
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    frequency = np.asarray(frequency, dtype=float)

    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if frequency.ndim != 1:
        raise ValueError("frequency should be one-dimensional")

    if center_data:
        y = y - y.mean()

    # Note: scipy input accepts angular frequencies
    p = signal.lombscargle(t, y, 2 * np.pi * frequency)

    if normalization == 'psd':
        pass
    elif normalization == 'standard':
        p *= 2 / (t.size * np.mean(y ** 2))
    elif normalization == 'log':
        p = -np.log(1 - 2 * p / (t.size * np.mean(y ** 2)))
    elif normalization == 'model':
        p /= 0.5 * t.size * np.mean(y ** 2) - p
    else:
        raise ValueError("normalization='{0}' "
                         "not recognized".format(normalization))
    return p
