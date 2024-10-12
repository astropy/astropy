import numpy as np

from astropy.utils.compat.optional_deps import HAS_SCIPY

from .utils import SCIPY_LT_1_15


def lombscargle_scipy(
        t, y, frequency, normalization="standard", fit_mean=False, center_data=True
):
    """Lomb-Scargle Periodogram.

    This is a wrapper of ``scipy.signal.lombscargle`` for computation of the
    Lomb-Scargle periodogram. This is a relatively fast version of the naive
    O[N^2] algorithm, but cannot handle heteroskedastic errors.

    Parameters
    ----------
    t, y : array-like
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
        case of incomplete phase coverage. Requires Scipy 1.15 and corresponds
        to the ``floating_mean`` argument in `scipy.signal.lombscargle`.
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data.

    Returns
    -------
    power : array-like
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    References
    ----------
    .. [1] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    .. [2] W. Press et al, Numerical Recipes in C (2002)
    .. [3] Scargle, J.D. 1982, ApJ 263:835-853
    """
    if not HAS_SCIPY:
        raise ModuleNotFoundError("scipy must be installed to use lombscargle_scipy")

    from scipy import signal

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

    if SCIPY_LT_1_15:
        if fit_mean:
            raise NotImplementedError("`fit_mean=True` requires Scipy 1.15+")
        else:
            kwargs = {}
    else:
        kwargs = {"floating_mean": fit_mean}

    # Note: scipy `freqs` input is in angular frequencies
    p = signal.lombscargle(t, y, 2 * np.pi * frequency, **kwargs)

    if normalization == "psd":
        pass
    elif normalization == "standard":
        p *= 2 / (t.size * np.mean(y**2))
    elif normalization == "log":
        p = -np.log(1 - 2 * p / (t.size * np.mean(y**2)))
    elif normalization == "model":
        p /= 0.5 * t.size * np.mean(y**2) - p
    else:
        raise ValueError(f"normalization='{normalization}' not recognized")
    return p
