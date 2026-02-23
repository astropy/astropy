import numpy as np

from astropy.utils import minversion
from astropy.utils.compat.optional_deps import HAS_SCIPY

SCIPY_LT_1_15 = not minversion("scipy", "1.15.0") if HAS_SCIPY else True


def lombscargle_scipy(
    t, y, frequency, normalization="standard", center_data=True, *, fit_mean=True
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
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data.
    fit_mean : bool (default: True)
        if True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
        Ignored for scipy versions older than 1.15

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
        kwargs = {}
    else:
        kwargs = {"floating_mean": fit_mean}

    # Note: scipy input accepts angular frequencies
    p = signal.lombscargle(t, y, 2 * np.pi * frequency, **kwargs)

    if normalization not in ("psd", "standard", "log", "model", "psd"):
        raise ValueError(f"{normalization=!r} not recognized")

    if normalization == "psd":
        return p

    if center_data or not kwargs.get("floating_mean", False):
        chi2_ref = 0.5 * t.size * np.mean(y**2)
    else:
        chi2_ref = 0.5 * t.size * np.mean((y - y.mean()) ** 2)

    if normalization == "standard":
        p /= chi2_ref
    elif normalization == "log":
        p = -np.log(1 - p / chi2_ref)
    elif normalization == "model":
        p /= chi2_ref - p

    return p
