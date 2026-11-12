from functools import cache

import numpy as np

from astropy.timeseries.periodograms.lombscargle.utils import (
    compute_chi2_ref,
    convert_normalization,
)
from astropy.utils import minversion
from astropy.utils.compat.optional_deps import HAS_SCIPY


@cache
def _scipy_supports_floating_mean():
    # The ``floating_mean`` and ``weights`` arguments were both added to
    # scipy.signal.lombscargle in scipy 1.15.
    # Deferred so that scipy is only imported and its version only probed when
    # the answer is needed; passing the module object lets minversion fall
    # back to scipy.__version__ when the distribution metadata is missing
    # (e.g. in some frozen applications).
    import scipy

    return minversion(scipy, "1.15.0.dev0")


def lombscargle_scipy(
    t,
    y,
    frequency,
    normalization="standard",
    center_data=True,
    *,
    fit_mean=False,
    dy=None,
):
    """Lomb-Scargle Periodogram.

    This is a wrapper of ``scipy.signal.lombscargle`` for computation of the
    Lomb-Scargle periodogram. This is a relatively fast version of the naive
    O[N^2] algorithm.

    Parameters
    ----------
    t, y : array-like
        Times and values of the data points. These should be
        broadcastable to the same shape. None should be `~astropy.units.Quantity`.
    frequency : array-like
        Frequencies (not angular frequencies) at which to calculate periodogram.
    normalization : str, optional
        Normalization to use for the periodogram.
        Options are 'standard', 'model', 'log', or 'psd'.
    center_data : bool, optional
        If True, pre-center the data by subtracting the weighted mean
        of the input data.
    fit_mean : bool, optional, keyword-only
        If True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage. Requires Scipy 1.15 and corresponds
        to the ``floating_mean`` argument in `scipy.signal.lombscargle`.
    dy : float or array-like, optional, keyword-only
        Errors of the data points, broadcastable to the shape of ``t`` and
        ``y``. Should not be a `~astropy.units.Quantity`. Errors that differ
        between data points require Scipy 1.15 and are passed as inverse
        variances to the ``weights`` argument in `scipy.signal.lombscargle`.

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

    if fit_mean and not _scipy_supports_floating_mean():
        raise ValueError("fit_mean=True requires scipy 1.15 or later")

    t, y, dy = np.broadcast_arrays(t, y, 1 if dy is None else dy)

    # Scipy requires floating-point input
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    dy = np.asarray(dy, dtype=float)
    frequency = np.asarray(frequency, dtype=float)

    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if frequency.ndim != 1:
        raise ValueError("frequency should be one-dimensional")

    w = dy**-2.0

    if center_data:
        y = y - np.dot(w, y) / w.sum()

    kwargs = {"floating_mean": True} if fit_mean else {}
    if np.any(dy != dy[0]):
        if not _scipy_supports_floating_mean():
            raise ValueError(
                "scipy method requires scipy 1.15 or later for non-uniform "
                "uncertainties dy"
            )
        kwargs["weights"] = w

    # Note: scipy `freqs` input is in angular frequencies, and scipy
    # normalizes the weights to have a mean of one
    p = signal.lombscargle(t, y, 2 * np.pi * frequency, **kwargs) * w.mean()

    if normalization == "psd":
        return p

    chi2_ref = compute_chi2_ref(y, dy, center_data=center_data, fit_mean=fit_mean)
    return convert_normalization(p, t.size, "psd", normalization, chi2_ref=chi2_ref)
