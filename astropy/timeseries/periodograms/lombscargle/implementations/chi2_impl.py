
import numpy as np

from .mle import design_matrix


def lombscargle_chi2(t, y, dy, frequency, normalization='standard',
                     fit_mean=True, center_data=True, nterms=1):
    """Lomb-Scargle Periodogram

    This implements a chi-squared-based periodogram, which is relatively slow
    but useful for validating the faster algorithms in the package.

    Parameters
    ----------
    t, y, dy : array_like (NOT astropy.Quantities)
        times, values, and errors of the data points. These should be
        broadcastable to the same shape.
    frequency : array_like
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
    nterms : int, optional
        Number of Fourier terms in the fit

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
    if dy is None:
        dy = 1

    t, y, dy = np.broadcast_arrays(t, y, dy)
    frequency = np.asarray(frequency)

    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if frequency.ndim != 1:
        raise ValueError("frequency should be one-dimensional")

    w = dy ** -2.0
    w /= w.sum()

    # if fit_mean is true, centering the data now simplifies the math below.
    if center_data or fit_mean:
        yw = (y - np.dot(w, y)) / dy
    else:
        yw = y / dy
    chi2_ref = np.dot(yw, yw)

    # compute the unnormalized model chi2 at each frequency
    def compute_power(f):
        X = design_matrix(t, f, dy=dy, bias=fit_mean, nterms=nterms)
        XTX = np.dot(X.T, X)
        XTy = np.dot(X.T, yw)
        return np.dot(XTy.T, np.linalg.solve(XTX, XTy))

    p = np.array([compute_power(f) for f in frequency])

    if normalization == 'psd':
        p *= 0.5
    elif normalization == 'model':
        p /= (chi2_ref - p)
    elif normalization == 'log':
        p = -np.log(1 - p / chi2_ref)
    elif normalization == 'standard':
        p /= chi2_ref
    else:
        raise ValueError("normalization='{}' "
                         "not recognized".format(normalization))
    return p
