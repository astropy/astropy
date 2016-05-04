from __future__ import print_function, division

import numpy as np

from .mle import design_matrix


def lombscargle_chi2(t, y, dy, frequency, normalization='normalized',
                     fit_bias=True, center_data=True, nterms=1):
    """Lomb-Scargle Periodogram

    This implements a chi-squared-based periodogram, which is relatively slow
    but useful for validating the faster algorithms in the package.

    Parameters
    ----------
    t, y, dy : array_like
        times, values, and errors of the data points. These should be
        broadcastable to the same shape.
    frequency : array_like
        frequencies (not angular frequencies) at which to calculate periodogram
    normalization : string (optional, default='normalized')
        Normalization to use for the periodogram
        TODO: figure out what options to use
    fit_bias : bool (optional, default=True)
        if True, include a constant offet as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool (optional, default=True)
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if ``fit_bias = False``
    nterms : int (optional, default=1)
        Number of Fourier terms in the fit

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
    if dy is None:
        dy = 1

    t, y, dy = np.broadcast_arrays(t, y, dy)
    frequency = np.asarray(frequency)
    assert t.ndim == 1
    assert frequency.ndim == 1

    w = dy ** -2.0
    w /= w.sum()

    # if fit_bias is true, centering the data now simplifies the math below.
    if center_data or fit_bias:
        yw = (y - np.dot(w, y)) / dy
    else:
        yw = y / dy
    chi2_ref = np.dot(yw, yw)

    # compute the unnormalized model chi2 at each frequency
    def compute_power(f):
        X = design_matrix(t, f, dy=dy, bias=fit_bias, nterms=nterms)
        XTX = np.dot(X.T, X)
        XTy = np.dot(X.T, yw)
        return np.dot(XTy.T, np.linalg.solve(XTX, XTy))

    p = np.array([compute_power(f) for f in frequency])

    if normalization == 'unnormalized':
        p *= 0.5 * t.size / (dy ** -2).sum()
    elif normalization == 'normalized':
        p /= chi2_ref
    else:
        raise ValueError("normalization='{0}' "
                         "not recognized".format(normalization))
    return p
