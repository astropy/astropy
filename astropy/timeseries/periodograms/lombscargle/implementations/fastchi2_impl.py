from itertools import chain

import numpy as np

from .utils import trig_sum


def lombscargle_fastchi2(
    t,
    y,
    dy,
    f0,
    df,
    Nf,
    normalization="standard",
    fit_mean=True,
    center_data=True,
    nterms=1,
    use_fft=True,
    trig_sum_kwds=None,
):
    """Lomb-Scargle Periodogram.

    This implements a fast chi-squared periodogram using the algorithm
    outlined in [4]_. The result is identical to the standard Lomb-Scargle
    periodogram. The advantage of this algorithm is the
    ability to compute multiterm periodograms relatively quickly.

    Parameters
    ----------
    t, y, dy : array-like
        times, values, and errors of the data points. These should be
        broadcastable to the same shape. None should be `~astropy.units.Quantity`.
    f0, df, Nf : (float, float, int)
        parameters describing the frequency grid, f = f0 + df * arange(Nf).
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
    power : array-like
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    References
    ----------
    .. [1] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    .. [2] W. Press et al, Numerical Recipes in C (2002)
    .. [3] Scargle, J.D. ApJ 263:835-853 (1982)
    .. [4] Palmer, J. ApJ 695:496-502 (2009)
    """
    if nterms == 0 and not fit_mean:
        raise ValueError("Cannot have nterms = 0 without fitting bias")

    if dy is None:
        dy = 1

    # Validate and setup input data
    t, y, dy = np.broadcast_arrays(t, y, dy)
    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")

    # Validate and setup frequency grid
    if f0 < 0:
        raise ValueError("Frequencies must be positive")
    if df <= 0:
        raise ValueError("Frequency steps must be positive")
    if Nf <= 0:
        raise ValueError("Number of frequencies must be positive")

    w = dy**-2.0
    ws = np.sum(w)

    # if fit_mean is true, centering the data now simplifies the math below.
    if center_data or fit_mean:
        y = y - np.dot(w, y) / ws

    yw = y / dy
    chi2_ref = np.dot(yw, yw)

    kwargs = dict.copy(trig_sum_kwds or {})
    kwargs.update(f0=f0, df=df, use_fft=use_fft, N=Nf)

    # Here we build-up the matrices XTX and XTy using pre-computed
    # sums. The relevant identities are
    # 2 sin(mx) sin(nx) = cos(m-n)x - cos(m+n)x
    # 2 cos(mx) cos(nx) = cos(m-n)x + cos(m+n)x
    # 2 sin(mx) cos(nx) = sin(m-n)x + sin(m+n)x

    yws = np.sum(y * w)

    SCw = [(np.zeros(Nf), ws * np.ones(Nf))]
    SCw.extend(
        [trig_sum(t, w, freq_factor=i, **kwargs) for i in range(1, 2 * nterms + 1)]
    )
    Sw, Cw = zip(*SCw)

    SCyw = [(np.zeros(Nf), yws * np.ones(Nf))]
    SCyw.extend(
        [trig_sum(t, w * y, freq_factor=i, **kwargs) for i in range(1, nterms + 1)]
    )
    Syw, Cyw = zip(*SCyw)

    # Now create an indexing scheme so we can quickly
    # build-up matrices at each frequency
    order = [("C", 0)] if fit_mean else []
    order.extend(chain(*([("S", i), ("C", i)] for i in range(1, nterms + 1))))

    funcs = dict(
        S=lambda m, i: Syw[m][i],
        C=lambda m, i: Cyw[m][i],
        SS=lambda m, n, i: 0.5 * (Cw[abs(m - n)][i] - Cw[m + n][i]),
        CC=lambda m, n, i: 0.5 * (Cw[abs(m - n)][i] + Cw[m + n][i]),
        SC=lambda m, n, i: 0.5 * (np.sign(m - n) * Sw[abs(m - n)][i] + Sw[m + n][i]),
        CS=lambda m, n, i: 0.5 * (np.sign(n - m) * Sw[abs(n - m)][i] + Sw[n + m][i]),
    )

    def compute_power(i):
        XTX = np.array(
            [[funcs[A[0] + B[0]](A[1], B[1], i) for A in order] for B in order]
        )
        XTy = np.array([funcs[A[0]](A[1], i) for A in order])
        return np.dot(XTy.T, np.linalg.solve(XTX, XTy))

    p = np.array([compute_power(i) for i in range(Nf)])

    if normalization == "psd":
        p *= 0.5
    elif normalization == "standard":
        p /= chi2_ref
    elif normalization == "log":
        p = -np.log(1 - p / chi2_ref)
    elif normalization == "model":
        p /= chi2_ref - p
    else:
        raise ValueError(f"normalization='{normalization}' not recognized")
    return p
