from itertools import chain

import numpy as np

from .mle import design_matrix
from .utils import trig_sum


def lombscargle_fastchi2(  # ruff: ignore[PLR0917]
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
    *,
    algorithm="lra",
):
    """Lomb-Scargle Periodogram.

    This implements a fast chi-squared periodogram using the algorithm
    outlined in [4]_. The result is identical to the standard Lomb-Scargle
    periodogram. The advantage of this algorithm is the
    ability to compute multiterm periodograms relatively quickly.

    Parameters
    ----------
    t, y, dy : array-like
        Times, values, and errors of the data points. These should be
        broadcastable to the same shape. None should be `~astropy.units.Quantity`.
    f0, df, Nf : (float, float, int)
        Parameters describing the frequency grid, f = f0 + df * arange(Nf).
    normalization : str, optional
        Normalization to use for the periodogram.
        Options are 'standard', 'model', 'log', or 'psd'.
    fit_mean : bool, optional
        If True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool, optional
        If True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if ``fit_mean = False``.
    nterms : int, optional
        Number of Fourier terms in the fit.
    use_fft : bool (default=True)
        If True, then use an approximate FFT-based O[NlogN] algorithm to
        compute the trigonometric sums. Otherwise, use a slower O[N^2]
        algorithm.
    trig_sum_kwds : dict or None, optional
        Extra keyword arguments to pass to the ``trig_sum`` utility.
        Options are ``oversampling``, ``Mfft`` and ``eps``. See documentation
        of ``trig_sum`` for details.
    algorithm : 'lra' (default), or 'fasper'
        This option is ignored if use_fft is False.
        Specify the approximation used to approximate the NUDFT of type 1. If the value is not valid falls back to the default option.
        Supported options are:

        - 'fasper': use Press & Rybicki's piecewise Lagrange polynomial extirpolation.
        - 'lra': Use the more accurate (but slower) Low Rank Approximation by Ruiz-Antolin and Townsend. This is the default option.

    Returns
    -------
    power : array-like
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    Notes
    -----
    The FFT trigonometric sums are approximate. At frequencies where the
    normal matrix is ill-conditioned, the power is instead computed by a
    direct singular-value decomposition of the weighted design matrix.
    Numerically unresolved singular vectors are discarded using the same
    default rank threshold as `numpy.linalg.lstsq`. This can increase the
    cost for sampling patterns with poor phase coverage.

    References
    ----------
    .. [1] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    .. [2] W. Press et al, Numerical Recipes in C (2002)
    .. [3] Scargle, J.D. ApJ 263:835-853 (1982)
    .. [4] Palmer, J. ApJ 695:496-502 (2009)
    .. [5] Ruiz-Antolin, D. and Townsend, A. "A nonuniform fast Fourier transform based on low rank approximation". SIAM 40.1 (2018)

    .. versionchanged:: 8.0
        The default algorithm has been changed from 'fasper' to 'lra'
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
    kwargs.update(f0=f0, df=df, use_fft=use_fft, N=Nf, algorithm=algorithm)

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

    # Assemble and solve several frequencies at a time. Limit the temporary
    # normal matrices to 2**16 entries, independently of the frequency grid.
    n_parameters = len(order)
    block_size = max(1, 2**16 // n_parameters**2)
    # Avoid amplifying the sum errors enough to lose half their precision.
    tolerance = np.finfo(float).eps
    if use_fft and algorithm == "lra":
        tolerance = max(tolerance, kwargs.get("eps", 5e-13))
    tolerance = np.sqrt(tolerance)
    p = np.empty(Nf)
    for start in range(0, Nf, block_size):
        stop = min(start + block_size, Nf)
        i = slice(start, stop)
        XTX = np.empty((stop - start, n_parameters, n_parameters))
        for a, A in enumerate(order):
            for b, B in enumerate(order[: a + 1]):
                XTX[:, a, b] = XTX[:, b, a] = funcs[A[0] + B[0]](A[1], B[1], i)
        XTy = np.stack([funcs[A[0]](A[1], i) for A in order], axis=-1)

        # Incomplete phase coverage can make the normal equations nearly
        # singular. Their construction subtracts trigonometric sums, so even
        # small FFT errors can overwhelm the smallest eigenvalue. First use
        # Gershgorin bounds, which cheaply certify well-conditioned matrices;
        # compute eigenvalues only when those bounds are inconclusive.
        diagonal = XTX.diagonal(axis1=1, axis2=2)
        radius = np.abs(XTX).sum(axis=2) - np.abs(diagonal)
        lower = (diagonal - radius).min(axis=1)
        upper = (diagonal + radius).max(axis=1)
        uncertain = lower <= tolerance * upper
        unstable = uncertain.copy()
        if np.any(uncertain):
            eigenvalues = np.linalg.eigvalsh(XTX[uncertain])
            unstable[uncertain] = eigenvalues[:, 0] <= tolerance * eigenvalues[:, -1]

        stable = ~unstable
        if np.any(stable):
            rhs = XTy[stable]
            solution = np.linalg.solve(XTX[stable], rhs[..., np.newaxis])[..., 0]
            p[start:stop][stable] = np.einsum("ij,ij->i", rhs, solution)

        # Recompute only the affected frequencies from the weighted design
        # matrix. Projection onto its resolved left singular vectors avoids
        # both squaring the condition number and cancellation in X @ theta.
        for j in np.flatnonzero(unstable):
            frequency = f0 + (start + j) * df
            X = design_matrix(t, frequency, dy=dy, bias=fit_mean, nterms=nterms)
            U, singular_values, _ = np.linalg.svd(X, full_matrices=False)
            cutoff = singular_values[0] * np.finfo(float).eps * max(X.shape)
            projected = U[:, singular_values > cutoff].T @ yw
            p[start + j] = np.dot(projected, projected)

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
