import numpy as np

from .utils import DEFAULT_EPS, LRA_EPS_FLOOR, trig_sum

# The power is divided by CC and SS, and both are obtained by subtracting
# terms of order one. The error of the trigonometric sums therefore enters the
# periodogram divided by the size of that difference, so wherever the
# difference is small the sums have to be evaluated more accurately than the
# periodogram itself. This is the absolute error of the power that the
# approximation of the sums is allowed to contribute.
SUM_ERROR_TARGET = 1e-8


def _tau_terms(t, w, y, kwargs, fit_mean):
    """Trigonometric sums needed for the power and for the time-shift tau.

    Returns
    -------
    Sh, Ch : ndarray
        Sine and cosine sums of the data.
    S2w, C2w, Cw, Sw : ndarray
        Sine and cosine of ``2 * omega * tau`` and of ``omega * tau``.
    CC, SS : ndarray
        Diagonal terms of the normal equations of the fit.
    """
    Sh, Ch = trig_sum(t, w * y, **kwargs)
    S2, C2 = trig_sum(t, w, freq_factor=2, **kwargs)

    if fit_mean:
        S, C = trig_sum(t, w, **kwargs)
        tan_2omega_tau = (S2 - 2 * S * C) / (C2 - (C * C - S * S))
    else:
        tan_2omega_tau = S2 / C2

    # This is what we're computing below; the straightforward way is slower
    # and less stable, so we use trig identities instead
    #
    # omega_tau = 0.5 * np.arctan(tan_2omega_tau)
    # S2w, C2w = np.sin(2 * omega_tau), np.cos(2 * omega_tau)
    # Sw, Cw = np.sin(omega_tau), np.cos(omega_tau)

    S2w = tan_2omega_tau / np.sqrt(1 + tan_2omega_tau * tan_2omega_tau)
    C2w = 1 / np.sqrt(1 + tan_2omega_tau * tan_2omega_tau)
    Cw = np.sqrt(0.5) * np.sqrt(1 + C2w)
    Sw = np.sqrt(0.5) * np.sign(S2w) * np.sqrt(1 - C2w)

    CC = 0.5 * (1 + C2 * C2w + S2 * S2w)
    SS = 0.5 * (1 - C2 * C2w - S2 * S2w)

    if fit_mean:
        CC = CC - (C * Cw + S * Sw) ** 2
        SS = SS - (S * Cw - C * Sw) ** 2

    return Sh, Ch, S2w, C2w, Cw, Sw, CC, SS


def lombscargle_fast(
    t,
    y,
    dy,
    f0,
    df,
    Nf,
    center_data=True,
    fit_mean=True,
    normalization="standard",
    use_fft=True,
    trig_sum_kwds=None,
    *,
    algorithm="lra",
):
    """Fast Lomb-Scargle Periodogram.

    This implements the Press & Rybicki [1]_ or Low Rank Approximation [4]_ method for fast O[N log(N)]
    Lomb-Scargle periodograms.

    Parameters
    ----------
    t, y, dy : array-like
        times, values, and errors of the data points. These should be
        broadcastable to the same shape. None should be `~astropy.units.Quantity`.
    f0, df, Nf : (float, float, int)
        parameters describing the frequency grid, f = f0 + df * arange(Nf).
    center_data : bool (default=True)
        Specify whether to subtract the mean of the data before the fit
    fit_mean : bool (default=True)
        If True, then compute the floating-mean periodogram; i.e. let the mean
        vary with the fit.
    normalization : str, optional
        Normalization to use for the periodogram.
        Options are 'standard', 'model', 'log', or 'psd'.
    use_fft : bool (default=True)
        If True, then use the Press & Rybicki O[NlogN] algorithm to compute
        the result. Otherwise, use a slower O[N^2] algorithm
    trig_sum_kwds : dict or None, optional
        extra keyword arguments to pass to the ``trig_sum`` utility.
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
    power : ndarray
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    Notes
    -----
    Note that the ``use_fft=True`` algorithm is an approximation to the true
    Lomb-Scargle periodogram, and as the number of points grows this
    approximation improves. On the other hand, for very small datasets
    (<~50 points or so) this approximation may not be useful.

    References
    ----------
    .. [1] Press W.H. and Rybicki, G.B, "Fast algorithm for spectral analysis
        of unevenly sampled data". ApJ 1:338, p277, 1989
    .. [2] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    .. [3] W. Press et al, Numerical Recipes in C (2002)
    .. [4] Ruiz-Antolin, D. and Townsend, A. "A nonuniform fast Fourier transform based on low rank approximation". SIAM 40.1 (2018)

    .. versionchanged:: 8.0
        The default algorithm has been changed from 'fasper' to 'lra'
    """
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
    w /= w.sum()

    # Center the data. Even if we're fitting the offset,
    # this step makes the expressions below more succinct
    if center_data or fit_mean:
        y = y - np.dot(w, y)

    # set up arguments to trig_sum
    kwargs = dict.copy(trig_sum_kwds or {})
    kwargs.update(f0=f0, df=df, use_fft=use_fft, N=Nf, algorithm=algorithm)

    # ----------------------------------------------------------------------
    # 1. compute functions of the time-shift tau at each frequency
    Sh, Ch, S2w, C2w, Cw, Sw, CC, SS = _tau_terms(t, w, y, kwargs, fit_mean)

    # CC and SS lose significant digits wherever they are much smaller than
    # the sums they are formed from, which happens when the fitted sinusoid is
    # nearly degenerate with the offset, as it is for tightly clustered
    # sampling. The approximation error of the sums is then amplified by the
    # division below and dominates the periodogram, so ask for a tolerance
    # matched to the conditioning. This only ever tightens the tolerance, hence
    # it is a single extra pass at worst and leaves well conditioned data
    # untouched.
    if use_fft and algorithm == "lra":
        wanted = SUM_ERROR_TARGET * min(np.min(np.abs(CC)), np.min(np.abs(SS)))
        eps = (trig_sum_kwds or {}).get("eps", DEFAULT_EPS)
        if wanted < eps:
            kwargs["eps"] = max(wanted, LRA_EPS_FLOOR)
            Sh, Ch, S2w, C2w, Cw, Sw, CC, SS = _tau_terms(t, w, y, kwargs, fit_mean)

    # ----------------------------------------------------------------------
    # 2. Compute the periodogram, following Zechmeister & Kurster
    #    and using tricks from Press & Rybicki.
    YY = np.dot(w, y**2)
    YC = Ch * Cw + Sh * Sw
    YS = Sh * Cw - Ch * Sw

    power = YC * YC / CC + YS * YS / SS

    if normalization == "standard":
        power /= YY
    elif normalization == "model":
        power /= YY - power
    elif normalization == "log":
        power = -np.log(1 - power / YY)
    elif normalization == "psd":
        power *= 0.5 * (dy**-2.0).sum()
    else:
        raise ValueError(f"normalization='{normalization}' not recognized")

    return power
