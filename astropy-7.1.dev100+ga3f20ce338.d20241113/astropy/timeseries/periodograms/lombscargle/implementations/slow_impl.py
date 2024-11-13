import numpy as np


def lombscargle_slow(
    t, y, dy, frequency, normalization="standard", fit_mean=True, center_data=True
):
    """Lomb-Scargle Periodogram.

    This is a pure-python implementation of the original Lomb-Scargle formalism
    (e.g. [1]_, [2]_), with the addition of the floating mean (e.g. [3]_)

    Parameters
    ----------
    t, y, dy : array-like
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
        case of incomplete phase coverage.
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if ``fit_mean = False``

    Returns
    -------
    power : array-like
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    References
    ----------
    .. [1] W. Press et al, Numerical Recipes in C (2002)
    .. [2] Scargle, J.D. 1982, ApJ 263:835-853
    .. [3] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    """
    if dy is None:
        dy = 1

    t, y, dy = np.broadcast_arrays(t, y, dy)
    frequency = np.asarray(frequency)

    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if frequency.ndim != 1:
        raise ValueError("frequency should be one-dimensional")

    w = dy**-2.0
    w /= w.sum()

    # if fit_mean is true, centering the data now simplifies the math below.
    if fit_mean or center_data:
        y = y - np.dot(w, y)

    omega = 2 * np.pi * frequency
    omega = omega.ravel()[np.newaxis, :]

    # make following arrays into column vectors
    t, y, dy, w = (x[:, np.newaxis] for x in (t, y, dy, w))

    sin_omega_t = np.sin(omega * t)
    cos_omega_t = np.cos(omega * t)

    # compute time-shift tau
    # S2 = np.dot(w.T, np.sin(2 * omega * t)
    S2 = 2 * np.dot(w.T, sin_omega_t * cos_omega_t)
    # C2 = np.dot(w.T, np.cos(2 * omega * t)
    C2 = 2 * np.dot(w.T, 0.5 - sin_omega_t**2)

    if fit_mean:
        S = np.dot(w.T, sin_omega_t)
        C = np.dot(w.T, cos_omega_t)

        S2 -= 2 * S * C
        C2 -= C * C - S * S

    # compute components needed for the fit
    omega_t_tau = omega * t - 0.5 * np.arctan2(S2, C2)

    sin_omega_t_tau = np.sin(omega_t_tau)
    cos_omega_t_tau = np.cos(omega_t_tau)

    Y = np.dot(w.T, y)

    wy = w * y

    YCtau = np.dot(wy.T, cos_omega_t_tau)
    YStau = np.dot(wy.T, sin_omega_t_tau)
    CCtau = np.dot(w.T, cos_omega_t_tau * cos_omega_t_tau)
    SStau = np.dot(w.T, sin_omega_t_tau * sin_omega_t_tau)

    if fit_mean:
        Ctau = np.dot(w.T, cos_omega_t_tau)
        Stau = np.dot(w.T, sin_omega_t_tau)

        YCtau -= Y * Ctau
        YStau -= Y * Stau
        CCtau -= Ctau * Ctau
        SStau -= Stau * Stau

    p = YCtau * YCtau / CCtau + YStau * YStau / SStau
    YY = np.dot(w.T, y * y)

    if normalization == "standard":
        p /= YY
    elif normalization == "model":
        p /= YY - p
    elif normalization == "log":
        p = -np.log(1 - p / YY)
    elif normalization == "psd":
        p *= 0.5 * (dy**-2.0).sum()
    else:
        raise ValueError(f"normalization='{normalization}' not recognized")
    return p.ravel()
