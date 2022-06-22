import numpy as np

"""
Modification of slow_impl.py to allow fitting a shifted frequency (Df) template.
Adrian Bayer abayer@berkeley.edu
"""

def lombscargle_slow_Df(t, y, dy, frequency, Df=None, normalization='standard',
                     fit_mean=True, center_data=True):
    """Lomb-Scargle Periodogram

    This is a pure-python implementation of the original Lomb-Scargle formalism
    (e.g. [1]_, [2]_), with the addition of the floating mean (e.g. [3]_)

    Parameters
    ----------
    t, y, dy, Df : array-like
        times, values, errors, and frequency shift of the data points. 
        These should be broadcastable to the same shape. 
        None should be `~astropy.units.Quantity`.
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
    if Df is None:
        Df = 0

    t, y, dy, Df = np.broadcast_arrays(t, y, dy, Df)
    frequency = np.asarray(frequency)

    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if frequency.ndim != 1:
        raise ValueError("frequency should be one-dimensional")

    w = dy ** -2.0
    w /= w.sum()

    # if fit_mean is true, centering the data now simplifies the math below.
    if fit_mean or center_data:
        y = y - np.dot(w, y)

    omega = 2 * np.pi * frequency
    omega = omega.ravel()[np.newaxis, :]

    # make following arrays into column vectors
    t, y, dy, w, Df = map(lambda x: x[:, np.newaxis], (t, y, dy, w, Df))

    # shift omega by 2 pi Df
    omega = omega + 2 * np.pi * Df
    #omega = omega * Df

    """
    sin_omega_t = np.sin(omega * t)
    cos_omega_t = np.cos(omega * t)

    # S2 = np.dot(w.T, np.sin(2 * omega * t)
    S2 = 2 * np.dot(w.T, sin_omega_t * cos_omega_t)
    # C2 = np.dot(w.T, np.cos(2 * omega * t)
    C2 = 2 * np.dot(w.T, 0.5 - sin_omega_t ** 2)
    
    CS = np.dot(w.T, )

    if fit_mean:
        S = np.dot(w.T, sin_omega_t)
        C = np.dot(w.T, cos_omega_t)

        S2 -= (2 * S * C)
        C2 -= (C * C - S * S)
    """
    if fit_mean:
        raise Exception("NotImplementedError")
    
    # compute components needed for the fit
    omega_t = omega * t

    sin_omega_t = np.sin(omega_t)
    cos_omega_t = np.cos(omega_t)

    Y = np.dot(w.T, y)

    wy = w * y

    YC = np.dot(wy.T, cos_omega_t)
    YS = np.dot(wy.T, sin_omega_t)
    CC = np.dot(w.T, cos_omega_t * cos_omega_t)
    SS = np.dot(w.T, sin_omega_t * sin_omega_t)
    CS = np.dot(w.T, cos_omega_t * sin_omega_t)
    
    D = CC * SS - CS**2

    if fit_mean:
        raise Exception("NotImplementedError")
        #Ctau = np.dot(w.T, cos_omega_t_tau)
        #Stau = np.dot(w.T, sin_omega_t_tau)

        #YCtau -= Y * Ctau
        #YStau -= Y * Stau
        #CCtau -= Ctau * Ctau
        #SStau -= Stau * Stau
        
    term1 = SS/D * YC - CS/D * YS
    term2 = -CS/D * YC + CC/D * YS
    
    #p = (YCtau * YCtau / CCtau + YStau * YStau / SStau)
    p = YC * term1 + YS * term2
    
    YY = np.dot(w.T, y * y)
    
    # FIXME, unsure if any normalization other than psd makes sense
    if normalization == 'standard':
        p /= YY
    elif normalization == 'model':
        p /= YY - p
    elif normalization == 'log':
        p = -np.log(1 - p / YY)
    elif normalization == 'psd':
        p *= 0.5 * (dy ** -2.0).sum()
    else:
        raise ValueError(f"normalization='{normalization}' not recognized")
    return p.ravel()
