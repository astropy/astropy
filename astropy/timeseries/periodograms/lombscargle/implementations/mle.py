
import numpy as np


def design_matrix(t, frequency, dy=None, bias=True, nterms=1):
    """Compute the Lomb-Scargle design matrix at the given frequency

    This is the matrix X such that the periodic model at the given frequency
    can be expressed :math:`\\hat{y} = X \\theta`.

    Parameters
    ----------
    t : array_like, shape=(n_times,)
        times at which to compute the design matrix
    frequency : float
        frequency for the design matrix
    dy : float or array_like (optional)
        data uncertainties: should be broadcastable with `t`
    bias : bool (default=True)
        If true, include a bias column in the matrix
    nterms : int (default=1)
        Number of Fourier terms to include in the model

    Returns
    -------
    X : ndarray, shape=(n_times, n_parameters)
        The design matrix, where n_parameters = bool(bias) + 2 * nterms
    """
    t = np.asarray(t)
    frequency = np.asarray(frequency)

    if t.ndim != 1:
        raise ValueError("t should be one dimensional")
    if frequency.ndim != 0:
        raise ValueError("frequency must be a scalar")

    if nterms == 0 and not bias:
        raise ValueError("cannot have nterms=0 and no bias")

    if bias:
        cols = [np.ones_like(t)]
    else:
        cols = []

    for i in range(1, nterms + 1):
        cols.append(np.sin(2 * np.pi * i * frequency * t))
        cols.append(np.cos(2 * np.pi * i * frequency * t))
    XT = np.vstack(cols)

    if dy is not None:
        XT /= dy

    return np.transpose(XT)


def periodic_fit(t, y, dy, frequency, t_fit,
                 center_data=True, fit_mean=True, nterms=1):
    """Compute the Lomb-Scargle model fit at a given frequency

    Parameters
    ----------
    t, y, dy : float or array_like
        The times, observations, and uncertainties to fit
    frequency : float
        The frequency at which to compute the model
    t_fit : float or array_like
        The times at which the fit should be computed
    center_data : bool (default=True)
        If True, center the input data before applying the fit
    fit_mean : bool (default=True)
        If True, include the bias as part of the model
    nterms : int (default=1)
        The number of Fourier terms to include in the fit

    Returns
    -------
    y_fit : ndarray
        The model fit evaluated at each value of t_fit
    """
    t, y, frequency = map(np.asarray, (t, y, frequency))
    if dy is None:
        dy = np.ones_like(y)
    else:
        dy = np.asarray(dy)

    t_fit = np.asarray(t_fit)

    if t.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if t_fit.ndim != 1:
        raise ValueError("t_fit should be one dimensional")
    if frequency.ndim != 0:
        raise ValueError("frequency should be a scalar")

    if center_data:
        w = dy ** -2.0
        y_mean = np.dot(y, w) / w.sum()
        y = (y - y_mean)
    else:
        y_mean = 0

    X = design_matrix(t, frequency, dy=dy, bias=fit_mean, nterms=nterms)
    theta_MLE = np.linalg.solve(np.dot(X.T, X),
                                np.dot(X.T, y / dy))

    X_fit = design_matrix(t_fit, frequency, bias=fit_mean, nterms=nterms)

    return y_mean + np.dot(X_fit, theta_MLE)
