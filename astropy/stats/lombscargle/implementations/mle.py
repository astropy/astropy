"""Tools for maximum likelihood estimation associated with Lomb-Scargle"""
import numpy as np

from .utils import validate_inputs


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
    assert t.ndim == 1
    assert frequency.ndim == 0

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
                 center_data=True, fit_bias=True, nterms=1):
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
    fit_bias : bool (default=True)
        If True, include the bias as part of the model
    nterms : int (default=1)
        The number of Fourier terms to include in the fit

    Returns
    -------
    y_fit : ndarray
        The model fit evaluated at each value of t_fit
    """
    if t_fit is None:
        raise ValueError('t_fit must not be None')
    if dy is None:
        dy = 1

    t, y, dy, frequency, t_fit, unit_dict = validate_inputs(t, y, dy,
                                                            frequency,
                                                            t_fit=t_fit)

    t_fit = np.asarray(t_fit)
    assert t.ndim == 1
    assert t_fit.ndim == 1
    assert frequency.ndim == 0

    if center_data:
        w = dy ** -2.0
        y_mean = np.dot(y, w) / w.sum()
        y = (y - y_mean)
    else:
        y_mean = 0

    X = design_matrix(t, frequency, dy=dy, bias=fit_bias, nterms=nterms)
    theta_MLE = np.linalg.solve(np.dot(X.T, X),
                                np.dot(X.T, y / dy))

    X_fit = design_matrix(t_fit, frequency, bias=fit_bias, nterms=nterms)

    return (y_mean + np.dot(X_fit, theta_MLE)) * unit_dict['y']
