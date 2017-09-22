# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Statistic functions used in `~astropy.modeling.fitting`.
"""

import numpy as np

__all__ = ['leastsquare']


def leastsquare(measured_vals, updated_model, weights, x, y=None):
    """
    Least square statistic with optional weights.

    Parameters
    ----------
    measured_vals : `~numpy.ndarray`
        Measured data values.
    updated_model : `~astropy.modeling.Model`
        Model with parameters set by the current iteration of the optimizer.
    weights : `~numpy.ndarray`
        Array of weights to apply to each residual.
    x : `~numpy.ndarray`
        Independent variable "x" to evaluate the model on.
    y : `~numpy.ndarray`, optional
        Independent variable "y" to evaluate the model on, for 2D models.

    Returns
    -------
    res : float
        The sum of least squares.
    """

    if y is None:
        model_vals = updated_model(x)
    else:
        model_vals = updated_model(x, y)
    if weights is None:
        return np.sum((model_vals - measured_vals) ** 2)
    else:
        return np.sum((weights * (model_vals - measured_vals)) ** 2)
