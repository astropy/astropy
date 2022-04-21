# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Statistic functions used in `~astropy.modeling.fitting`.
"""
# pylint: disable=invalid-name
import numpy as np

__all__ = ["leastsquare", "leastsquare_1d", "leastsquare_2d", "leastsquare_3d"]


def leastsquare(measured_vals, updated_model, weights, *x):
    """Least square statistic, with optional weights, in N-dimensions.

    Parameters
    ----------
    measured_vals : ndarray or sequence
        Measured data values. Will be cast to array whose
        shape must match the array-cast of the evaluated model.
    updated_model : :class:`~astropy.modeling.Model` instance
        Model with parameters set by the current iteration of the optimizer.
        when evaluated on "x", must return array of shape "measured_vals"
    weights : ndarray or None
        Array of weights to apply to each residual.
    *x : ndarray
        Independent variables on which to evaluate the model.

    Returns
    -------
    res : float
        The sum of least squares.

    See Also
    --------
    :func:`~astropy.modeling.statistic.leastsquare_1d`
    :func:`~astropy.modeling.statistic.leastsquare_2d`
    :func:`~astropy.modeling.statistic.leastsquare_3d`

    Notes
    -----
    Models in :mod:`~astropy.modeling` have broadcasting rules that try to
    match inputs with outputs with Model shapes. Numpy arrays have flexible
    broadcasting rules, so mismatched shapes can often be made compatible. To
    ensure data matches the model we must perform shape comparison and leverage
    the Numpy arithmetic functions. This can obfuscate arithmetic computation
    overrides, like with Quantities. Implement a custom statistic for more
    direct control.

    """

    model_vals = updated_model(*x)

    if np.shape(model_vals) != np.shape(measured_vals):
        raise ValueError(f"Shape mismatch between model ({np.shape(model_vals)}) "
                         f"and measured ({np.shape(measured_vals)})")

    if weights is None:
        weights = 1.0

    return np.sum(np.square(weights * np.subtract(model_vals, measured_vals)))


# -------------------------------------------------------------------

def leastsquare_1d(measured_vals, updated_model, weights, x):
    """
    Least square statistic with optional weights.
    Safer than the general :func:`~astropy.modeling.statistic.leastsquare`
    for 1D models by avoiding numpy methods that support broadcasting.

    Parameters
    ----------
    measured_vals : ndarray
        Measured data values.
    updated_model : `~astropy.modeling.Model`
        Model with parameters set by the current iteration of the optimizer.
    weights : ndarray or None
        Array of weights to apply to each residual.
    x : ndarray
        Independent variable "x" on which to evaluate the model.

    Returns
    -------
    res : float
        The sum of least squares.

    See Also
    --------
    :func:`~astropy.modeling.statistic.leastsquare`

    """
    model_vals = updated_model(x)

    if weights is None:
        return np.sum((model_vals - measured_vals) ** 2)
    return np.sum((weights * (model_vals - measured_vals)) ** 2)


def leastsquare_2d(measured_vals, updated_model, weights, x, y):
    """
    Least square statistic with optional weights.
    Safer than the general :func:`~astropy.modeling.statistic.leastsquare`
    for 2D models by avoiding numpy methods that support broadcasting.

    Parameters
    ----------
    measured_vals : ndarray
        Measured data values.
    updated_model : `~astropy.modeling.Model`
        Model with parameters set by the current iteration of the optimizer.
    weights : ndarray or None
        Array of weights to apply to each residual.
    x : ndarray
        Independent variable "x" on which to evaluate the model.
    y : ndarray
        Independent variable "y" on which to evaluate the model.

    Returns
    -------
    res : float
        The sum of least squares.

    See Also
    --------
    :func:`~astropy.modeling.statistic.leastsquare`

    """
    model_vals = updated_model(x, y)

    if weights is None:
        return np.sum((model_vals - measured_vals) ** 2)
    return np.sum((weights * (model_vals - measured_vals)) ** 2)


def leastsquare_3d(measured_vals, updated_model, weights, x, y, z):
    """
    Least square statistic with optional weights.
    Safer than the general :func:`~astropy.modeling.statistic.leastsquare`
    for 3D models by avoiding numpy methods that support broadcasting.

    Parameters
    ----------
    measured_vals : ndarray
        Measured data values.
    updated_model : `~astropy.modeling.Model`
        Model with parameters set by the current iteration of the optimizer.
    weights : ndarray or None
        Array of weights to apply to each residual.
    x : ndarray
        Independent variable "x" on which to evaluate the model.
    y : ndarray
        Independent variable "y" on which to evaluate the model.
    z : ndarray
        Independent variable "z" on which to evaluate the model.

    Returns
    -------
    res : float
        The sum of least squares.

    See Also
    --------
    :func:`~astropy.modeling.statistic.leastsquare`

    """
    model_vals = updated_model(x, y, z)

    if weights is None:
        return np.sum((model_vals - measured_vals) ** 2)
    return np.sum((weights * (model_vals - measured_vals)) ** 2)
