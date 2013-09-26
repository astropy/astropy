# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import warnings

import numpy as np

from ..modeling.core import Parametric1DModel, Parametric2DModel


__all__ = ['discretize_model']


class DiscretizationError(Exception):
    """
    Called when discretization of models goes wrong.
    """


class KernelSizeError(Exception):
    """
    Called when size of kernels is even.
    """


def add_kernel_arrays_1D(array_1, array_2):
    """
    Add two 1D kernel arrays of different size.

    The arrays are added with the centers lying upon each other.
    """
    if array_1.size > array_2.size:
        center = array_1.size // 2
        slice_ = slice(center - array_2.size // 2,
                       center + array_2.size // 2 + 1)
        array_1[slice_] += array_2
        return array_1
    elif array_2.size > array_1.size:
        center = array_2.size // 2
        slice_ = slice(center - array_1.size // 2,
                       center + array_1.size // 2 + 1)
        array_2[slice_] += array_1
        return array_2
    return array_2 + array_1


def add_kernel_arrays_2D(array_1, array_2):
    """
    Add two 1D kernel arrays of different size.

    The arrays are added with the centers lying upon each other.
    """
    if array_1.size > array_2.size:
        center = [axes_size // 2 for axes_size in array_1.shape]
        slice_x = slice(center[1] - array_2.shape[1] // 2,
                        center[1] + array_2.shape[1] // 2 + 1)
        slice_y = slice(center[0] - array_2.shape[0] // 2,
                        center[0] + array_2.shape[0] // 2 + 1)
        array_1[slice_y, slice_x] += array_2
        return array_1
    elif array_2.size > array_1.size:
        center = [axes_size // 2 for axes_size in array_2.shape]
        slice_x = slice(center[1] - array_1.shape[1] // 2,
                        center[1] + array_1.shape[1] // 2 + 1)
        slice_y = slice(center[0] - array_1.shape[0] // 2,
                        center[0] + array_1.shape[0] // 2 + 1)
        array_2[slice_y, slice_x] += array_1
        return array_2
    return array_2 + array_1


def discretize_model(model, x_range, y_range=None, mode='center', factor=10):
    """
    Function to evaluate analytical models on a grid.

    Parameters
    ----------
    model : :class:`~astropy.modeling.core.ParametricModel` instance
        Instance of a :class:`~astropy.modeling.core.ParametricModel` to be evaluated.
    x_range : tuple
        x range in which the model is evaluated.
    y_range : tuple, optional
        y range in which the model is evaluated.
        Necessary only for 2D models.
    mode : str, optional
        One of the following modes:
            * ``'center'`` (default)
                Discretize model by taking the value
                at the center of the bin.
            * ``'linear_interp'``
                Discretize model by linearly interpolating
                between the values at the corners of the bin.
                For 2D models interpolation is bilinear.
            * ``'oversample'``
                Discretize model by taking the average
                on an oversampled grid.
            * ``'integrate'``
                Discretize model by integrating the model 
                over the bin using `scipy.integrate.quad`.
                Very slow.
    factor : float or int
        Factor of oversampling. Default = 10.

    Notes
    -----

    The ``oversample`` mode allows to conserve the integral on a subpixel
    scale. Here is the example of a normalized Gaussian1DModel:

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import numpy as np
        from astropy.modeling.models import Gaussian1DModel
        from astropy.convolution.utils import discretize_model
        gauss_1D = Gaussian1DModel(1 / (0.5 * np.sqrt(2 * np.pi)), 0, 0.5)
        y_center = discretize_model(gauss_1D, (-2, 3), mode='center')
        y_corner = discretize_model(gauss_1D, (-2, 3), mode='linear_interp')
        y_oversample = discretize_model(gauss_1D, (-2, 3), mode='oversample')
        plt.plot(y_center, label='center sum = {0:3f}'.format(y_center.sum()))
        plt.plot(y_corner, label='linear_interp sum = {0:3f}'.format(y_corner.sum()))
        plt.plot(y_oversample, label='oversample sum = {0:3f}'.format(y_oversample.sum()))
        plt.xlabel('pixels')
        plt.ylabel('value')
        plt.legend()
        plt.show()


    """
    if isinstance(model, Parametric2DModel) and y_range is None:
        raise Exception("Please specify y range.")
    if mode == "center":
        if isinstance(model, Parametric1DModel):
            return discretize_center_1D(model, x_range)
        if isinstance(model, Parametric2DModel):
            return discretize_center_2D(model, x_range, y_range)
    elif mode == "linear_interp":
        if isinstance(model, Parametric1DModel):
            return discretize_linear_1D(model, x_range)
        if isinstance(model, Parametric2DModel):
            return discretize_bilinear_2D(model, x_range, y_range)
    elif mode == "oversample":
        if y_range is not None:
            N = factor * (x_range[1] - x_range[0]) * (y_range[1] - y_range[0])
        else:
            N = factor * (x_range[1] - x_range[0])
        if N >= 10000000:
            warnings.warn("Large oversample factor or data, computing can be very slow.")
        if isinstance(model, Parametric1DModel):
            return discretize_oversample_1D(model, x_range, factor)
        if isinstance(model, Parametric2DModel):
            return discretize_oversample_2D(model, x_range, y_range, factor)
    elif mode == "integrate":
        warnings.warn("Mode 'integrate' is very slow. Use only if highest "
                        "accuracy is required.")
        if isinstance(model, Parametric1DModel):
            return discretize_integrate_1D(model, x_range)
        if isinstance(model, Parametric2DModel):
            return discretize_integrate_2D(model, x_range, y_range)
    else:
        raise DiscretizationError('Invalid mode.')


def discretize_center_1D(model, x_range):
    """
    Discretize model by taking the value at the center of the bin.
    """
    x = np.arange(*x_range)
    return model(x)


def discretize_center_2D(model, x_range, y_range):
    """
    Discretize model by taking the value at the center of the pixel.
    """
    x = np.arange(*x_range)
    y = np.arange(*y_range)
    x, y = np.meshgrid(x, y)
    return model(x, y)


def discretize_linear_1D(model, x_range):
    """
    Discretize model by performing a linear interpolation.
    """
    # Evaluate model 0.5 pixel outside the boundaries
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    values_intermediate_grid = model(x)
    return 0.5 * (values_intermediate_grid[1:] + values_intermediate_grid[:-1])


def discretize_bilinear_2D(model, x_range, y_range):
    """
    Discretize model by performing a bilinear interpolation.
    """
    # Evaluate model 0.5 pixel outside the boundaries
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    y = np.arange(y_range[0] - 0.5, y_range[1] + 0.5)
    x, y = np.meshgrid(x, y)
    values_intermediate_grid = model(x, y)

    # Mean in y direction
    values = 0.5 * (values_intermediate_grid[1:, :]
                    + values_intermediate_grid[:-1, :])
    # Mean in x direction
    values = 0.5 * (values[:, 1:]
                    + values[:, :-1])
    return values


def discretize_oversample_1D(model, x_range, factor=10):
    """
    Discretize model by taking the average on an oversampled grid.
    """
    # Evaluate model on oversampled grid
    x = np.arange(x_range[0] - 0.5 * (1 - 1 / factor),
                  x_range[1] + 0.5 * (1 + 1 / factor), 1. / factor)

    values = model(x)

    # Reshape and compute mean
    values = np.reshape(values, (x.size / factor, factor))
    return values.mean(axis=1)[:-1]


def discretize_oversample_2D(model, x_range, y_range, factor=10):
    """
    Discretize model by taking the average on an oversampled grid.
    """
    # Evaluate model on oversampled grid
    x = np.arange(x_range[0] - 0.5 * (1 - 1 / factor),
                  x_range[1] + 0.5 * (1 + 1 / factor), 1. / factor)

    y = np.arange(y_range[0] - 0.5 * (1 - 1 / factor),
                  y_range[1] + 0.5 * (1 + 1 / factor), 1. / factor)
    x_grid, y_grid = np.meshgrid(x, y)
    values = model(x_grid, y_grid)

    # Reshape and compute mean
    shape = (y.size / factor, factor, x.size / factor, factor)
    values = np.reshape(values, shape)
    return values.mean(axis=3).mean(axis=1)[:-1, :-1]


def discretize_integrate_1D(model, x_range):
    """
    Discretize model by integrating numerically the model over the bin.
    """
    try:
        from scipy.integrate import quad
    except ImportError:
        raise Exception("Mode 'integrate' requires scipy.")

    # Set up grid
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    values = np.array([])

    # Integrate over all bins
    for i in range(x.size - 1):
        values = np.append(values, quad(model, x[i], x[i + 1])[0])
    return values


def discretize_integrate_2D(model, x_range, y_range):
    """
    Discretize model by integrating the model over the pixel.
    """
    try:
        from scipy.integrate import dblquad
    except ImportError:
        raise Exception("Mode 'integrate' requires scipy.")

    # Set up grid
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    y = np.arange(y_range[0] - 0.5, y_range[1] + 0.5)
    values = np.empty((y.size, x.size))

    # Integrate over all pixels
    for i in range(x.size - 1):
        for j in range(y.size - 1):
            values[j, i] = dblquad(model, x[i], x[i + 1],
                                   lambda x: y[j], lambda x: y[j + 1])[0]
    return values
