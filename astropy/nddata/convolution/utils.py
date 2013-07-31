from __future__ import division
import numpy as np


class NormalizationError(Exception):
    """
    Called when normalization of kernels is wrong.
    """
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message


class DiscretizationError(Exception):
    """
    Called when discretization of models goes wrong.
    """
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message


class KernelSizeError(Exception):
    """
    Called when size of kernels is wrong.
    """
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message


def next_larger_odd(number):
    """
    Round int to next larger odd int.
    """
    number = int(number + 0.5)
    if number % 2 == 0:
        return number + 1
    else:
        return number



default_sizes = {
                 }


def discretize_model(model, range=None, mode='center'):
    """
    Function to evaluate analytical models on a grid.

    Parameters
    ----------
    model : Instance of ParametricModel
        Instance of a astropy.ParametricModel to be evaluated.
    range : tuple
        Range in which the model is evaluated.
    mode: string
        One of the following modes:
            * 'center'
                Discretize model by taking the value
                at the center of the bin.
            * 'corner'
                Discretize model by taking average of
                the values at the corners of the bin.
            * 'oversample'
                Discretize model by taking the average
                on an oversampled grid.
            * 'integrate'
                Discretize model by integrating the
                model over the bin.

    Examples
    --------
    """
    if range == None:
        pass

    if mode == "center":
        pass
    elif mode == "corner":
        pass
    elif mode == "oversample":
        pass
    elif mode == "integrate":
        pass
    return 0


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


def discretize_corner_1D(model, x_range):
    """
    Discretize model by taking average of the values at the corners of the bin.
    """
    # Evaluate model 0.5 pixel outside the boundaries
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    values_intermediate_grid = model(x)

    # Convolve array to calculate the mean
    return np.convolve(values_intermediate_grid, [0.5, 0.5])[1:-1]


def discretize_corner_2D(model, x_range, y_range):
    """
    Discretize model by taking average of the values at the corners of the pixel.
    """
    # Evaluate model 0.5 pixel outside the boundaries
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    y = np.arange(y_range[0] - 0.5, y_range[1] + 0.5)
    x, y = np.meshgrid(x, y)
    values_intermediate_grid = model(x, y)

    # Save shape
    shape = values_intermediate_grid.shape

    # Convolve in x direction to calculate the mean
    convolved = np.convolve(np.ravel(values_intermediate_grid), [0.5, 0.5])

    # Reshape and transpose
    convolved = np.reshape(convolved[:-1], shape).T

    # Convolve in y direction to calculate the mean
    convolved = np.convolve(np.ravel(convolved), [0.5, 0.5])
    convolved = np.reshape(convolved[:-1], shape)

    # Cut border
    return convolved[1:, 1:]


def discretize_oversample_1D(model, x_range, factor=10):
    """
    Discretize model by taking the average on an oversampled grid.
    """
    # Evaluate model on oversampled grid
    x = np.arange(x_range[0] - 0.5 * (1 - 1 / factor),
                  x_range[1] - 0.5 * (1 - 1 / factor), 1. / factor)

    values = model(x)

    # Reshape and compute mean
    values = np.reshape(values, (x.size, factor))
    return values.mean(axis=1)


def discretize_oversample_2D(model, x_range, y_range, factor=10):
    """
    Discretize model by taking the average on an oversampled grid.
    """
    # Evaluate model on oversampled grid
    x = np.arange(x_range[0] - 0.5 * (1 - 1 / factor),
                  x_range[1] - 0.5 * (1 - 1 / factor), 1. / factor)

    y = np.arange(y_range[0] - 0.5 * (1 - 1 / factor),
                  y_range[1] - 0.5 * (1 - 1 / factor), 1. / factor)
    x_grid, y_grid = np.meshgrid(x, y)
    values = model(x_grid, y_grid)

    # Reshape and compute mean
    shape = (x.size / factor, y.size / factor, factor, factor)
    values = np.reshape(values, shape)
    values.transpose((0, 2, 1, 3))
    return values.mean(axis=3).mean(axis=2)


def discretize_integrate(model, range_, mode='analytical'):
    """
    Discretize model by integrating the model over the bin.
    """
    if getattr(model, "integral", False):
        raise NotImplementedError("Currently not supported")
    else:
        raise DiscretizationError("Model needs integrate function.")

