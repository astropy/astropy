# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from astropy.modeling.core import Model, custom_model

__all__ = [
    "discretize_model",
    "KernelError",
    "KernelSizeError",
    "KernelArithmeticError",
]


class KernelError(Exception):
    """
    Base error class for kernel errors.
    """


class KernelSizeError(KernelError):
    """
    Called when size of kernels is even.
    """


class KernelArithmeticError(KernelError):
    """Called when doing invalid arithmetic with a kernel."""


def has_even_axis(array):
    if isinstance(array, (list, tuple)):
        return not len(array) % 2
    else:
        return any(not axes_size % 2 for axes_size in array.shape)


def add_kernel_arrays_1D(array_1, array_2):
    """
    Add two 1D kernel arrays of different size.

    The arrays are added with the centers lying upon each other.
    """
    if array_1.size > array_2.size:
        new_array = array_1.copy()
        center = array_1.size // 2
        slice_ = slice(center - array_2.size // 2, center + array_2.size // 2 + 1)
        new_array[slice_] += array_2
        return new_array
    elif array_2.size > array_1.size:
        new_array = array_2.copy()
        center = array_2.size // 2
        slice_ = slice(center - array_1.size // 2, center + array_1.size // 2 + 1)
        new_array[slice_] += array_1
        return new_array
    return array_2 + array_1


def add_kernel_arrays_2D(array_1, array_2):
    """
    Add two 2D kernel arrays of different size.

    The arrays are added with the centers lying upon each other.
    """
    if array_1.size > array_2.size:
        new_array = array_1.copy()
        center = [axes_size // 2 for axes_size in array_1.shape]
        slice_x = slice(
            center[1] - array_2.shape[1] // 2, center[1] + array_2.shape[1] // 2 + 1
        )
        slice_y = slice(
            center[0] - array_2.shape[0] // 2, center[0] + array_2.shape[0] // 2 + 1
        )
        new_array[slice_y, slice_x] += array_2
        return new_array
    elif array_2.size > array_1.size:
        new_array = array_2.copy()
        center = [axes_size // 2 for axes_size in array_2.shape]
        slice_x = slice(
            center[1] - array_1.shape[1] // 2, center[1] + array_1.shape[1] // 2 + 1
        )
        slice_y = slice(
            center[0] - array_1.shape[0] // 2, center[0] + array_1.shape[0] // 2 + 1
        )
        new_array[slice_y, slice_x] += array_1
        return new_array
    return array_2 + array_1


def discretize_model(model, x_range, y_range=None, mode="center", factor=10):
    """
    Evaluate an analytical model function on a pixel grid.

    Parameters
    ----------
    model : `~astropy.modeling.Model` or callable.
        Analytical model function to be discretized. A callable that is
        not a `~astropy.modeling.Model` instance is converted to a model
        using `~astropy.modeling.custom_model`.
    x_range : 2-tuple
        Lower and upper bounds of x pixel values at which the model is
        evaluated. The upper bound is non-inclusive. A ``x_range`` of
        ``(0, 3)`` means the model will be evaluated at x pixels 0, 1,
        and 2. The difference between the upper and lower bound must be
        a whole number so that the output array size is well defined.
    y_range : 2-tuple or `None`, optional
        Lower and upper bounds of y pixel values at which the model is
        evaluated. The upper bound is non-inclusive. A ``y_range`` of
        ``(0, 3)`` means the model will be evaluated at y pixels of 0,
        1, and 2. The difference between the upper and lower bound must
        be a whole number so that the output array size is well defined.
        ``y_range`` is necessary only for 2D models.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
        One of the following modes:
            * ``'center'`` (default)
                Discretize model by taking the value at the center of
                the pixel bins.
            * ``'linear_interp'``
                Discretize model by linearly interpolating between the
                values at the edges (1D) or corners (2D) of the pixel
                bins. For 2D models, the interpolation is bilinear.
            * ``'oversample'``
                Discretize model by taking the average of model values
                in the pixel bins on an oversampled grid. Use the
                ``factor`` keyword to set the integer oversampling
                factor.
            * ``'integrate'``
                Discretize model by integrating the model over the pixel
                bins using `scipy.integrate.quad`. This mode conserves
                the model integral on a subpixel scale, but is very
                slow.
    factor : int, optional
        The integer oversampling factor used when ``mode='oversample'``.
        Ignored otherwise.

    Returns
    -------
    array : `numpy.ndarray`
        The discretized model array.

    Examples
    --------
    In this example, we define a
    `~astropy.modeling.functional_models.Gaussian1D` model that has been
    normalized so that it sums to 1.0. We then discretize this model
    using the ``'center'``, ``'linear_interp'``, and ``'oversample'``
    (with ``factor=10``) modes.

    .. plot::
        :show-source-link:

        import matplotlib.pyplot as plt
        import numpy as np
        from astropy.convolution.utils import discretize_model
        from astropy.modeling.models import Gaussian1D

        gauss_1D = Gaussian1D(1 / (0.5 * np.sqrt(2 * np.pi)), 0, 0.5)
        x_range = (-2, 3)
        x = np.arange(*x_range)
        y_center = discretize_model(gauss_1D, x_range, mode='center')
        y_edge = discretize_model(gauss_1D, x_range, mode='linear_interp')
        y_oversample = discretize_model(gauss_1D, x_range, mode='oversample')

        fig, ax = plt.subplots(figsize=(8, 6))
        label = f'center (sum={y_center.sum():.3f})'
        ax.plot(x, y_center, '.-', label=label)
        label = f'linear_interp (sum={y_edge.sum():.3f})'
        ax.plot(x, y_edge, '.-', label=label)
        label = f'oversample (sum={y_oversample.sum():.3f})'
        ax.plot(x, y_oversample, '.-', label=label)
        ax.set_xlabel('x')
        ax.set_ylabel('Value')
        plt.legend()
    """
    if not callable(model):
        raise TypeError("Model must be callable.")
    if not isinstance(model, Model):
        model = custom_model(model)()
    ndim = model.n_inputs
    if ndim > 2:
        raise ValueError("discretize_model supports only 1D and 2D models.")

    dxrange = np.diff(x_range)[0]
    if dxrange != int(dxrange):
        raise ValueError(
            "The difference between the upper and lower limit of"
            " 'x_range' must be a whole number."
        )

    if y_range:
        dyrange = np.diff(y_range)[0]
        if dyrange != int(dyrange):
            raise ValueError(
                "The difference between the upper and lower limit of"
                " 'y_range' must be a whole number."
            )

    if factor != int(factor):
        raise ValueError("factor must have an integer value")
    factor = int(factor)

    if ndim == 2 and y_range is None:
        raise ValueError("y_range must be specified for a 2D model")
    if ndim == 1 and y_range is not None:
        raise ValueError("y_range should not be input for a 1D model")
    if mode == "center":
        if ndim == 1:
            return discretize_center_1D(model, x_range)
        elif ndim == 2:
            return discretize_center_2D(model, x_range, y_range)
    elif mode == "linear_interp":
        if ndim == 1:
            return discretize_linear_1D(model, x_range)
        if ndim == 2:
            return discretize_bilinear_2D(model, x_range, y_range)
    elif mode == "oversample":
        if ndim == 1:
            return discretize_oversample_1D(model, x_range, factor)
        if ndim == 2:
            return discretize_oversample_2D(model, x_range, y_range, factor)
    elif mode == "integrate":
        if ndim == 1:
            return discretize_integrate_1D(model, x_range)
        if ndim == 2:
            return discretize_integrate_2D(model, x_range, y_range)
    else:
        raise ValueError("Invalid mode for discretize_model.")


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
    values = 0.5 * (values_intermediate_grid[1:, :] + values_intermediate_grid[:-1, :])
    # Mean in x direction
    values = 0.5 * (values[:, 1:] + values[:, :-1])
    return values


def discretize_oversample_1D(model, x_range, factor=10):
    """
    Discretize model by taking the average on an oversampled grid.
    """
    # Evaluate model on oversampled grid
    x = np.linspace(
        x_range[0] - 0.5 * (1 - 1 / factor),
        x_range[1] - 0.5 * (1 + 1 / factor),
        num=int((x_range[1] - x_range[0]) * factor),
    )

    values = model(x)

    # Reshape and compute mean
    values = np.reshape(values, (x.size // factor, factor))
    return values.mean(axis=1)


def discretize_oversample_2D(model, x_range, y_range, factor=10):
    """
    Discretize model by taking the average on an oversampled grid.
    """
    # Evaluate model on oversampled grid
    x = np.linspace(
        x_range[0] - 0.5 * (1 - 1 / factor),
        x_range[1] - 0.5 * (1 + 1 / factor),
        num=int((x_range[1] - x_range[0]) * factor),
    )
    y = np.linspace(
        y_range[0] - 0.5 * (1 - 1 / factor),
        y_range[1] - 0.5 * (1 + 1 / factor),
        num=int((y_range[1] - y_range[0]) * factor),
    )

    x_grid, y_grid = np.meshgrid(x, y)
    values = model(x_grid, y_grid)

    # Reshape and compute mean
    shape = (y.size // factor, factor, x.size // factor, factor)
    values = np.reshape(values, shape)
    return values.mean(axis=3).mean(axis=1)


def discretize_integrate_1D(model, x_range):
    """
    Discretize model by integrating numerically the model over the bin.
    """
    from scipy.integrate import quad

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
    from scipy.integrate import dblquad

    # Set up grid
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    y = np.arange(y_range[0] - 0.5, y_range[1] + 0.5)
    values = np.empty((y.size - 1, x.size - 1))

    # Integrate over all pixels
    for i in range(x.size - 1):
        for j in range(y.size - 1):
            values[j, i] = dblquad(
                func=lambda y, x: model(x, y),
                a=x[i],
                b=x[i + 1],
                gfun=lambda x: y[j],
                hfun=lambda x: y[j + 1],
            )[0]
    return values
