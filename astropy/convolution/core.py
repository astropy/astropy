# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains the convolution and filter functionalities of astropy.

A few conceptual notes:
A filter kernel is mainly characterized by its response function. In the 1D
case we speak of "impulse response function", in the 2D case we call it "point
spread function". This response function is given for every kernel by an
astropy `FittableModel`, which is evaluated on a grid to obtain a filter array,
which can then be applied to binned data.

The model is centered on the array and should have an amplitude such that the array
integrates to one per default.

Currently only symmetric 2D kernels are supported.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import warnings
import copy

import numpy as np
from ..utils.exceptions import AstropyUserWarning
from .utils import (discretize_model, add_kernel_arrays_1D,
                    add_kernel_arrays_2D)

MAX_NORMALIZATION = 100

__all__ = ['Kernel', 'Kernel1D', 'Kernel2D', 'kernel_arithmetics']


class Kernel(object):
    """
    Convolution kernel base class.

    Parameters
    ----------
    array : `~numpy.ndarray`
        Kernel array.
    """
    _separable = False
    _is_bool = True
    _model = None

    def __init__(self, array):
        self._array = np.asanyarray(array)
        self.calculate_normalization()

    @property
    def truncation(self):
        """
        Deviation from the normalization to one.
        """
        return self._truncation

    @property
    def is_bool(self):
        """
        Indicates if kernel is bool.

        If the kernel is bool the multiplication in the convolution could
        be omitted, to increase the performance.
        """
        return self._is_bool

    @property
    def model(self):
        """
        Kernel response model.
        """
        return self._model

    @property
    def dimension(self):
        """
        Kernel dimension.
        """
        return self.array.ndim

    @property
    def center(self):
        """
        Index of the kernel center.
        """
        return [axes_size // 2 for axes_size in self._array.shape]

    def calculate_normalization(self, mode='integral'):
        """
        Calculate the kernel normalization factor.

        Parameters
        ----------
        mode : {'integral', 'peak'}
            One of the following modes:
                * 'integral' (default)
                    Kernel is normalized such that its integral = 1.
                * 'peak'
                    Kernel is normalized such that its peak = 1.
        """

        if mode == 'integral':
            if self._array.sum() == 0:
                self._normalization = np.inf
            else:
                self._normalization = 1. / self._array.sum()
        elif mode == 'peak':
            self._normalization = 1. / self._array.max()
        else:
            raise ValueError("invalid mode, must be 'integral' or 'peak'")
        return self._normalization

    @property
    def normalization(self):
        """
        The kernel normalization factor.
        """
        return self._normalization

    def normalize(self, mode='integral'):
        """
        Normalize the filter kernel.

        Parameters
        ----------
        mode : {'integral', 'peak'}
            One of the following modes:
                * 'integral' (default)
                    Kernel is normalized such that its integral = 1.
                * 'peak'
                    Kernel is normalized such that its peak = 1.
        """

        np.multiply(self._array, self.calculate_normalization(mode=mode),
                    self.array)

        # Warn the user for kernels that sum to zero
        if np.isinf(self._normalization):
            warnings.warn('Kernel cannot be normalized because the '
                          'normalization factor is infinite.',
                          AstropyUserWarning)
            return

        if np.abs(self._normalization) > MAX_NORMALIZATION:
            warnings.warn('The kernel normalization factor is exceptionally '
                          'large, > {0}.'.format(MAX_NORMALIZATION),
                          AstropyUserWarning)

    @property
    def shape(self):
        """
        Shape of the kernel array.
        """
        return self._array.shape

    @property
    def separable(self):
        """
        Indicates if the filter kernel is separable.

        A 2D filter is separable, when its filter array can be written as the
        outer product of two 1D arrays.

        If a filter kernel is separable, higher dimension convolutions will be
        performed by applying the 1D filter array consecutively on every dimension.
        This is significantly faster, than using a filter array with the same
        dimension.
        """
        return self._separable

    @property
    def array(self):
        """
        Filter kernel array.
        """
        return self._array

    def __add__(self, kernel):
        """
        Add two filter kernels.
        """
        return kernel_arithmetics(self, kernel, 'add')

    def __sub__(self,  kernel):
        """
        Subtract two filter kernels.
        """
        return kernel_arithmetics(self, kernel, 'sub')

    def __mul__(self, value):
        """
        Multiply kernel with number or convolve two kernels.
        """
        return kernel_arithmetics(self, value, "mul")

    def __rmul__(self, value):
        """
        Multiply kernel with number or convolve two kernels.
        """
        return kernel_arithmetics(self, value, "mul")

    def __array__(self):
        """
        Array representation of the kernel.
        """
        return self._array

    def __array_wrap__(self, array, context=None):
        """
        Wrapper for multiplication with numpy arrays.
        """
        if type(context[0]) == np.ufunc:
            return NotImplemented
        else:
            return array


class Kernel1D(Kernel):
    """
    Base class for 1D filter kernels.

    Parameters
    ----------
    model : `~astropy.modeling.FittableModel`
        Model to be evaluated.
    x_size : odd int, optional
        Size of the kernel array. Default = 8 * width.
    array : `~numpy.ndarray`
        Kernel array.
    width : number
        Width of the filter kernel.
    mode : str, optional
        One of the following discretization modes:
            * 'center' (default)
                Discretize model by taking the value
                at the center of the bin.
            * 'linear_interp'
                Discretize model by linearly interpolating
                between the values at the corners of the bin.
            * 'oversample'
                Discretize model by taking the average
                on an oversampled grid.
            * 'integrate'
                Discretize model by integrating the
                model over the bin.
    factor : number, optional
        Factor of oversampling. Default factor = 10.
    """
    def __init__(self, model=None, x_size=None, array=None, **kwargs):
        # Initialize from model
        if array is None:
            if self._model is None:
                raise TypeError("Must specify either array or model.")

            if x_size is None:
                x_size = self._default_size
            elif x_size != int(x_size):
                raise TypeError("x_size should be an integer")

            # Set ranges where to evaluate the model

            if x_size % 2 == 0:  # even kernel
                x_range = (-(int(x_size)) // 2 + 0.5, (int(x_size)) // 2 + 0.5)
            else:  # odd kernel
                x_range = (-(int(x_size) - 1) // 2, (int(x_size) - 1) // 2 + 1)

            array = discretize_model(self._model, x_range, **kwargs)

        # Initialize from array
        elif array is not None:
            self._model = None

        super(Kernel1D, self).__init__(array)


class Kernel2D(Kernel):
    """
    Base class for 2D filter kernels.

    Parameters
    ----------
    model : `~astropy.modeling.FittableModel`
        Model to be evaluated.
    x_size : odd int, optional
        Size in x direction of the kernel array. Default = 8 * width.
    y_size : odd int, optional
        Size in y direction of the kernel array. Default = 8 * width.
    array : `~numpy.ndarray`
        Kernel array.
    mode : str, optional
        One of the following discretization modes:
            * 'center' (default)
                Discretize model by taking the value
                at the center of the bin.
            * 'linear_interp'
                Discretize model by performing a bilinear interpolation
                between the values at the corners of the bin.
            * 'oversample'
                Discretize model by taking the average
                on an oversampled grid.
            * 'integrate'
                Discretize model by integrating the
                model over the bin.
    width : number
        Width of the filter kernel.
    factor : number, optional
        Factor of oversampling. Default factor = 10.
    """
    def __init__(self, model=None, x_size=None, y_size=None, array=None, **kwargs):

        # Initialize from model
        if array is None:
            if self._model is None:
                raise TypeError("Must specify either array or model.")

            if x_size is None:
                x_size = self._default_size
            elif x_size != int(x_size):
                raise TypeError("x_size should be an integer")

            if y_size is None:
                y_size = x_size
            elif y_size != int(y_size):
                raise TypeError("y_size should be an integer")

            # Set ranges where to evaluate the model

            if x_size % 2 == 0:  # even kernel
                x_range = (-(int(x_size)) // 2 + 0.5, (int(x_size)) // 2 + 0.5)
            else:  # odd kernel
                x_range = (-(int(x_size) - 1) // 2, (int(x_size) - 1) // 2 + 1)

            if y_size % 2 == 0:  # even kernel
                y_range = (-(int(y_size)) // 2 + 0.5, (int(y_size)) // 2 + 0.5)
            else:  # odd kernel
                y_range = (-(int(y_size) - 1) // 2, (int(y_size) - 1) // 2 + 1)

            array = discretize_model(self._model, x_range, y_range, **kwargs)

        # Initialize from array
        elif array is not None:
            self._model = None

        super(Kernel2D, self).__init__(array)


def kernel_arithmetics(kernel, value, operation):
    """
    Add, subtract or multiply two kernels.

    Parameters
    ----------
    kernel : `astropy.convolution.Kernel`
        Kernel instance
    value : kernel, float or int
        Value to operate with
    operation : {'add', 'sub', 'mul'}
        One of the following operations:
            * 'add'
                Add two kernels
            * 'sub'
                Subtract two kernels
            * 'mul'
                Multiply kernel with number or convolve two kernels.
    """
    # 1D kernels
    if isinstance(kernel, Kernel1D) and isinstance(value, Kernel1D):
        if operation == "add":
            new_array = add_kernel_arrays_1D(kernel.array, value.array)
        if operation == "sub":
            new_array = add_kernel_arrays_1D(kernel.array, -value.array)
        if operation == "mul":
            raise Exception("Kernel operation not supported. Maybe you want "
                            "to use convolve(kernel1, kernel2) instead.")
        new_kernel = Kernel1D(array=new_array)
        new_kernel._separable = kernel._separable and value._separable
        new_kernel._is_bool = kernel._is_bool or value._is_bool

    # 2D kernels
    elif isinstance(kernel, Kernel2D) and isinstance(value, Kernel2D):
        if operation == "add":
            new_array = add_kernel_arrays_2D(kernel.array, value.array)
        if operation == "sub":
            new_array = add_kernel_arrays_2D(kernel.array, -value.array)
        if operation == "mul":
            raise Exception("Kernel operation not supported. Maybe you want "
                            "to use convolve(kernel1, kernel2) instead.")
        new_kernel = Kernel2D(array=new_array)
        new_kernel._separable = kernel._separable and value._separable
        new_kernel._is_bool = kernel._is_bool or value._is_bool

    # kernel and number
    elif ((isinstance(kernel, Kernel1D) or isinstance(kernel, Kernel2D))
        and np.isscalar(value)):
        if operation == "mul":
            new_kernel = copy.copy(kernel)
            new_kernel._array *= value
            new_kernel._normalization /= value
        else:
            raise Exception("Kernel operation not supported.")
    else:
        raise Exception("Kernel operation not supported.")
    return new_kernel
