# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains the convolution and filter functionalities of astropy.

A few conceptual notes:
A filter kernel is mainly characterized by its response function. In the 1D
case we speak of "impulse response function", in the 2D case we call it
"point spread function". This response function is given for every kernel
by an astropy ParametricModel, which is evaluated on a grid to obtain a filter
array, which can then be applied to binned data.

The model is centered on the array and should have an amplitude such that the array
integrates to one per default.

Currently only symmetric 2D kernels are supported.
"""
from __future__ import division
import warnings
import copy

import numpy as np

from .utils import (discretize_model, add_kernel_arrays_1D,
                    add_kernel_arrays_2D)

MAX_NORMALIZATION = 100


class Kernel(object):
    """
    Convolution kernel base class.

    Parameters
    ----------
    array : ndarray
        Kernel array.
    """
    _separable = False
    _is_bool = True
    _model = None

    def __init__(self, array):
        self._array = array
        self._normalization = 1. / self._array.sum()

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

    @property
    def normalization(self):
        """
        Kernel normalization factor
        """
        return self._normalization

    def normalize(self, mode='integral'):
        """
        Force normalization of filter kernel.

        Parameters
        ----------
        mode : string
            One of the following modes:
                * 'integral' (default)
                    Kernel normalized such that its integral = 1.
                * 'peak'
                    Kernel normalized such that its peak = 1.
        """
        # There are kernel that sum to zero and
        # the user should be warned in this case
        if np.abs(self._normalization) > MAX_NORMALIZATION:
            warnings.warn("Normalization factor of kernel is "
                          "exceptionally large > {0}.".format(MAX_NORMALIZATION))
        if mode == 'integral':
            self._array *= self._normalization
        if mode == 'peak':
            self._array /= self._array.max()
            self._normalization = 1. / self._array.sum()

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


class Kernel1D(Kernel):
    """
    Base class for 1D filter kernels.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    x_size : odd int, optional
        Size of the kernel array. Default = 8 * width.
    mode: string, optional
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
        if array == None:
            if model != None:
                self._model = model
            if x_size == None:
                x_size = self._default_size
            x_range = (-np.rint(x_size / 2), np.rint(x_size / 2) + 1)
            array = discretize_model(self._model, x_range, **kwargs)

        # Initialize from array
        elif array != None:
            self._model = None
        else:
            raise TypeError("Must specify either array or model.")
        super(Kernel1D, self).__init__(array)


class Kernel2D(Kernel):
    """
    Base class for 2D filter kernels.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    x_size : odd int, optional
        Size in x direction of the kernel array. Default = 8 * width.
    y_size : odd int, optional
        Size in y direction of the kernel array. Default = 8 * width.
    mode: string, optional
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
    factor : number, optional
        Factor of oversampling. Default factor = 10.
    """
    def __init__(self, model=None, x_size=None, y_size=None, array=None, **kwargs):
        # Initialize from model
        if array == None:
            if x_size == None:
                x_size = self._default_size
            if y_size == None:
                y_size = x_size
            # Set ranges where to evaluate the model
            x_range = (-np.rint(x_size / 2), np.rint(x_size / 2) + 1)
            y_range = (-np.rint(y_size / 2), np.rint(y_size / 2) + 1)
            array = discretize_model(self._model, x_range, y_range, **kwargs)

        # Initialize from array
        elif array != None:
            self._model = None
        else:
            raise TypeError("Must specify either array or model.")
        super(Kernel2D, self).__init__(array)


def kernel_arithmetics(kernel, value, operation):
    """
    Add, subtract or multiply two kernels.

    Parameters
    ----------
    kernel : astropy.nddata.convolution.kernel
        Kernel instance
    values : kernel, float or int
        Value to operate with
    operation : string
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
            raise Exception("Kernel operation not supported. Maybe you want to"
                             "use convolve(kernel1, kernel2) instead.")
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
            raise Exception("Kernel operation not supported. Maybe you want to"
                            "use convolve(kernel1, kernel2) instead.")
        new_kernel = Kernel2D(array=new_array)
        new_kernel._separable = kernel._separable and value._separable
        new_kernel._is_bool = kernel._is_bool or value._is_bool

    # kernel and number
    elif ((isinstance(kernel, Kernel1D) or isinstance(kernel, Kernel2D))
        and isinstance(value, (int, float))):
        if operation == "mul":
            new_kernel = copy.copy(kernel)
            new_kernel._array *= value
            new_kernel._normalization /= value
        else:
            raise Exception("Kernel operation not supported.")
    else:
        raise Exception("Kernel operation not supported.")
    return new_kernel
