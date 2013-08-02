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
import copy
import warnings

import numpy as np

from .utils import (discretize_model, add_kernel_arrays_1D,
                    add_kernel_arrays_2D)


class Kernel(object):
    """
    Convolution kernel base class.
    """
    _odd = True
    _separable = False
    _weighted = False
    _model = None

    def __init__(self, array):
        self._array = array
        self._normalization = 1. / self._array.sum()
        # The value of 100 is kind of arbitrary
        # there are kernel that sum to zero and
        # the user should be warned in this case
        if np.abs(self._normalization) > 100:
            warnings.warn("Normalization factor of kernel is" +
                                        "exceptionally large > 100.")

    @property
    def truncation(self):
        """
        Deviation from the normalization to one.
        """
        return self._truncation

    @property
    def weighted(self):
        """
        Indicates if kernel is weighted.

        If the kernel is not weighted the multiplication in the convolution will
        be omitted, to increase the performance.
        """
        return self._weighted

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
    def odd(self):
        """
        Indicates if kernel size is odd in all axes.
        """
        return self._odd

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
        if mode == 'integral':
            self._array *= self._normalization
            self._normalization = 1.
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
        if isinstance(self, Kernel1D) and isinstance(kernel, Kernel1D):
            # As convolution is linear we can add two kernels
            add_array = add_kernel_arrays_1D(self.array, kernel.array)
            add_kernel = Kernel1D(array=add_array)
        elif isinstance(self, Kernel2D) and isinstance(kernel, Kernel2D):
            # As convolution is linear we can add two kernels
            add_array = add_kernel_arrays_2D(self.array, kernel.array)
            add_kernel = Kernel2D(array=add_array)
        else:
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(kernel)))
        add_kernel._separable = self._separable and kernel._separable
        add_kernel._weighted = self._weighted or kernel._weighted
        return add_kernel

    def __sub__(self,  kernel):
        """
        Subtract two filter kernels.
        """
        if isinstance(self, Kernel1D) and isinstance(kernel, Kernel1D):
            # As convolution is linear we can add two kernels
            sub_array = add_kernel_arrays_1D(self.array, -kernel.array)
            sub_kernel = Kernel1D(array=sub_array)
        elif isinstance(self, Kernel2D) and isinstance(kernel, Kernel2D):
            # As convolution is linear we can add two kernels
            sub_array = add_kernel_arrays_2D(self.array, -kernel.array)
            sub_kernel = Kernel2D(array=sub_array)
        else:
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(kernel)))

        # As convolution is linear we can subtract two kernels
        sub_kernel._separable = self._separable and kernel._separable
        sub_kernel._weighted = self._weighted or kernel._weighted
        return sub_kernel

    def __mul__(self, value):
        """
        Multiply kernel with number or convolve two kernels.
        """
        # As convolution is linear we can multiply with a scalar
        if isinstance(value, (int, float)):
            kernel = copy.copy(self)
            kernel._array *= value
            kernel._normalization /= value
            return kernel
        else:
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(value)))

    def __rmul__(self, value):
        """
        Multiply kernel with number or convolve two kernels.
        """
        # As convolution is linear we can multiply with a scalar
        if isinstance(value, (int, float)):
            kernel = copy.copy(self)
            kernel._array *= value
            kernel._normalization /= value
            return kernel
        else:
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(value)))


class Kernel1D(Kernel):
    """
    Base class for 1D filter kernels

    Parameters
    ----------
    width : float
        Width of the kernel model.
    size : odd int
        Size of the kernel array.
    """
    def __init__(self, model=None, x_size=None, mode='center', array=None):
        if array == None:
            if model != None:
                self._model = model
            if x_size == None:
                x_size = self._default_size
            x_range = (-np.rint(x_size / 2), np.rint(x_size / 2) + 1)
            array = discretize_model(self._model, x_range, mode)
        else:
            self._model = None
        super(Kernel1D, self).__init__(array)


class Kernel2D(Kernel):
    """
    Abstract base class for 1D filter kernels
    """
    def __init__(self, x_size=None, y_size=None, mode='center', array=None):
        if array == None:
            if x_size == None:
                x_size = self._default_size
            if y_size == None:
                y_size = x_size

            # Set ranges where to evaluate the model
            x_range = (-np.rint(x_size / 2), np.rint(x_size / 2) + 1)
            y_range = (-np.rint(y_size / 2), np.rint(y_size / 2) + 1)
            array = discretize_model(self._model, x_range, y_range, mode)
        else:
            self._model = None
        super(Kernel2D, self).__init__(array)
