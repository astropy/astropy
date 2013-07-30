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
"""

import numpy as np
import copy
from astropy.nddata.convolution.utils import KernelSizeError


class Kernel(object):
    """
    Abstract convolution kernel class
    """
    _odd = True
    _separable = False
    _weighted = False
    _model = None

    def __init__(self, array):
        self._array = array
        self._normalization = 1. / self._array.sum()
        if self._normalization > 0.001:
            self._truncation = 1. - 1. / self._normalization
        else:
            # There kernels that are normalized to zero
            self._truncation = self._array.sum()

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
        return [axes_size / 2 for axes_size in self._array.shape]

    @property
    def normalization(self):
        """
        Kernel normalization factor
        """
        return self._normalization

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
    def __init__(self, size):
        self.axes = 0
        if not isinstance(size, int):
            raise TypeError("Size must be integer.")
        if size % 2 == 0:
            raise KernelSizeError("Kernel size must be odd.")
        self._odd = True
        array = self._init_array(size)
        super(Kernel1D, self).__init__(array)

    def _init_array(self, size):
        """
        Evaluate kernel model on grid.

        Parameters
        ----------
        size : odd int
            Total size of the array.
        """
        x = np.arange(-(size / 2), (size / 2) + 1)
        return self._model(x)

    def __add__(self, kernel):
        """
        Add two filter kernels and do a renormalization.
        """
        # As convolution is linear we can add two kernels
        # This will fail, if the kernel shapes don't match
        add_kernel = Kernel1D(self._array + kernel._array)
        add_kernel._separable = self._separable and kernel._separable
        add_kernel._weighted = self._weighted or kernel._weighted
        return add_kernel

    def __sub__(self,  kernel):
        """
        Subtract two filter kernels and do a renormalization.
        """
        # As convolution is linear we can subtract two kernels
        # This will fail, if the kernel shapes don't match
        sub_kernel = Kernel1D(self._array - kernel._array)
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
            kernel._normalization *= value
            return kernel
        elif isinstance(value, Kernel):
            #Convolve the two kernels with each other
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(value)))
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
            kernel._normalization *= value
            return kernel
        elif isinstance(value, Kernel):
            #Convolve the two kernels with each other
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(value)))
        else:
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(value)))


class Kernel2D(Kernel):
    """
    Abstract base class for 1D filter kernels
    """
    def __init__(self, shape):
        if not np.all([isinstance(number, int) for number in shape]):
            raise TypeError("Size must be integer.")
        if np.any([number % 2 == 0 for number in shape]):
            raise KernelSizeError("Kernel size must be odd.")
        self._odd = True
        array = self._init_array(shape)
        super(Kernel2D, self).__init__(array)

    def _init_array(self, shape):
        """
        Evaluate kernel model on grid.

        Parameters
        ----------
        shape : tuple
            Shape of the array. Sizes must be odd.
        """
        x = np.arange(-(shape[1] / 2), (shape[1] / 2) + 1)
        y = np.arange(-(shape[0] / 2), (shape[0] / 2) + 1)
        x, y = np.meshgrid(x, y)
        return self._model(x, y)


