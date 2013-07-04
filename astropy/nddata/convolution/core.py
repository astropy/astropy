# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains the convolution and filter functionalities of astropy.
It is in preliminary state.

A few conceptual notes:
A filter kernel is mainly characterized by its response function. In the 1D
case we speak of "impulse response function", in the 2D case we call it
"point spread function". This response function is given for every kernel
by an astropy ParametricModel, which is evaluated on a grid to obtain a filter
mask, which can then be applied to binned data.

The model is centered on the mask and should have an amplitude such that the mask
integrates to one per default.
"""

import abc

import numpy as np
import copy


class Kernel(object):
    """
    Abstract convolution kernel class
    """
    __metaclass__ = abc.ABCMeta
    _odd = True
    _separable = False
    _weighted = False
    _model = None

    def __init__(self, mask):
        self._mask = mask
        self._normalization = 1. / self._mask.sum()
        self._truncation = 1. - 1. / self._normalization

    @property
    def truncation(self):
        """
        Difference of sum(mask) to 1.
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
        return len(self.mask.shape)

    @property
    def odd(self):
        """
        Check if kernel size is odd in all axes.
        """
        return self._odd

    @property
    def center(self):
        """
        Index of the kernels center.
        """
        return [axes_size / 2 for axes_size in self._mask.shape]

    @property
    def normalization(self):
        """
        Kernel normalization factor
        """
        return self._normalization

    @property
    def shape(self):
        """
        Shape of the kernel mask.
        """
        return self._mask.shape

    @property
    def separable(self):
        """
        Indicates if the filter kernel is separable.

        A 2D filter is separable, when its filter mask can be written as the
        outer product of two 1D masks.

        If a filter kernel is separable, higher dimension convolutions will be
        performed by applying the 1D filter mask consecutively on every dimension.
        This is significantly faster, than using a filter mask with the same
        dimension.
        """
        return self._separable

    @property
    def mask(self):
        """
        Filter kernel mask.
        """
        return self._mask

    def __add__(self, kernel):
        """
        Add two filter kernels and do a renormalization.
        """
        # As convolution is linear we can add two kernels
        # This will fail, if the kernel shapes don't match
        add_kernel = Kernel1D(self._mask + kernel._mask)
        add_kernel._separable = self._separable and kernel._separable
        add_kernel._weighted = self._weighted or kernel._weighted
        return add_kernel

    def __sub__(self,  kernel):
        """
        Subtract two filter kernels and do a renormalization.
        """
        # As convolution is linear we can subtract two kernels
        # This will fail, if the kernel shapes don't match
        sub_kernel = Kernel1D(self._mask - kernel._mask)
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
            pass
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
            pass
        else:
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(value)))


class Kernel1D(Kernel):
    """
    Base class for 1D filter kernels
    """
    def __init__(self, mask):
        self.axes = 0
        super(Kernel1D, self).__init__(mask)


class Kernel2D(Kernel):
    """
    Abstract base class for 1D filter kernels
    """
    def __init__(self):
        super(Kernel2D, self).__init__()

