# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains the convolution and filter functionalities of astropy.
It is in preliminary state. It tries to combine the concepts of astropy models
and PSF classes.

A few conceptual notes:
A filter kernel is mainly characterized by its response function. In the 1D
case we speak of "impulse response function", in the 2D case we call it
"point spread function". This response function is given for every kernel
by an astropy ParametricModel, which is evaluated on a grid to obtain a filter
mask, which can then be applied to binned data.

The model is centered on the mask and has an amplitude such that the mask
integrate to ones per default.
"""

import abc
import numpy as np

from .utils import NormalizationError
from ...modeling.core import Parametric1DModel, Parametric2DModel


class Kernel(object):
    """
    Abstract convolution kernel class
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self._mask
        self._normalized
        self._normalization
        self._separable
        self._model
        self._maskfft
        self.truncation = 1E-5  # Lower limit at which the mask is truncated
        self.pixsize
        # This attributes and functions could be useful
        # self.transfer_function
        # self.symmetric
        # self.antisymmetric

    @property
    def truncation(self):
        """
        Lower limit at which the mask is truncated
        """
        return self._truncation

    @truncation.setter
    def truncation(self, value):
        """
        Lower limit at which the mask is truncated
        """
        self._truncation = value
        # Note: This is preliminary because it can create "holes" in the mask,
        # if the mask contains smaller values than truncation, not only at the edges.
        # A correct way would be to measure the containment fraction and assure
        # that it is > 1 - truncation.
        self._greater_truncation = np.abs(self._mask) > self._truncation

    @property
    def weighted(self):
        """
        Indicates if kernel is weighted
        """
        if getattr(self, "_weighted", False):
            return self._weighted
        else:
            # Check whether mask contains only ones and zeros
            ones = self.mask == 1
            zeros = self.mask == 0
            weighted = not np.all(np.logical_or(ones, zeros))
            setattr(self, "_weighted", weighted)
            return weighted

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
        odd = [True for axes_size in self.mask.shape if axes_size % 2 != 0]
        return np.all(odd)

    @property
    def center(self):
        """
        Index of the kernels center.
        """
        return [axes_size / 2 for axes_size in self.mask.shape]

    @property
    def normalization(self):
        """
        Kernel normalization factor
        """
        self._normalization = self._mask.sum()
        return self._normalization

    @property
    def normalized(self):
        """
        Check if the kernel mask is normalized.
        """
        if self._normalization > (1 - self._truncation):
            self._normalized = True
        else:
            self._normalized = False
        return self._normalized

    @normalized.setter
    def normalized(self, value):
        """
        Set kernel normalization.
        """
        if abs(self._normalization) < self._truncation:
            raise NormalizationError("Can't normalize {0} filter kernel"
                                     .format(__class__.__name__))
        if value and not self._normalized:
            self.mask = self.mask / self._normalization
            self._normalized = value
        elif value and self._normalized:
            raise NormalizationError("Kernel is already normalized")
        elif not value and self._normalized:
            raise NormalizationError("Can't undo normalization")

    @property
    def shape(self):
        """
        Shape of the kernel mask.
        """
        return self.mask.shape

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
        return self._mask[self._greater_truncat]

    @property
    def model(self):
        """
        Kernel response model.
        """
        return self._model

    @abc.abstractmethod
    @model.setter
    def model(self):
        """
        Evaluate kernel response model.

        I a kernel has no model it should be None.
        """
        return NotImplementedError("Subclasses should implement this")

    @mask.setter
    def mask(self, array):
        """
        Filter kernel mask.
        """
        if isinstance(array, np.array):
            self._mask = array
        elif isinstance(array, list):
            self._mask = np.array(list)
        else:
            raise TypeError("Must be list or array.")

    def _init_mask(self):
        """
        Initialize mask from response function.
        """
        pass

    def __add__(self, kernel):
        """
        Add two filter kernels and do a renormalization.
        """
        # As convolution is linear we can add two kernels
        if self.mask.shape == kernel.mask.shape:
            self.mask = self.mask + kernel.mask
        else:
            pass  # new_shape = np.maximum(self.mask.shape, kernel.mask.shape)
        self.normalized = True

    def __sub__(self,  kernel):
        """
        Subtract two filter kernels and do a renormalization.
        """
        # As convolution is linear we can subtract two kernels
        if self.mask.shape == kernel.mask.shape:
            self.mask = self.mask - kernel.mask
        else:
            pass
        self.normalized = True

    def __mul__(self, value):
        """
        Multiply kernel with number or convolve two kernels.
        """
        # As convolution is linear we can multiply with a scalar
        if isinstance(value, (int, float)):
            self._mask = self._normalisation * value
        elif isinstance(value, Kernel):
            #Convolve the two kernels with each other
            pass
        else:
            raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'"
                            .format(self.__class__, type(value)))


class Kernel1D(Kernel):
    """
    Abstract base class for 1D filter kernels
    """
    def __init__(self):
        self.axes = 0

    @model.setter
    def model(self, model):
        """
        """
        if isinstance(model, Parametric1DModel):
            self._model = model
        else:
            TypeError("Must be Parametric1DModel")


class Kernel2D(Kernel):
    """
    Abstract base class for 1D filter kernels
    """
    def __init__(self):
        pass

    @model.setter
    def model(self, model):
        if isinstance(model, Parametric2DModel):
            self._model = model
        else:
            TypeError("Must be Parametric2DModel")
