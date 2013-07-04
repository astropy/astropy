from .core import Kernel1D, Kernel2D, Kernel
from .utils import KernelSizeError, NormalizationError

from ...modeling.models import Gaussian1DModel
from ...modeling.core import ParametricModel  # Box1DModel
#from ...modeling.models import Disk2DModel, Box2DModel

import numpy as np

__all__ = ['GaussianKernel', 'CustomKernel']


class GaussianKernel(Kernel1D):
    """
    Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifact.
    """
    _separable = True
    _weighted = True

    def __init__(self, width):
        amplitude = 1. / (np.sqrt(2 * np.pi) * width)
        self._model = Gaussian1DModel(amplitude=amplitude, mean=0.,
                                                 stddev=width)
        width = np.ceil(width)
        x = np.arange(-4 * width, 4 * width + 1)
        super(GaussianKernel, self).__init__(self.model(x))


class BoxKernel(Kernel1D):
    """
    Box filter kernel.

    The Box filter or running mean is a smoothing filter. It is not isotropic 
    and can produce artifact, when applied repeatedly to the same data. It is
    faster than a Gaussian smoothing filter.
    """
    _separable = True

    def __init__(self, width):
        amplitude = 1. / width
        #self._model = Box1DModel(amplitude=amplitude, mean=0.,
        #                                       stddev=width)


class Tophat2DKernel(Kernel2D):
    """
    Tophat filter kernel.

    The Tophat filter is an isotropic smoothing filter. It can produce artifact,
    when applied repeatedly on the same data.
    """

    def __init__(self, radius):
        amplitude = 1. / (np.pi * radius ** 2)
        #self.model = Disk2DModel(amplitude=amplitude, mean=0., radius=radius)


class Ring2DKernel(Kernel2D):
    """
    Ring filter kernel.

    The Ring filter kernel is the difference between two Tophat kernels of
    different width. It is useful for e.g background estimation.
    """
    def __init__(self, radius_in, radius_out):
        pass
        #self.model = (Disk2DModel(amplitude=1., mean=0., radius=radius_out)
        #                     - Disk2DModel(amplitude=1., mean=0., radius=radius_in))


class Trapezoid1DKernel(Kernel1D):
    """
    1D trapezoid kernel.
    """
    _weighted = True

    def __init__(self):
        pass
        #self.model = Trapezoid1DModel()


class MexicanHat1DKernel(Kernel1D):
    """
    1D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a background free smoothing
    filter. It does not conserve the mean. It is useful for peak or multi-scale
    detection.
    """
    _weighted = True

    def __init__(self):
        pass
        #self.model = MexicanHat1DModel()


class MexicanHat2DKernel(Kernel2D):
    """
    2D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a background free smoothing
    filter. It does not conserve the mean. It is useful for peak or multi-scale
    detection.
    """
    def __init__(self):
        pass
        #self.model = MexicanHat2DModel()


class Airy2DKernel(Kernel2D):
    def __init__(self):
        self._separable = False
        #self.model = Airy2DModel()


class Model1DKernel(Kernel1D):
    """
    Initialize kernel from astropy Parametric1DModel
    """
    def __init__(self, model):
        self._separable = False
        if isinstance(model, ParametricModel):
            self._model = model
        else:
            raise TypeError("Must be Parametric1DModel")


class Model2DKernel(Kernel2D):
    """
    Initialize kernel from astropy Parametric2DModel
    """
    def __init__(self, model):
        self._separable = False
        if isinstance(model, ParametricModel):
            self._model = model
        else:
            raise TypeError("Must be Parametric2DModel")


class PSFKernel(Kernel2D):
    """
    Initialize filter kernel from astropy PSF instance.
    """
    _separable = False

    def __init__(self):
        raise NotImplementedError('Not yet implemented')


class CustomKernel(Kernel):
    """
    Initialize filter custom kernel from mask.
    """
    
    def __init__(self, mask):
        # Pass 'None' because mask is overridden in the next line
        super(GaussianKernel, self).__init__(None)
        self.mask = mask

    @property
    def mask(self):
        """
        Filter kernel mask.
        """
        return self._mask

    @mask.setter
    def mask(self, mask):
        """
        Filter kernel mask setter
        """
        if isinstance(mask, np.ndarray):
            self._mask = mask
        elif isinstance(mask, list):
            self._mask = np.array(mask)
        else:
            raise TypeError("Must be list or array.")

        #Check if mask is odd in all axis
        self._odd = np.all([axes_size % 2 != 0 for axes_size in self.shape])

        if not self.odd:
            raise KernelSizeError("Kernel size must be odd.")

        # Check if mask is weighted
        ones = self._mask == 1.
        zeros = self._mask == 0
        self._weighted = not np.all(np.logical_or(ones, zeros))

        # Set normalization
        self._normalization = 1. / self._mask.sum()
