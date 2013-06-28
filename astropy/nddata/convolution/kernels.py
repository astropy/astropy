from .core import Kernel1D, Kernel2D
from .utils import KernelSizeError

from ...modeling.models import Gaussian1DModel, Box1DModel
from ...modeling.models import Disk2DModel, Box2DModel

import numpy as np


class Laplace1DKernel(Kernel1D):
    """
    Laplace filter kernel.

    The Laplace filter is a derivative filter. It is useful for e.g. detection
    of fast changing structures in the data e.g edges.
    """
    def __init__(self):
        self._separable = False


class Laplace2DKernel(Kernel2D):
    """
    Laplace filter kernel.

    The Laplace filter is a derivative filter. It is useful for detection
    of fast changing structures in the data e.g edges.
    """
    def __init__(self):
        self._separable = False


class GaussianKernel(Kernel1D):
    """
    Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifact.
    """
    def __init__(self, width):
        self._separable = True
        self._normalized = True
        self._weighted = True
        self.model = Gaussian1DModel(amplitude=1., mean=0.,
                                                 stddev=width)


class Binomial1DKernel(Kernel1D):
    """
    Binomial filter kernel
    """
    def __init__(self, size):
        self._separable = True
        self._normalized = True
        self._weighted = True
        self.model = None
        self.mask = np.array([1, 2, 1])

    def _init_mask(self, size):
        """
        Init binomial filter mask
        """
        if size % 2 == 0:
            from scipy.special import binom
            mask = np.empty(size)
            for i in range(size):
                mask[i] = binom(size, i)
        else:
            raise KernelSizeError("Size must be odd.")


class BoxKernel(Kernel1D):
    """
    Box filter kernel.

    The Box filter or running mean is a smoothing filter. It is not isotropic 
    and can produce artifact, when applied repeatedly to the same data. It is
    faster than a Gaussian smoothing filter.
    """
    def __init__(self, width):
        self._separable = True
        self._normalized = True
        self._weighted = False
        self.model = Box1DModel(amplitude=1., mean=0.,
                                                 stddev=width)


class Tophat2DKernel(Kernel2D):
    """
    Tophat filter kernel.

    The Tophat filter is an isotropic smoothing filter. It can produce artifact,
    when applied repeatedly on the same data.
    """
    def __init__(self, radius):
        self._separable = False
        self._normalized = True
        self._weighted = False
        self.model = Disk2DModel(amplitude=1., mean=0., radius=radius)


class Ring2DKernel(Kernel2D):
    """
    Ring filter kernel.

    The Ring filter kernel is the difference between two Tophat kernels of
    different width. It is useful for e.g background estimation.
    """
    def __init__(self, radius_in, radius_out):
        self._separable = False
        self._normalized = True
        self._weighted = False
        self.model = (Disk2DModel(amplitude=1., mean=0., radius=radius_out)
                                 - Disk2DModel(amplitude=1., mean=0., radius=radius_in))


class Trapezoid1DKernel(Kernel1D):
    """
    1D trapezoid kernel.
    """
    def __init__(self):
        self._separable = False
        self._weighted = True
        self.model = Trapezoid1DModel()


class MexicanHat1DKernel(Kernel1D):
    """
    1D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a background free smoothing
    filter. It does not conserve the mean. It is useful for peak or multi-scale
    detection.
    """
    def __init__(self):
        self._separable = False
        self._normalized = False
        self._weighted = True
        self.model = MexicanHat1DModel()


class MexicanHat2DKernel(Kernel2D):
    """
    2D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a background free smoothing
    filter. It does not conserve the mean. It is useful for peak or multi-scale
    detection.
    """
    def __init__(self):
        self._separable = False
        self._normalized = False
        self.model = MexicanHat2DModel()

class Airy2DKernel(Kernel2D):
    def __init__(self):
        self._separable = False
        self.model = Airy2DModel()


class Model1DKernel(Kernel1D):
    """
    Initialize kernel from astropy ParametricModel
    """
    def __init__(self):
        self._separable = False


class PSFKernel(Kernel2D):
    """
    Initialize filter kernel from astropy PSF model.
    """
    def __init__(self):
        self._separable = False


class Custom1DKernel(Kernel1D):
    """
    Initialize filter kernel from mask.
    """
    def __init__(self, mask):
        self._separable = False
        self.model = None
        self._mask = np.array(mask)
        if not self.odd:
            raise KernelSizeError("Mask must be of odd size.")
