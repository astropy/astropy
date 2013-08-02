from __future__ import division
import warnings

import numpy as np

from .core import Kernel1D, Kernel2D, Kernel
from .utils import KernelSizeError
from ...modeling.functional_models import *
from ...modeling.core import Parametric1DModel, Parametric2DModel

__all__ = sorted(['Gaussian1DKernel', 'Gaussian2DKernel', 'CustomKernel',
                  'Box1DKernel', 'Box2DKernel', 'Tophat2DKernel',
                  'Trapezoid1DKernel', 'MexicanHat1DKernel',
                  'MexicanHat2DKernel', 'AiryDisk2DKernel',
                  'Model1DKernel', 'Model2DKernel',
                  'TrapezoidDisk2DKernel'])


class Gaussian1DKernel(Kernel1D):
    """
    Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifact.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    x_size : odd int
        Size of the kernel array. Default = 8 * width.

    See Also
    --------
    Box1DKernel, Trapezoid1DKernel, MexicanHat1DKernel

    """
    _separable = True
    _weighted = True

    def __init__(self, width, **kwargs):
        self._model = Gaussian1DModel(1. / (np.sqrt(2 * np.pi) * width), 0, width)
        self._default_size = 8 * width
        super(Gaussian1DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(1. - 1 / self._normalization)


class Gaussian2DKernel(Kernel2D):
    """
    Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifact.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    x_size : odd int
        Size in x direction of the kernel array. Default = 8 * width.
    y_size : odd int
        Size in y direction of the kernel array. Default = 8 * width.

    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel, 
    TrapezoidDisk2DKernel, AiryDisk2DKernel

    """
    _separable = True
    _weighted = True

    def __init__(self, width, **kwargs):
        self._model = Gaussian2DModel(1. / (2 * np.pi * width ** 2), 0, 0, width, width)
        self._default_size = 8 * width
        super(Gaussian2DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(1. - 1 / self._normalization)


class Box1DKernel(Kernel1D):
    """
    Box filter kernel.

    The Box filter or running mean is a smoothing filter. It is not isotropic
    and can produce artifact, when applied repeatedly to the same data. It is
    faster than a Gaussian smoothing filter.

    Parameters
    ----------
    width : number
        Width of the filter kernel.

    See Also
    --------
    Gaussian1DKernel, Trapezoid1DKernel, MexicanHat1DKernel
    """
    _separable = True
    _weighted = False

    def __init__(self, width, **kwargs):
        self._model = Box1DModel(1. / width, 0, width)
        self._default_size = width
        super(Box1DKernel, self).__init__(**kwargs)
        self._truncation = 0
        self.normalize()


class Box2DKernel(Kernel2D):
    """
    Box filter kernel.

    The Box filter or running mean is a smoothing filter. It is not isotropic
    and can produce artifact, when applied repeatedly to the same data. It is
    faster than a Gaussian smoothing filter.

    Parameters
    ----------
    width : number
        Width of the filter kernel.

    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel, 
    TrapezoidDisk2DKernel, AiryDisk2DKernel
    """
    _separable = True
    _weighted = False

    def __init__(self, width, **kwargs):
        self._model = Box2DModel(1. / width ** 2, 0, 0, width, width)
        self._default_size = width
        super(Box2DKernel, self).__init__(**kwargs)
        self._truncation = 0


class Tophat2DKernel(Kernel2D):
    """
    Tophat filter kernel.

    The Tophat filter is an isotropic smoothing filter. It can produce artifact,
    when applied repeatedly on the same data.

    Parameters
    ----------
    radius : int
        Radius of the filter kernel.

    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel, 
    TrapezoidDisk2DKernel, AiryDisk2DKernel
    """

    def __init__(self, radius, **kwargs):
        self._model = Disk2DModel(1. / (np.pi * radius ** 2), 0, 0, radius)
        self._default_size = 2 * radius
        super(Tophat2DKernel, self).__init__(**kwargs)
        self._truncation = 0


class Ring2DKernel(Kernel2D):
    """
    Ring filter kernel.

    The Ring filter kernel is the difference between two Tophat kernels of
    different width. It is useful for e.g background estimation.
    """
    def __init__(self, radius_in, radius_out):
        raise NotImplementedError


class Trapezoid1DKernel(Kernel1D):
    """
    1D trapezoid kernel.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    slope : number
        Slope of the filter kernels tails


    See Also
    --------
    Box1DKernel, Gaussian1DKernel, MexicanHat1DKernel
    """
    _weighted = True

    def __init__(self, width, slope=1., **kwargs):
        self._model = Trapezoid1DModel(1, 0, width, slope)
        self._default_size = width + 2. / slope
        super(Trapezoid1DKernel, self).__init__(**kwargs)
        self._truncation = 0


class TrapezoidDisk2DKernel(Kernel2D):
    """
    2D trapezoid kernel.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    slope : number
        Slope of the filter kernels tails


    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel, 
    TrapezoidDisk2DKernel, AiryDisk2DKernel
    """
    _weighted = True

    def __init__(self, width, slope=1., **kwargs):
        self._model = TrapezoidDisk2DModel(1, 0, 0, width, slope)
        self._default_size = width + 2. / slope
        super(TrapezoidDisk2DKernel, self).__init__(**kwargs)
        self._truncation = 0


class MexicanHat1DKernel(Kernel1D):
    """
    1D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a background free smoothing
    filter. It does not conserve the mean. It is useful for peak or multi-scale
    detection.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    x_size : odd int
        Size in x direction of the kernel array. Default = 8 * width.

    See Also
    --------
    Box1DKernel, Gaussian1DKernel, Trapezoid1DKernel
    """
    _weighted = True

    def __init__(self, width, **kwargs):
        self._default_size = 8 * width
        self._model = MexicanHat1DModel(1, 0, width)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            super(MexicanHat1DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(self._array.sum() / self._array.size)
        self._normalization = 0


class MexicanHat2DKernel(Kernel2D):
    """
    2D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a background free smoothing
    filter. It does not conserve the mean. It is useful for peak or multi-scale
    detection.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    x_size : odd int
        Size in x direction of the kernel array. Default = 8 * width.
    y_size : odd int
        Size in y direction of the kernel array. Default = 8 * width.

    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel
    """
    _weighted = True

    def __init__(self, width, **kwargs):
        self._default_size = 8 * width
        self._model = MexicanHat2DModel(1, 0, 0, width)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            super(MexicanHat2DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(self._array.sum() / self._array.size)
        self._normalization = 0


class AiryDisk2DKernel(Kernel2D):
    """
    Airy disk 2D kernel.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    x_size : odd int
        Size in x direction of the kernel array. Default = 8 * width.
    y_size : odd int
        Size in y direction of the kernel array. Default = 8 * width.

    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel, 
    TrapezoidDisk2DKernel, AiryDisk2DKernel
    """
    _weighted = True

    def __init__(self, width, **kwargs):
        self._default_size = 8 * width
        self._model = AiryDisk2DModel(1, 0, 0, width)
        super(AiryDisk2DKernel, self).__init__(**kwargs)


class Model1DKernel(Kernel1D):
    """
    Create kernel from astropy.models.Parametric1DModel.

    Parameters
    ----------
    model : Parametric1DModel
        Kernel response function model

    Raises
    ------
    TypeError
        If model is not an instance of astropy.models.Parametric1DModel

    See also
    --------
    Model2DKernel : Create kernel from astropy.models.Parametric2DModel
    CustomKernel : Create kernel from list or array
    """
    _separable = False
    _weighted = True

    def __init__(self, model, **kwargs):
        if isinstance(model, Parametric1DModel):
            self._model = model
        else:
            raise TypeError("Must be Parametric1DModel")
        super(Model1DKernel, self).__init__(**kwargs)


class Model2DKernel(Kernel2D):
    """
    Create kernel from astropy.models.Parametric2DModel.

    Parameters
    ----------
    model : Parametric2DModel
        Kernel response function model

    Raises
    ------
    TypeError
        If model is not an instance of astropy.models.Parametric2DModel

    See also
    --------
    Model1DKernel : Create kernel from astropy.models.Parametric1DModel
    CustomKernel : Create kernel from list or array
    """
    _weighted = True
    _separable = False

    def __init__(self, model, **kwargs):
        self._separable = False
        if isinstance(model, Parametric2DModel):
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
    Create filter kernel from list or array.

    Parameters
    ----------
    array : list or array
        Filter kernel array. Size must be odd.

    Raises
    ------
    TypeError
        If array is not a list or array.
    KernelSizeError
        If array size is even.

    See also
    --------
    Model2DKernel, Model1DKernel

    Examples
    --------
    Define one dimensional array:

        >>> from astropy.nddata.convolution.kernels import CustomKernel
        >>> import numpy as np
        >>> array = np.array([1, 2, 3, 2, 1])
        >>> kernel = CustomKernel(array)
        >>> kernel.dimension
        1

    Define two dimensional array:

        >>> array = np.array([[1, 1, 1], [1, 2, 1], [1, 1, 1]])
        >>> kernel = CustomKernel(array)
        >>> kernel.dimension
        2
    """

    def __init__(self, array):
        self.array = array
        super(CustomKernel, self).__init__(self._array)

    @property
    def array(self):
        """
        Filter kernel array.
        """
        return self._array

    @array.setter
    def array(self, array):
        """
        Filter kernel array setter
        """
        if isinstance(array, np.ndarray):
            self._array = array
        elif isinstance(array, list):
            self._array = np.array(array)
        else:
            raise TypeError("Must be list or array.")

        #Check if array is odd in all axis
        self._odd = np.all([axes_size % 2 != 0 for axes_size in self.shape])

        if not self.odd:
            raise KernelSizeError("Kernel size must be odd in all axes.")

        # Check if array is weighted
        ones = self._array == 1.
        zeros = self._array == 0
        self._weighted = not np.all(np.logical_or(ones, zeros))

        # Set normalization
        self._normalization = 1. / self._array.sum()
