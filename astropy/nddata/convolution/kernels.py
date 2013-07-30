import numpy as np

from .core import Kernel1D, Kernel2D, Kernel
from .utils import KernelSizeError

from ...modeling.models import *
from ...modeling.core import Parametric1DModel, Parametric2DModel
from astropy.modeling.functional_models import Trapezoid1DModel


__all__ = ['GaussianKernel', 'CustomKernel', 'BoxKernel', 'Tophat2DKernel',
           'Trapezoid1DKernel', 'MexicanHat1DKernel', 'MexicanHat2DKernel',
           'AiryDisk2DKernel']


class GaussianKernel(Kernel1D):
    """
    Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifact.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    size : odd int
        Size of the kernel mask. Default = 8 * width + 1

    See Also
    --------
    BoxKernel : Box filter kernel.

    """
    _separable = True
    _weighted = True

    def __init__(self, width, size=None):
        if size == None:
            #Default Kernel size for GaussianKernel
            size = 8 * int(width) + 1
        amplitude = 1. / (np.sqrt(2 * np.pi) * width)
        self._model = Gaussian1DModel(amplitude=amplitude, mean=0.,
                                                 stddev=width)
        super(GaussianKernel, self).__init__(size)


class BoxKernel(Kernel1D):
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
    GaussianKernel : Gaussian filter kernel.
    """
    _separable = True
    _weighted = False

    def __init__(self, width):
        self._model = Box1DModel(amplitude=1., x_0=0., width=width)
        super(BoxKernel, self).__init__(width)
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
    """

    def __init__(self, radius):
        self._model = Disk2DModel(amplitude=1., x_0=0., y_0=0., R_0=radius)
        super(Tophat2DKernel, self).__init__([2 * int(radius) + 1, 2 * int(radius) + 1])
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
    """
    _weighted = True

    def __init__(self, width, slope=1):
        self._model = Trapezoid1DModel(amplitude=1, x_0=0., width=width, slope=slope)
        size = 2 * int(width / 2. + 1. / slope) + 1 
        super(Trapezoid1DKernel, self).__init__(size)
        self._truncation = 0


class MexicanHat1DKernel(Kernel1D):
    """
    1D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a background free smoothing
    filter. It does not conserve the mean. It is useful for peak or multi-scale
    detection.
    """
    _weighted = True

    def __init__(self, width, size=None):
        if size == None:
            #Default Kernel size for GaussianKernel
            size = 8 * int(width) + 1
        self._model = MexicanHat1DModel(amplitude=1, x_0=0., sigma=width)
        super(MexicanHat1DKernel, self).__init__(size)


class MexicanHat2DKernel(Kernel2D):
    """
    2D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a background free smoothing
    filter. It does not conserve the mean. It is useful for peak or multi-scale
    detection.
    """
    def __init__(self, width, shape=None):
        if shape == None:
            #Default Kernel size for GaussianKernel
            shape = (8 * int(width) + 1, 8 * int(width) + 1)
        self._model = MexicanHat2DModel(amplitude=1, x_0=0., y_0=0, sigma=width)
        super(MexicanHat2DKernel, self).__init__(shape)


class AiryDisk2DKernel(Kernel2D):
    """
    Airy 2D kernel.
    """
    def __init__(self, width, shape=None):
        if shape == None:
            #Default Kernel size for GaussianKernel
            shape = (8 * int(width) + 1, 8 * int(width) + 1)
        self._model = AiryDisk2DModel(amplitude=1, x_0=0., y_0=0, width=width)
        super(AiryDisk2DKernel, self).__init__(shape)


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

    def __init__(self, model):
        if isinstance(model, Parametric1DModel):
            self._model = model
        else:
            raise TypeError("Must be Parametric1DModel")


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
    def __init__(self, model):
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
    mask : list or array
        Filter kernel mask. Size must be odd.

    Raises
    ------
    TypeError
        If mask is not a list or array.
    KernelSizeError
        If mask size is even.

    See also
    --------
    Model2DKernel : Create kernel from astropy.models.Parametric2DModel
    Model1DKernel : Create kernel from astropy.models.Parametric1DModel

    Examples
    --------
    Define one dimensional mask:

        >>> from astropy.nddata.convolution.kernels import CustomKernel
        >>> import numy as np
        >>> mask = np.array([1, 2, 3, 2, 1])
        >>> kernel = CustomKernel(mask)
        >>> kernel.dimension
        1

    Define two dimensional mask:

        >>> mask = np.array([[1, 1, 1], [1, 2, 1], [1, 1, 1]])
        >>> kernel = CustomKernel(mask)
        >>> kernel.dimension
        2
    """

    def __init__(self, mask):
        # Pass 'None' because mask is overridden in the next line
        self.mask = mask
        super(CustomKernel, self).__init__(self._mask)

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
