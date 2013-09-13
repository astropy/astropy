from __future__ import division

import numpy as np

from .core import Kernel1D, Kernel2D, Kernel
from .utils import KernelSizeError
from ...modeling import models
from ...modeling.core import Parametric1DModel, Parametric2DModel

__all__ = sorted(['Gaussian1DKernel', 'Gaussian2DKernel', 'CustomKernel',
                  'Box1DKernel', 'Box2DKernel', 'Tophat2DKernel',
                  'Trapezoid1DKernel', 'MexicanHat1DKernel',
                  'MexicanHat2DKernel', 'AiryDisk2DKernel',
                  'Model1DKernel', 'Model2DKernel',
                  'TrapezoidDisk2DKernel', 'Ring2DKernel'])


class Gaussian1DKernel(Kernel1D):
    """
    1D Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifacts.

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


    See Also
    --------
    Box1DKernel, Trapezoid1DKernel, MexicanHat1DKernel


    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import Gaussian1DKernel
        gauss_1D_kernel = Gaussian1DKernel(10)
        plt.plot(gauss_1D_kernel, drawstyle='steps')
        plt.xlabel('x [pixels]')
        plt.ylabel('value')
        plt.show()
    """
    _separable = True
    _is_bool = False

    def __init__(self, width, **kwargs):
        self._model = models.Gaussian1DModel(1. / (np.sqrt(2 * np.pi) * width), 0, width)
        self._default_size = 8 * width
        super(Gaussian1DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(1. - 1 / self._normalization)


class Gaussian2DKernel(Kernel2D):
    """
    2D Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifacts.

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


    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import Gaussian2DKernel
        gaussian_2D_kernel = Gaussian2DKernel(10)
        plt.imshow(gaussian_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()

    """
    _separable = True
    _is_bool = False

    def __init__(self, width, **kwargs):
        self._model = models.Gaussian2DModel(1. / (2 * np.pi * width ** 2), 0, 0, width, width)
        self._default_size = 8 * width
        super(Gaussian2DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(1. - 1 / self._normalization)


class Box1DKernel(Kernel1D):
    """
    1D Box filter kernel.

    The Box filter or running mean is a smoothing filter. It is not isotropic
    and can produce artifacts, when applied repeatedly to the same data.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
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

    See Also
    --------
    Gaussian1DKernel, Trapezoid1DKernel, MexicanHat1DKernel


    Examples
    --------
    Kernel response function:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import Box1DKernel
        box_1D_kernel = Box1DKernel(9)
        plt.plot(box_1D_kernel, drawstyle='steps')
        plt.xlim(-1, 9)
        plt.xlabel('x [pixels]')
        plt.ylabel('value')
        plt.show()

    """
    _separable = True
    _is_bool = True

    def __init__(self, width, **kwargs):
        self._model = models.Box1DModel(1. / width, 0, width)
        self._default_size = width
        super(Box1DKernel, self).__init__(**kwargs)
        self._truncation = 0
        self.normalize()


class Box2DKernel(Kernel2D):
    """
    2D Box filter kernel.

    The Box filter or running mean is a smoothing filter. It is not isotropic
    and can produce artifact, when applied repeatedly to the same data.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
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


    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import Box2DKernel
        box_2D_kernel = Box2DKernel(9)
        plt.imshow(box_2D_kernel, interpolation='none', origin='lower')
        plt.xlim(-1, 9)
        plt.ylim(-1, 9)
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """
    _separable = True
    _is_bool = True

    def __init__(self, width, **kwargs):
        self._model = models.Box2DModel(1. / width ** 2, 0, 0, width, width)
        self._default_size = width
        super(Box2DKernel, self).__init__(**kwargs)
        self._truncation = 0
        self.normalize()


class Tophat2DKernel(Kernel2D):
    """
    2D Tophat filter kernel.

    The Tophat filter is an isotropic smoothing filter. It can produce artifact,
    when applied repeatedly on the same data.

    Parameters
    ----------
    radius : int
        Radius of the filter kernel.
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


    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel, 
    TrapezoidDisk2DKernel, AiryDisk2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import Tophat2DKernel
        tophat_2D_kernel = Tophat2DKernel(40)
        plt.imshow(tophat_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()

    """
    def __init__(self, radius, **kwargs):
        self._model = models.Disk2DModel(1. / (np.pi * radius ** 2), 0, 0, radius)
        self._default_size = 2 * radius
        super(Tophat2DKernel, self).__init__(**kwargs)
        self._truncation = 0


class Ring2DKernel(Kernel2D):
    """
    2D Ring filter kernel.

    The Ring filter kernel is the difference between two Tophat kernels of
    different width. This kernel is useful for, e.g., background estimation.

    Parameters
    ----------
    radius_in : number
        Inner radius of the ring kernel.
    radius_out : number
        Outer radius of the ring kernel.
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

    See Also
    --------
    Box2DKernel, Gaussian2DKernel, MexicanHat2DKernel, Tophat2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import Ring2DKernel
        ring_2D_kernel = Ring2DKernel(9, 17)
        plt.imshow(ring_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """
    def __init__(self, radius_in, radius_out, **kwargs):
        self._model = models.Ring2DModel(1. / (np.pi * (radius_out ** 2 - radius_in ** 2)), 
                                        0, 0, radius_in, radius_out)
        self._default_size = 2 * radius_out
        super(Ring2DKernel, self).__init__(**kwargs)
        self._truncation = 0


class Trapezoid1DKernel(Kernel1D):
    """
    1D trapezoid kernel.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    slope : number
        Slope of the filter kernel's tails
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

    See Also
    --------
    Box1DKernel, Gaussian1DKernel, MexicanHat1DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import Trapezoid1DKernel
        trapezoid_1D_kernel = Trapezoid1DKernel(17, slope=0.2)
        plt.plot(trapezoid_1D_kernel, drawstyle='steps')
        plt.xlabel('x [pixels]')
        plt.ylabel('amplitude')
        plt.xlim(-1, 28)
        plt.show()
    """
    _is_bool = False

    def __init__(self, width, slope=1., **kwargs):
        self._model = models.Trapezoid1DModel(1, 0, width, slope)
        self._default_size = 2 * (width / 2 + 1. / slope) + 1 
        super(Trapezoid1DKernel, self).__init__(**kwargs)
        self._truncation = 0
        self.normalize()


class TrapezoidDisk2DKernel(Kernel2D):
    """
    2D trapezoid kernel.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    slope : number
        Slope of the filter kernel's tails
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

    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import TrapezoidDisk2DKernel
        trapezoid_2D_kernel = TrapezoidDisk2DKernel(20, slope=0.2)
        plt.imshow(trapezoid_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()

    """
    _is_bool = False

    def __init__(self, radius, slope=1., **kwargs):
        self._model = models.TrapezoidDisk2DModel(1, 0, 0, radius, slope)
        self._default_size = 2 * (radius + 1. / slope) - 1
        super(TrapezoidDisk2DKernel, self).__init__(**kwargs)
        self._truncation = 0
        self.normalize()


class MexicanHat1DKernel(Kernel1D):
    """
    1D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a bandpass filter. It
    smoothes the data and removes slowly varying or constant structures
    (e.g. Background). It is useful for peak or multi-scale detection.
    This kernel is derived from a normalized Gaussian.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    x_size : odd int, optional
        Size in x direction of the kernel array. Default = 8 * width.
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


    See Also
    --------
    Box1DKernel, Gaussian1DKernel, Trapezoid1DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import MexicanHat1DKernel
        mexicanhat_1D_kernel = MexicanHat1DKernel(10)
        plt.plot(mexicanhat_1D_kernel, drawstyle='steps')
        plt.xlabel('x [pixels]')
        plt.ylabel('value')
        plt.show()

    """
    _is_bool = True

    def __init__(self, width, **kwargs):
        self._default_size = 8 * width
        self._model = models.MexicanHat1DModel(-1. / (np.sqrt(2 * np.pi) * width ** 3), 0, width)
        super(MexicanHat1DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(self._array.sum() / self._array.size)
        self._normalization = 0


class MexicanHat2DKernel(Kernel2D):
    """
    2D Mexican hat filter kernel.

    The Mexican Hat or Gaussian-Laplace filter is a bandpass filter. It
    smoothes the data and removes slowly varying or constant structures
    (e.g. Background). It is useful for peak or multi-scale detection.
    This kernel is derived from a normalized Gaussian.

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


    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import MexicanHat2DKernel
        mexicanhat_2D_kernel = MexicanHat2DKernel(10)
        plt.imshow(mexicanhat_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """
    _is_bool = False

    def __init__(self, width, **kwargs):
        self._default_size = 8 * width
        self._model = models.MexicanHat2DModel(-1. / (np.pi * width ** 4), 0, 0, width)
        super(MexicanHat2DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(self._array.sum() / self._array.size)
        self._normalization = 0


class AiryDisk2DKernel(Kernel2D):
    """
    2D Airy disk kernel.

    This kernel models the diffraction pattern of a circular aperture. This 
    kernel is normalized to a peak value of 1.

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

    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.nddata.convolution import AiryDisk2DKernel
        airydisk_2D_kernel = AiryDisk2DKernel(10)
        plt.imshow(airydisk_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """
    _is_bool = False

    def __init__(self, width, **kwargs):
        self._default_size = 8 * width
        self._model = models.AiryDisk2DModel(1, 0, 0, width)
        super(AiryDisk2DKernel, self).__init__(**kwargs)
        self.normalize()
        self._truncation = None


class Model1DKernel(Kernel1D):
    """
    Create kernel from astropy.models.Parametric1DModel.

    The model has to be centered on x = 0.

    Parameters
    ----------
    model : Parametric1DModel
        Kernel response function model
    x_size : odd int, optional
        Size in x direction of the kernel array. Default = 8 * width.
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

    Raises
    ------
    TypeError
        If model is not an instance of astropy.models.Parametric1DModel

    See also
    --------
    Model2DKernel : Create kernel from astropy.models.Parametric2DModel
    CustomKernel : Create kernel from list or array

    Examples
    --------
    Define a Gaussian1D model:

        >>> from astropy.modeling.models import Gaussian1DModel
        >>> from astropy.nddata.convolution.kernels import Model1DKernel
        >>> gauss = Gaussian1DModel(1, 0, 2)

    And create a custom one dimensional kernel from it:

        >>> gauss_kernel = Model1DKernel(gauss, x_size=9)

    This kernel can now be used like a usual astropy kernel.
    """
    _separable = False
    _is_bool = False

    def __init__(self, model, **kwargs):
        if isinstance(model, Parametric1DModel):
            self._model = model
        else:
            raise TypeError("Must be Parametric1DModel")
        super(Model1DKernel, self).__init__(**kwargs)


class Model2DKernel(Kernel2D):
    """
    Create kernel from astropy.models.Parametric2DModel.

    The model has to be centered on x = 0 and y= 0.

    Parameters
    ----------
    model : Parametric2DModel
        Kernel response function model
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

    Raises
    ------
    TypeError
        If model is not an instance of astropy.models.Parametric2DModel

    See also
    --------
    Model1DKernel : Create kernel from astropy.models.Parametric1DModel
    CustomKernel : Create kernel from list or array

    Examples
    --------
    Define a Gaussian2D model:

        >>> from astropy.modeling.models import Gaussian2DModel
        >>> from astropy.nddata.convolution.kernels import Model2DKernel
        >>> gauss = Gaussian2DModel(1, 0, 0, 2, 2)

    And create a custom two dimensional kernel from it:

        >>> gauss_kernel = Model2DKernel(gauss, x_size=9)

    This kernel can now be used like a usual astropy kernel.

    """
    _is_bool = False
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
        odd = np.all([axes_size % 2 != 0 for axes_size in self.shape])

        if not odd:
            raise KernelSizeError("Kernel size must be odd in all axes.")

        # Check if array is bool
        ones = self._array == 1.
        zeros = self._array == 0
        self._is_bool = np.all(np.logical_or(ones, zeros))

        # Set normalization
        self._normalization = 1. / self._array.sum()
