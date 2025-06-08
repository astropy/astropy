# Licensed under a 3-clause BSD style license - see LICENSE.rst

import math

import numpy as np

from astropy.modeling import models
from astropy.modeling.core import Fittable1DModel, Fittable2DModel

from .core import Kernel, Kernel1D, Kernel2D
from .utils import KernelSizeError, has_even_axis

__all__ = [
    "AiryDisk2DKernel",
    "Box1DKernel",
    "Box2DKernel",
    "CustomKernel",
    "Gaussian1DKernel",
    "Gaussian2DKernel",
    "Model1DKernel",
    "Model2DKernel",
    "Moffat2DKernel",
    "RickerWavelet1DKernel",
    "RickerWavelet2DKernel",
    "Ring2DKernel",
    "Tophat2DKernel",
    "Trapezoid1DKernel",
    "TrapezoidDisk2DKernel",
]


def _round_up_to_odd_integer(value):
    i = math.ceil(value)
    return i + 1 if i % 2 == 0 else i


class Gaussian1DKernel(Kernel1D):
    """
    1D Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifacts.

    The generated kernel is normalized so that it integrates to 1.

    Parameters
    ----------
    stddev : number
        Standard deviation of the Gaussian kernel.
    x_size : int, optional
        Size of the kernel array. Default = ⌊8*stddev+1⌋.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
                model over the bin. Very slow.
    factor : number, optional
        Factor of oversampling. Default factor = 10. If the factor
        is too large, evaluation can be very slow.

    See Also
    --------
    Box1DKernel, Trapezoid1DKernel, RickerWavelet1DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import Gaussian1DKernel
        gauss_1D_kernel = Gaussian1DKernel(10)
        plt.plot(gauss_1D_kernel, drawstyle='steps')
        plt.xlabel('x [pixels]')
        plt.ylabel('value')
        plt.show()
    """

    _separable = True
    _is_bool = False

    def __init__(self, stddev, **kwargs):
        self._model = models.Gaussian1D(1.0 / (np.sqrt(2 * np.pi) * stddev), 0, stddev)
        self._default_size = _round_up_to_odd_integer(8 * stddev)
        super().__init__(**kwargs)
        self.normalize()


class Gaussian2DKernel(Kernel2D):
    """
    2D Gaussian filter kernel.

    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifacts.

    The generated kernel is normalized so that it integrates to 1.

    Parameters
    ----------
    x_stddev : float
        Standard deviation of the Gaussian in x before rotating by theta.
    y_stddev : float
        Standard deviation of the Gaussian in y before rotating by theta.
    theta : float or `~astropy.units.Quantity` ['angle']
        Rotation angle. If passed as a float, it is assumed to be in radians.
        The rotation angle increases counterclockwise.
    x_size : int, optional
        Size in x direction of the kernel array. Default = ⌊8*stddev + 1⌋.
    y_size : int, optional
        Size in y direction of the kernel array. Default = ⌊8*stddev + 1⌋.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
    Box2DKernel, Tophat2DKernel, RickerWavelet2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel, Moffat2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import Gaussian2DKernel
        gaussian_2D_kernel = Gaussian2DKernel(10)
        plt.imshow(gaussian_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """

    _separable = True
    _is_bool = False

    def __init__(self, x_stddev, y_stddev=None, theta=0.0, **kwargs):
        if y_stddev is None:
            y_stddev = x_stddev
        self._model = models.Gaussian2D(
            amplitude=1.0 / (2 * np.pi * x_stddev * y_stddev),
            x_mean=0,
            y_mean=0,
            x_stddev=x_stddev,
            y_stddev=y_stddev,
            theta=theta,
        )
        self._default_size = _round_up_to_odd_integer(8 * np.max([x_stddev, y_stddev]))
        super().__init__(**kwargs)
        self.normalize()


class Box1DKernel(Kernel1D):
    """
    1D Box filter kernel.

    The Box filter or running mean is a smoothing filter. It is not isotropic
    and can produce artifacts when applied repeatedly to the same data.

    The generated kernel is normalized so that it integrates to 1.

    By default the Box kernel uses the ``linear_interp`` discretization mode,
    which allows non-shifting, even-sized kernels.  This is achieved by
    weighting the edge pixels with 1/2. E.g a Box kernel with an effective
    smoothing of 4 pixel would have the following array: [0.5, 1, 1, 1, 0.5].

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    mode : {'linear_interp', 'center', 'oversample', 'integrate'}, optional
        One of the following discretization modes:
            * 'linear_interp' (default)
                Discretize model by linearly interpolating
                between the values at the corners of the bin.
            * 'center'
                Discretize model by taking the value
                at the center of the bin.
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
    Gaussian1DKernel, Trapezoid1DKernel, RickerWavelet1DKernel

    Examples
    --------
    Kernel response function:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import Box1DKernel
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
        self._model = models.Box1D(1.0 / width, 0, width)
        self._default_size = _round_up_to_odd_integer(width)
        kwargs["mode"] = "linear_interp"
        super().__init__(**kwargs)
        self.normalize()


class Box2DKernel(Kernel2D):
    """
    2D Box filter kernel.

    The Box filter or running mean is a smoothing filter. It is not isotropic
    and can produce artifacts when applied repeatedly to the same data.

    The generated kernel is normalized so that it integrates to 1.

    By default the Box kernel uses the ``linear_interp`` discretization mode,
    which allows non-shifting, even-sized kernels.  This is achieved by
    weighting the edge pixels with 1/2.

    Parameters
    ----------
    width : number
        Width of the filter kernel.
    mode : {'linear_interp', 'center', 'oversample', 'integrate'}, optional
        One of the following discretization modes:
            * 'linear_interp' (default)
                Discretize model by performing a bilinear interpolation
                between the values at the corners of the bin.
            * 'center'
                Discretize model by taking the value
                at the center of the bin.
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
    Gaussian2DKernel, Tophat2DKernel, RickerWavelet2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel, Moffat2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import Box2DKernel
        box_2D_kernel = Box2DKernel(9)
        plt.imshow(box_2D_kernel, interpolation='none', origin='lower',
                   vmin=0.0, vmax=0.015)
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
        self._model = models.Box2D(1.0 / width**2, 0, 0, width, width)
        self._default_size = _round_up_to_odd_integer(width)
        kwargs["mode"] = "linear_interp"
        super().__init__(**kwargs)
        self.normalize()


class Tophat2DKernel(Kernel2D):
    """
    2D Tophat filter kernel.

    The Tophat filter is an isotropic smoothing filter. It can produce
    artifacts when applied repeatedly on the same data.

    The generated kernel is normalized so that it integrates to 1.

    Parameters
    ----------
    radius : int
        Radius of the filter kernel.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
    Gaussian2DKernel, Box2DKernel, RickerWavelet2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel, Moffat2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import Tophat2DKernel
        tophat_2D_kernel = Tophat2DKernel(40)
        plt.imshow(tophat_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """

    def __init__(self, radius, **kwargs):
        self._model = models.Disk2D(1.0 / (np.pi * radius**2), 0, 0, radius)
        self._default_size = _round_up_to_odd_integer(2 * radius)
        super().__init__(**kwargs)
        self.normalize()


class Ring2DKernel(Kernel2D):
    """
    2D Ring filter kernel.

    The Ring filter kernel is the difference between two Tophat kernels of
    different width. This kernel is useful for, e.g., background estimation.

    The generated kernel is normalized so that it integrates to 1.

    Parameters
    ----------
    radius_in : number
        Inner radius of the ring kernel.
    width : number
        Width of the ring kernel.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
    Gaussian2DKernel, Box2DKernel, Tophat2DKernel, RickerWavelet2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel, Moffat2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import Ring2DKernel
        ring_2D_kernel = Ring2DKernel(9, 8)
        plt.imshow(ring_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """

    def __init__(self, radius_in, width, **kwargs):
        radius_out = radius_in + width
        self._model = models.Ring2D(
            1.0 / (np.pi * (radius_out**2 - radius_in**2)), 0, 0, radius_in, width
        )
        self._default_size = _round_up_to_odd_integer(2 * radius_out)
        super().__init__(**kwargs)
        self.normalize()


class Trapezoid1DKernel(Kernel1D):
    """
    1D trapezoid kernel.

    The generated kernel is normalized so that it integrates to 1.

    Parameters
    ----------
    width : number
        Width of the filter kernel, defined as the width of the constant part,
        before it begins to slope down.
    slope : number
        Slope of the filter kernel's tails
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
    Box1DKernel, Gaussian1DKernel, RickerWavelet1DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import Trapezoid1DKernel
        trapezoid_1D_kernel = Trapezoid1DKernel(17, slope=0.2)
        plt.plot(trapezoid_1D_kernel, drawstyle='steps')
        plt.xlabel('x [pixels]')
        plt.ylabel('amplitude')
        plt.xlim(-1, 28)
        plt.show()
    """

    _is_bool = False

    def __init__(self, width, slope=1.0, **kwargs):
        self._model = models.Trapezoid1D(1, 0, width, slope)
        self._default_size = _round_up_to_odd_integer(width + 2.0 / slope)
        super().__init__(**kwargs)
        self.normalize()


class TrapezoidDisk2DKernel(Kernel2D):
    """
    2D trapezoid kernel.

    The generated kernel is normalized so that it integrates to 1.

    Parameters
    ----------
    radius : number
        Width of the filter kernel, defined as the width of the constant part,
        before it begins to slope down.
    slope : number
        Slope of the filter kernel's tails
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
    Gaussian2DKernel, Box2DKernel, Tophat2DKernel, RickerWavelet2DKernel,
    Ring2DKernel, AiryDisk2DKernel, Moffat2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import TrapezoidDisk2DKernel
        trapezoid_2D_kernel = TrapezoidDisk2DKernel(20, slope=0.2)
        plt.imshow(trapezoid_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """

    _is_bool = False

    def __init__(self, radius, slope=1.0, **kwargs):
        self._model = models.TrapezoidDisk2D(1, 0, 0, radius, slope)
        self._default_size = _round_up_to_odd_integer(2 * radius + 2.0 / slope)
        super().__init__(**kwargs)
        self.normalize()


class RickerWavelet1DKernel(Kernel1D):
    """
    1D Ricker wavelet filter kernel (sometimes known as a "Mexican Hat"
    kernel).

    The Ricker wavelet, or inverted Gaussian-Laplace filter, is a
    bandpass filter. It smooths the data and removes slowly varying
    or constant structures (e.g. Background). It is useful for peak or
    multi-scale detection.

    This kernel is derived from a normalized Gaussian function, by
    computing the second derivative. This results in an amplitude
    at the kernels center of 1. / (sqrt(2 * pi) * width ** 3). The
    normalization is the same as for `scipy.ndimage.gaussian_laplace`,
    except for a minus sign.

    .. note::

        See https://github.com/astropy/astropy/pull/9445 for discussions
        related to renaming of this kernel.

    Parameters
    ----------
    width : number
        Width of the filter kernel, defined as the standard deviation
        of the Gaussian function from which it is derived.
    x_size : int, optional
        Size in x direction of the kernel array. Default = ⌊8*width +1⌋.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
        from astropy.convolution import RickerWavelet1DKernel
        ricker_1d_kernel = RickerWavelet1DKernel(10)
        plt.plot(ricker_1d_kernel, drawstyle='steps')
        plt.xlabel('x [pixels]')
        plt.ylabel('value')
        plt.show()
    """

    _is_bool = True

    def __init__(self, width, **kwargs):
        amplitude = 1.0 / (np.sqrt(2 * np.pi) * width**3)
        self._model = models.RickerWavelet1D(amplitude, 0, width)
        self._default_size = _round_up_to_odd_integer(8 * width)
        super().__init__(**kwargs)


class RickerWavelet2DKernel(Kernel2D):
    """
    2D Ricker wavelet filter kernel (sometimes known as a "Mexican Hat"
    kernel).

    The Ricker wavelet, or inverted Gaussian-Laplace filter, is a
    bandpass filter. It smooths the data and removes slowly varying
    or constant structures (e.g. Background). It is useful for peak or
    multi-scale detection.

    This kernel is derived from a normalized Gaussian function, by
    computing the second derivative. This results in an amplitude
    at the kernels center of 1. / (pi * width ** 4). The normalization
    is the same as for `scipy.ndimage.gaussian_laplace`, except
    for a minus sign.

    .. note::

        See https://github.com/astropy/astropy/pull/9445 for discussions
        related to renaming of this kernel.

    Parameters
    ----------
    width : number
        Width of the filter kernel, defined as the standard deviation
        of the Gaussian function from which it is derived.
    x_size : int, optional
        Size in x direction of the kernel array. Default = ⌊8*width +1⌋.
    y_size : int, optional
        Size in y direction of the kernel array. Default = ⌊8*width +1⌋.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
    Gaussian2DKernel, Box2DKernel, Tophat2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel, Moffat2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import RickerWavelet2DKernel
        ricker_2d_kernel = RickerWavelet2DKernel(10)
        plt.imshow(ricker_2d_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """

    _is_bool = False

    def __init__(self, width, **kwargs):
        amplitude = 1.0 / (np.pi * width**4)
        self._model = models.RickerWavelet2D(amplitude, 0, 0, width)
        self._default_size = _round_up_to_odd_integer(8 * width)
        super().__init__(**kwargs)


class AiryDisk2DKernel(Kernel2D):
    """
    2D Airy disk kernel.

    This kernel models the diffraction pattern of a circular aperture.

    The generated kernel is normalized so that it integrates to 1.

    Parameters
    ----------
    radius : float
        The radius of the Airy disk kernel (radius of the first zero).
    x_size : int, optional
        Size in x direction of the kernel array. Default = ⌊8*radius + 1⌋.
    y_size : int, optional
        Size in y direction of the kernel array. Default = ⌊8*radius + 1⌋.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
    Gaussian2DKernel, Box2DKernel, Tophat2DKernel, RickerWavelet2DKernel,
    Ring2DKernel, TrapezoidDisk2DKernel, Moffat2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import AiryDisk2DKernel
        airydisk_2D_kernel = AiryDisk2DKernel(10)
        plt.imshow(airydisk_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """

    _is_bool = False

    def __init__(self, radius, **kwargs):
        self._model = models.AiryDisk2D(1, 0, 0, radius)
        self._default_size = _round_up_to_odd_integer(8 * radius)
        super().__init__(**kwargs)
        self.normalize()


class Moffat2DKernel(Kernel2D):
    """
    2D Moffat kernel.

    This kernel is a typical model for a seeing limited PSF.

    The generated kernel is normalized so that it integrates to 1.

    Parameters
    ----------
    gamma : float
        Core width of the Moffat model.
    alpha : float
        Power index of the Moffat model.
    x_size : int, optional
        Size in x direction of the kernel array. Default = ⌊8*radius + 1⌋.
    y_size : int, optional
        Size in y direction of the kernel array. Default = ⌊8*radius + 1⌋.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
    Gaussian2DKernel, Box2DKernel, Tophat2DKernel, RickerWavelet2DKernel,
    Ring2DKernel, TrapezoidDisk2DKernel, AiryDisk2DKernel

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.convolution import Moffat2DKernel
        moffat_2D_kernel = Moffat2DKernel(3, 2)
        plt.imshow(moffat_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """

    _is_bool = False

    def __init__(self, gamma, alpha, **kwargs):
        # Compute amplitude, from
        # https://en.wikipedia.org/wiki/Moffat_distribution
        amplitude = (alpha - 1.0) / (np.pi * gamma * gamma)
        self._model = models.Moffat2D(amplitude, 0, 0, gamma, alpha)
        self._default_size = _round_up_to_odd_integer(4.0 * self._model.fwhm)
        super().__init__(**kwargs)
        self.normalize()


class Model1DKernel(Kernel1D):
    """
    Create kernel from 1D model.

    The model has to be centered on x = 0.

    Parameters
    ----------
    model : `~astropy.modeling.Fittable1DModel`
        Kernel response function model
    x_size : int, optional
        Size in x direction of the kernel array. Default = ⌊8*width +1⌋.
        Must be odd.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
        If model is not an instance of `~astropy.modeling.Fittable1DModel`

    See Also
    --------
    Model2DKernel : Create kernel from `~astropy.modeling.Fittable2DModel`
    CustomKernel : Create kernel from list or array

    Examples
    --------
    Define a Gaussian1D model:

        >>> from astropy.modeling.models import Gaussian1D
        >>> from astropy.convolution.kernels import Model1DKernel
        >>> gauss = Gaussian1D(1, 0, 2)

    And create a custom one dimensional kernel from it:

        >>> gauss_kernel = Model1DKernel(gauss, x_size=9)

    This kernel can now be used like a usual Astropy kernel.
    """

    _separable = False
    _is_bool = False

    def __init__(self, model, **kwargs):
        if isinstance(model, Fittable1DModel):
            self._model = model
        else:
            raise TypeError("Must be Fittable1DModel")
        super().__init__(**kwargs)


class Model2DKernel(Kernel2D):
    """
    Create kernel from 2D model.

    The model has to be centered on x = 0 and y = 0.

    Parameters
    ----------
    model : `~astropy.modeling.Fittable2DModel`
        Kernel response function model
    x_size : int, optional
        Size in x direction of the kernel array. Default = ⌊8*width +1⌋.
        Must be odd.
    y_size : int, optional
        Size in y direction of the kernel array. Default = ⌊8*width +1⌋.
    mode : {'center', 'linear_interp', 'oversample', 'integrate'}, optional
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
        If model is not an instance of `~astropy.modeling.Fittable2DModel`

    See Also
    --------
    Model1DKernel : Create kernel from `~astropy.modeling.Fittable1DModel`
    CustomKernel : Create kernel from list or array

    Examples
    --------
    Define a Gaussian2D model:

        >>> from astropy.modeling.models import Gaussian2D
        >>> from astropy.convolution.kernels import Model2DKernel
        >>> gauss = Gaussian2D(1, 0, 0, 2, 2)

    And create a custom two dimensional kernel from it:

        >>> gauss_kernel = Model2DKernel(gauss, x_size=9)

    This kernel can now be used like a usual astropy kernel.
    """

    _is_bool = False
    _separable = False

    def __init__(self, model, **kwargs):
        self._separable = False
        if isinstance(model, Fittable2DModel):
            self._model = model
        else:
            raise TypeError("Must be Fittable2DModel")
        super().__init__(**kwargs)


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
    `~astropy.convolution.KernelSizeError`
        If array size is even.

    See Also
    --------
    Model2DKernel, Model1DKernel

    Examples
    --------
    Define one dimensional array:

        >>> from astropy.convolution.kernels import CustomKernel
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
        super().__init__(self._array)

    @property
    def array(self):
        """
        Filter kernel array.
        """
        return self._array

    @array.setter
    def array(self, array):
        """
        Filter kernel array setter.
        """
        if isinstance(array, np.ndarray):
            self._array = array.astype(np.float64)
        elif isinstance(array, list):
            self._array = np.array(array, dtype=np.float64)
        else:
            raise TypeError("Must be list or array.")

        # Check if array is odd in all axes
        if has_even_axis(self):
            raise KernelSizeError("Kernel size must be odd in all axes.")

        # Check if array is bool
        ones = self._array == 1.0
        zeros = self._array == 0
        self._is_bool = bool(np.all(np.logical_or(ones, zeros)))
