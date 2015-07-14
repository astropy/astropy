# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Mathematical models."""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import inspect

import numpy as np

from .core import (Fittable1DModel, Fittable2DModel, Model,
                   ModelDefinitionError, custom_model)
from .parameters import Parameter, InputParameterError
from ..utils import deprecated
from ..extern import six

__all__ = sorted([
    'AiryDisk2D', 'Moffat1D', 'Moffat2D', 'Box1D',
    'Box2D', 'Const1D', 'Const2D', 'Ellipse2D', 'Disk2D',
    'Gaussian1D', 'GaussianAbsorption1D', 'Gaussian2D', 'Linear1D',
    'Lorentz1D', 'MexicanHat1D', 'MexicanHat2D', 'Scale', 'Redshift', 'Shift',
    'Sine1D', 'Trapezoid1D', 'TrapezoidDisk2D', 'Ring2D',
    'custom_model_1d'
])


class Gaussian1D(Fittable1DModel):
    """
    One dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the Gaussian.
    mean : float
        Mean of the Gaussian.
    stddev : float
        Standard deviation of the Gaussian.

    Notes
    -----

    Model formula:

        .. math:: f(x) = A e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}

    Examples
    --------
    >>> from astropy.modeling import models
    >>> def tie_center(model):
    ...         mean = 50 * model.stddev
    ...         return mean
    >>> tied_parameters = {'mean': tie_center}

    Specify that 'mean' is a tied parameter in one of two ways:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                             tied=tied_parameters)

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.mean.tied
    False
    >>> g1.mean.tied = tie_center
    >>> g1.mean.tied
    <function tie_center at 0x...>

    Fixed parameters:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                        fixed={'stddev': True})
    >>> g1.stddev.fixed
    True

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.stddev.fixed
    False
    >>> g1.stddev.fixed = True
    >>> g1.stddev.fixed
    True

    See Also
    --------
    Gaussian2D, Box1D, Moffat1D, Lorentz1D
    """

    amplitude = Parameter(default=1)
    mean = Parameter(default=0)
    stddev = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, mean, stddev):
        """
        Gaussian1D model function.
        """
        return amplitude * np.exp(- 0.5 * (x - mean) ** 2 / stddev ** 2)

    @staticmethod
    def fit_deriv(x, amplitude, mean, stddev):
        """
        Gaussian1D model function derivatives.
        """

        d_amplitude = np.exp(-0.5 / stddev ** 2 * (x - mean) ** 2)
        d_mean = amplitude * d_amplitude * (x - mean) / stddev ** 2
        d_stddev = amplitude * d_amplitude * (x - mean) ** 2 / stddev ** 3
        return [d_amplitude, d_mean, d_stddev]


class GaussianAbsorption1D(Fittable1DModel):
    """
    One dimensional Gaussian absorption line model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the gaussian absorption.
    mean : float
        Mean of the gaussian.
    stddev : float
        Standard deviation of the gaussian.

    Notes
    -----

    Model formula:

        .. math:: f(x) = 1 - A e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}

    See Also
    --------
    Gaussian1D
    """

    amplitude = Parameter(default=1)
    mean = Parameter(default=0)
    stddev = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, mean, stddev):
        """
        GaussianAbsorption1D model function.
        """
        return 1.0 - Gaussian1D.evaluate(x, amplitude, mean, stddev)

    @staticmethod
    def fit_deriv(x, amplitude, mean, stddev):
        """
        GaussianAbsorption1D model function derivatives.
        """
        import operator
        return list(six.moves.map(
            operator.neg, Gaussian1D.fit_deriv(x, amplitude, mean, stddev)))


class Gaussian2D(Fittable2DModel):
    """
    Two dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the Gaussian.
    x_mean : float
        Mean of the Gaussian in x.
    y_mean : float
        Mean of the Gaussian in y.
    x_stddev : float
        Standard deviation of the Gaussian in x before rotating by theta.
        ``x_stddev`` and ``y_stddev`` must be specified unless a covariance
        matrix (``cov_matrix``) is input.
    y_stddev : float
        Standard deviation of the Gaussian in y before rotating by theta.
        ``x_stddev`` and ``y_stddev`` must be specified unless a covariance
        matrix (``cov_matrix``) is input.
    theta : float, optional
        Rotation angle in radians. The rotation angle increases
        counterclockwise, from the positive x-axis.
    cov_matrix : ndarray, optional
        A 2x2 covariance matrix. If specified, overrides the ``x_stddev``,
        ``y_stddev``, and ``theta`` specification.

    Notes
    -----
    Model formula:

        .. math::

            f(x, y) = A e^{-a\\left(x - x_{0}\\right)^{2}  -b\\left(x - x_{0}\\right)
            \\left(y - y_{0}\\right)  -c\\left(y - y_{0}\\right)^{2}}

    Using the following definitions:

        .. math::
            a = \\left(\\frac{\\cos^{2}{\\left (\\theta \\right )}}{2 \\sigma_{x}^{2}} +
            \\frac{\\sin^{2}{\\left (\\theta \\right )}}{2 \\sigma_{y}^{2}}\\right)

            b = \\left(\\frac{\\sin{\\left (2 \\theta \\right )}}{2 \\sigma_{x}^{2}} -
            \\frac{\\sin{\\left (2 \\theta \\right )}}{2 \\sigma_{y}^{2}}\\right)

            c = \\left(\\frac{\\sin^{2}{\\left (\\theta \\right )}}{2 \\sigma_{x}^{2}} +
            \\frac{\\cos^{2}{\\left (\\theta \\right )}}{2 \\sigma_{y}^{2}}\\right)

    If using a ``cov_matrix``, the model is of the form:
        .. math::
            f(x, y) = A e^{-0.5 \\left(\\vec{x} - \\vec{x}_{0}\\right)^{T} \\Sigma^{-1} \\left(\\vec{x} - \\vec{x}_{0}\\right)}

    where :math:`\\vec{x} = [x, y]`, :math:`\\vec{x}_{0} = [x_{0}, y_{0}]`,
    and :math:`\\Sigma` is the covariance matrix:

        .. math::
            \\Sigma = \\left(\\begin{array}{ccc}
            \\sigma_x^2               & \\rho \\sigma_x \\sigma_y \\\\
            \\rho \\sigma_x \\sigma_y & \\sigma_y^2
            \end{array}\\right)

    :math:`\\rho` is the correlation between ``x`` and ``y``, which should
    be between -1 and +1.  Positive correlation corresponds to a
    ``theta`` in the range 0 to 90 degrees.  Negative correlation
    corresponds to a ``theta`` in the range of 0 to -90 degrees.

    See [1]_ for more details about the 2D Gaussian function.

    See Also
    --------
    Gaussian1D, Box2D, Moffat2D

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Gaussian_function
    """

    amplitude = Parameter(default=1)
    x_mean = Parameter(default=0)
    y_mean = Parameter(default=0)
    x_stddev = Parameter(default=1)
    y_stddev = Parameter(default=1)
    theta = Parameter(default=0)


    def __init__(self, amplitude=amplitude.default, x_mean=x_mean.default,
                 y_mean=y_mean.default, x_stddev=None, y_stddev=None,
                 theta=0.0, cov_matrix=None, **kwargs):
        if y_stddev is None and cov_matrix is None:
            raise InputParameterError(
                "Either x/y_stddev must be specified, or a "
                "covariance matrix.")
        elif x_stddev is None and cov_matrix is None:
            raise InputParameterError(
                "Either x/y_stddev must be specified, or a "
                "covariance matrix.")
        elif cov_matrix is not None and (x_stddev is not None or
                                         y_stddev is not None):
            raise InputParameterError(
                "Cannot specify both cov_matrix and x/y_stddev")

        # Compute principle coordinate system transformation
        elif cov_matrix is not None:
            if (x_stddev is not None or y_stddev is not None):
                raise InputParameterError(
                    "Cannot specify both cov_matrix and x/y_stddev")
            # Compute principle coordinate system transformation
            cov_matrix = np.array(cov_matrix)

            if cov_matrix.shape != (2, 2):
                # TODO: Maybe it should be possible for the covariance matrix
                # to be some (x, y, ..., z, 2, 2) array to be broadcast with
                # other parameters of shape (x, y, ..., z)
                # But that's maybe a special case to work out if/when needed
                raise ValueError("Covariance matrix must be 2x2")

            eig_vals, eig_vecs = np.linalg.eig(cov_matrix)
            x_stddev, y_stddev = np.sqrt(eig_vals)
            y_vec = eig_vecs[:, 0]
            theta = np.arctan2(y_vec[1], y_vec[0])

        super(Gaussian2D, self).__init__(
            amplitude=amplitude, x_mean=x_mean, y_mean=y_mean,
            x_stddev=x_stddev, y_stddev=y_stddev, theta=theta, **kwargs)

    @staticmethod
    def evaluate(x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta):
        """Two dimensional Gaussian function"""

        cost2 = np.cos(theta) ** 2
        sint2 = np.sin(theta) ** 2
        sin2t = np.sin(2. * theta)
        xstd2 = x_stddev ** 2
        ystd2 = y_stddev ** 2
        xdiff = x - x_mean
        ydiff = y - y_mean
        a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2))
        b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2))
        c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2))
        return amplitude * np.exp(-((a * xdiff ** 2) + (b * xdiff * ydiff) +
                                    (c * ydiff ** 2)))

    @staticmethod
    def fit_deriv(x, y, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta):
        """Two dimensional Gaussian function derivative with respect to parameters"""

        cost = np.cos(theta)
        sint = np.sin(theta)
        cost2 = np.cos(theta) ** 2
        sint2 = np.sin(theta) ** 2
        cos2t = np.cos(2. * theta)
        sin2t = np.sin(2. * theta)
        xstd2 = x_stddev ** 2
        ystd2 = y_stddev ** 2
        xstd3 = x_stddev ** 3
        ystd3 = y_stddev ** 3
        xdiff = x - x_mean
        ydiff = y - y_mean
        xdiff2 = xdiff ** 2
        ydiff2 = ydiff ** 2
        a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2))
        b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2))
        c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2))
        g = amplitude * np.exp(-((a * xdiff2) + (b * xdiff * ydiff) +
                                 (c * ydiff2)))
        da_dtheta = (sint * cost * ((1. / ystd2) - (1. / xstd2)))
        da_dx_stddev = -cost2 / xstd3
        da_dy_stddev = -sint2 / ystd3
        db_dtheta = (cos2t / xstd2) - (cos2t / ystd2)
        db_dx_stddev = -sin2t / xstd3
        db_dy_stddev = sin2t / ystd3
        dc_dtheta = -da_dtheta
        dc_dx_stddev = -sint2 / xstd3
        dc_dy_stddev = -cost2 / ystd3
        dg_dA = g / amplitude
        dg_dx_mean = g * ((2. * a * xdiff) + (b * ydiff))
        dg_dy_mean = g * ((b * xdiff) + (2. * c * ydiff))
        dg_dx_stddev = g * (-(da_dx_stddev * xdiff2 +
                              db_dx_stddev * xdiff * ydiff +
                              dc_dx_stddev * ydiff2))
        dg_dy_stddev = g * (-(da_dy_stddev * xdiff2 +
                              db_dy_stddev * xdiff * ydiff +
                              dc_dy_stddev * ydiff2))
        dg_dtheta = g * (-(da_dtheta * xdiff2 +
                           db_dtheta * xdiff * ydiff +
                           dc_dtheta * ydiff2))
        return [dg_dA, dg_dx_mean, dg_dy_mean, dg_dx_stddev, dg_dy_stddev,
                dg_dtheta]


class Shift(Model):
    """
    Shift a coordinate.

    Parameters
    ----------
    offset : float
        Offset to add to a coordinate.
    """

    inputs = ('x',)
    outputs = ('x',)

    offset = Parameter(default=0)

    @property
    def inverse(self):
        inv = self.copy()
        inv.offset *= -1
        return inv

    @staticmethod
    def evaluate(x, offset):
        return x + offset


class Scale(Model):
    """
    Multiply a model by a factor.

    Parameters
    ----------
    factor : float
        Factor by which to scale a coordinate.
    """

    inputs = ('x',)
    outputs = ('x',)

    factor = Parameter(default=1)
    linear = True

    @property
    def inverse(self):
        inv = self.copy()
        inv.factor = 1 / self.factor
        return inv

    @staticmethod
    def evaluate(x, factor):
        return factor * x


class Redshift(Fittable1DModel):
    """
    One dimensional redshift model.

    Parameters
    ----------
    z : float or a list of floats
        Redshift value(s).

    Notes
    -----
    Model formula:

        .. math:: \\lambda_{obs} = (1 + z) \\lambda_{rest}

    """

    z = Parameter(description='redshift', default=0)

    @staticmethod
    def evaluate(x, z):
        """One dimensional Redshift model function"""

        return (1 + z) * x

    @staticmethod
    def fit_deriv(x, z):
        """One dimensional Redshift model derivative"""
        d_z = x
        return [d_z]

    @property
    def inverse(self):
        """Inverse Redshift model"""

        inv = self.copy()
        inv.z = 1.0 / (1.0 + self.z) - 1.0
        return inv


class Sine1D(Fittable1DModel):
    """
    One dimensional Sine model.

    Parameters
    ----------
    amplitude : float
        Oscillation amplitude
    frequency : float
        Oscillation frequency

    See Also
    --------
    Const1D, Linear1D


    Notes
    -----
    Model formula:

        .. math:: f(x) = A \\sin(2 \\pi f x)
    """

    amplitude = Parameter(default=1)
    frequency = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, frequency):
        """One dimensional Sine model function"""

        return amplitude * np.sin(2 * np.pi * frequency * x)

    @staticmethod
    def fit_deriv(x, amplitude, frequency):
        """One dimensional Sine model derivative"""

        d_amplitude = np.sin(2 * np.pi * frequency * x)
        d_frequency = (2 * np.pi * x * amplitude *
                       np.cos(2 * np.pi * frequency * x))
        return [d_amplitude, d_frequency]


class Linear1D(Fittable1DModel):
    """
    One dimensional Line model.

    Parameters
    ----------
    slope : float
        Slope of the straight line

    intercept : float
        Intercept of the straight line

    See Also
    --------
    Const1D

    Notes
    -----
    Model formula:

        .. math:: f(x) = a x + b
    """

    slope = Parameter(default=1)
    intercept = Parameter(default=0)
    linear = True

    @staticmethod
    def evaluate(x, slope, intercept):
        """One dimensional Line model function"""

        return slope * x + intercept

    @staticmethod
    def fit_deriv(x, slope, intercept):
        """One dimensional Line model derivative with respect to parameters"""

        d_slope = x
        d_intercept = np.ones_like(x)
        return [d_slope, d_intercept]


class Lorentz1D(Fittable1DModel):
    """
    One dimensional Lorentzian model.

    Parameters
    ----------
    amplitude : float
        Peak value
    x_0 : float
        Position of the peak
    fwhm : float
        Full width at half maximum

    See Also
    --------
    Gaussian1D, Box1D, MexicanHat1D

    Notes
    -----
    Model formula:

    .. math::

        f(x) = \\frac{A \\gamma^{2}}{\\gamma^{2} + \\left(x - x_{0}\\right)^{2}}
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    fwhm = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, fwhm):
        """One dimensional Lorentzian model function"""

        return (amplitude * ((fwhm / 2.) ** 2) / ((x - x_0) ** 2 +
                (fwhm / 2.) ** 2))

    @staticmethod
    def fit_deriv(x, amplitude, x_0, fwhm):
        """One dimensional Lorentzian model derivative with respect to parameters"""

        d_amplitude = fwhm ** 2 / (fwhm ** 2 + (x - x_0) ** 2)
        d_x_0 = (amplitude * d_amplitude * (2 * x - 2 * x_0) /
                 (fwhm ** 2 + (x - x_0) ** 2))
        d_fwhm = 2 * amplitude * d_amplitude / fwhm * (1 - d_amplitude)
        return [d_amplitude, d_x_0, d_fwhm]


class Const1D(Fittable1DModel):
    """
    One dimensional Constant model.

    Parameters
    ----------
    amplitude : float
        Value of the constant function

    See Also
    --------
    Const2D

    Notes
    -----
    Model formula:

        .. math:: f(x) = A
    """

    amplitude = Parameter(default=1)
    linear = True

    @staticmethod
    def evaluate(x, amplitude):
        """One dimensional Constant model function"""

        if amplitude.size == 1:
            # This is slightly faster than using ones_like and multiplying
            x = np.empty_like(x)
            x.fill(amplitude.item())
        else:
            # This case is less likely but could occur if the amplitude
            # parameter is given an array-like value
            x = amplitude * np.ones_like(x)

        return x

    @staticmethod
    def fit_deriv(x, amplitude):
        """One dimensional Constant model derivative with respect to parameters"""

        d_amplitude = np.ones_like(x)
        return [d_amplitude]


class Const2D(Fittable2DModel):
    """
    Two dimensional Constant model.

    Parameters
    ----------
    amplitude : float
        Value of the constant function

    See Also
    --------
    Const1D

    Notes
    -----
    Model formula:

        .. math:: f(x, y) = A
    """

    amplitude = Parameter(default=1)
    linear = True

    @staticmethod
    def evaluate(x, y, amplitude):
        """Two dimensional Constant model function"""

        if amplitude.size == 1:
            # This is slightly faster than using ones_like and multiplying
            x = np.empty_like(x)
            x.fill(amplitude.item())
        else:
            # This case is less likely but could occur if the amplitude
            # parameter is given an array-like value
            x = amplitude * np.ones_like(x)

        return x


class Ellipse2D(Fittable2DModel):
    """
    A 2D Ellipse model.

    Parameters
    ----------
    amplitude : float
        Value of the ellipse.

    x_0 : float
        x position of the center of the disk.

    y_0 : float
        y position of the center of the disk.

    a : float
        The length of the semimajor axis.

    b : float
        The length of the semiminor axis.

    theta : float
        The rotation angle in radians of the semimajor axis.  The
        rotation angle increases counterclockwise from the positive x
        axis.

    See Also
    --------
    Disk2D, Box2D

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = \\left \\{
                    \\begin{array}{ll}
                      \\mathrm{amplitude} & : \\left[\\frac{(x - x_0) \\cos
                        \\theta + (y - y_0) \\sin \\theta}{a}\\right]^2 +
                        \\left[\\frac{-(x - x_0) \\sin \\theta + (y - y_0)
                        \\cos \\theta}{b}\\right]^2  \\leq 1 \\\\
                      0 & : \\mathrm{otherwise}
                    \\end{array}
                  \\right.

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        from astropy.modeling.models import Ellipse2D
        from astropy.coordinates import Angle
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        x0, y0 = 25, 25
        a, b = 20, 10
        theta = Angle(30, 'deg')
        e = Ellipse2D(amplitude=100., x_0=x0, y_0=y0, a=a, b=b,
                      theta=theta.radian)
        y, x = np.mgrid[0:50, 0:50]
        fig, ax = plt.subplots(1, 1)
        ax.imshow(e(x, y), origin='lower', interpolation='none', cmap='Greys_r')
        e2 = mpatches.Ellipse((x0, y0), 2*a, 2*b, theta.degree,
                              edgecolor='red', facecolor='none')
        ax.add_patch(e2)
        plt.show()
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    a = Parameter(default=1)
    b = Parameter(default=1)
    theta = Parameter(default=0)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, a, b, theta):
        """Two dimensional Ellipse model function."""

        xx = x - x_0
        yy = y - y_0
        cost = np.cos(theta)
        sint = np.sin(theta)
        numerator1 = (xx * cost) + (yy * sint)
        numerator2 = -(xx * sint) + (yy * cost)
        in_ellipse = (((numerator1 / a) ** 2 + (numerator2 / b) ** 2) <= 1.)
        return np.select([in_ellipse], [amplitude])


class Disk2D(Fittable2DModel):
    """
    Two dimensional radial symmetric Disk model.

    Parameters
    ----------
    amplitude : float
        Value of the disk function
    x_0 : float
        x position center of the disk
    y_0 : float
        y position center of the disk
    R_0 : float
        Radius of the disk

    See Also
    --------
    Box2D, TrapezoidDisk2D

    Notes
    -----
    Model formula:

        .. math::

            f(r) = \\left \\{
                     \\begin{array}{ll}
                       A & : r \\leq R_0 \\\\
                       0 & : r > R_0
                     \\end{array}
                   \\right.
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    R_0 = Parameter(default=1)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, R_0):
        """Two dimensional Disk model function"""

        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        return np.select([rr <= R_0 ** 2], [amplitude])


class Ring2D(Fittable2DModel):

    """
    Two dimensional radial symmetric Ring model.

    Parameters
    ----------
    amplitude : float
        Value of the disk function
    x_0 : float
        x position center of the disk
    y_0 : float
        y position center of the disk
    r_in : float
        Inner radius of the ring
    width : float
        Width of the ring.
    r_out : float
        Outer Radius of the ring. Can be specified instead of width.

    See Also
    --------
    Disk2D, TrapezoidDisk2D

    Notes
    -----
    Model formula:

        .. math::

            f(r) = \\left \\{
                     \\begin{array}{ll}
                       A & : r_{in} \\leq r \\leq r_{out} \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.

    Where :math:`r_{out} = r_{in} + r_{width}`.
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    r_in = Parameter(default=1)
    width = Parameter(default=1)

    def __init__(self, amplitude=amplitude.default, x_0=x_0.default,
                 y_0=y_0.default, r_in=r_in.default, width=width.default,
                 r_out=None, **kwargs):
        if r_out is not None:
            if width is not None:
                raise InputParameterError(
                    "Cannot specify both width and outer radius separately.")
            width = r_out - r_in
        elif width is None:
            width = self.width.default

        super(Ring2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, r_in=r_in, width=width,
            **kwargs)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, r_in, width):
        """Two dimensional Ring model function."""

        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        r_range = np.logical_and(rr >= r_in ** 2, rr <= (r_in + width) ** 2)
        return np.select([r_range], [amplitude])


class Delta1D(Fittable1DModel):
    """One dimensional Dirac delta function."""

    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Delta2D(Fittable2DModel):
    """Two dimensional Dirac delta function."""

    def __init__(self):
        raise ModelDefinitionError("Not implemented")


class Box1D(Fittable1DModel):
    """
    One dimensional Box model.

    Parameters
    ----------
    amplitude : float
        Amplitude A
    x_0 : float
        Position of the center of the box function
    width : float
        Width of the box

    See Also
    --------
    Box2D, TrapezoidDisk2D

    Notes
    -----
    Model formula:

      .. math::

            f(x) = \\left \\{
                     \\begin{array}{ll}
                       A & : x_0 - w/2 \\geq x \\geq x_0 + w/2 \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    width = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, width):
        """One dimensional Box model function"""

        return np.select([np.logical_and(x >= x_0 - width / 2.,
                                         x <= x_0 + width / 2.)],
                         [amplitude], 0)

    @classmethod
    def fit_deriv(cls, x, amplitude, x_0, width):
        """One dimensional Box model derivative with respect to parameters"""

        d_amplitude = cls.evaluate(x, 1, x_0, width)
        d_x_0 = np.zeros_like(x)
        d_width = np.zeros_like(x)
        return [d_amplitude, d_x_0, d_width]


class Box2D(Fittable2DModel):
    """
    Two dimensional Box model.

    Parameters
    ----------
    amplitude : float
        Amplitude A
    x_0 : float
        x position of the center of the box function
    x_width : float
        Width in x direction of the box
    y_0 : float
        y position of the center of the box function
    y_width : float
        Width in y direction of the box

    See Also
    --------
    Box1D, Gaussian2D, Moffat2D

    Notes
    -----
    Model formula:

      .. math::

            f(x, y) = \\left \\{
                     \\begin{array}{ll}
                       A & : x_0 - w_x/2 \\geq x \\geq x_0 + w_x/2 \\\\
                       A & : y_0 - w_y/2 \\geq y \\geq y_0 + w_y/2 \\\\
                       0 & : \\textnormal{else}
                     \\end{array}
                   \\right.

    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    x_width = Parameter(default=1)
    y_width = Parameter(default=1)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, x_width, y_width):
        """Two dimensional Box model function"""

        x_range = np.logical_and(x >= x_0 - x_width / 2.,
                                 x <= x_0 + x_width / 2.)
        y_range = np.logical_and(y >= y_0 - y_width / 2.,
                                 y <= y_0 + y_width / 2.)
        return np.select([np.logical_and(x_range, y_range)], [amplitude], 0)


class Trapezoid1D(Fittable1DModel):
    """
    One dimensional Trapezoid model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the trapezoid
    x_0 : float
        Center position of the trapezoid
    width : float
        Width of the constant part of the trapezoid.
    slope : float
        Slope of the tails of the trapezoid

    See Also
    --------
    Box1D, Gaussian1D, Moffat1D
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    width = Parameter(default=1)
    slope = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, width, slope):
        """One dimensional Trapezoid model function"""

        # Compute the four points where the trapezoid changes slope
        # x1 <= x2 <= x3 <= x4
        x2 = x_0 - width / 2.
        x3 = x_0 + width / 2.
        x1 = x2 - amplitude / slope
        x4 = x3 + amplitude / slope

        # Compute model values in pieces between the change points
        range_a = np.logical_and(x >= x1, x < x2)
        range_b = np.logical_and(x >= x2, x < x3)
        range_c = np.logical_and(x >= x3, x < x4)
        val_a = slope * (x - x1)
        val_b = amplitude
        val_c = slope * (x4 - x)
        return np.select([range_a, range_b, range_c], [val_a, val_b, val_c])


class TrapezoidDisk2D(Fittable2DModel):
    """
    Two dimensional circular Trapezoid model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the trapezoid
    x_0 : float
        x position of the center of the trapezoid
    y_0 : float
        y position of the center of the trapezoid
    R_0 : float
        Radius of the constant part of the trapezoid.
    slope : float
        Slope of the tails of the trapezoid in x direction.

    See Also
    --------
    Disk2D, Box2D
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    R_0 = Parameter(default=1)
    slope = Parameter(default=1)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, R_0, slope):
        """Two dimensional Trapezoid Disk model function"""

        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2)
        range_1 = r <= R_0
        range_2 = np.logical_and(r > R_0,  r <= R_0 + amplitude / slope)
        val_1 = amplitude
        val_2 = amplitude + slope * (R_0 - r)
        return np.select([range_1, range_2], [val_1, val_2])


class MexicanHat1D(Fittable1DModel):
    """
    One dimensional Mexican Hat model.

    Parameters
    ----------
    amplitude : float
        Amplitude
    x_0 : float
        Position of the peak
    sigma : float
        Width of the Mexican hat

    See Also
    --------
    MexicanHat2D, Box1D, Gaussian1D, Trapezoid1D

    Notes
    -----
    Model formula:

    .. math::

        f(x) = {A \\left(1 - \\frac{\\left(x - x_{0}\\right)^{2}}{\\sigma^{2}}\\right)
        e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}}

    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    sigma = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, sigma):
        """One dimensional Mexican Hat model function"""

        xx_ww = (x - x_0) ** 2 / (2 * sigma ** 2)
        return amplitude * (1 - 2 * xx_ww) * np.exp(-xx_ww)


class MexicanHat2D(Fittable2DModel):
    """
    Two dimensional symmetric Mexican Hat model.

    Parameters
    ----------
    amplitude : float
        Amplitude
    x_0 : float
        x position of the peak
    y_0 : float
        y position of the peak
    sigma : float
        Width of the Mexican hat

    See Also
    --------
    MexicanHat1D, Gaussian2D

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 - \\frac{\\left(x - x_{0}\\right)^{2}
        + \\left(y - y_{0}\\right)^{2}}{\\sigma^{2}}\\right)
        e^{\\frac{- \\left(x - x_{0}\\right)^{2}
        - \\left(y - y_{0}\\right)^{2}}{2 \\sigma^{2}}}
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    sigma = Parameter(default=1)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, sigma):
        """Two dimensional Mexican Hat model function"""

        rr_ww = ((x - x_0) ** 2 + (y - y_0) ** 2) / (2 * sigma ** 2)
        return amplitude * (1 - rr_ww) * np.exp(- rr_ww)


class AiryDisk2D(Fittable2DModel):
    """
    Two dimensional Airy disk model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the Airy function.
    x_0 : float
        x position of the maximum of the Airy function.
    y_0 : float
        y position of the maximum of the Airy function.
    radius : float
        The radius of the Airy disk (radius of the first zero).

    See Also
    --------
    Box2D, TrapezoidDisk2D, Gaussian2D

    Notes
    -----
    Model formula:

        .. math:: f(r) = A \\left[\\frac{2 J_1(\\frac{\\pi r}{R/R_z})}{\\frac{\\pi r}{R/R_z}}\\right]^2

    Where :math:`J_1` is the first order Bessel function of the first
    kind, :math:`r` is radial distance from the maximum of the Airy
    function (:math:`r = \\sqrt{(x - x_0)^2 + (y - y_0)^2}`), :math:`R`
    is the input ``radius`` parameter, and :math:`R_z =
    1.2196698912665045`).

    For an optical system, the radius of the first zero represents the
    limiting angular resolution and is approximately 1.22 * lambda / D,
    where lambda is the wavelength of the light and D is the diameter of
    the aperture.

    See [1]_ for more details about the Airy disk.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Airy_disk
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    radius = Parameter(default=1)
    _j1 = None

    def __init__(self, amplitude=amplitude.default, x_0=x_0.default,
                 y_0=y_0.default, radius=radius.default, **kwargs):
        if self._j1 is None:
            try:
                from scipy.special import j1, jn_zeros
                self.__class__._j1 = j1
                self.__class__._rz = jn_zeros(1, 1)[0] / np.pi
            # add a ValueError here for python3 + scipy < 0.12
            except ValueError:
                raise ImportError("AiryDisk2D model requires scipy > 0.11.")

        super(AiryDisk2D, self).__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, radius=radius, **kwargs)

    # TODO: Why does this particular model have its own special __deepcopy__
    # and __copy__?  If it has anything to do with the use of the j_1 function
    # that should be reworked.
    def __deepcopy__(self, memo):
        new_model = self.__class__(self.amplitude.value, self.x_0.value,
                                   self.y_0.value, self.radius.value)
        return new_model

    def __copy__(self):
        new_model = self.__class__(self.amplitude.value, self.x_0.value,
                                   self.y_0.value, self.radius.value)
        return new_model

    @classmethod
    def evaluate(cls, x, y, amplitude, x_0, y_0, radius):
        """Two dimensional Airy model function"""

        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2) / (radius / cls._rz)
        # Since r can be zero, we have to take care to treat that case
        # separately so as not to raise a numpy warning
        z = np.ones(r.shape)
        rt = np.pi * r[r > 0]
        z[r > 0] = (2.0 * cls._j1(rt) / rt)**2
        z *= amplitude
        return z


class Moffat1D(Fittable1DModel):
    """
    One dimensional Moffat model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the model.
    x_0 : float
        x position of the maximum of the Moffat model.
    gamma : float
        Core width of the Moffat model.
    alpha : float
        Power index of the Moffat model.

    See Also
    --------
    Gaussian1D, Box1D

    Notes
    -----
    Model formula:

    .. math::

        f(x) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    gamma = Parameter(default=1)
    alpha = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, gamma, alpha):
        """One dimensional Moffat model function"""

        return amplitude * (1 + ((x - x_0) / gamma) ** 2) ** (-alpha)

    @staticmethod
    def fit_deriv(x, amplitude, x_0, gamma, alpha):
        """One dimensional Moffat model derivative with respect to parameters"""

        d_A = (1 + (x - x_0) ** 2 / gamma ** 2) ** (-alpha)
        d_x_0 = (-amplitude * alpha * d_A * (-2 * x + 2 * x_0) /
                 (gamma ** 2 * d_A ** alpha))
        d_gamma = (2 * amplitude * alpha * d_A * (x - x_0) ** 2 /
                   (gamma ** 3 * d_A ** alpha))
        d_alpha = -amplitude * d_A * np.log(1 + (x - x_0) ** 2 / gamma ** 2)
        return [d_A, d_x_0, d_gamma, d_alpha]


class Moffat2D(Fittable2DModel):
    """
    Two dimensional Moffat model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the model.
    x_0 : float
        x position of the maximum of the Moffat model.
    y_0 : float
        y position of the maximum of the Moffat model.
    gamma : float
        Core width of the Moffat model.
    alpha : float
        Power index of the Moffat model.

    See Also
    --------
    Gaussian2D, Box2D

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2} +
        \\left(y - y_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    gamma = Parameter(default=1)
    alpha = Parameter(default=1)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, gamma, alpha):
        """Two dimensional Moffat model function"""

        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        return amplitude * (1 + rr_gg) ** (-alpha)

    @staticmethod
    def fit_deriv(x, y, amplitude, x_0, y_0, gamma, alpha):
        """Two dimensional Moffat model derivative with respect to parameters"""

        rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        d_A = (1 + rr_gg) ** (-alpha)
        d_x_0 = (-amplitude * alpha * d_A * (-2 * x + 2 * x_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_y_0 = (-amplitude * alpha * d_A * (-2 * y + 2 * y_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_alpha = -amplitude * d_A * np.log(1 + rr_gg)
        d_gamma = 2 * amplitude * alpha * d_A * (rr_gg / (gamma * (1 + rr_gg)))
        return [d_A, d_x_0, d_y_0, d_gamma, d_alpha]


@deprecated('1.0', alternative='astropy.modeling.models.custom_model',
            pending=True)
def custom_model_1d(func, func_fit_deriv=None):
    argspec = inspect.getargspec(func)
    param_values = argspec.defaults
    nparams = len(param_values)

    if len(argspec.args) == nparams + 1:
        param_names = argspec.args[-nparams:]
    else:
        raise ModelDefinitionError(
            "All parameters must be keyword arguments")

    return custom_model(func, fit_deriv=func_fit_deriv)
