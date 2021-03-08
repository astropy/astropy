# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Mathematical models."""
# pylint: disable=line-too-long, too-many-lines, too-many-arguments, invalid-name
import numpy as np

from astropy import units as u
from astropy.units import Quantity, UnitsError
from astropy.utils.decorators import deprecated
from .core import (Fittable1DModel, Fittable2DModel)

from .parameters import Parameter, InputParameterError
from .utils import ellipse_extent


__all__ = ['AiryDisk2D', 'Moffat1D', 'Moffat2D', 'Box1D', 'Box2D', 'Const1D',
           'Const2D', 'Ellipse2D', 'Disk2D', 'Gaussian1D', 'Gaussian2D',
           'Linear1D', 'Lorentz1D', 'RickerWavelet1D', 'RickerWavelet2D',
           'RedshiftScaleFactor', 'Multiply', 'Planar2D', 'Scale',
           'Sersic1D', 'Sersic2D', 'Shift', 'Sine1D', 'Trapezoid1D',
           'TrapezoidDisk2D', 'Ring2D', 'Voigt1D', 'KingProjectedAnalytic1D',
           'Exponential1D', 'Logarithmic1D']

TWOPI = 2 * np.pi
FLOAT_EPSILON = float(np.finfo(np.float32).tiny)

# Note that we define this here rather than using the value defined in
# astropy.stats to avoid importing astropy.stats every time astropy.modeling
# is loaded.
GAUSSIAN_SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


class Gaussian1D(Fittable1DModel):
    """
    One dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float or `~astropy.units.Quantity`.
        Amplitude (peak value) of the Gaussian - for a normalized profile
        (integrating to 1), set amplitude = 1 / (stddev * np.sqrt(2 * np.pi))
    mean : float or `~astropy.units.Quantity`.
        Mean of the Gaussian.
    stddev : float or `~astropy.units.Quantity`.
        Standard deviation of the Gaussian with FWHM = 2 * stddev * np.sqrt(2 * np.log(2)).

    Notes
    -----
    Either all or none of input ``x``, ``mean`` and ``stddev`` must be provided
    consistently with compatible units or as unitless numbers.

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

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Gaussian1D

        plt.figure()
        s1 = Gaussian1D()
        r = np.arange(-5, 5, .01)

        for factor in range(1, 4):
            s1.amplitude = factor
            plt.plot(r, s1(r), color=str(0.25 * factor), lw=2)

        plt.axis([-5, 5, -1, 4])
        plt.show()

    See Also
    --------
    Gaussian2D, Box1D, Moffat1D, Lorentz1D
    """

    amplitude = Parameter(default=1)
    mean = Parameter(default=0)

    # Ensure stddev makes sense if its bounds are not explicitly set.
    # stddev must be non-zero and positive.
    stddev = Parameter(default=1, bounds=(FLOAT_EPSILON, None))

    def bounding_box(self, factor=5.5):
        """
        Tuple defining the default ``bounding_box`` limits,
        ``(x_low, x_high)``

        Parameters
        ----------
        factor : float
            The multiple of `stddev` used to define the limits.
            The default is 5.5, corresponding to a relative error < 1e-7.

        Examples
        --------
        >>> from astropy.modeling.models import Gaussian1D
        >>> model = Gaussian1D(mean=0, stddev=2)
        >>> model.bounding_box
        (-11.0, 11.0)

        This range can be set directly (see: `Model.bounding_box
        <astropy.modeling.Model.bounding_box>`) or by using a different factor,
        like:

        >>> model.bounding_box = model.bounding_box(factor=2)
        >>> model.bounding_box
        (-4.0, 4.0)
        """

        x0 = self.mean
        dx = factor * self.stddev

        return (x0 - dx, x0 + dx)

    @property
    def fwhm(self):
        """Gaussian full width at half maximum."""
        return self.stddev * GAUSSIAN_SIGMA_TO_FWHM

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

    @property
    def input_units(self):
        if self.mean.unit is None:
            return None
        return {self.inputs[0]: self.mean.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'mean': inputs_unit[self.inputs[0]],
                'stddev': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


class Gaussian2D(Fittable2DModel):
    r"""
    Two dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float or `~astropy.units.Quantity`.
        Amplitude (peak value) of the Gaussian.
    x_mean : float or `~astropy.units.Quantity`.
        Mean of the Gaussian in x.
    y_mean : float or `~astropy.units.Quantity`.
        Mean of the Gaussian in y.
    x_stddev : float or `~astropy.units.Quantity` or None.
        Standard deviation of the Gaussian in x before rotating by theta. Must
        be None if a covariance matrix (``cov_matrix``) is provided. If no
        ``cov_matrix`` is given, ``None`` means the default value (1).
    y_stddev : float or `~astropy.units.Quantity` or None.
        Standard deviation of the Gaussian in y before rotating by theta. Must
        be None if a covariance matrix (``cov_matrix``) is provided. If no
        ``cov_matrix`` is given, ``None`` means the default value (1).
    theta : float or `~astropy.units.Quantity`, optional.
        Rotation angle (value in radians). The rotation angle increases
        counterclockwise.  Must be None if a covariance matrix (``cov_matrix``)
        is provided. If no ``cov_matrix`` is given, ``None`` means the default
        value (0).
    cov_matrix : ndarray, optional
        A 2x2 covariance matrix. If specified, overrides the ``x_stddev``,
        ``y_stddev``, and ``theta`` defaults.

    Notes
    -----
    Either all or none of input ``x, y``, ``[x,y]_mean`` and ``[x,y]_stddev``
    must be provided consistently with compatible units or as unitless numbers.

    Model formula:

        .. math::

            f(x, y) = A e^{-a\left(x - x_{0}\right)^{2}  -b\left(x - x_{0}\right)
            \left(y - y_{0}\right)  -c\left(y - y_{0}\right)^{2}}

    Using the following definitions:

        .. math::
            a = \left(\frac{\cos^{2}{\left (\theta \right )}}{2 \sigma_{x}^{2}} +
            \frac{\sin^{2}{\left (\theta \right )}}{2 \sigma_{y}^{2}}\right)

            b = \left(\frac{\sin{\left (2 \theta \right )}}{2 \sigma_{x}^{2}} -
            \frac{\sin{\left (2 \theta \right )}}{2 \sigma_{y}^{2}}\right)

            c = \left(\frac{\sin^{2}{\left (\theta \right )}}{2 \sigma_{x}^{2}} +
            \frac{\cos^{2}{\left (\theta \right )}}{2 \sigma_{y}^{2}}\right)

    If using a ``cov_matrix``, the model is of the form:
        .. math::
            f(x, y) = A e^{-0.5 \left(\vec{x} - \vec{x}_{0}\right)^{T} \Sigma^{-1} \left(\vec{x} - \vec{x}_{0}\right)}

    where :math:`\vec{x} = [x, y]`, :math:`\vec{x}_{0} = [x_{0}, y_{0}]`,
    and :math:`\Sigma` is the covariance matrix:

        .. math::
            \Sigma = \left(\begin{array}{ccc}
            \sigma_x^2               & \rho \sigma_x \sigma_y \\
            \rho \sigma_x \sigma_y   & \sigma_y^2
            \end{array}\right)

    :math:`\rho` is the correlation between ``x`` and ``y``, which should
    be between -1 and +1.  Positive correlation corresponds to a
    ``theta`` in the range 0 to 90 degrees.  Negative correlation
    corresponds to a ``theta`` in the range of 0 to -90 degrees.

    See [1]_ for more details about the 2D Gaussian function.

    See Also
    --------
    Gaussian1D, Box2D, Moffat2D

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Gaussian_function
    """

    amplitude = Parameter(default=1)
    x_mean = Parameter(default=0)
    y_mean = Parameter(default=0)
    x_stddev = Parameter(default=1)
    y_stddev = Parameter(default=1)
    theta = Parameter(default=0.0)

    def __init__(self, amplitude=amplitude.default, x_mean=x_mean.default,
                 y_mean=y_mean.default, x_stddev=None, y_stddev=None,
                 theta=None, cov_matrix=None, **kwargs):
        if cov_matrix is None:
            if x_stddev is None:
                x_stddev = self.__class__.x_stddev.default
            if y_stddev is None:
                y_stddev = self.__class__.y_stddev.default
            if theta is None:
                theta = self.__class__.theta.default
        else:
            if x_stddev is not None or y_stddev is not None or theta is not None:
                raise InputParameterError("Cannot specify both cov_matrix and "
                                          "x/y_stddev/theta")
            # Compute principle coordinate system transformation
            cov_matrix = np.array(cov_matrix)

            if cov_matrix.shape != (2, 2):
                raise ValueError("Covariance matrix must be 2x2")

            eig_vals, eig_vecs = np.linalg.eig(cov_matrix)
            x_stddev, y_stddev = np.sqrt(eig_vals)
            y_vec = eig_vecs[:, 0]
            theta = np.arctan2(y_vec[1], y_vec[0])

        # Ensure stddev makes sense if its bounds are not explicitly set.
        # stddev must be non-zero and positive.
        # TODO: Investigate why setting this in Parameter above causes
        #       convolution tests to hang.
        kwargs.setdefault('bounds', {})
        kwargs['bounds'].setdefault('x_stddev', (FLOAT_EPSILON, None))
        kwargs['bounds'].setdefault('y_stddev', (FLOAT_EPSILON, None))

        super().__init__(
            amplitude=amplitude, x_mean=x_mean, y_mean=y_mean,
            x_stddev=x_stddev, y_stddev=y_stddev, theta=theta, **kwargs)

    @property
    def x_fwhm(self):
        """Gaussian full width at half maximum in X."""
        return self.x_stddev * GAUSSIAN_SIGMA_TO_FWHM

    @property
    def y_fwhm(self):
        """Gaussian full width at half maximum in Y."""
        return self.y_stddev * GAUSSIAN_SIGMA_TO_FWHM

    def bounding_box(self, factor=5.5):
        """
        Tuple defining the default ``bounding_box`` limits in each dimension,
        ``((y_low, y_high), (x_low, x_high))``

        The default offset from the mean is 5.5-sigma, corresponding
        to a relative error < 1e-7. The limits are adjusted for rotation.

        Parameters
        ----------
        factor : float, optional
            The multiple of `x_stddev` and `y_stddev` used to define the limits.
            The default is 5.5.

        Examples
        --------
        >>> from astropy.modeling.models import Gaussian2D
        >>> model = Gaussian2D(x_mean=0, y_mean=0, x_stddev=1, y_stddev=2)
        >>> model.bounding_box
        ((-11.0, 11.0), (-5.5, 5.5))

        This range can be set directly (see: `Model.bounding_box
        <astropy.modeling.Model.bounding_box>`) or by using a different factor
        like:

        >>> model.bounding_box = model.bounding_box(factor=2)
        >>> model.bounding_box
        ((-4.0, 4.0), (-2.0, 2.0))
        """

        a = factor * self.x_stddev
        b = factor * self.y_stddev
        theta = self.theta.value
        dx, dy = ellipse_extent(a, b, theta)

        return ((self.y_mean - dy, self.y_mean + dy),
                (self.x_mean - dx, self.x_mean + dx))

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

    @property
    def input_units(self):
        if self.x_mean.unit is None and self.y_mean.unit is None:
            return None
        return {self.inputs[0]: self.x_mean.unit,
                self.inputs[1]: self.y_mean.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_mean': inputs_unit[self.inputs[0]],
                'y_mean': inputs_unit[self.inputs[0]],
                'x_stddev': inputs_unit[self.inputs[0]],
                'y_stddev': inputs_unit[self.inputs[0]],
                'theta': u.rad,
                'amplitude': outputs_unit[self.outputs[0]]}


class Shift(Fittable1DModel):
    """
    Shift a coordinate.

    Parameters
    ----------
    offset : float
        Offset to add to a coordinate.
    """

    offset = Parameter(default=0)
    linear = True

    @property
    def input_units(self):
        if self.offset.unit is None:
            return None
        return {self.inputs[0]: self.offset.unit}

    @property
    def inverse(self):
        """One dimensional inverse Shift model function"""
        inv = self.copy()
        inv.offset *= -1
        return inv

    @staticmethod
    def evaluate(x, offset):
        """One dimensional Shift model function"""
        return x + offset

    @staticmethod
    def sum_of_implicit_terms(x):
        """Evaluate the implicit term (x) of one dimensional Shift model"""
        return x

    @staticmethod
    def fit_deriv(x, *params):
        """One dimensional Shift model derivative with respect to parameter"""

        d_offset = np.ones_like(x)
        return [d_offset]

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'offset': outputs_unit[self.outputs[0]]}


class Scale(Fittable1DModel):
    """
    Multiply a model by a dimensionless factor.

    Parameters
    ----------
    factor : float
        Factor by which to scale a coordinate.

    Notes
    -----

    If ``factor`` is a `~astropy.units.Quantity` then the units will be
    stripped before the scaling operation.

    """

    factor = Parameter(default=1)
    linear = True
    fittable = True

    _input_units_strict = True
    _input_units_allow_dimensionless = True

    @property
    def input_units(self):
        if self.factor.unit is None:
            return None
        return {self.inputs[0]: self.factor.unit}

    @property
    def inverse(self):
        """One dimensional inverse Scale model function"""
        inv = self.copy()
        inv.factor = 1 / self.factor
        return inv

    @staticmethod
    def evaluate(x, factor):
        """One dimensional Scale model function"""
        if isinstance(factor, u.Quantity):
            factor = factor.value

        return factor * x

    @staticmethod
    def fit_deriv(x, *params):
        """One dimensional Scale model derivative with respect to parameter"""

        d_factor = x
        return [d_factor]

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'factor': outputs_unit[self.outputs[0]]}


class Multiply(Fittable1DModel):
    """
    Multiply a model by a quantity or number.

    Parameters
    ----------
    factor : float
        Factor by which to multiply a coordinate.
    """

    factor = Parameter(default=1)
    linear = True
    fittable = True

    @property
    def inverse(self):
        """One dimensional inverse multiply model function"""
        inv = self.copy()
        inv.factor = 1 / self.factor
        return inv

    @staticmethod
    def evaluate(x, factor):
        """One dimensional multiply model function"""
        return factor * x

    @staticmethod
    def fit_deriv(x, *params):
        """One dimensional multiply model derivative with respect to parameter"""

        d_factor = x
        return [d_factor]

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'factor': outputs_unit[self.outputs[0]]}


class RedshiftScaleFactor(Fittable1DModel):
    """
    One dimensional redshift scale factor model.

    Parameters
    ----------
    z : float
        Redshift value.

    Notes
    -----
    Model formula:

        .. math:: f(x) = x (1 + z)
    """

    z = Parameter(description='redshift', default=0)

    @staticmethod
    def evaluate(x, z):
        """One dimensional RedshiftScaleFactor model function"""

        return (1 + z) * x

    @staticmethod
    def fit_deriv(x, z):
        """One dimensional RedshiftScaleFactor model derivative"""

        d_z = x
        return [d_z]

    @property
    def inverse(self):
        """Inverse RedshiftScaleFactor model"""

        inv = self.copy()
        inv.z = 1.0 / (1.0 + self.z) - 1.0
        return inv


class Sersic1D(Fittable1DModel):
    r"""
    One dimensional Sersic surface brightness profile.

    Parameters
    ----------
    amplitude : float
        Surface brightness at r_eff.
    r_eff : float
        Effective (half-light) radius
    n : float
        Sersic Index.

    See Also
    --------
    Gaussian1D, Moffat1D, Lorentz1D

    Notes
    -----
    Model formula:

    .. math::

        I(r)=I_e\exp\left\{-b_n\left[\left(\frac{r}{r_{e}}\right)^{(1/n)}-1\right]\right\}

    The constant :math:`b_n` is defined such that :math:`r_e` contains half the total
    luminosity, and can be solved for numerically.

    .. math::

        \Gamma(2n) = 2\gamma (b_n,2n)

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        from astropy.modeling.models import Sersic1D
        import matplotlib.pyplot as plt

        plt.figure()
        plt.subplot(111, xscale='log', yscale='log')
        s1 = Sersic1D(amplitude=1, r_eff=5)
        r=np.arange(0, 100, .01)

        for n in range(1, 10):
             s1.n = n
             plt.plot(r, s1(r), color=str(float(n) / 15))

        plt.axis([1e-1, 30, 1e-2, 1e3])
        plt.xlabel('log Radius')
        plt.ylabel('log Surface Brightness')
        plt.text(.25, 1.5, 'n=1')
        plt.text(.25, 300, 'n=10')
        plt.xticks([])
        plt.yticks([])
        plt.show()

    References
    ----------
    .. [1] http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
    """

    amplitude = Parameter(default=1)
    r_eff = Parameter(default=1)
    n = Parameter(default=4)
    _gammaincinv = None

    @classmethod
    def evaluate(cls, r, amplitude, r_eff, n):
        """One dimensional Sersic profile function."""

        if cls._gammaincinv is None:
            try:
                from scipy.special import gammaincinv
                cls._gammaincinv = gammaincinv
            except ValueError:
                raise ImportError('Sersic1D model requires scipy.')

        return (amplitude * np.exp(
            -cls._gammaincinv(2 * n, 0.5) * ((r / r_eff) ** (1 / n) - 1)))

    @property
    def input_units(self):
        if self.r_eff.unit is None:
            return None
        return {self.inputs[0]: self.r_eff.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'r_eff': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


class Sine1D(Fittable1DModel):
    """
    One dimensional Sine model.

    Parameters
    ----------
    amplitude : float
        Oscillation amplitude
    frequency : float
        Oscillation frequency
    phase : float
        Oscillation phase

    See Also
    --------
    Const1D, Linear1D


    Notes
    -----
    Model formula:

        .. math:: f(x) = A \\sin(2 \\pi f x + 2 \\pi p)

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Sine1D

        plt.figure()
        s1 = Sine1D(amplitude=1, frequency=.25)
        r=np.arange(0, 10, .01)

        for amplitude in range(1,4):
             s1.amplitude = amplitude
             plt.plot(r, s1(r), color=str(0.25 * amplitude), lw=2)

        plt.axis([0, 10, -5, 5])
        plt.show()
    """

    amplitude = Parameter(default=1)
    frequency = Parameter(default=1)
    phase = Parameter(default=0)

    @staticmethod
    def evaluate(x, amplitude, frequency, phase):
        """One dimensional Sine model function"""
        # Note: If frequency and x are quantities, they should normally have
        # inverse units, so that argument ends up being dimensionless. However,
        # np.sin of a dimensionless quantity will crash, so we remove the
        # quantity-ness from argument in this case (another option would be to
        # multiply by * u.rad but this would be slower overall).
        argument = TWOPI * (frequency * x + phase)
        if isinstance(argument, Quantity):
            argument = argument.value
        return amplitude * np.sin(argument)

    @staticmethod
    def fit_deriv(x, amplitude, frequency, phase):
        """One dimensional Sine model derivative"""

        d_amplitude = np.sin(TWOPI * frequency * x + TWOPI * phase)
        d_frequency = (TWOPI * x * amplitude *
                       np.cos(TWOPI * frequency * x + TWOPI * phase))
        d_phase = (TWOPI * amplitude *
                   np.cos(TWOPI * frequency * x + TWOPI * phase))
        return [d_amplitude, d_frequency, d_phase]

    @property
    def input_units(self):
        if self.frequency.unit is None:
            return None
        return {self.inputs[0]: 1. / self.frequency.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'frequency': inputs_unit[self.inputs[0]] ** -1,
                'amplitude': outputs_unit[self.outputs[0]]}


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
    def fit_deriv(x, *params):
        """One dimensional Line model derivative with respect to parameters"""

        d_slope = x
        d_intercept = np.ones_like(x)
        return [d_slope, d_intercept]

    @property
    def inverse(self):
        new_slope = self.slope ** -1
        new_intercept = -self.intercept / self.slope
        return self.__class__(slope=new_slope, intercept=new_intercept)

    @property
    def input_units(self):
        if self.intercept.unit is None and self.slope.unit is None:
            return None
        return {self.inputs[0]: self.intercept.unit / self.slope.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'intercept': outputs_unit[self.outputs[0]],
                'slope': outputs_unit[self.outputs[0]] / inputs_unit[self.inputs[0]]}


class Planar2D(Fittable2DModel):
    """
    Two dimensional Plane model.

    Parameters
    ----------
    slope_x : float
        Slope of the straight line in X

    slope_y : float
        Slope of the straight line in Y

    intercept : float
        Z-intercept of the straight line

    Notes
    -----
    Model formula:

        .. math:: f(x, y) = a x + b y + c
    """

    slope_x = Parameter(default=1)
    slope_y = Parameter(default=1)
    intercept = Parameter(default=0)
    linear = True

    @staticmethod
    def evaluate(x, y, slope_x, slope_y, intercept):
        """Two dimensional Plane model function"""

        return slope_x * x + slope_y * y + intercept

    @staticmethod
    def fit_deriv(x, y, *params):
        """Two dimensional Plane model derivative with respect to parameters"""

        d_slope_x = x
        d_slope_y = y
        d_intercept = np.ones_like(x)
        return [d_slope_x, d_slope_y, d_intercept]

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'intercept': outputs_unit['z'],
                'slope_x': outputs_unit['z'] / inputs_unit['x'],
                'slope_y': outputs_unit['z'] / inputs_unit['y']}


class Lorentz1D(Fittable1DModel):
    """
    One dimensional Lorentzian model.

    Parameters
    ----------
    amplitude : float or `~astropy.units.Quantity`.
        Peak value - for a normalized profile (integrating to 1),
        set amplitude = 2 / (np.pi * fwhm)
    x_0 : float or `~astropy.units.Quantity`.
        Position of the peak
    fwhm : float or `~astropy.units.Quantity`.
        Full width at half maximum (FWHM)

    See Also
    --------
    Gaussian1D, Box1D, RickerWavelet1D

    Notes
    -----
    Either all or none of input ``x``, position ``x_0`` and ``fwhm`` must be provided
    consistently with compatible units or as unitless numbers.

    Model formula:

    .. math::

        f(x) = \\frac{A \\gamma^{2}}{\\gamma^{2} + \\left(x - x_{0}\\right)^{2}}

    where :math:`\\gamma` is half of given FWHM.

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Lorentz1D

        plt.figure()
        s1 = Lorentz1D()
        r = np.arange(-5, 5, .01)

        for factor in range(1, 4):
            s1.amplitude = factor
            plt.plot(r, s1(r), color=str(0.25 * factor), lw=2)

        plt.axis([-5, 5, -1, 4])
        plt.show()
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

    def bounding_box(self, factor=25):
        """Tuple defining the default ``bounding_box`` limits,
        ``(x_low, x_high)``.

        Parameters
        ----------
        factor : float
            The multiple of FWHM used to define the limits.
            Default is chosen to include most (99%) of the
            area under the curve, while still showing the
            central feature of interest.

        """
        x0 = self.x_0
        dx = factor * self.fwhm

        return (x0 - dx, x0 + dx)

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'x_0': inputs_unit[self.inputs[0]],
                'fwhm': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


class Voigt1D(Fittable1DModel):
    """
    One dimensional model for the Voigt profile.

    Parameters
    ----------
    x_0 : float or `~astropy.units.Quantity`
        Position of the peak
    amplitude_L : float or `~astropy.units.Quantity`.
        The Lorentzian amplitude (peak of the associated Lorentz function)
        - for a normalized profile (integrating to 1), set
        amplitude_L = 2 / (np.pi * fwhm_L)
    fwhm_L : float or `~astropy.units.Quantity`
        The Lorentzian full width at half maximum
    fwhm_G : float or `~astropy.units.Quantity`.
        The Gaussian full width at half maximum
    method : str, optional
        Algorithm for computing the complex error function; one of
        'Humlicek2' (default, fast and generally more accurate than ``rtol=3.e-5``) or
        'Scipy', alternatively 'wofz' (requires `scipy`, almost as fast and reference in accuracy).

    See Also
    --------
    Gaussian1D, Lorentz1D

    Notes
    -----
    Either all or none of input ``x``, position ``x_0`` and the ``fwhm_*`` must be provided
    consistently with compatible units or as unitless numbers.
    Voigt function is calculated as real part of the complex error function computed from either
    Humlicek's rational approximations (JQSRT 21:309, 1979; 27:437, 1982) following
    Schreier 2018 (MNRAS 479, 3068; and ``hum2zpf16m`` from his cpfX.py module); or
    `~scipy.special.wofz` (implementing 'Faddeeva.cc').

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        from astropy.modeling.models import Voigt1D
        import matplotlib.pyplot as plt

        plt.figure()
        x = np.arange(0, 10, 0.01)
        v1 = Voigt1D(x_0=5, amplitude_L=10, fwhm_L=0.5, fwhm_G=0.9)
        plt.plot(x, v1(x))
        plt.show()
    """

    x_0 = Parameter(default=0)
    amplitude_L = Parameter(default=1)     # noqa: N815
    fwhm_L = Parameter(default=2/np.pi)    # noqa: N815
    fwhm_G = Parameter(default=np.log(2))  # noqa: N815

    sqrt_pi = np.sqrt(np.pi)
    sqrt_ln2 = np.sqrt(np.log(2))
    sqrt_ln2pi = np.sqrt(np.log(2) * np.pi)
    _last_z = np.zeros(1, dtype=complex)
    _last_w = np.zeros(1, dtype=float)
    _faddeeva = None

    def __init__(self, x_0=x_0.default, amplitude_L=amplitude_L.default,            # noqa: N803
                 fwhm_L=fwhm_L.default, fwhm_G=fwhm_G.default, method='humlicek2',  # noqa: N803
                 **kwargs):
        if str(method).lower() in ('wofz', 'scipy'):
            try:
                from scipy.special import wofz
            except (ValueError, ImportError) as err:
                raise ImportError(f'Voigt1D method {method} requires scipy: {err}.') from err
            self._faddeeva = wofz
        elif str(method).lower() == 'humlicek2':
            self._faddeeva = self._hum2zpf16c
        else:
            raise ValueError(f'Not a valid method for Voigt1D Faddeeva function: {method}.')
        self.method = self._faddeeva.__name__

        super().__init__(x_0=x_0, amplitude_L=amplitude_L, fwhm_L=fwhm_L, fwhm_G=fwhm_G, **kwargs)

    def _wrap_wofz(self, z):
        """Call complex error (Faddeeva) function w(z) implemented by algorithm `method`;
        cache results for consecutive calls from `evaluate`, `fit_deriv`."""

        if (z.shape == self._last_z.shape and
                np.allclose(z, self._last_z, rtol=1.e-14, atol=1.e-15)):
            return self._last_w

        self._last_w = self._faddeeva(z)
        self._last_z = z
        return self._last_w

    def evaluate(self, x, x_0, amplitude_L, fwhm_L, fwhm_G):  # noqa: N803
        """One dimensional Voigt function scaled to Lorentz peak amplitude."""

        z = np.atleast_1d(2 * (x - x_0) + 1j * fwhm_L) * self.sqrt_ln2 / fwhm_G
        # The normalised Voigt profile is w.real * self.sqrt_ln2 / (self.sqrt_pi * fwhm_G) * 2 ;
        # for the legacy definition we multiply with np.pi * fwhm_L / 2 * amplitude_L
        return self._wrap_wofz(z).real * self.sqrt_ln2pi / fwhm_G * fwhm_L * amplitude_L

    def fit_deriv(self, x, x_0, amplitude_L, fwhm_L, fwhm_G):  # noqa: N803
        """Derivative of the one dimensional Voigt function with respect to parameters."""

        s = self.sqrt_ln2 / fwhm_G
        z = np.atleast_1d(2 * (x - x_0) + 1j * fwhm_L) * s
        # V * constant from McLean implementation (== their Voigt function)
        w = self._wrap_wofz(z) * s * fwhm_L * amplitude_L * self.sqrt_pi

        # Schreier (2018) Eq. 6 == (dvdx + 1j * dvdy) / (sqrt(pi) * fwhm_L * amplitude_L)
        dwdz = -2 * z * w + 2j * s * fwhm_L * amplitude_L

        return [-dwdz.real * 2 * s,
                w.real / amplitude_L,
                w.real / fwhm_L - dwdz.imag * s,
                (-w.real - s * (2 * (x - x_0) * dwdz.real - fwhm_L * dwdz.imag)) / fwhm_G]

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'x_0': inputs_unit[self.inputs[0]],
                'fwhm_L': inputs_unit[self.inputs[0]],
                'fwhm_G': inputs_unit[self.inputs[0]],
                'amplitude_L': outputs_unit[self.outputs[0]]}

    @staticmethod
    def _hum2zpf16c(z, s=10.0):
        """Complex error function w(z) for z = x + iy combining Humlicek's rational approximations:

        |x| + y > 10:  Humlicek (JQSRT, 1982) rational approximation for region II;
        else:          Humlicek (JQSRT, 1979) rational approximation with n=16 and delta=y0=1.35

        Version using a mask and np.place;
        single complex argument version of Franz Schreier's cpfX.hum2zpf16m.
        Originally licensed under a 3-clause BSD style license - see
        https://atmos.eoc.dlr.de/tools/lbl4IR/cpfX.py
        """

        # Optimized (single fraction) Humlicek region I rational approximation for n=16, delta=1.35

        AA = np.array([+46236.3358828121,   -147726.58393079657j,   # noqa: N806
                       -206562.80451354137,  281369.1590631087j,
                       +183092.74968253175, -184787.96830696272j,
                       -66155.39578477248,   57778.05827983565j,
                       +11682.770904216826, -9442.402767960672j,
                       -1052.8438624933142,  814.0996198624186j,
                       +45.94499030751872,  -34.59751573708725j,
                       -0.7616559377907136,  0.5641895835476449j])  # 1j/sqrt(pi) to the 12. digit

        bb = np.array([+7918.06640624997, 0.0,
                       -126689.0625,      0.0,
                       +295607.8125,      0.0,
                       -236486.25,        0.0,
                       +84459.375,        0.0,
                       -15015.0,          0.0,
                       +1365.0,           0.0,
                       -60.0,             0.0,
                       +1.0])

        sqrt_piinv = 1.0 / np.sqrt(np.pi)

        zz = z * z
        w  = 1j * (z * (zz * sqrt_piinv - 1.410474)) / (0.75 + zz*(zz - 3.0))

        if np.any(z.imag < s):
            mask  = abs(z.real) + z.imag < s  # returns true for interior points
            # returns small complex array covering only the interior region
            Z     = z[np.where(mask)] + 1.35j
            ZZ    = Z * Z
            numer = (((((((((((((((AA[15]*Z + AA[14])*Z + AA[13])*Z + AA[12])*Z + AA[11])*Z +
                               AA[10])*Z + AA[9])*Z + AA[8])*Z + AA[7])*Z + AA[6])*Z +
                          AA[5])*Z + AA[4])*Z+AA[3])*Z + AA[2])*Z + AA[1])*Z + AA[0])
            denom = (((((((ZZ + bb[14])*ZZ + bb[12])*ZZ + bb[10])*ZZ+bb[8])*ZZ + bb[6])*ZZ +
                      bb[4])*ZZ + bb[2])*ZZ + bb[0]
            np.place(w, mask, numer / denom)

        return w


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

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Const1D

        plt.figure()
        s1 = Const1D()
        r = np.arange(-5, 5, .01)

        for factor in range(1, 4):
            s1.amplitude = factor
            plt.plot(r, s1(r), color=str(0.25 * factor), lw=2)

        plt.axis([-5, 5, -1, 4])
        plt.show()
    """

    amplitude = Parameter(default=1)
    linear = True

    @staticmethod
    def evaluate(x, amplitude):
        """One dimensional Constant model function"""

        if amplitude.size == 1:
            # This is slightly faster than using ones_like and multiplying
            x = np.empty_like(x, subok=False)
            x.fill(amplitude.item())
        else:
            # This case is less likely but could occur if the amplitude
            # parameter is given an array-like value
            x = amplitude * np.ones_like(x, subok=False)

        if isinstance(amplitude, Quantity):
            return Quantity(x, unit=amplitude.unit, copy=False)
        return x

    @staticmethod
    def fit_deriv(x, amplitude):
        """One dimensional Constant model derivative with respect to parameters"""

        d_amplitude = np.ones_like(x)
        return [d_amplitude]

    @property
    def input_units(self):
        return None

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'amplitude': outputs_unit[self.outputs[0]]}


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
            x = np.empty_like(x, subok=False)
            x.fill(amplitude.item())
        else:
            # This case is less likely but could occur if the amplitude
            # parameter is given an array-like value
            x = amplitude * np.ones_like(x, subok=False)

        if isinstance(amplitude, Quantity):
            return Quantity(x, unit=amplitude.unit, copy=False)
        return x

    @property
    def input_units(self):
        return None

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'amplitude': outputs_unit[self.outputs[0]]}


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
        e2 = mpatches.Ellipse((x0, y0), 2*a, 2*b, theta.degree, edgecolor='red',
                              facecolor='none')
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
        result = np.select([in_ellipse], [amplitude])

        if isinstance(amplitude, Quantity):
            return Quantity(result, unit=amplitude.unit, copy=False)
        return result

    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box`` limits.

        ``((y_low, y_high), (x_low, x_high))``
        """

        a = self.a
        b = self.b
        theta = self.theta.value
        dx, dy = ellipse_extent(a, b, theta)

        return ((self.y_0 - dy, self.y_0 + dy),
                (self.x_0 - dx, self.x_0 + dx))

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit,
                self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[0]],
                'a': inputs_unit[self.inputs[0]],
                'b': inputs_unit[self.inputs[0]],
                'theta': u.rad,
                'amplitude': outputs_unit[self.outputs[0]]}


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
        result = np.select([rr <= R_0 ** 2], [amplitude])

        if isinstance(amplitude, Quantity):
            return Quantity(result, unit=amplitude.unit, copy=False)
        return result

    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box`` limits.

        ``((y_low, y_high), (x_low, x_high))``
        """

        return ((self.y_0 - self.R_0, self.y_0 + self.R_0),
                (self.x_0 - self.R_0, self.x_0 + self.R_0))

    @property
    def input_units(self):
        if self.x_0.unit is None and self.y_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit,
                self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[0]],
                'R_0': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


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
                       0 & : \\text{else}
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
        # If outer radius explicitly given, it overrides default width.
        if r_out is not None:
            if width != self.width.default:
                raise InputParameterError(
                    "Cannot specify both width and outer radius separately.")
            width = r_out - r_in
        elif width is None:
            width = self.width.default

        super().__init__(
            amplitude=amplitude, x_0=x_0, y_0=y_0, r_in=r_in, width=width,
            **kwargs)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, r_in, width):
        """Two dimensional Ring model function."""

        rr = (x - x_0) ** 2 + (y - y_0) ** 2
        r_range = np.logical_and(rr >= r_in ** 2, rr <= (r_in + width) ** 2)
        result = np.select([r_range], [amplitude])

        if isinstance(amplitude, Quantity):
            return Quantity(result, unit=amplitude.unit, copy=False)
        return result

    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box``.

        ``((y_low, y_high), (x_low, x_high))``
        """

        dr = self.r_in + self.width

        return ((self.y_0 - dr, self.y_0 + dr),
                (self.x_0 - dr, self.x_0 + dr))

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit,
                self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[0]],
                'r_in': inputs_unit[self.inputs[0]],
                'width': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


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
                       A & : x_0 - w/2 \\leq x \\leq x_0 + w/2 \\\\
                       0 & : \\text{else}
                     \\end{array}
                   \\right.

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Box1D

        plt.figure()
        s1 = Box1D()
        r = np.arange(-5, 5, .01)

        for factor in range(1, 4):
            s1.amplitude = factor
            s1.width = factor
            plt.plot(r, s1(r), color=str(0.25 * factor), lw=2)

        plt.axis([-5, 5, -1, 4])
        plt.show()
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    width = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, width):
        """One dimensional Box model function"""

        inside = np.logical_and(x >= x_0 - width / 2., x <= x_0 + width / 2.)
        return np.select([inside], [amplitude], 0)

    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box`` limits.

        ``(x_low, x_high))``
        """

        dx = self.width / 2

        return (self.x_0 - dx, self.x_0 + dx)

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit}

    @property
    def return_units(self):
        if self.amplitude.unit is None:
            return None
        return {self.outputs[0]: self.amplitude.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'x_0': inputs_unit[self.inputs[0]],
                'width': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


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
            A : & x_0 - w_x/2 \\leq x \\leq x_0 + w_x/2 \\text{ and} \\\\
                & y_0 - w_y/2 \\leq y \\leq y_0 + w_y/2 \\\\
            0 : & \\text{else}
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

        result = np.select([np.logical_and(x_range, y_range)], [amplitude], 0)

        if isinstance(amplitude, Quantity):
            return Quantity(result, unit=amplitude.unit, copy=False)
        return result

    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box``.

        ``((y_low, y_high), (x_low, x_high))``
        """

        dx = self.x_width / 2
        dy = self.y_width / 2

        return ((self.y_0 - dy, self.y_0 + dy),
                (self.x_0 - dx, self.x_0 + dx))

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit,
                self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[1]],
                'x_width': inputs_unit[self.inputs[0]],
                'y_width': inputs_unit[self.inputs[1]],
                'amplitude': outputs_unit[self.outputs[0]]}


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

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Trapezoid1D

        plt.figure()
        s1 = Trapezoid1D()
        r = np.arange(-5, 5, .01)

        for factor in range(1, 4):
            s1.amplitude = factor
            s1.width = factor
            plt.plot(r, s1(r), color=str(0.25 * factor), lw=2)

        plt.axis([-5, 5, -1, 4])
        plt.show()
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
        result = np.select([range_a, range_b, range_c], [val_a, val_b, val_c])

        if isinstance(amplitude, Quantity):
            return Quantity(result, unit=amplitude.unit, copy=False)
        return result

    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box`` limits.

        ``(x_low, x_high))``
        """

        dx = self.width / 2 + self.amplitude / self.slope

        return (self.x_0 - dx, self.x_0 + dx)

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'x_0': inputs_unit[self.inputs[0]],
                'width': inputs_unit[self.inputs[0]],
                'slope': outputs_unit[self.outputs[0]] / inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


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
        range_2 = np.logical_and(r > R_0, r <= R_0 + amplitude / slope)
        val_1 = amplitude
        val_2 = amplitude + slope * (R_0 - r)
        result = np.select([range_1, range_2], [val_1, val_2])

        if isinstance(amplitude, Quantity):
            return Quantity(result, unit=amplitude.unit, copy=False)
        return result

    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box``.

        ``((y_low, y_high), (x_low, x_high))``
        """

        dr = self.R_0 + self.amplitude / self.slope

        return ((self.y_0 - dr, self.y_0 + dr),
                (self.x_0 - dr, self.x_0 + dr))

    @property
    def input_units(self):
        if self.x_0.unit is None and self.y_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit,
                self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit['x'] != inputs_unit['y']:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[0]],
                'R_0': inputs_unit[self.inputs[0]],
                'slope': outputs_unit[self.outputs[0]] / inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


class RickerWavelet1D(Fittable1DModel):
    """
    One dimensional Ricker Wavelet model (sometimes known as a "Mexican Hat"
    model).

    .. note::

        See https://github.com/astropy/astropy/pull/9445 for discussions
        related to renaming of this model.

    Parameters
    ----------
    amplitude : float
        Amplitude
    x_0 : float
        Position of the peak
    sigma : float
        Width of the Ricker wavelet

    See Also
    --------
    RickerWavelet2D, Box1D, Gaussian1D, Trapezoid1D

    Notes
    -----
    Model formula:

    .. math::

        f(x) = {A \\left(1 - \\frac{\\left(x - x_{0}\\right)^{2}}{\\sigma^{2}}\\right)
        e^{- \\frac{\\left(x - x_{0}\\right)^{2}}{2 \\sigma^{2}}}}

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import RickerWavelet1D

        plt.figure()
        s1 = RickerWavelet1D()
        r = np.arange(-5, 5, .01)

        for factor in range(1, 4):
            s1.amplitude = factor
            s1.width = factor
            plt.plot(r, s1(r), color=str(0.25 * factor), lw=2)

        plt.axis([-5, 5, -2, 4])
        plt.show()
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    sigma = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, x_0, sigma):
        """One dimensional Ricker Wavelet model function"""

        xx_ww = (x - x_0) ** 2 / (2 * sigma ** 2)
        return amplitude * (1 - 2 * xx_ww) * np.exp(-xx_ww)

    def bounding_box(self, factor=10.0):
        """Tuple defining the default ``bounding_box`` limits,
        ``(x_low, x_high)``.

        Parameters
        ----------
        factor : float
            The multiple of sigma used to define the limits.

        """
        x0 = self.x_0
        dx = factor * self.sigma

        return (x0 - dx, x0 + dx)

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'x_0': inputs_unit[self.inputs[0]],
                'sigma': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


class RickerWavelet2D(Fittable2DModel):
    """
    Two dimensional Ricker Wavelet model (sometimes known as a "Mexican Hat"
    model).

    .. note::

        See https://github.com/astropy/astropy/pull/9445 for discussions
        related to renaming of this model.

    Parameters
    ----------
    amplitude : float
        Amplitude
    x_0 : float
        x position of the peak
    y_0 : float
        y position of the peak
    sigma : float
        Width of the Ricker wavelet

    See Also
    --------
    RickerWavelet1D, Gaussian2D

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
        """Two dimensional Ricker Wavelet model function"""

        rr_ww = ((x - x_0) ** 2 + (y - y_0) ** 2) / (2 * sigma ** 2)
        return amplitude * (1 - rr_ww) * np.exp(- rr_ww)

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit,
                self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[0]],
                'sigma': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


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
    .. [1] https://en.wikipedia.org/wiki/Airy_disk
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    radius = Parameter(default=1)
    _rz = None
    _j1 = None

    @classmethod
    def evaluate(cls, x, y, amplitude, x_0, y_0, radius):
        """Two dimensional Airy model function"""

        if cls._rz is None:
            try:
                from scipy.special import j1, jn_zeros
                cls._rz = jn_zeros(1, 1)[0] / np.pi
                cls._j1 = j1
            except ValueError:
                raise ImportError('AiryDisk2D model requires scipy.')

        r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2) / (radius / cls._rz)

        if isinstance(r, Quantity):
            # scipy function cannot handle Quantity, so turn into array.
            r = r.to_value(u.dimensionless_unscaled)

        # Since r can be zero, we have to take care to treat that case
        # separately so as not to raise a numpy warning
        z = np.ones(r.shape)
        rt = np.pi * r[r > 0]
        z[r > 0] = (2.0 * cls._j1(rt) / rt) ** 2

        if isinstance(amplitude, Quantity):
            # make z quantity too, otherwise in-place multiplication fails.
            z = Quantity(z, u.dimensionless_unscaled, copy=False)

        z *= amplitude
        return z

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit,
                self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[0]],
                'radius': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


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

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt

        from astropy.modeling.models import Moffat1D

        plt.figure()
        s1 = Moffat1D()
        r = np.arange(-5, 5, .01)

        for factor in range(1, 4):
            s1.amplitude = factor
            s1.width = factor
            plt.plot(r, s1(r), color=str(0.25 * factor), lw=2)

        plt.axis([-5, 5, -1, 4])
        plt.show()
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    gamma = Parameter(default=1)
    alpha = Parameter(default=1)

    @property
    def fwhm(self):
        """
        Moffat full width at half maximum.
        Derivation of the formula is available in
        `this notebook by Yoonsoo Bach <https://nbviewer.jupyter.org/github/ysbach/AO_2017/blob/master/04_Ground_Based_Concept.ipynb#1.2.-Moffat>`_.
        """
        return 2.0 * np.abs(self.gamma) * np.sqrt(2.0 ** (1.0 / self.alpha) - 1.0)

    @staticmethod
    def evaluate(x, amplitude, x_0, gamma, alpha):
        """One dimensional Moffat model function"""

        return amplitude * (1 + ((x - x_0) / gamma) ** 2) ** (-alpha)

    @staticmethod
    def fit_deriv(x, amplitude, x_0, gamma, alpha):
        """One dimensional Moffat model derivative with respect to parameters"""

        fac = (1 + (x - x_0) ** 2 / gamma ** 2)
        d_A = fac ** (-alpha)
        d_x_0 = (2 * amplitude * alpha * (x - x_0) * d_A / (fac * gamma ** 2))
        d_gamma = (2 * amplitude * alpha * (x - x_0) ** 2 * d_A /
                   (fac * gamma ** 3))
        d_alpha = -amplitude * d_A * np.log(fac)
        return [d_A, d_x_0, d_gamma, d_alpha]

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'x_0': inputs_unit[self.inputs[0]],
                'gamma': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


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

    @property
    def fwhm(self):
        """
        Moffat full width at half maximum.
        Derivation of the formula is available in
        `this notebook by Yoonsoo Bach <https://nbviewer.jupyter.org/github/ysbach/AO_2017/blob/master/04_Ground_Based_Concept.ipynb#1.2.-Moffat>`_.
        """
        return 2.0 * np.abs(self.gamma) * np.sqrt(2.0 ** (1.0 / self.alpha) - 1.0)

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
        d_x_0 = (2 * amplitude * alpha * d_A * (x - x_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_y_0 = (2 * amplitude * alpha * d_A * (y - y_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_alpha = -amplitude * d_A * np.log(1 + rr_gg)
        d_gamma = (2 * amplitude * alpha * d_A * rr_gg /
                   (gamma * (1 + rr_gg)))
        return [d_A, d_x_0, d_y_0, d_gamma, d_alpha]

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        else:
            return {self.inputs[0]: self.x_0.unit,
                    self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[0]],
                'gamma': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


class Sersic2D(Fittable2DModel):
    r"""
    Two dimensional Sersic surface brightness profile.

    Parameters
    ----------
    amplitude : float
        Surface brightness at r_eff.
    r_eff : float
        Effective (half-light) radius
    n : float
        Sersic Index.
    x_0 : float, optional
        x position of the center.
    y_0 : float, optional
        y position of the center.
    ellip : float, optional
        Ellipticity.
    theta : float, optional
        Rotation angle in radians, counterclockwise from
        the positive x-axis.

    See Also
    --------
    Gaussian2D, Moffat2D

    Notes
    -----
    Model formula:

    .. math::

        I(x,y) = I(r) = I_e\exp\left\{-b_n\left[\left(\frac{r}{r_{e}}\right)^{(1/n)}-1\right]\right\}

    The constant :math:`b_n` is defined such that :math:`r_e` contains half the total
    luminosity, and can be solved for numerically.

    .. math::

        \Gamma(2n) = 2\gamma (2n,b_n)

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        from astropy.modeling.models import Sersic2D
        import matplotlib.pyplot as plt

        x,y = np.meshgrid(np.arange(100), np.arange(100))

        mod = Sersic2D(amplitude = 1, r_eff = 25, n=4, x_0=50, y_0=50,
                       ellip=.5, theta=-1)
        img = mod(x, y)
        log_img = np.log10(img)


        plt.figure()
        plt.imshow(log_img, origin='lower', interpolation='nearest',
                   vmin=-1, vmax=2)
        plt.xlabel('x')
        plt.ylabel('y')
        cbar = plt.colorbar()
        cbar.set_label('Log Brightness', rotation=270, labelpad=25)
        cbar.set_ticks([-1, 0, 1, 2], update_ticks=True)
        plt.show()

    References
    ----------
    .. [1] http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
    """

    amplitude = Parameter(default=1)
    r_eff = Parameter(default=1)
    n = Parameter(default=4)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    ellip = Parameter(default=0)
    theta = Parameter(default=0)
    _gammaincinv = None

    @classmethod
    def evaluate(cls, x, y, amplitude, r_eff, n, x_0, y_0, ellip, theta):
        """Two dimensional Sersic profile function."""

        if cls._gammaincinv is None:
            try:
                from scipy.special import gammaincinv
                cls._gammaincinv = gammaincinv
            except ValueError:
                raise ImportError('Sersic2D model requires scipy.')

        bn = cls._gammaincinv(2. * n, 0.5)
        a, b = r_eff, (1 - ellip) * r_eff
        cos_theta, sin_theta = np.cos(theta), np.sin(theta)
        x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
        x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta
        z = np.sqrt((x_maj / a) ** 2 + (x_min / b) ** 2)

        return amplitude * np.exp(-bn * (z ** (1 / n) - 1))

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {self.inputs[0]: self.x_0.unit,
                self.inputs[1]: self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit[self.inputs[0]],
                'y_0': inputs_unit[self.inputs[0]],
                'r_eff': inputs_unit[self.inputs[0]],
                'theta': u.rad,
                'amplitude': outputs_unit[self.outputs[0]]}


class KingProjectedAnalytic1D(Fittable1DModel):
    """
    Projected (surface density) analytic King Model.


    Parameters
    ----------
    amplitude : float
        Amplitude or scaling factor.
    r_core : float
        Core radius (f(r_c) ~ 0.5 f_0)
    r_tide : float
        Tidal radius.


    Notes
    -----

    This model approximates a King model with an analytic function. The derivation of this
    equation can be found in King '62 (equation 14). This is just an approximation of the
    full model and the parameters derived from this model should be taken with caution.
    It usually works for models with a concentration (c = log10(r_t/r_c) paramter < 2.

    Model formula:

    .. math::

        f(x) = A r_c^2  \\left(\\frac{1}{\\sqrt{(x^2 + r_c^2)}} -
        \\frac{1}{\\sqrt{(r_t^2 + r_c^2)}}\\right)^2

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        from astropy.modeling.models import KingProjectedAnalytic1D
        import matplotlib.pyplot as plt

        plt.figure()
        rt_list = [1, 2, 5, 10, 20]
        for rt in rt_list:
            r = np.linspace(0.1, rt, 100)

            mod = KingProjectedAnalytic1D(amplitude = 1, r_core = 1., r_tide = rt)
            sig = mod(r)


            plt.loglog(r, sig/sig[0], label='c ~ {:0.2f}'.format(mod.concentration))

        plt.xlabel("r")
        plt.ylabel(r"$\\sigma/\\sigma_0$")
        plt.legend()
        plt.show()

    References
    ----------
    .. [1] http://articles.adsabs.harvard.edu/pdf/1962AJ.....67..471K
    """

    amplitude = Parameter(default=1, bounds=(FLOAT_EPSILON, None))
    r_core = Parameter(default=1, bounds=(FLOAT_EPSILON, None))
    r_tide = Parameter(default=2, bounds=(FLOAT_EPSILON, None))

    @property
    def concentration(self):
        """Concentration parameter of the king model"""
        return np.log10(np.abs(self.r_tide/self.r_core))

    @staticmethod
    def evaluate(x, amplitude, r_core, r_tide):
        """
        Analytic King model function.
        """

        result = amplitude * r_core ** 2 * (1/np.sqrt(x ** 2 + r_core ** 2) -
                                            1/np.sqrt(r_tide ** 2 + r_core ** 2)) ** 2

        # Set invalid r values to 0
        bounds = (x >= r_tide) | (x < 0)
        result[bounds] = result[bounds] * 0.

        return result

    @staticmethod
    def fit_deriv(x, amplitude, r_core, r_tide):
        """
        Analytic King model function derivatives.
        """
        d_amplitude = r_core ** 2 * (1/np.sqrt(x ** 2 + r_core ** 2) -
                                     1/np.sqrt(r_tide ** 2 + r_core ** 2)) ** 2

        d_r_core = 2 * amplitude * r_core ** 2 * (r_core/(r_core ** 2 + r_tide ** 2) ** (3/2) -
                                                  r_core/(r_core ** 2 + x ** 2) ** (3/2)) * \
                   (1./np.sqrt(r_core ** 2 + x ** 2) - 1./np.sqrt(r_core ** 2 + r_tide ** 2)) + \
                   2 * amplitude * r_core * (1./np.sqrt(r_core ** 2 + x ** 2) -
                                             1./np.sqrt(r_core ** 2 + r_tide ** 2)) ** 2

        d_r_tide = (2 * amplitude * r_core ** 2 * r_tide *
                    (1./np.sqrt(r_core ** 2 + x ** 2) -
                     1./np.sqrt(r_core ** 2 + r_tide ** 2)))/(r_core ** 2 + r_tide ** 2) ** (3/2)

        # Set invalid r values to 0
        bounds = (x >= r_tide) | (x < 0)
        d_amplitude[bounds] = d_amplitude[bounds]*0
        d_r_core[bounds] = d_r_core[bounds]*0
        d_r_tide[bounds] = d_r_tide[bounds]*0

        return [d_amplitude, d_r_core, d_r_tide]

    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box`` limits.

        The model is not defined for r > r_tide.

        ``(r_low, r_high)``
        """

        return (0 * self.r_tide, 1 * self.r_tide)

    @property
    def input_units(self):
        if self.r_core.unit is None:
            return None
        return {self.inputs[0]: self.r_core.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'r_core': inputs_unit[self.inputs[0]],
                'r_tide': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


class Logarithmic1D(Fittable1DModel):
    """
    One dimensional logarithmic model.

    Parameters
    ----------
    amplitude : float, optional
    tau : float, optional

    See Also
    --------
    Exponential1D, Gaussian1D
    """

    amplitude = Parameter(default=1)
    tau = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, tau):
        return amplitude * np.log(x / tau)

    @staticmethod
    def fit_deriv(x, amplitude, tau):
        d_amplitude = np.log(x / tau)
        d_tau = np.zeros(x.shape) - (amplitude / tau)
        return [d_amplitude, d_tau]

    @property
    def inverse(self):
        new_amplitude = self.tau
        new_tau = self.amplitude
        return Exponential1D(amplitude=new_amplitude, tau=new_tau)

    @tau.validator
    def tau(self, val):
        if val == 0:
            raise ValueError("0 is not an allowed value for tau")

    @property
    def input_units(self):
        if self.tau.unit is None:
            return None
        return {self.inputs[0]: self.tau.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'tau': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


class Exponential1D(Fittable1DModel):
    """
    One dimensional exponential model.

    Parameters
    ----------
    amplitude : float, optional
    tau : float, optional

    See Also
    --------
    Logarithmic1D, Gaussian1D
    """
    amplitude = Parameter(default=1)
    tau = Parameter(default=1)

    @staticmethod
    def evaluate(x, amplitude, tau):
        return amplitude * np.exp(x / tau)

    @staticmethod
    def fit_deriv(x, amplitude, tau):
        ''' Derivative with respect to parameters'''
        d_amplitude = np.exp(x / tau)
        d_tau = -amplitude * (x / tau**2) * np.exp(x / tau)
        return [d_amplitude, d_tau]

    @property
    def inverse(self):
        new_amplitude = self.tau
        new_tau = self.amplitude
        return Logarithmic1D(amplitude=new_amplitude, tau=new_tau)

    @tau.validator
    def tau(self, val):
        ''' tau cannot be 0'''
        if val == 0:
            raise ValueError("0 is not an allowed value for tau")

    @property
    def input_units(self):
        if self.tau.unit is None:
            return None
        return {self.inputs[0]: self.tau.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'tau': inputs_unit[self.inputs[0]],
                'amplitude': outputs_unit[self.outputs[0]]}


@deprecated('4.0', alternative='RickerWavelet1D')
class MexicanHat1D(RickerWavelet1D):
    """ Deprecated."""


@deprecated('4.0', alternative='RickerWavelet2D')
class MexicanHat2D(RickerWavelet2D):
    """  Deprecated."""
