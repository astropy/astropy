Defining New Model Classes
==========================

This document describes how to add a model to the package or to create a
user-defined model. In short, one needs to define all model parameters and
write a function which evaluates the model, that is, computes the mathematical
function that implements the model.  If the model is fittable, a function to
compute the derivatives with respect to parameters is required if a linear
fitting algorithm is to be used and optional if a non-linear fitter is to be
used.


Custom 1-D models
-----------------

For 1-D models, the `~astropy.modeling.functional_models.custom_model_1d`
decorator is provided to make it very easy to define new models. The following
example demonstrates how to set up a model consisting of two Gaussians:

.. plot::
   :include-source:

    import numpy as np
    from astropy.modeling.models import custom_model_1d
    from astropy.modeling.fitting import LevMarLSQFitter

    # Define model
    @custom_model_1d
    def sum_of_gaussians(x, amplitude1=1., mean1=-1., sigma1=1.,
                            amplitude2=1., mean2=1., sigma2=1.):
        return (amplitude1 * np.exp(-0.5 * ((x - mean1) / sigma1)**2) +
                amplitude2 * np.exp(-0.5 * ((x - mean2) / sigma2)**2))

    # Generate fake data
    np.random.seed(0)
    x = np.linspace(-5., 5., 200)
    m_ref = sum_of_gaussians(amplitude1=2., mean1=-0.5, sigma1=0.4,
                             amplitude2=0.5, mean2=2., sigma2=1.0)
    y = m_ref(x) + np.random.normal(0., 0.1, x.shape)

    # Fit model to data
    m_init = sum_of_gaussians()
    fit = LevMarLSQFitter()
    m = fit(m_init, x, y)

    # Plot the data and the best fit
    plt.plot(x, y, 'o', color='k')
    plt.plot(x, m(x), color='r', lw=2)

.. note::

    Currently this shortcut for model definition only works for 1-D models, but
    it is being expanded to support 2 or greater dimension models.


A step by step definition of a 1-D Gaussian model
-------------------------------------------------

The example described in `Custom 1-D models`_ can be used for most 1-D cases,
but the following section described how to construct model classes in general.
The details are explained below with a 1-D Gaussian model as an example.  There
are two base classes for models. If the model is fittable, it should inherit
from `~astropy.modeling.FittableModel`; if not it should subclass
`~astropy.modeling.Model`.

If the model takes parameters they should be specified as class attributes in
the model's class definition using the `~astropy.modeling.Parameter`
descriptor.  All arguments to the Parameter constructor are optional, and may
include a default value for that parameter, a text description of the parameter
(useful for `help` and documentation generation), as well default constraints
and custom getters/setters for the parameter value.

If the first argument ``name`` is specified it must be identical to the class
attribute being assigned that Parameter.  As such, Parameters take their name
from this attribute by default.  In other words, ``amplitude =
Parameter('amplitude')`` is equivalent to ``amplitude = Parameter()``.  This
differs from Astropy v0.3.x, where it was necessary to provide the name twice.

::

    from astropy.modeling import FittableModel, Parameter

    class Gaussian1D(FittableModel):
        amplitude = Parameter()
        mean = Parameter()
        stddev = Parameter()

At a minimum, the ``__init__`` method takes all parameters and a few keyword
arguments such as values for constraints::

    def __init__(self, amplitude, mean, stddev, **kwargs):
        # Note that this __init__ does nothing different from the base class's
        # __init__.  The main point of defining it is so that the function
        # signature is more informative.
        super(Gaussian1D, self).__init__(
            amplitude=amplitude, mean=mean, stddev=stddev, **kwargs)

.. note::

    If a parameter is defined with a default value you may make the argument
    for that parameter in the ``__init__`` optional.  Otherwise it is
    recommended to make it a required argument.  In the above example none of
    the parameters have default values.

Fittable models can be linear or nonlinear in a regression sense. The default
value of the `~astropy.modeling.Model.linear` attribute is ``False``.  Linear
models should define the ``linear`` class attribute as ``True``.
The `~astropy.modeling.Model.n_inputs` attribute stores the number of
input variables the model expects.  The
`~astropy.modeling.Model.n_outputs` attribute stores the number of output
variables returned after evaluating the model.  These two attributes are used
with composite models.

Next, provide methods called ``evaluate`` to evaluate the model and
``fit_deriv``, to compute its derivatives with respect to parameters.  These
may be normal methods, `classmethod`, or `staticmethod`, though the convention
is to use `staticmethod` when the function does not depend on any of the
object's other attributes (i.e., it does not reference ``self``).  The
evaluation method takes all input coordinates as separate arguments and all of
the model's parameters in the same order they would be listed by
`~astropy.modeling.Model.param_names`.

For this example::

    @staticmethod
    def evaluate(x, amplitude, mean, stddev):
        return amplitude * np.exp((-(1 / (2. * stddev**2)) * (x - mean)**2))

The ``fit_deriv`` method takes as input all coordinates as separate arguments.
There is an option to compute numerical derivatives for nonlinear models in
which case the ``fit_deriv`` method should be ``None``::

    @staticmethod
    def fit_deriv(x, amplitude, mean, stddev):
        d_amplitude = np.exp((-(1 / (stddev**2)) * (x - mean)**2))
        d_mean = (2 * amplitude *
                  np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
                  (x - mean) / (stddev**2))
        d_stddev = (2 * amplitude *
                    np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
                    ((x - mean)**2) / (stddev**3))
        return [d_amplitude, d_mean, d_stddev]


Finally, the ``__call__`` method takes the input coordinates as arguments,
checks that the input arrays are correctly formatted and broadcastable, and
then passes them along with the the model's parameter values to the
``evaluate`` method.

In most cases it is not necessary to reimplement ``__call__``, though it may be
desirable to override it so that its signature reflects the desired names of
the input coordinates.  When overriding ``__call__`` it should simply call the
super-class's ``__call__`` and not do anything else.  The only reason currently
to override it is to provide a specific argument signature.  It should also
support ``model_set_axis`` keyword argument, at a minimum::

    def __call__(self, x, model_set_axis=None):
        return super(Gaussian1D, self).__call__(
            x, model_set_axis=model_set_axis)


Full example
^^^^^^^^^^^^

.. code-block:: python

    from astropy.modeling import FittableModel, Parameter

    class Gaussian1D(FittableModel):
        amplitude = Parameter()
        mean = Parameter()
        stddev = Parameter()

        def __init__(self, amplitude, mean, stddev, **kwargs):
            # Note that this __init__ does nothing different from the base
            # class's __init__.  The main point of defining it is so that the
            # function signature is more informative.
            super(Gaussian1D, self).__init__(
                amplitude=amplitude, mean=mean, stddev=stddev, **kwargs)

        @staticmethod
        def evaluate(x, amplitude, mean, stddev):
            return amplitude * np.exp((-(1 / (2. * stddev**2)) * (x - mean)**2))

        @staticmethod
        def fit_deriv(x, amplitude, mean, stddev):
            d_amplitude = np.exp((-(1 / (stddev**2)) * (x - mean)**2))
            d_mean = (2 * amplitude *
                      np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
                      (x - mean) / (stddev**2))
            d_stddev = (2 * amplitude *
                        np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
                        ((x - mean)**2) / (stddev**3))
            return [d_amplitude, d_mean, d_stddev]

        def __call__(self, x, model_set_axis=None):
            return super(Gaussian1D, self).__call__(
                x, model_set_axis=model_set_axis)


A full example of a LineModel
-----------------------------

.. code-block:: python

    from astropy.modeling import FittableModel, Parameter
    import numpy as np

    class LineModel(FittableModel):
        slope = Parameter()
        intercept = Parameter()
        linear = True

        def __init__(self, slope, intercept, **kwargs):
            super(LineModel, self).__init__(slope=slope, intercept=intercept,
                                            **kwargs)

        @staticmethod
        def evaluate(x, slope, intercept):
            return slope * x + intercept

        @staticmethod
        def fit_deriv(x, slope, intercept):
            d_slope = x
            d_intercept = np.ones_like(x)
            return [d_slope, d_intercept]

        def __call__(self, x, model_set_axis=None):
            return super(LineModel, self).__call__(
                x, model_set_axis=model_set_axis)


Defining New Fitter Classes
===========================

This section describes how to add a new nonlinear fitting algorithm to this
package or write a user-defined fitter.  In short, one needs to define an error
function and a ``__call__`` method and define the types of constraints which
work with this fitter (if any).

The details are described below using scipy's SLSQP algorithm as an example.
The base class for all fitters is `~astropy.modeling.fitting.Fitter`::

    class SLSQPFitter(Fitter):
        supported_constraints = ['bounds', 'eqcons', 'ineqcons', 'fixed',
                                 'tied']

        def __init__(self):
            # Most currently defined fitters take no arguments in their
            # __init__, but the option certainly exists for custom fitters
            super(SLSQPFitter, self).__init__()

All fitters take a model (their ``__call__`` method modifies the model's
parameters) as their first argument.

Next, the error function takes a list of parameters returned by an iteration of
the fitting algorithm and input coordinates, evaluates the model with them and
returns some type of a measure for the fit.  In the example the sum of the
squared residuals is used as a measure of fitting.::

    def objective_function(self, fps, *args):
        model = args[0]
        meas = args[-1]
        model.fitparams(fps)
        res = self.model(*args[1:-1]) - meas
        return np.sum(res**2)

The ``__call__`` method performs the fitting. As a minimum it takes all
coordinates as separate arguments. Additional arguments are passed as
necessary.::

    def __call__(self, model, x, y , maxiter=MAXITER, epsilon=EPS):
        if model.linear:
                raise ModelLinearityException(
                    'Model is linear in parameters; '
                    'non-linear fitting methods should not be used.')
        model_copy = model.copy()
        init_values, _ = _model_to_fit_params(model_copy)
        self.fitparams = optimize.fmin_slsqp(self.errorfunc, p0=init_values,
                                             args=(y, x),
                                             bounds=self.bounds,
                                             eqcons=self.eqcons,
                                             ineqcons=self.ineqcons)
        return model_copy


Using a Custom Statistic Function
=================================

This section describes how to write a new fitter with a user-defined statistic
function.  The example below shows a specialized class which fits a straight
line with uncertainties in both variables.

The following import statements are needed.::

    import numpy as np
    from astropy.modeling.fitting import (_validate_model,
                                          _fitter_to_model_params,
                                          _model_to_fit_params, Fitter,
                                          _convert_input)
    from astropy.modeling.optimizers import Simplex

First one needs to define a statistic. This can be a function or a callable
class.::

    def chi_line(measured_vals, updated_model, x_sigma, y_sigma, x):
        """
        Chi^2 statistic for fitting a straight line with uncertainties in x and
        y.

        Parameters
        ----------
        measured_vals : array
        updated_model : `~astropy.modeling.ParametricModel`
            model with parameters set by the current iteration of the optimizer
        x_sigma : array
            uncertainties in x
        y_sigma : array
            uncertainties in y

        """
        model_vals = updated_model(x)
        if x_sigma is None and y_sigma is None:
            return np.sum((model_vals - measured_vals) ** 2)
        elif x_sigma is not None and y_sigma is not None:
            weights = 1 / (y_sigma ** 2 + updated_model.parameters[1] ** 2 *
                           x_sigma ** 2)
            return np.sum((weights * (model_vals - measured_vals)) ** 2)
        else:
            if x_sigma is not None:
                weights = 1 / x_sigma ** 2
            else:
                weights = 1 / y_sigma ** 2
            return np.sum((weights * (model_vals - measured_vals)) ** 2)

In general, to define a new fitter, all one needs to do is provide a statistic
function and an optimizer. In this example we will let the optimizer be an
optional argument to the fitter and will set the statistic to ``chi_line``
above.::

    class LineFitter(Fitter):
        """
        Fit a straight line with uncertainties in both variables

        Parameters
        ----------
        optimizer : class or callable
            one of the classes in optimizers.py (default: Simplex)
        """

        def __init__(self, optimizer=Simplex):
            self.statistic = chi_line
            super(LineFitter, self).__init__(optimizer,
                                             statistic=self.statistic)

The last thing to define is the ``__call__`` method.::

    def __call__(self, model, x, y, x_sigma=None, y_sigma=None, **kwargs):
        """
        Fit data to this model.

        Parameters
        ----------
        model : `~astropy.modeling.core.ParametricModel`
            model to fit to x, y
        x : array
            input coordinates
        y : array
            input coordinates
        x_sigma : array
            uncertainties in x
        y_sigma : array
            uncertainties in y
        kwargs : dict
            optional keyword arguments to be passed to the optimizer

        Returns
        ------
        model_copy : `~astropy.modeling.core.ParametricModel`
            a copy of the input model with parameters set by the fitter

        """
        model_copy = _validate_model(model,
                                     self._opt_method.supported_constraints)

        farg = _convert_input(x, y)
        farg = (model_copy, x_sigma, y_sigma) + farg
        p0, _ = _model_to_fit_params(model_copy)

        fitparams, self.fit_info = self._opt_method(
            self.objective_function, p0, farg, **kwargs)
        _fitter_to_model_params(model_copy, fitparams)

        return model_copy
