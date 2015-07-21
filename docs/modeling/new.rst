.. _modeling-new-classes:

Defining New Model Classes
==========================

This document describes how to add a model to the package or to create a
user-defined model. In short, one needs to define all model parameters and
write a function which evaluates the model, that is, computes the mathematical
function that implements the model.  If the model is fittable, a function to
compute the derivatives with respect to parameters is required if a linear
fitting algorithm is to be used and optional if a non-linear fitter is to be
used.


Basic custom models
-------------------

For most cases, the `~astropy.modeling.custom_model` decorator provides an
easy way to make a new `~astropy.modeling.Model` class from an existing Python
callable. The following example demonstrates how to set up a model consisting
of two Gaussians:

.. plot::
   :include-source:

    import numpy as np
    from astropy.modeling.models import custom_model
    from astropy.modeling.fitting import LevMarLSQFitter

    # Define model
    @custom_model
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


This decorator also supports setting a model's
`~astropy.modeling.FittableModel.fit_deriv` as well as creating models with
more than one inputs.  It can also be used as a normal factory function (for
example ``SumOfGaussians = custom_model(sum_of_gaussians)``) rather than as a
decorator.  See the `~astropy.modeling.custom_model` documentation for more
examples.


A step by step definition of a 1-D Gaussian model
-------------------------------------------------

The example described in `Basic custom models`_ can be used for most simple
cases, but the following section describes how to construct model classes in
general.  Defining a full model class may be desirable, for example, to
provide more specialized parameters, or to implement special functionality not
supported by the basic `~astropy.modeling.custom_model` factory function.

The details are explained below with a 1-D Gaussian model as an example.  There
are two base classes for models. If the model is fittable, it should inherit
from `~astropy.modeling.FittableModel`; if not it should subclass
`~astropy.modeling.Model`.

If the model takes parameters they should be specified as class attributes in
the model's class definition using the `~astropy.modeling.parameters.Parameter`
descriptor.  All arguments to the Parameter constructor are optional, and may
include a default value for that parameter, a text description of the parameter
(useful for `help` and documentation generation), as well default constraints
and custom getters/setters for the parameter value.  It is also possible to
define a "validator" method for each parameter, enabling custom code to check
whether that parameter's value is valid according to the model definition (for
example if it must be non-negative).  See the example in
`Parameter.validator <astropy.modeling.Parameter.validator>` for more details.

::

    from astropy.modeling import FittableModel, Parameter

    class Gaussian1D(FittableModel):
        inputs = ('x',)
        outputs = ('y',)

        amplitude = Parameter()
        mean = Parameter()
        stddev = Parameter()

The ``inputs`` and ``outputs`` class attributes must be tuples of strings
indicating the number of independent variables that are input to evaluate the
model, and the number of outputs it returns.  The labels of the inputs and
outputs (in this case ``'x'`` and ``'y'`` respectively) are currently used for
informational purposes only and have no requirements on them other than that
they do not conflict with parameter names.  Outputs may have the same labels as
inputs (eg. ``inputs = ('x', 'y')`` and ``outputs = ('x', 'y')``).  However,
inputs must not conflict with each other (eg. ``inputs = ('x', 'x')`` is
incorrect) and likewise for outputs.  The lengths of these tuples are
important for specifying the correct number of inputs and outputs.  These
attributes supersede the ``n_inputs`` and ``n_outputs`` attributes in older
versions of this package.

There are two helpful base classes in the modeling package that can be used to
avoid specifying ``inputs`` and ``outputs`` for most common models.  These are
`~astropy.modeling.Fittable1DModel` and `~astropy.modeling.Fittable2DModel`.
For example, the real `~astropy.modeling.functional_models.Gaussian1D` model is
actually a subclass of `~astropy.modeling.Fittable1DModel`.  This helps cut
down on boilerplate by not having to specify ``inputs`` and ``outputs`` for
many models (follow the link to Gaussian1D to see its source code, for
example).

Fittable models can be linear or nonlinear in a regression sense. The default
value of the `~astropy.modeling.Model.linear` attribute is ``False``.  Linear
models should define the ``linear`` class attribute as ``True``.  Because this
model is non-linear we can stick with the default.

Next, provide methods called ``evaluate`` to evaluate the model and
``fit_deriv``, to compute its derivatives with respect to parameters.  These
may be normal methods, `classmethod`, or `staticmethod`, though the convention
is to use `staticmethod` when the function does not depend on any of the
object's other attributes (i.e., it does not reference ``self``) or any of the
class's other attributes as in the case of `classmethod`.  The evaluation
method takes all input coordinates as separate arguments and all of the model's
parameters in the same order they would be listed by
`~astropy.modeling.Model.param_names`.

For this example::

    @staticmethod
    def evaluate(x, amplitude, mean, stddev):
        return amplitude * np.exp((-(1 / (2. * stddev**2)) * (x - mean)**2))

It should be made clear that the ``evaluate`` method must be designed to take
the model's parameter values as arguments.  This may seem at odds with the fact
that the parameter values are already available via attribute of the model
(eg. ``model.amplitude``).  However, passing the parameter values directly to
``evaluate`` is a more efficient way to use it in many cases, such as fitting.

Users of your model would not generally use ``evaluate`` directly.  Instead
they create an instance of the model and call it on some input.  The
``__call__`` method of models uses ``evaluate`` internally, but users do not
need to be aware of it.  The default ``__call__`` implementation also handles
details such as checking that the inputs are correctly formatted and follow
Numpy's broadcasting rules before attempting to evaluate the model.

Like ``evaluate``, the ``fit_deriv`` method takes as input all coordinates and
all parameter values as arguments.  There is an option to compute numerical
derivatives for nonlinear models in which case the ``fit_deriv`` method should
be ``None``::

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


Note that we did *not* have to define an ``__init__`` method or a ``__call__``
method for our model (this contrasts with Astropy versions 0.4.x and earlier).
For most models the ``__init__`` follows the same pattern, taking the parameter
values as positional arguments, followed by several optional keyword arguments
(constraints, etc.).  The modeling framework automatically generates an
``__init__`` for your class that has the correct calling signature (see for
yourself by calling ``help(Gaussian1D.__init__)`` on the example model we just
defined).

There are cases where it might be desirable to define a custom ``__init__``.
For example, the `~astropy.modeling.functional_models.Gaussian2D` model takes
an optional ``cov_matrix`` argument which can be used as an alternative way to
specify the x/y_stddev and theta parameters.  This is perfectly valid so long
as the ``__init__`` determines appropriate values for the actual parameters and
then calls the super ``__init__`` with the standard arguments.  Schematically
this looks something like:

.. code-block:: python

    def __init__(self, amplitude, x_mean, y_mean, x_stddev=None,
                 y_stddev=None, theta=None, cov_matrix=None, **kwargs):
        # The **kwargs here should be understood as other keyword arguments
        # accepted by the basic Model.__init__ (such as constraints)
        if cov_matrix is not None:
            # Set x/y_stddev and theta from the covariance matrix
            x_stddev = ...
            y_stddev = ...
            theta = ...

        # Don't pass on cov_matrix since it doesn't mean anything to the base
        # class
        super(Gaussian2D, self).__init__(amplitude, x_mean, y_mean, x_stddev,
                                         y_stddev, theta, **kwargs)


Full example
^^^^^^^^^^^^

.. code-block:: python

    from astropy.modeling import FittableModel, Parameter

    class Gaussian1D(FittableModel):
        amplitude = Parameter()
        mean = Parameter()
        stddev = Parameter()

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


A full example of a LineModel
-----------------------------

This example demonstrates one other optional feature for model classes, which
is an *inverse*.  An `~astropy.modeling.Model.inverse` implementation should be
a `property` that returns a new model instance (not necessarily of the same
class as the model being inverted) that computes the inverse of that model, so
that for some model instance with an inverse, ``model.inverse(model(*input)) ==
input``.

.. code-block:: python

    from astropy.modeling import FittableModel, Parameter
    import numpy as np

    class LineModel(FittableModel):
        slope = Parameter()
        intercept = Parameter()
        linear = True

        @staticmethod
        def evaluate(x, slope, intercept):
            return slope * x + intercept

        @staticmethod
        def fit_deriv(x, slope, intercept):
            d_slope = x
            d_intercept = np.ones_like(x)
            return [d_slope, d_intercept]

        @property
        def inverse(self):
            new_slope = self.slope ** -1
            new_intercept = -self.intercept / self.slope
            return LineModel(slope=new_slope, intercept=new_intercept)


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
necessary::

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

The following import statements are needed::

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
above::

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

The last thing to define is the ``__call__`` method::

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
