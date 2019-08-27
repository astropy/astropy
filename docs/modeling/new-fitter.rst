.. _new_fitter:

Defining New Fitter Classes
***************************

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
            super().__init__()

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

Defining a Plugin Fitter
========================

`astropy.modeling` includes a plugin mechanism which allows fitters
defined outside of astropy's core to be inserted into the
`astropy.modeling.fitting` namespace through the use of entry points.
Entry points are references to importable objects. A tutorial on defining
entry points can be found in `setuptools' documentation <https://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins>`_.
Plugin fitters must to extend from the `~astropy.modeling.fitting.Fitter`
base class. For the fitter to be discovered and inserted into
`astropy.modeling.fitting` the entry points must be inserted into
the `astropy.modeling` entry point group

.. doctest-skip::

    setup(
          # ...
          entry_points = {'astropy.modeling': 'PluginFitterName = fitter_module:PlugFitterClass'}
    )

This would allow users to import the ``PlugFitterName`` through `astropy.modeling.fitting` by

.. doctest-skip::

    from astropy.modeling.fitting import PlugFitterName

One project which uses this functionality is `Saba <https://saba.readthedocs.io/>`_
and be can be used as a reference.

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
            super().__init__(optimizer, statistic=self.statistic)

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
