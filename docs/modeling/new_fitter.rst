***************************
Defining New Fitter Classes
***************************

This section describes how to add a new nonlinear fitting algorithm to this
package or write a user-defined fitter.  In short, one needs to define an error
function and a ``__call__`` method and define the types of constraints which
work with this fitter (if any).

The details are described below using scipy's SLSQP algorithm as an example.
The base class for all fitters is `~astropy.modeling.fitting.Fitter`::

    class SLSQPFitter(Fitter):
        supported_constraints = ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']

        def __init__(self):
            super(SLSQPFitter,self).__init__()

All fitters take a model (their ``__call__`` method modifies the model's
parameters) as their first argument.

Next, the error function takes a list of parameters returned by an iteration of
the fitting algorithm and input coordinates, evaluates the model with them and
returns some type of a measure for the fit.  In the example the sum of the
squared residuals is used as a measure of fitting.::

    def errorfunc(self, fps, *args):
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
        init_values, _ = model_copy._model_to_fit_params()
        self.fitparams = optimize.fmin_slsqp(self.errorfunc, p0=init_values,
                                             args=(y, x),
                                             bounds=self.bounds,
                                             eqcons=self.eqcons,
                                             ineqcons=self.ineqcons)
        return model_copy
