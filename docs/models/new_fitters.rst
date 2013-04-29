*****************************
Creating a New Type of Fitter
*****************************

This document describes how to add a new nonlinear fitting algorithm
to this package. In short, one needs to define an error function and a ``__call__``
method and define the types of constraints which work with this fitter (if any).

The details are described below using scipy's SLSQP algorithm as an example.
The base class for all fitters is `~astropy.models.fitting.Fitter`. ::

    class SLSQPFitter(Fitter):
        def __init__(self, model, fixed=None, tied=None, bounds=None,
                            eqcons=None, ineqcons=None):
            Fitter.__init__(self, model, fixed=fixed, tied=tied, bounds=bounds, 
                                      eqcons=eqcons, ineqcons=ineqcons)
            if self.model.linear:
                raise ModelLinearityException('Model is linear in parameters, '
                             'non-linear fitting methods should not be used.')

All fitters take a model (their ``__call__`` method modifies the model's parameters).
If the fitter does not support constraint fitting, this may be the only argument 
passed to the constructor. In our example the rest of the arguments represent 
different types of constraints.

Next, the error function takes a list of parameters returned by an iteration of the 
fitting algorithm and input coordinates, evaluates the model with them and 
returns some type of a measure for the fit. In the example the sum of the 
squared residuals is used as a measure of fitting. ::

    def errorfunc(self, fps, *args):
        meas = args[0]
        self.fitpars = fps
        res = self.model(*args[1:]) - meas
        return np.sum(res**2)
    
The ``__call__`` method performs the fitting. As a minimum it takes all coordinates 
as separate arguments. Additional arguments are passed as necessary.::

    def __call__(self, x, y , maxiter=MAXITER, epsilon=EPS):
        self.fitpars = optimize.fmin_slsqp(self.errorfunc, p0=self.model.parameters[:], args=(y, x), 
            bounds=self.constraints._bounds, eqcons=self.constraints.eqcons, 
            ieqcons=self.constraints.ineqcons)

A Full Example of a new fitter
------------------------------

.. literalinclude:: new_fitter.py
   :language: python
   :linenos:
