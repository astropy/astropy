.. _new:

********************
Creating a New Model
********************

This document describes how to add a model to the package. 
In short, one needs to define all model parameters and write an eval function
which evaluates the model. If the model is fittable, a function to compute the 
derivatives may be provided as well.

The details are explained below with a 1D gaussian model as an example.
There are two base classes for models. If the model is fittable, it 
should inherit from `models.ParametricModel`,
if not it should subclass `models.Model`. If the model takes parameters, 
their names are stored in a list as a class attribute named `parnames`.
Pass the list of parameter names and the number of parameter sets to the base 
class. Note, that if the method which evaluates the model cannot work
with multiple parameter sets, `paramdim` should not be given as an argument
in the __init__ method. The default for `paramdim` is set in the base class to 1.

::

    class Gauss1DModel(ParametricModel):
        parnames = ['amplitude', 'xcen', 'xsigma']


As a minimum the __init__ method takes all parameters and the number of parameter sets, `paramdim`:

::

    def __init__(self, amplitude, xcen, xsigma, paramdim=1):
        self.linear = False
        self._amplitude = parameters._Parameter(name='amplitude', val=amplitude, mclass=self, paramdim=paramdim)
        self._xsigma = parameters._Parameter(name='xsigma', val=xsigma, mclass=self, paramdim=paramdim)
        self._xcen = parameters._Parameter(name='xcen', val=xcen, mclass=self, paramdim=paramdim)
        ParametricModel.__init__(self, self.parnames, ndim=1, outdim=1, paramdim=paramdim)
    
Parametric models can be linear or nonlinear in a regression sense. The default 
value of the linear attribute is True. 
The `ndim` attribute stores the number of input coordinates.
The `outdim` attribute stores the number of output coordinates.
These two attributes are used with composite models.
Each parameter must be defined as a private attribute of the model class. 
Parameters are instances of `parameters._Parameter` class which takes as
arguments the name of the parameter, its value, the instance of the class 
and the number of parameter sets.

Next, provide a method, called `eval,  to evaluate the model and a method,
called `deriv`,  to compute its derivatives. The evaluation method takes all
input coordinates as separate arguments and a parameter set. For this example:

::

    def eval(self, x, params):
        return params[0] * np.exp((-(1/(params[2]**2)) * (x-params[1])**2))
                                                

The `deriv` method takes as input all coordinates as separate arguments.
There is an option to compute numerical derivatives for nonlinear models
in which case the `deriv` method should return None.

Finally, the __call__ method takes input coordinates as separate arguments.
It reformats them (if necessary) and calls the `eval` method to perform the 
model evaluation using model.psets as parameters. 
The reason there is a separate `eval` method is to allow fitters to call the `eval`
method with different parameters which is necessary for fitting with constraints.

::

    def __call__(self, x):
        x, format = _convert_input(x, self.paramdim)
        result = self.eval(x, self.psets)
        return _convert_output(result, format)
    
*********************
Creating a New Fitter
*********************

This document describes how to add a new nonlinear fitting algorithm
to this package. In short, one needs to define an error function and a __call__
method and define the types of constraints which work with this fitter (if any).

The details are described below using scipy's SLSQP algorithm as an example.
The base class for all fitters is `fitting.Fitter`. 

::

    class SLSQPFitter(Fitter):
        def __init__(self, model, fixed=None, tied=None, bounds=None,
                            eqcons=None, ineqcons=None):
            Fitter.__init__(self, model, fixed=fixed, tied=tied, bounds=bounds, 
                                      eqcons=eqcons, ineqcons=ineqcons)
            if self.model.linear:
                raise ModelLinearityException('Model is linear in parameters, '
                             'non-linear fitting methods should not be used.')

All fitters take a model (their __call__ method modifies the model's parameters).
If the fitter does not support constraint fitting, this may be the only argument 
passed to the constructor. In our example the rest of the arguments represent 
different types of constraints.

Next, the error function takes a list of parameters returned by an iteration of the 
fitting algorithm and input coordinates, evaluates the model with them and 
returns some type of a measure for the fit. In the example the sum of the 
squared residuals is used as a measure of fitting.

::

    def errorfunc(self, fps, *args):
        meas = args[0]
        self.fitpars = fps
        res = self.model(*args[1:]) - meas
        return np.sum(res**2)
    
The __call__ method performs the fitting. As a minimum it takes all coordinates 
as separate arguments. Additional arguments are passed as necessary.

::

    def __call__(self, x, y , maxiter=MAXITER, epsilon=EPS):
        self.fitpars = optimize.fmin_slsqp(self.errorfunc, p0=self.model.parameters[:], args=(y, x), 
            bounds=self.constraints._bounds, eqcons=self.constraints.eqcons, 
            ieqcons=self.constraints.ineqcons)