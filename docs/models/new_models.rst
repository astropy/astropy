****************************
Defining a New Type of Model
****************************

This document describes how to add a model to the package or to define a
user model. In short, one needs to define all model parameters and write
an eval function which evaluates the model. If the model is fittable,
a function to compute the derivatives is required if linear fitting
algorithm is to be used and optional if a non-liner fitter is to be used.

A Step by Step Definition of a 1D Gaussian Model
------------------------------------------------

The details are explained below with a 1D gaussian model as an example.
There are two base classes for models. If the model is fittable, it 
should inherit from `~astropy.models.core.ParametricModel`,
if not it should subclass `~astropy.models.core.Model`. If the model
takes parameters, their names are stored in a list as a class attribute
named `~astropy.models.core.Model.parnames`. Pass the list of parameter
names and the number of parameter sets to the base class. Note, that if
the method which evaluates the model cannot work with multiple parameter sets,
`~astropy.models.core.Model.paramdim` should not be given
as an argument in the __init__ method. The default for
`~astropy.models.core.Model.paramdim` is set in the base class to 1.::

    from astropy.models import *
    
    class Gauss1DModel(ParametricModel):
        parnames = ['amplitude', 'xcen', 'xsigma']


As a minimum the __init__ method takes all parameters and the number of
parameter sets, `~astropy.models.Model.paramdim`::

    def __init__(self, amplitude, xcen, xsigma, paramdim=1):
        self.linear = False
        self._amplitude = Parameter(name='amplitude', val=amplitude, mclass=self, paramdim=paramdim)
        self._xsigma = Parameter(name='xsigma', val=xsigma, mclass=self, paramdim=paramdim)
        self._xcen = Parameter(name='xcen', val=xcen, mclass=self, paramdim=paramdim)
        ParametricModel.__init__(self, self.parnames, ndim=1, outdim=1, paramdim=paramdim)
    
Parametric models can be linear or nonlinear in a regression sense. The default 
value of the `~astropy.models.core.Model.linear` attribute is True. 
The `~astropy.models.core.Model.ndim` attribute stores the number of input
variables the model expects.. The `~astropy.models.core.Model.outdim` attribute
stores the number of output variables returned after evaluating the model.
These two attributes are used with composite models.
Each parameter must be defined as a private attribute of the model class. 
Parameters are instances of `~astropy.models.parameters.Parameter` class which takes as
arguments the name of the parameter, its value, the instance of the class 
and the number of parameter sets.

Next, provide a method, called "eval" to evaluate the model and a method,
called "deriv",  to compute its derivatives. The evaluation method takes all
input coordinates as separate arguments and a parameter set. For this example::

    def eval(self, x, params):
        return params[0] * np.exp((-(1/(params[2]**2)) * (x-params[1])**2))
                                                

The "deriv" method takes as input all coordinates as separate arguments.
There is an option to compute numerical derivatives for nonlinear models
in which case the "deriv" method should return None.

Finally, the __call__ method takes input coordinates as separate arguments.
It reformats them (if necessary) and calls the eval method to perform the 
model evaluation using model.psets as parameters. 
The reason there is a separate eval method is to allow fitters to call the eval
method with different parameters which is necessary for fitting with constraints.::

    def __call__(self, x):
        x, format = _convert_input(x, self.paramdim)
        result = self.eval(x, self.psets)
        return _convert_output(result, format)
    
A Full Example of a line model
------------------------------

.. literalinclude:: new_model.py
   :language: python
   :linenos:
