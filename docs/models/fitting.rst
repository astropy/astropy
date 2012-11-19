.. _fitting:

*******
Fitting
*******

This module provides wrappers, called Fitters, around some Numpy and Scipy 
fitting functions. All Fitters take an instance of
`models.ParametricModel` as input and define a __call__ method
which fits the model to the data and changes the model's parameters 
attribute. The idea is to make this extensible and allow users to easily add 
other fitters.

Linear fitting is done using Numpy's linalg.lstsq function.
There are currently two non-linear fitters which use leastsq and slsqp functions
in scipy.optimize.

The rules for passing input to fitters are:

* Non-linear fitters work only with single data sets

* The linear fitter can fit single input to multiple data sets creating multiple 
  parameter sets. For example fitting a 2D model with input x, y arrays 
  of shape (n, m) to a z array of shape (p, n, m), will set 
  model.parameters.ndim to p, even if it was 1 when the model was created.

* Attempting to fit a model with multiple parameter sets to a single 
  data set results in an error.

Fitters support constraint fitting through `fitting.Constraints`.

Fitting Examples
----------------

- Fitting simultaneously a polynomial model to multiple data sets


>>> p1 = models.Poly1DModel(3)
>>> p1.c0=1
>>> p1.c1=2
>>> p1.parameters
[1.0, 2.0, 0.0, 0.0]
>>> x=np.arange(10)
>>> y=p1(x)
>>> yy=np.array([y,y]).T
>>> p2=models.Poly1DModel(3, paramdim=2)
>>> pfit=fitting.LinearLSQFitter(p2)
>>> pfit(x,yy)
>>> print p2.psets
array([[  1.00000000e+00,   1.00000000e+00],
       [  2.00000000e+00,   2.00000000e+00],
       [  3.91115939e-16,   3.91115939e-16],
       [ -2.99676984e-17,  -2.99676984e-17]])

- All fitters support fixed parameters. 

For linear fitters fixing a polynomial coefficient means that a 
polynomial without that term will be fitted to the data. For example the 
fitter, LinearLSQFitter(p1, fixed=['c0']),  will fit a polynomial 
with the zero-th order term missing.

- Print a list of available fitting constraints

>>> fitting.Constraints.fitters
{'LinearLSQFitter': ['fixed'],
 'NonLinearLSQFitter': ['fixed', 'tied'],
 'SLSQPFitter': ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']}


    
