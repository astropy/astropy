.. _fitting:

*******
Fitting
*******

This module provides wrappers, called Fitters, around some Numpy and Scipy 
fitting functions. All Fitters take an instance of
`models.ParametricModel` as input and define a `__call__` method
which fits the model to the data and changes `model.parameters
attribute. The idea is to make this extensible and allow users to easily add 
other fitters.

Linear fitting is done using Numpy's ~numpy.linalg.lstsq` function.
There are currently two non-linear fitters which use `~scipy.optimize.leastsq`
and `~scipy.optimize.slsqp`.

The rules for passing input to fitters are:

* Non-linear fitters work only with single data sets

* The linear fitter can fit single input to multiple data sets creating multiple 
  parameter sets. For example fitting a 2D model with input x, y arrays 
  of shape (n, m) to a z array of shape (p, n, m), will set 
  model.parameters.ndim to p, even if it was 1 when the model was created.

* Attempting to fit a model with multiple parameter sets to a single 
  data set results in an error.



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

Fitters support constraint fitting through `~models.constraints.Constraints`.

- All fitters support fixed (frozen) parameters through the **fixed** argument to models or setting the fixed attribute directly on a parameter.

For linear fitters freezing a polynomial coefficient means that a 
polynomial without that term will be fitted to the data. For example, fixing
c0 in a polynomial model will fit a polynomial with the zero-th order term missing.

>>> x=np.arange(1,10,.1)
>>> p1= models.Poly1DModel(2, paramdim=2)
>>> p1.parameters=[1,1,2,2,3,3]
>>> p1.psets
array([[ 1.,  1.],
       [ 2.,  2.],
       [ 3.,  3.]])
>>> y = p1(x)
>>> p1.c1.fixed = True
>>> pfit=fitting.LinearLSQFitter(p1)
>>> pfit(x, y)
>>> p1.psets
array([[ 5.50225913,  5.50225913],
       [ 2.        ,  2.        ],
       [ 3.17551299,  3.17551299]])

       
- Parameters can be tied. This can be done in two ways:

>>> def tiedfunc(g1):
    ...    xcen = 3*g1.xsigma[0]
    ...    return xcen
>>> g1 = models.Gauss1D(amplitude=10., xcen=3, xsigma=.5, tied={'xcen':tiedfunc})

or

>>> g1 = models.Gauss1D(amplitude=10., xcen=3, xsigma=.5)
>>> g1.xcen.tied = tiedfunc
>>> gfit = fitting.NonLinearLSQFitter(g1)


- Print a list of available fitting constraints

>>> fitting.Constraints.fitters
{'LinearLSQFitter': ['fixed'],
 'NonLinearLSQFitter': ['fixed', 'tied', 'bounds'],
 'SLSQPFitter': ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']}


    
