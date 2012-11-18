.. _introduction:

************
Introduction
************

The `~fitting.models` and `~fitting.fitting` modules described here are designed to work as 
peers with each other. The goal is to be able to add models without 
explicit reference to fitting algorithms (though exceptions are 
sometimes necessary) and likewise, add different fitting algorithms 
without changing the existing models. The mechanism that allows this 
is the special `~fitting.parameters` module that both models and fitters use to 
interact with each other. Nevertheless, most users won't need to 
interact with this module unless they wish to add new models or 
fitters (the term used hereafter for specific fitting algorithms) to 
the existing suites of models and fitters.

Furthermore, the models are designed to be combined in many ways. It 
is possible, for example, to combine models serially `~fitting.models.SCompositeModel`, so that the 
output values of one model are used as input values to another. It is 
also possible to form a new model as a sum of two other models (in 
effect, a parallel combination of models), `~fitting.models.PCompositeModel`. Since models may have 
multiple input values, machinery is provided that allows assigning 
outputs from one model into the appropriate input of another in 
flexible way, `~fitting.models.LabeledInput`. Finally, it is permitted to combine any number of models 
using all of these mechanisms simultaneously. A composite model can be 
used to make further composite models. The goal is to eventually 
provide a rich toolset of models and fitters such that most users will 
not need to define new model classes, nor special purpose fitting 
routines (but not making that hard to do if it is necessary).

The following documentation will begin with a few examples 
illustrating some of these features, before going into more detail 
about each of the three important components of this system.

While this system is initially being developed to support WCS models, 
its generality extends well beyond WCS cases, and thus the initial 
examples are not specifically geared to WCS problems.

Examples
--------

All examples assume the following modules have been imported

>>> import numpy as np
>>> from fitting import models, fitting

- Working with 1D models

Create data using 1D Chebyshev model

>>> x = np.arange(1,10,.1)
>>> ch1 = models.ChebyshevModel(3)
>>> ch1.parameters
[0.0, 0.0, 0.0, 0.0]
>>> ch1.parameters = [1,2,3,4]
>>> ch1.parameters
[1.0, 2.0, 3.0, 4.0]
>>> print ch1
Model: ChebyshevModel
Dim:   1
Degree: 3
Parameter sets: 1
Parameters: 
           c0:  [1.0]
           c1:  [2.0]
           c2:  [3.0]
           c3:  [4.0]
>>> y = ch1(x)
>>> n=np.random.randn(90)
>>> ny = y + n 


Fit a Chebyshev polynomial to the data

>>> ch2 = models.ChebyshevModel(3)
>>> chfit = fitting.LinearLSQFitter(ch2)
>>> chfit(x,y)
>>> ch2.parameters
[0.7863153, 1.7473515, 2.8038203, 3.9106717]

.. figure:: images/cheb_fit.png
   :scale: 75 %

   Fit a Chebyshev polynomial of order 3 to some noisy data
   
Fix one of the parameters and fit again. This is done using a dictionary 
with {key: value} pairs where key is the name of the parameter to be 
held fixed and value is set to False. The fixed coefficient will not be
used in the fitting and a polynomial without this term will be fitted. 
It will be used when evaluating the polynomial, so if the entire term is 
to be excluded from the polynomial, the coefficient should be set to zero.

>>> pmask = {'c2': False}
>>> ch2.parnames
['c0', 'c1', 'c2', 'c3']
>>> ch2.c2 = 0
>>> chfit = fitting.LinearLSQFitter(ch2, pmask=pmask)
>>> ch2.parameters
[1.00000000, 2.00000000, 0.0, 4.0]

- Working with 2D models

First create some data to be fitted with a 2D polynomial

>>> x, y = np.mgrid[:10, :10]
>>> def poly2(x, y):
        return 1+2*x+3*x**2+4*y+5*y**2+6*x*y
z = poly2(x, y)

Fit a 2D polynomial to the data

>>> p2 = models.Poly2DModel(2, xdomain = [0,9], ydomain = [0,9], xwindow = [0,9], ywindow = [0,9])
>>> print p2
Model: Poly2DModel
Dim:   2
Degree: 2
Parameter sets: 1
Parameters: 
           c0_0:  [0.0]
           c1_0:  [0.0]
           c2_0:  [0.0]
           c0_1:  [0.0]
           c0_2:  [0.0]
           c1_1:  [0.0]
>>>pfit = fitting.LinearLSQFitter(p2)
>>>n = np.random.randn(200)
>>>pfit(x, y, z+n)
>>>fitter = fitting.LinearLSQFitter(self.model)
>>> p2.parameters
[1.1022004, 2.0094546, 2.9953080, 4.0000000, 4.9999999, 6.0000000]

.. figure:: images/poly2d_fit.png
   :scale: 75 %
   
   A 2D polynomial model fit to Z and some noise
   
- Create and evaluate a parallel composite model

>>> x = np.arange(1,10,.1)
>>> p1 = models.Poly1DModel(1)

Parameters have a dual interface, explained :doc:`later <parameters>`,
which allows them to be changed as a whole or using their names.

>>> p1.parameters
[0.0, 0.0
>>> p1.parameters = [3, 0.]
>>> p1.parameters
[3.0, 0.0]
>>> p1.parnames
['c0', 'c1']
>>> p1.c1 = 6.7
>>> p1.parameters
[3.0, 6.7]
>>> g1 = models.Gauss1DModel(10., xsigma=2.1, xcen=4.2)
>>> pcomptr = models.PCompositeModel([g1, p1])
>>> y = pcomptr(x)

This is equivalent to applying the two models in parallel:

>>> y = x + (g1(x) - x) + (p1(x) - x)

In more complex cases the input and output may be mapped:

>>> x, y = np.mgrid[:10, :10]
>>> off = models.ShiftModel(-3.2)
>>> poly2 = models.Poly2DModel(2)
>>> scomptr = models.SCompositeModel([off, poly2], inmap=[['x'], ['x', 'y']], outmap=[['x'], ['z']])

The above composite transform will apply an inplace shift to x, followed by a 2D 
polynomial and will save the result in an array, labeled 'z'.
To evaluate this model use a LabeledInput object

>>> ado = models.LabeledInput([x, y], ['x', 'y'])
>>> result = scomptr(ado)

The output is also a LabeledInput object and the result is stored in label 'z'.

>>> print result
{'x': array([[-3.2, -3.2, -3.2, -3.2, -3.2, -3.2, -3.2, -3.2, -3.2, -3.2],
       [-2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2],
       [-1.2, -1.2, -1.2, -1.2, -1.2, -1.2, -1.2, -1.2, -1.2, -1.2],
       [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2],
       [ 0.8,  0.8,  0.8,  0.8,  0.8,  0.8,  0.8,  0.8,  0.8,  0.8],
       [ 1.8,  1.8,  1.8,  1.8,  1.8,  1.8,  1.8,  1.8,  1.8,  1.8],
       [ 2.8,  2.8,  2.8,  2.8,  2.8,  2.8,  2.8,  2.8,  2.8,  2.8],
       [ 3.8,  3.8,  3.8,  3.8,  3.8,  3.8,  3.8,  3.8,  3.8,  3.8],
       [ 4.8,  4.8,  4.8,  4.8,  4.8,  4.8,  4.8,  4.8,  4.8,  4.8],
       [ 5.8,  5.8,  5.8,  5.8,  5.8,  5.8,  5.8,  5.8,  5.8,  5.8]]),
 'y': array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]),
 'z': array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])}



Fit a data set with a gaussian model.

>>> print g1
Model: Gauss1DModel
Dim:   1
Degree: N/A
Parameter sets: 1
Parameters: 
           amplitude:  [10.0]
           xcen:  [4.2000000000000002]
           xsigma:  [2.1000000000000001]
>>> y = g1(x)
>>> n = np.random.randn(90)
>>> ny = y + n
>>> gfit = fitting.NonLinearLSQFitter(g1)
>>> gfit(x, ny)
>>> print g1
Model: Gauss1DModel
Dim:   1
Degree: N/A
Parameter sets: 1
Parameters: 
           amplitude:  [10.141697966089579]
           xcen:  [4.2140429078454309]
           xsigma:  [2.0780002458907352]

* Freezing parameters in nonlinear models

Hold the sigma of the gaussian fixed when fitting.

>>> g1 = models.Gauss1DModel(10., xsigma=2.1, xcen=4.2)
>>> mask = {'xsigma': False}
>>> gfit = fitting.NonLinearLSQFitter(g1, pmask=mask)
>>> gfit(x, ny)
>>> print g1
Model: Gauss1DModel
Dim:   1
Degree: N/A
Parameter sets: 1
Parameters: 
           amplitude:  [10.088672698995229]
           xcen:  [4.2136244532915494]
           xsigma:  [2.1000000000000001]

- Fitting two models with shared parameters

Create 2 gaussian models and 2 data sets.
Create fitters for the two models keeping the sigma common for the two models.

>>> g1 = models.Gauss1DModel(10., xsigma=2.1, xcen=4.2
>>> g2 = models.Gauss1DModel(15., xsigma=2.8, xcen=2.2)
>>> y1 = g1(x)
>>> y2 = g2(x)
>>> n = np.random.randn(90)
>>> ny1 = y1 + n
>>> ny2 = y2 + n
>>> jf = fitting.JointFitter([g1, g2], jointparameters={g1:['xsigma'], g2:['xsigma']}, initvals=[2.5])
>>> jf(x, ny1, x, ny2)
>>> g1.parameters
[9.34640124, 4.20004213, 2.42646818]
g2.parameters
[15.44834757, 2.40199975, 2.42646818]
