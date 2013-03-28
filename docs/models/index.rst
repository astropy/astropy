**************************************
Models  And Fitting (`astropy.models`)
**************************************

.. _designed: design.rst

Introduction
============
`~astropy.models` provides a framework for representing models and 
performing model evaluation and fitting. It supports 1D and 2D models 
parameter constraints.

It is `designed`_ to be easily extensible and flexible.
Models do not reference fitting algorithms explicitely
(though exceptions are sometimes necessary) and new fitting 
algorithms may be added without changing the existing models.
In addition models can be combined in different ways using a machinery 
that allows assigning outputs from one model into the appropriate input 
of another in a flexible way, `~astropy.models.models.LabeledInput`.
The goal is to eventually provide a rich toolset of models and fitters
such that most users will not need to define new model classes, nor
special purpose fitting routines (but not making that hard to do if it is necessary).

This is a work in progress but the main infrastructure is in place and usable 
in many ways now.

Getting Started
===============

All examples assume the following modules have been imported

>>> import numpy as np
>>> from astropy.models import models, fitting

Working with 1D models
======================

Fit a data set with a gaussian model.

>>> g1 = models.Gauss1DModel(10., xsigma=2.1, xcen=4.2)
>>> g1
<Gauss1DModel(amplitude= [10.0],xcen= [4.2000000000000002],xsigma= [2.1000000000000001],paramdim=1)>
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
         
Create data using 1D Chebyshev model

>>> x = np.arange(1,10,.1)
>>> ch1 = models.ChebyshevModel(3, domain=[x.min(), x.max()])
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

Add some noise

>>> n=np.random.randn(90)
>>> ny = y + n 

Fit a Chebyshev polynomial to the data

>>> ch2 = models.ChebyshevModel(3)
>>> chfit = fitting.LinearLSQFitter(ch2)
>>> chfit(x, ny)
>>> ch2.parameters
[0.957, 1.931, 3.029, 4.305]

.. figure:: images/cheb_fit.png
    :scale: 25 %

Working with 2D models
======================

First create some data to be fitted with a 2D polynomial

>>> x, y = np.mgrid[:10, :10]
>>> def poly2(x, y):
        return 1+2*x+3*x**2+4*y+5*y**2+6*x*y
>>> z = poly2(x, y)

Fit a 2D polynomial to the data

>>> p2 = models.Poly2DModel(2)
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
>>>n = np.random.randn(100)
>>>n.shape = (10, 10)
>>>pfit(x, y, z+n)
>>> p2.parameters
[0.6354845, 2.016544, 3.0035796, 4.0907439, 4.989999, 6.000127]



Using `models`
==============


.. toctree::
    parameters
    models
    fitting
    util
    new
    algorithms




Reference/API
=============

.. automodapi:: astropy.models.models
.. automodapi:: astropy.models.projections
.. automodapi:: astropy.models.rotations
.. automodapi:: astropy.models.fitting
.. automodapi:: astropy.models.parameters

