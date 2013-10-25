.. _astropy-modeling:

***************************************
Models And Fitting (`astropy.modeling`)
***************************************

Introduction
============
`~astropy.modeling` provides a framework for representing models and
performing model evaluation and fitting. It supports 1D and 2D models
and fitting with parameter constraints.

It is :ref:`designed <modeling-design>` to be easily extensible and flexible.
Models do not reference fitting algorithms explicitly (though exceptions are
sometimes necessary) and new fitting algorithms may be added without changing
the existing models.  In addition models can be combined in different ways
using a machinery that allows assigning outputs from one model into the
appropriate input of another in a flexible way,
`~astropy.modeling.core.LabeledInput`.  The goal is to eventually provide a
rich toolset of models and fitters such that most users will not need to define
new model classes, nor special purpose fitting routines (but not making that
hard to do if it is necessary).

.. warning::
    `~astropy.modeling` is currently a work-in-progress, and thus it is
    likely there will be significant API changes in later versions of
    Astropy.


Getting started
===============

The examples here use the predefined models and assume the following modules
have been imported::

    >>> import numpy as np
    >>> from astropy.modeling import models, fitting


Working with 1D models
======================

Fit a data set with a Gaussian model::

    >>> x = np.arange(1, 10, .1)
    >>> g1 = models.Gaussian1DModel(10., stddev=2.1, mean=4.2)
    >>> g1
    <Gaussian1DModel(amplitude=Parameter('amplitude', value=10.0), mean=Parameter('mean', value=4.2000000000000002), stddev=Parameter('stddev', value=2.1000000000000001), param_dim=1)>
    >>> y = g1(x)
    >>> n = np.random.randn(90)
    >>> ny = y + n
    >>> gfit = fitting.NonLinearLSQFitter(g1)
    >>> gfit(x, ny)
    >>> print(g1)
    Model: Gaussian1DModel
    n_inputs:   1
    Degree: N/A
    Parameter sets: 1
    Parameters:
               amplitude: Parameter('amplitude', value=10.315151071186117)
               mean: Parameter('mean', value=4.1780718585869554)
               stddev:  Parameter('stddev', value=2.0585618180947494)

Create data using a 1D Chebyshev model::

    >>> ch1 = models.Chebyshev1DModel(3, domain=[x.min(), x.max()])
    >>> ch1.parameters
    array([0., 0., 0., 0.])
    >>> ch1.parameters = [1, 2, 3, 4]
    >>> ch1.parameters
    array([1., 2., 3., 4.])
    >>> print(ch1)
    Model: Chebyshev1DModel
    n_inputs:   1
    Degree: 3
    Parameter sets: 1
    Parameters:
               c0: Parameter('c0', value=1.0)
               c1: Parameter('c1', value=2.0)
               c2: Parameter('c2', value=3.0)
               c3: Parameter('c3', value=4.0)
    >>> y = ch1(x)

Add some noise::

    >>> n = np.random.randn(90)
    >>> ny = y + n

Fit a Chebyshev polynomial to the data::

    >>> ch2 = models.Chebyshev1DModel(3)
    >>> chfit = fitting.LinearLSQFitter(ch2)
    >>> chfit(x, ny)
    >>> ch2.parameters
    array([ 1.08612543,  1.79746444,  3.15233293,  4.06529137])

.. plot::

   import matplotlib.pyplot as plt
   import numpy as np
   from astropy.modeling import models, fitting
   x = np.arange(1, 10, .1)
   ch1 = models.Chebyshev1DModel(3, domain=[x.min(), x.max()])
   ch1.parameters = [1, 2, 3, 4]
   y = ch1(x)
   n = np.random.randn(90)
   ny = y + n
   ch2 = models.Chebyshev1DModel(3)
   chfit = fitting.LinearLSQFitter(ch2)
   chfit(x, ny)
   plt.plot(x, y, label='y - Chebyshev polynomial')
   plt.plot(x, ny, label='ny - Chebyshev polynomial with noise')
   plt.plot(x, ch2(x), label='ch2(x) - Fitted model')
   plt.legend()
   plt.show()


Working with 2D models
======================

First create some data to be fitted with a 2D polynomial::

    >>> x, y = np.mgrid[:10, :10]
    >>> def poly2(x, y):
            return 1+2*x+3*x**2+4*y+5*y**2+6*x*y
    >>> z = poly2(x, y)

Fit a 2D polynomial to the data::

    >>> p2 = models.Polynomial2DModel(2)
    >>> print(p2)
    Model: Polynomial2DModel
    n_inputs:   2
    Degree: 2
    Parameter sets: 1
    Parameters:
               c0_0: Parameter('c0_0', value=0.0)
               c1_0: Parameter('c1_0', value=0.0)
               c2_0: Parameter('c2_0', value=0.0)
               c0_1: Parameter('c0_1', value=0.0)
               c0_2: Parameter('c0_2', value=0.0)
               c1_1: Parameter('c1_1', value=0.0)
    >>> pfit = fitting.LinearLSQFitter(p2)
    >>> n = np.random.randn(100)
    >>> n.shape = (10, 10)
    >>> pfit(x, y, z+n)
    >>> p2.parameters
    array([ 0.97599264,  1.95050208,  3.00524297,  4.01663038,  5.00150801,
            5.999489  ])


Using `modeling`
================

.. toctree::
   :maxdepth: 1

   parameters
   models
   fitting
   new
   algorithms
   design


Reference/API
=============

.. automodapi:: astropy.modeling
.. automodapi:: astropy.modeling.fitting
.. automodapi:: astropy.modeling.functional_models
.. automodapi:: astropy.modeling.powerlaws
.. automodapi:: astropy.modeling.polynomial
.. automodapi:: astropy.modeling.projections
.. automodapi:: astropy.modeling.rotations
