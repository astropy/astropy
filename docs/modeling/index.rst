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
    Astropy. If you have specific ideas for how it might be improved,
    feel free to let us know on the `astropy-dev mailing list`_ or at
    http://feedback.astropy.org


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
    >>> g1 = models.Gaussian1D(10., stddev=2.1, mean=4.2)
    >>> g1
    <Gaussian1D(amplitude=Parameter('amplitude', value=10.0),
                mean=Parameter('mean', value=4.2000000000000002),
                stddev=Parameter('stddev', value=2.1000000000000001),
                param_dim=1)>
    >>> y = g1(x)
    >>> np.random.seed(0)
    >>> n = np.random.randn(90)
    >>> ny = y + n
    >>> gfit = fitting.NonLinearLSQFitter()
    >>> new_model = gfit(g1, x, ny)
    >>> print(new_model)
    Model: Gaussian1D
    n_inputs:   1
    Degree: N/A
    Parameter sets: 1
    Parameters:
               amplitude: Parameter('amplitude', value=9.8931826765510706)
               mean: Parameter('mean', value=4.0263781556737737)
               stddev: Parameter('stddev', value=2.152396119425859)

Create data using a 1D Chebyshev model::

    >>> ch1 = models.Chebyshev1D(3, domain=[x.min(), x.max()])
    >>> ch1.parameters
    array([ 0.,  0.,  0.,  0.])
    >>> ch1.parameters = [1, 2, 3, 4]
    >>> ch1.parameters
    array([ 1.,  2.,  3.,  4.])
    >>> print(ch1)
    Model: Chebyshev1D
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

    >>> np.random.seed(0)
    >>> n = np.random.randn(90)
    >>> ny = y + n

Fit a Chebyshev polynomial to the data::

    >>> ch2 = models.Chebyshev1D(3)
    >>> chfit = fitting.LinearLSQFitter()
    >>> new_model = chfit(ch2, x, ny)
    >>> new_model.parameters
    array([ 1.17789166,  1.67145195,  3.53825251,  4.05892813])

.. plot::

   import matplotlib.pyplot as plt
   import numpy as np
   from astropy.modeling import models, fitting
   x = np.arange(1, 10, .1)
   ch1 = models.Chebyshev1D(3, domain=[x.min(), x.max()])
   ch1.parameters = [1, 2, 3, 4]
   y = ch1(x)
   np.random.seed(0)
   n = np.random.randn(90)
   ny = y + n
   ch2 = models.Chebyshev1D(3)
   chfit = fitting.LinearLSQFitter()
   model = chfit(ch2, x, ny)
   plt.plot(x, y, label='y - Chebyshev polynomial')
   plt.plot(x, ny, label='ny - Chebyshev polynomial with noise')
   plt.plot(x, model(x), label='ch2(x) - Fitted model')
   plt.legend()
   plt.show()


Working with 2D models
======================

First create some data to be fitted with a 2D polynomial::

    >>> x, y = np.mgrid[:10, :10]
    >>> def poly2(x, y):
    ...     return 1+2*x+3*x**2+4*y+5*y**2+6*x*y
    >>> z = poly2(x, y)

Fit a 2D polynomial to the data::

    >>> p2 = models.Polynomial2D(2)
    >>> print(p2)
    Model: Polynomial2D
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

    >>> pfit = fitting.LinearLSQFitter()
    >>> np.random.seed(0)
    >>> n = np.random.randn(100)
    >>> n.shape = (10, 10)
    >>> new_model = pfit(p2, x, y, z+n)
    >>> new_model.parameters
    array([ 1.79964917,  1.44891526,  3.05358047,  4.08895144,  4.98756933,
            6.00824639])


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
