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


Using Models
------------

The `astropy.modeling` package defines a number of models that live inside
`astropy.modeling.models` and behave like parametrized functions::

    >>> from astropy.modeling import models
    >>> g = models.Gaussian1D(amplitude=1.2, mean=0.9, stddev=0.5)
    >>> print(g)
    Model: Gaussian1DModel
    n_inputs:   1
    Degree: N/A
    Parameter sets: 1
    Parameters:
               amplitude: Parameter('amplitude', value=1.2)
               mean: Parameter('mean', value=0.9000...)
               stddev: Parameter('stddev', value=0.5)

Model parameters can be accessed as attributes:

    >>> g.amplitude
    Parameter('amplitude', value=1.2)
    >>> g.mean
    Parameter('mean', value=0.9000...)
    >>> g.stddev
    Parameter('stddev', value=0.5)

and can also be set using the attributes::

    >>> g.amplitude = 0.8
    >>> g.amplitude
    Parameter('amplitude', value=0.8000...)

Models can be evaluated by calling them as functions::

    >>> g(0.1)
    0.22242984036255528
    >>> g(np.linspace(0.5, 1.5, 7))
    array([ 0.58091923,  0.71746405,  0.7929204 ,  0.78415894,  0.69394278,
            0.54952605,  0.3894018 ])

Models can therefore already be useful to evaluate common functions,
independent of the fitting part of the package.

Simple 1D model fitting
-----------------------

In this section, we look at a simple example of fitting a Gaussian to a
simulated dataset. We use the :class:`~astropy.modeling.functional_models.Gaussian1D`
and :class:`~astropy.modeling.functional_models.Trapezoid1D` models and the
:class:`~astropy.modeling.fitting.NonLinearLSQFitter` fitter to
fit the data:

.. plot::
   :include-source:

    import numpy as np
    from astropy.modeling import models, fitting

    # Generate fake data
    np.random.seed(0)
    x = np.linspace(-5., 5., 200)
    y = 3 * np.exp(-0.5 * (x - 1.3)**2 / 0.8**2)
    y += np.random.normal(0., 0.2, x.shape)

    # Fit the data using a box model
    t_init = models.Trapezoid1D(amplitude=1., x_0=0., width=1., slope=0.5)
    f1 = fitting.NonLinearLSQFitter()
    t = f1(t_init, x, y)

    # Fit the data using a Gaussian
    g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
    f2 = fitting.NonLinearLSQFitter()
    g = f2(g_init, x, y)

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,5))
    plt.plot(x, y, 'ko')
    plt.plot(x, t(x), 'b-', lw=2, label='Trapezoid')
    plt.plot(x, g(x), 'r-', lw=2, label='Gaussian')
    plt.xlabel('Position')
    plt.ylabel('Flux')
    plt.legend(loc=2)

As shown above, once instantiated, the fitter class can be used as a function
that takes the initial model (``t_init`` or ``g_init``) and the data values
(``x`` and ``y``), and returns a fitted model (``t`` or ``g``).

Simple 2D model fitting
-----------------------

Similarly to the 1-d example, we can create a simulated 2-d data dataset, and fit a polynomial model to it. This could be used for example to fit the background in an image.

.. plot::
   :include-source:

    import numpy as np
    from astropy.modeling import models, fitting

    # Generate fake data
    np.random.seed(0)
    x, y = np.mgrid[:128, :128]
    z = 2. * x ** 2 - 0.5 * x ** 2 + 1.5 * x * y - 1.
    z += np.random.normal(0., 0.1, z.shape) * 50000.

    # Fit the data using astropy.modeling
    p_init = models.Polynomial2D(degree=2)
    f = fitting.NonLinearLSQFitter()
    p = f(p_init, x, y, z)

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,3,1)
    plt.imshow(z, interpolation='nearest', vmin=-1e4, vmax=5e4)
    plt.title("Data")
    plt.subplot(1,3,2)
    plt.imshow(p(x, y), interpolation='nearest', vmin=-1e4, vmax=5e4)
    plt.title("Model")
    plt.subplot(1,3,3)
    plt.imshow(z - p(x, y), interpolation='nearest', vmin=-1e4, vmax=5e4)
    plt.title("Residual")

A list of models is provided in the `Reference/API`_ section. The fitting
framework includes many useful features that are not demonstrated here, such as
weighting of datapoints, fixing or linking parameters, and placing lower or
upper limits on parameters. For more information on these, take a look at the
:doc:`fitting` documentation.

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
