.. include:: links.inc

.. _astropy-modeling:

***************************************
Models and Fitting (`astropy.modeling`)
***************************************

Introduction
============

`astropy.modeling` provides a framework for representing models and performing
model evaluation and fitting. The major components of the package are models,
parameters and fitters. Models are behave like parametrized functions, with
Parameters set as attributes. Fitters are used to fit specific models to
datasets, and Parameters can easily be tied or constrained for the fitting
routine.

This module currently supports 1-D and 2-D models. It is designed to be
easily extensible and flexible.  Models do not reference fitting algorithms
explicitly and new fitting algorithms may be added without changing the
existing models (though not all models can be used with all fitting algorithms
due to constraints such as model linearity).

The goal is to eventually provide a rich toolset of models and fitters such
that most users will not need to define new model classes, nor special purpose
fitting routines (while making it reasonably easy to do when necessary).

.. note::

    `astropy.modeling` is currently a work-in-progress, and thus it is likely
    there will still be API changes in later versions of Astropy.  Backwards
    compatibility support between versions will still be maintained as much as
    possible, but new features and enhancements are coming in future versions.
    If you have specific ideas for how it might be improved, feel free to let
    us know on the `astropy-dev mailing list`_ or at
    http://feedback.astropy.org


Using `astropy.modeling`
========================

.. toctree::
   :maxdepth: 1

   models
   parameters
   fitting
   compound-models
   new
   algorithms
   api


.. _quick_tutorial:

Quick Tutorial
==============

The examples here use the predefined models and assume the following modules
have been imported::

    >>> import numpy as np
    >>> from astropy.modeling import models, fitting


Using Models
------------

The `astropy.modeling` package defines a number of models that are collected
under a single namespace as ``astropy.modeling.models``.  Models behave like
parametrized functions::

    >>> from astropy.modeling import models
    >>> g = models.Gaussian1D(amplitude=1.2, mean=0.9, stddev=0.5)
    >>> print(g)
    Model: Gaussian1D
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Parameters:
        amplitude mean stddev
        --------- ---- ------
              1.2  0.9    0.5

Model parameters can be accessed as attributes::

    >>> g.amplitude
    Parameter('amplitude', value=1.2)
    >>> g.mean
    Parameter('mean', value=0.9)
    >>> g.stddev
    Parameter('stddev', value=0.5)

and can also be updated via those attributes::

    >>> g.amplitude = 0.8
    >>> g.amplitude
    Parameter('amplitude', value=0.8)

Models can be evaluated by calling them as functions::

    >>> g(0.1)
    0.22242984036255528
    >>> g(np.linspace(0.5, 1.5, 7))
    array([ 0.58091923,  0.71746405,  0.7929204 ,  0.78415894,  0.69394278,
            0.54952605,  0.3894018 ])

As the above example demonstrates, in general most models evaluate array-like
inputs according to the standard `Numpy broadcasting rules`_ for arrays.

Models can therefore already be useful to evaluate common functions,
independently of the fitting features of the package.

.. _1D-fitting:

Simple 1-D model fitting
------------------------

In this section, we look at a simple example of fitting a Gaussian to a
simulated dataset. We use the `~astropy.modeling.functional_models.Gaussian1D`
and `~astropy.modeling.functional_models.Trapezoid1D` models and the
`~astropy.modeling.fitting.LevMarLSQFitter` fitter to fit the data:

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
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, x, y)

    # Fit the data using a Gaussian
    g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y)

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

The fitting framework includes many useful features that are not demonstrated
here, such as weighting of datapoints, fixing or linking parameters, and
placing lower or upper limits on parameters. For more information on these,
take a look at the :doc:`fitting` documentation. A list of available models is
provided :ref:`here <models_list>`.

.. _compound-models-intro:

Compound models
---------------
.. versionadded:: 1.0

    This feature is experimental and expected to see significant further
    development, but the basic usage is stable and expected to see wide use.

While the Astropy modeling package makes it very easy to define :doc:`new
models <new>` either from existing functions, or by writing a
`~astropy.modeling.Model` subclass, an additional way to create new models is
by combining them using arithmetic expressions.  This works with models built
into Astropy, and most user-defined models as well.  For example, it is
possible to create a superposition of two Gaussians like so::

    >>> from astropy.modeling import models
    >>> g1 = models.Gaussian1D(1, 0, 0.2)
    >>> g2 = models.Gaussian1D(2.5, 0.5, 0.1)
    >>> g1_plus_2 = g1 + g2

The resulting object ``g1_plus_2`` is itself a new model.  Evaluating, say,
``g1_plus_2(0.25)`` is the same as evaluating ``g1(0.25) + g2(0.25)``::

    >>> g1_plus_2(0.25)  # doctest: +FLOAT_CMP
    0.5676756958301329
    >>> g1_plus_2(0.25) == g1(0.25) + g2(0.25)
    True

This model can be further combined with other models in new expressions.  It is
also possible to define entire new model *classes* using arithmetic expressions
of other model classes.  This allows general compound models to be created
without specifying any parameter values up front.  This more advanced usage is
explained in more detail in the :ref:`compound model documentation
<compound-model-classes>`.

These new compound models can also be fitted to data, like most other models:

.. plot::
    :include-source:

    import numpy as np
    from astropy.modeling import models, fitting

    # Generate fake data
    np.random.seed(42)
    g1 = models.Gaussian1D(1, 0, 0.2)
    g2 = models.Gaussian1D(2.5, 0.5, 0.1)
    x = np.linspace(-1, 1, 200)
    y = g1(x) + g2(x) + np.random.normal(0., 0.2, x.shape)

    # Now to fit the data create a new superposition with initial
    # guesses for the parameters:
    gg_init = models.Gaussian1D(1, 0, 0.1) + models.Gaussian1D(2, 0.5, 0.1)
    fitter = fitting.SLSQPLSQFitter()
    gg_fit = fitter(gg_init, x, y)

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,5))
    plt.plot(x, y, 'ko')
    plt.plot(x, gg_fit(x), 'r-', lw=2)
    plt.xlabel('Position')
    plt.ylabel('Flux')

This works for 1-D models, 2-D models, and combinations thereof, though there
are some complexities involved in correctly matching up the inputs and outputs
of all models used to build a compound model.  You can learn more details in
the :doc:`compound-models` documentation.
