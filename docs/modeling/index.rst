.. include:: links.inc

.. _astropy-modeling:

***************************************
Models and Fitting (`astropy.modeling`)
***************************************

Introduction
============

`astropy.modeling` provides a framework for representing models and performing
model evaluation and fitting. It currently supports 1-D and 2-D models and
fitting with parameter constraints.

It is :ref:`designed <modeling-design>` to be easily extensible and flexible.
Models do not reference fitting algorithms explicitly and new fitting
algorithms may be added without changing the existing models (though not all
models can be used with all fitting algorithms due to constraints such as model
linearity).

The goal is to eventually provide a rich toolset of models and fitters such
that most users will not need to define new model classes, nor special purpose
fitting routines (while making it reasonably easy to do when necessary).

.. warning::

    `astropy.modeling` is currently a work-in-progress, and thus it is likely
    there will be significant API changes in later versions of Astropy. If you
    have specific ideas for how it might be improved, feel free to let us know
    on the `astropy-dev mailing list`_ or at http://feedback.astropy.org


Getting started
===============

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
    Inputs: 1
    Outputs: 1
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


Simple 2-D model fitting
------------------------

Similarly to the 1-D example, we can create a simulated 2-D data dataset, and
fit a polynomial model to it.  This could be used for example to fit the
background in an image.

.. plot::
   :include-source:

    import warnings
    import numpy as np
    from astropy.modeling import models, fitting

    # Generate fake data
    np.random.seed(0)
    y, x = np.mgrid[:128, :128]
    z = 2. * x ** 2 - 0.5 * x ** 2 + 1.5 * x * y - 1.
    z += np.random.normal(0., 0.1, z.shape) * 50000.

    # Fit the data using astropy.modeling
    p_init = models.Polynomial2D(degree=2)
    fit_p = fitting.LevMarLSQFitter()

    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        p = fit_p(p_init, x, y, z)

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


Model sets
----------

In some cases it is necessary to describe many models of the same type but with
different parameter values.  This could be done simply by instantiating as many
instances of a `~astropy.modeling.Model` as are needed.  But that can be
inefficient for a large number of models.  To that end, all model classes in
`astropy.modeling` can also be used to represent a model *set* which is a
collection of models of the same type, but with different values for their
parameters.

To instantiate a model set, use argument ``n_models=N`` where ``N`` is the
number of models in the set when constructing the model.  The value of each
parameter must be a list or array of length ``N``, such that each item in
the array corresponds to one model in the set::

    >>> g = models.Gaussian1D(amplitude=[1, 2], mean=[0, 0],
    ...                       stddev=[0.1, 0.2], n_models=2)
    >>> print(g)
    Model: Gaussian1D
    Inputs: 1
    Outputs: 1
    Model set size: 2
    Parameters:
        amplitude mean stddev
        --------- ---- ------
              1.0  0.0    0.1
              2.0  0.0    0.2

This is equivalent to two Gaussians with the parameters ``amplitude=1, mean=0,
stddev=0.1`` and ``amplitude=2, mean=0, stddev=0.2`` respectively.  When
printing the model the parameter values are displayed as a table, with each row
corresponding to a single model in the set.

The number of models in a model set can be determined using the `len` builtin::

    >>> len(g)
    2

Single models have a length of 1, and are not considered a model set as such.

When evaluating a model set, by default the input must be the same length as
the number of models, with one input per model::

    >>> g([0, 0.1])
    array([ 1.        ,  1.76499381])

The result is an array with one result per model in the set.  It is also
possible to broadcast a single value to all models in the set::

    >>> g(0)
    array([ 1.,  2.])

Model sets are used primarily for fitting, allowing a large number of models of
the same type to be fitted simultaneously (and independently from each other)
to some large set of inputs.  For example, fitting a polynomial to the time
response of each pixel in a data cube.  This can greatly speed up the fitting
process, especially for linear models.


Using `astropy.modeling`
========================

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
.. automodapi:: astropy.modeling.optimizers
.. automodapi:: astropy.modeling.powerlaws
.. automodapi:: astropy.modeling.polynomial
.. automodapi:: astropy.modeling.projections
.. automodapi:: astropy.modeling.statistic
.. automodapi:: astropy.modeling.rotations

