**********************
Fitting Models to Data
**********************

This module provides wrappers, called Fitters, around some Numpy and Scipy
fitting functions. All Fitters can be called as functions. They take an
instance of `~astropy.modeling.FittableModel` as input and modify its
``parameters`` attribute. The idea is to make this extensible and allow
users to easily add other fitters.

Linear fitting is done using Numpy's `numpy.linalg.lstsq` function. There are
currently two non-linear fitters which use `scipy.optimize.leastsq` and
`scipy.optimize.fmin_slsqp`.

The rules for passing input to fitters are:

* Non-linear fitters currently work only with single models (not model sets).

* The linear fitter can fit a single input to multiple model sets creating
  multiple fitted models. This may require specifying the ``model_set_axis``
  argument just as used when evaluating models; this may be required for the
  fitter to know how to broadcast the input data.

.. _constraining_parameters:

Constrained fitting
===================

Fitters support constrained fitting.

- All fitters support fixed (frozen) parameters through the ``fixed`` argument
  to models or setting the `~astropy.modeling.Parameter.fixed`
  attribute directly on a parameter.

  For linear fitters, freezing a polynomial coefficient means that a polynomial
  without that term will be fitted to the data. For example, fixing ``c0`` in a
  polynomial model will fit a polynomial with the zero-th order term missing.
  However, the fixed value of the coefficient is used when evaluating the
  model::

      >>> x = np.arange(1, 10, .1)
      >>> p1 = models.Polynomial1D(2, c0=[1, 1], c1=[2, 2], c2=[3, 3],
      ...                          n_models=2)
      >>> p1
      <Polynomial1D(2, c0=[ 1., 1.], c1=[ 2., 2.], c2=[ 3., 3.], n_models=2)>
      >>> y = p1(x, model_set_axis=False)
      >>> p1.c0.fixed = True
      >>> pfit = fitting.LinearLSQFitter()
      >>> new_model = pfit(p1, x, y)
      >>> print(new_model)  # doctest: +SKIP
      Model: Polynomial1D
      Inputs: 1
      Outputs: 1
      Model set size: 2
      Degree: 2
      Parameters:
           c0     c1         c2
          --- ------------- -------------
          1.0 2.38641216243 2.96827885742
          1.0 2.38641216243 2.96827885742

- A parameter can be `~astropy.modeling.Parameter.tied` (linked to
  another parameter). This can be done in two ways::

      >>> def tiedfunc(g1):
      ...    mean = 3 * g1.stddev
      ...    return mean
      >>> g1 = models.Gaussian1D(amplitude=10., mean=3, stddev=.5,
      ...                        tied={'mean': tiedfunc})

  or::

      >>> g1 = models.Gaussian1D(amplitude=10., mean=3, stddev=.5)
      >>> g1.mean.tied = tiedfunc

- Bounded fitting is supported through the ``bounds`` arguments to models or by
  setting `~astropy.modeling.Parameter.min` and `~astropy.modeling.Parameter.max`
  attributes on a parameter. Bounds for the
  `~astropy.modeling.fitting.LevMarLSQFitter` are always exactly satisfied--if
  the value of the parameter is outside the fitting interval, it will be reset
  to the value at the bounds. The `~astropy.modeling.fitting.SLSQPLSQFitter`
  handles bounds internally.

      >>> g1 = models.Gaussian1D(amplitude=10., mean=3, stddev=.5)
      >>> g1.mean.min = 2.5
      >>> g1.mean.max = 3.2

- Different fitters support different types of constraints::

    >>> fitting.LinearLSQFitter.supported_constraints
    ['fixed']
    >>> fitting.LevMarLSQFitter.supported_constraints
    ['fixed', 'tied', 'bounds']
    >>> fitting.SLSQPLSQFitter.supported_constraints
    ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']


Fitting examples
================

Simple 2-D model fitting
------------------------

Similarly to the :ref:`1-D example <1D-fitting>`, we can create a simulated 2-D data dataset, and
fit a polynomial model to it. This could be used for example to fit the
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
    plt.imshow(z, origin='lower', interpolation='nearest', vmin=-1e4, vmax=5e4)
    plt.title("Data")
    plt.subplot(1,3,2)
    plt.imshow(p(x, y), origin='lower', interpolation='nearest', vmin=-1e4,
               vmax=5e4)
    plt.title("Model")
    plt.subplot(1,3,3)
    plt.imshow(z - p(x, y), origin='lower', interpolation='nearest', vmin=-1e4,
               vmax=5e4)
    plt.title("Residual")

Fitting to multiple datasets
----------------------------

The package also allows fitting a polynomial model to multiple data sets
simultaneously::

      >>> from astropy.modeling import models, fitting
      >>> import numpy as np
      >>> p1 = models.Polynomial1D(3)
      >>> p1.c0 = 1
      >>> p1.c1 = 2
      >>> print(p1)
      Model: Polynomial1D
      Inputs: ('x',)
      Outputs: ('y',)
      Model set size: 1
      Degree: 3
      Parameters:
           c0  c1  c2  c3
          --- --- --- ---
          1.0 2.0 0.0 0.0
      >>> x = np.arange(10)
      >>> y = p1(x)
      >>> yy = np.array([y, y])
      >>> p2 = models.Polynomial1D(3, n_models=2)
      >>> pfit = fitting.LinearLSQFitter()
      >>> new_model = pfit(p2, x, yy)
      >>> print(new_model)  # doctest: +SKIP
      Model: Polynomial1D
      Inputs: 1
      Outputs: 1
      Model set size: 2
      Degree: 3
      Parameters:
           c0  c1         c2                 c3
          --- --- ------------------ -----------------
          1.0 2.0 -5.86673908219e-16 3.61636197841e-17
          1.0 2.0 -5.86673908219e-16 3.61636197841e-17
