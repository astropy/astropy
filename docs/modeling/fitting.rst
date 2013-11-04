**********************
Fitting Models to Data
**********************

This module provides wrappers, called Fitters, around some Numpy and Scipy
fitting functions. All Fitters can be called as functions. They take an
instance of `~astropy.modeling.core.ParametricModel` as input and modify
`~astropy.modeling.core.ParametricModel.parameters` attribute. The idea is to
make this extensible and allow users to easily add other fitters.

Linear fitting is done using Numpy's `~numpy.linalg.lstsq` function.  There are
currently two non-linear fitters which use `~scipy.optimize.leastsq` and
`~scipy.optimize.fmin_slsqp`.

The rules for passing input to fitters are:

* Non-linear fitters work only with single data sets.

* The linear fitter can fit a single input to multiple data sets creating
  multiple parameter sets. For example fitting a 2D model with input x, y
  arrays of shape (n, m) to a z array of shape (p, n, m), will set
  model.parameters.n_inputs to p, even if it was 1 when the model was created.

* Attempting to fit a model with multiple parameter sets to a single data set
  results in an error.


Fitting examples
----------------

- Fitting a polynomial model to multiple data sets simultaneously::

    >>> from astropy.modeling import models, fitting
    >>> import numpy as np
    >>> p1 = models.Polynomial1D(3)
    >>> p1.c0 = 1
    >>> p1.c1 = 2
    >>> p1.parameters
    array([ 1.,  2.,  0.,  0.])
    >>> x = np.arange(10)
    >>> y = p1(x)
    >>> yy = np.array([y, y]).T
    >>> p2 = models.Polynomial1D(3, param_dim=2)
    >>> pfit = fitting.LinearLSQFitter()
    >>> new_model = pfit(p2, x, yy)
    >>> print(new_model.param_sets)
    [[  1.00000000e+00   1.00000000e+00]
     [  2.00000000e+00   2.00000000e+00]
     [  3.88335494e-16   3.88335494e-16]
     [ -2.99749607e-17  -2.99749607e-17]]

Fitters support constrained fitting.

- All fitters support fixed (frozen) parameters through the ``fixed`` argument
  to models or setting the `~astropy.modeling.parameters.Parameter.fixed`
  attribute directly on a parameter.

  For linear fitters, freezing a polynomial coefficient means that a polynomial
  without that term will be fitted to the data. For example, fixing ``c0`` in a
  polynomial model will fit a polynomial with the zero-th order term missing.
  However, the fixed value of the coefficient is used when evaluating the
  model::

      >>> x = np.arange(1, 10, .1)
      >>> p1 = models.Polynomial1D(2, param_dim=2)
      >>> p1.parameters = [1, 1, 2, 2, 3, 3]
      >>> p1.param_sets
      array([[ 1.,  1.],
             [ 2.,  2.],
             [ 3.,  3.]])
      >>> y = p1(x)
      >>> p1.c0.fixed = True
      >>> pfit = fitting.LinearLSQFitter()
      >>> new_model = pfit(p1, x, y)
      ...
      >>> new_model.param_sets  # doctest: +SKIP
      array([[ 1.,          1.        ],
             [ 2.38641216,  2.38641216],
             [ 2.96827886,  2.96827886]])


- A parameter can be `~astropy.modeling.parameters.Parameter.tied` (linked to
  another parameter). This can be done in two ways::

      >>> def tiedfunc(g1):
      ...    mean = 3 * g1.stddev
      ...    return mean
      >>> g1 = models.Gaussian1D(amplitude=10., mean=3, stddev=.5,
      ...                        tied={'mean': tiedfunc})

  or::

      >>> g1 = models.Gaussian1D(amplitude=10., mean=3, stddev=.5)
      >>> g1.mean.tied = tiedfunc
      >>> gfit = fitting.NonLinearLSQFitter()

Bounded fitting is supported through the ``bounds`` arguments to models or by
setting `~astropy.modeling.parameters.Parameter.min` and
`~astropy.modeling.parameters.Parameter.max` attributes on a parameter.  Bounds
for the `~astropy.modeling.fitting.NonLinearLSQFitter` are always exactly
satisfied--if the value of the parameter is outside the fitting interval, it
will be reset to the value at the bounds. The
`~astropy.modeling.fitting.SLSQPFitter` handles bounds internally.

- Different fitters support different types of constraints::

    >>> fitting.LinearLSQFitter.supported_constraints
    ['fixed']
    >>> fitting.NonLinearLSQFitter.supported_constraints
    ['fixed', 'tied', 'bounds']
    >>> fitting.SLSQPFitter.supported_constraints
    ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']
