Fitting with constraints
========================

Fitters support constrained fitting.

- All fitters support fixed (frozen) parameters through the ``fixed`` argument
  to models or setting the `~astropy.modeling.Parameter.fixed`
  attribute directly on a parameter.

  For linear fitters, freezing a polynomial coefficient means that the
  corresponding term will be subtracted from the data before fitting a
  polynomial without that term to the result. For example, fixing ``c0`` in a
  polynomial model will fit a polynomial with the zero-th order term missing
  to the data minus that constant. However, the fixed coefficient value is
  restored when evaluating the model, to fit the original data values::

      >>> import numpy as np
      >>> from astropy.modeling import models, fitting
      >>> x = np.arange(1, 10, .1)
      >>> p1 = models.Polynomial1D(2, c0=[1, 1], c1=[2, 2], c2=[3, 3],
      ...                          n_models=2)
      >>> p1  # doctest: +FLOAT_CMP
      <Polynomial1D(2, c0=[1., 1.], c1=[2., 2.], c2=[3., 3.], n_models=2)>
      >>> y = p1(x, model_set_axis=False)
      >>> p1.c0.fixed = True
      >>> pfit = fitting.LinearLSQFitter()
      >>> new_model = pfit(p1, x, y)
      >>> print(new_model)  # doctest: +SKIP
      Model: Polynomial1D
      Inputs: ('x',)
      Outputs: ('y',)
      Model set size: 2
      Degree: 2
      Parameters:
           c0  c1  c2
          --- --- ---
          1.0 2.0 3.0
          1.0 2.0 3.0

  The syntax to fix the same parameter ``c0`` using an argument to the model
  instead of ``p1.c0.fixed = True`` would be::

      >>> p1 = models.Polynomial1D(2, c0=[1, 1], c1=[2, 2], c2=[3, 3],
      ...                          n_models=2, fixed={'c0': True})


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

Bounded fitting is supported through the ``bounds`` arguments to models or by
setting `~astropy.modeling.Parameter.min` and `~astropy.modeling.Parameter.max`
attributes on a parameter.  Bounds for the
`~astropy.modeling.fitting.LevMarLSQFitter` are always exactly satisfied--if
the value of the parameter is outside the fitting interval, it will be reset to
the value at the bounds. The `~astropy.modeling.fitting.SLSQPLSQFitter` handles
bounds internally.

- Different fitters support different types of constraints::

    >>> fitting.LinearLSQFitter.supported_constraints
    ['fixed']
    >>> fitting.LevMarLSQFitter.supported_constraints
    ['fixed', 'tied', 'bounds']
    >>> fitting.SLSQPLSQFitter.supported_constraints
    ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']

Note that there are two "constraints" (``prior`` and ``posterior``) that are
not currently used by any of the built-in fitters.  They are provided to allow
possible user code that might implement Bayesian fitters (e.g.,
https://gist.github.com/rkiman/5c5e6f80b455851084d112af2f8ed04f).
