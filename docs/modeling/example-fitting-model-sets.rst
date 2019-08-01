Fitting Model Sets
==================

Fitting a polynomial model to multiple data sets simultaneously::

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
