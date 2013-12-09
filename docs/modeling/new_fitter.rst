***************************
Defining a new Fitter class
***************************

If you would like to use a fit statistic or minimizer different from the ones built into `~astropy.modeling`,
you have to create a subclass of `~astropy.modeling.fitting.Fitter`. 

Here is an example how to define a least squares fitter for nonlinear models
by wrapping the `~scipy.optimize.fmin_slsqp`  minimisation function:

.. literalinclude:: new_fitter.py
   :linenos:

The user calls ``Chi2Fitter.__call__`` with a given ``model`` and data ``x`` and ``y``,
which checks and reformats the inputs a bit and then calls `~scipy.optimize.fmin_slsqp`,
which repeatedly calls ``Chi2Fitter.errorfunc`` in each fit iteration. 

Let's go through it line by line:

* Line TODO: Create ``Chi2Fitter`` as a sub-class of `~astropy.modeling.fitting.Fitter`.
* Line TODO: Declare the supported contraints of this fitter as a list of strings.
  TODO: link to complete list and explanation of supported constraints types.

The special method ``Chi2Fitter.errorfunc`` computes the chi-square fit statistic:

* The ``fitparams`` array contains the parameter values in the current fit iteration.
* The ``args`` parameter tuple can be used to pass other parameters needed to compute the
  fit statistic, in this case ``args[0] = model``, ``args[1] = x`` and ``args[2] = y``.
* Line TODO: The `~astropy.modeling.fitting.Fitter._fitter_to_model_params` helper function is used
  to set the model parameters.
* Line TODO: The model is evaluated and the fit statistic
  (in this case the sum of squared residuals) computed and returned.

The special method ``Chi2Fitter.__call__`` runs the minimisation:

* The main inputs are the ``model`` and the data ``x`` and ``y``.
  Extra arguments like ``maxiter`` to control the fitting can be added. 
* Line TODO: Import the `~scipy.optimize.fmin_slsqp` minimisation function.
  Lines TODO to TODO: Checks if the model is fittable and has valid constraints.
* Line TODO: Make a copy of the model, because fitters are supposed to leave the state of the
  input model unchanged and return a copy with optimised parameters.
* Line TODO: Call the `~astropy.modeling.fitting.Fitter._model_to_fit_params` helper function
  to extract the input parameter values as a numpy array ``x0``.
* Line TODO: Call `~scipy.optimize.fmin_slsqp`, passing the fit statistic function ``self.errorfunc``,
  the initial parameter values ``x0``, the ``model`` and data ``x`` and ``y`` as well a the ``iter`` parameter
  that defines the maximum allowed number of ``self.errorfunc`` evaluations. 
* Line TODO: Call the `~astropy.modeling.fitting.Fitter._fitter_to_model_params` helper function
  to copy the ``fitparams`` returned by ``fmin_slsqp`` into the ``model_copy``.

Usage
-----

You can use this user-defined `SLSQPFitter` just like the built-in astropy Fitters:: 

   import numpy as np
   from astropy.modeling import models
   
   g = models.Gaussian1D(amplitude=1.2, mean=0.9, stddev=0.5)
   
   # Generate fake data
   np.random.seed(0)
   x = np.linspace(-5., 5., 100)
   y = g(x) + np.random.normal(0., 0.2, x.shape)
   
   # Fit the data using a Gaussian
   fitter = SLSQPFitter()
   g = fitter(g, x, y, maxiter=50)
