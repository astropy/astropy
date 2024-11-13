
.. _astropy-modeling-performance:

Performance Tips
****************

Initializing a compound model with many constituent models can be time consuming.
If your code uses the same compound model repeatedly consider initializing it
once and reusing the model.

Consider the :ref:`performance tips <astropy-units-performance>` that apply to
quantities when initializing and evaluating models with quantities.

When fitting models with one of the fitter classes, by default a copy of the
model is returned, with parameters set to those determined by the fitting. If
you do not need to preserve the initial model used in the fitting, you can
optionally pass ``inplace=True`` when calling the fitter, and the parameters
will be updated on the model you supply rather than returning a copy of the
model - this can improve performance in some cases.
