
.. _astropy-modeling-performance:

Performance Tips
****************

Initializing a compound model with many constituent models can be time consuming.
If your code uses the same compound model repeatedly consider initializing it
once and reusing the model.

Consider the :ref:`performance tips <astropy-units-performance>` that apply to
quantities when initializing and evaluating models with quantities.
