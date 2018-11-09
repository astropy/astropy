.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-modeling-performance:

.. Performance Tips
.. ================
..
.. Here we provide some tips and tricks for how to optimize performance of code
using `astropy.modeling`.

Initializing a compound model with many constituent models can be time consuming. If your code uses the same compound model repeatedly consider initializing it once and reusing the model.

Consider the performance tips that apply to quantities when initializing and evaluating models with quantities.
