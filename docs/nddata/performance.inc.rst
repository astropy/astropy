.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-nddata-performance:

Performance Tips
================

+ Using the uncertainty class `~astropy.nddata.VarianceUncertainty` will
  be somewhat more efficient than the other two uncertainty classes,
  `~astropy.nddata.InverseVariance` and `~astropy.nddata.StdDevUncertainty`.
  The latter two are converted to variance for the purposes of error
  propagation and then converted from variance back to the original
  uncertainty type. The performance difference should be small.
+ When possible, mask values by setting them to ``np.nan`` and use the
  ``numpy`` functions and methods that automatically exclude ``np.nan``,
  like ``np.nanmedian`` and ``np.nanstd``. This will typically be much
  faster than using `numpy.ma.MaskedArray`.
