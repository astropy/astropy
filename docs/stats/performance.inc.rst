.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-stats-performance:

Performance Tips
================

If you are finding sigma clipping to be slow, and if you have not already done
so, consider installing the `bottleneck <https://pypi.org/project/Bottleneck/>`_
package, which will speed up some of the internal computations. In addition, if
you are using standard functions for ``cenfunc`` and/or ``stdfunc``, make sure
you specify these as strings rather than passing a NumPy function — that is,
use::

    >>> sigma_clip(array, cenfunc='median')  # doctest: +SKIP

instead of::

    >>> sigma_clip(array, cenfunc=np.nanmedian)  # doctest: +SKIP

Using strings will allow the sigma-clipping algorithm to pick the fastest
implementation available for finding the median.
