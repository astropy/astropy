.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the sub-package toctree

.. _astropy-convolution-performance:

Performance Tips
================

The :func:`~astropy.convolution.convolve` function is best suited to small
kernels, and can become very slow for larger kernels. In this case, consider
using :func:`~astropy.convolution.convolve_fft` (though note that this function
uses more memory).
