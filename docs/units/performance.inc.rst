.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-units-performance:

Performance Tips
================

Here we provide some tips and tricks for how to optimize performance of code
using `astropy.units`.

Initialization without a copy::

  >>> import numpy as np
  >>> import astropy.units as u
  >>> a = np.arange(5.)
  >>> q1 = u.Quantity(a, u.m, copy=False)
  >>> np.may_share_memory(a, q1)
  True
  >>> q2 = a << u.m
  >>> np.may_share_memory(a, q2)
  True

Conversion that copies only when necessary::

  >>> q3 = q2 << u.m
  >>> np.may_share_memory(q2, q3)
  True
  >>> q3 = q2 << u.cm
  >>> np.may_share_memory(q2, q3)
  False

In-place conversion to a different unit::

  >>> q2 <<= u.cm
  >>> q2  # doctest: +FLOAT_CMP
  <Quantity [  0., 100., 200., 300., 400.] cm>
