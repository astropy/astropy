.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-time-performance:

Performance Tips
================

Here we provide some tips and tricks for how to optimize performance of code
using `astropy.time`.

Broadcasting
------------

Like most of Astropy's classes, |Time| can be array-valued and fully supports
NumPy's `numpy broadcasting <https://numpy.org/doc/stable/user/basics.broadcasting.html>`_.
The best performance is generally achieved by making full use of broadcasting.
For example, when calculating light travel times for many sources, it is much
faster to group all coordinates into a single |SkyCoord| array and call
:meth:`~astropy.time.Time.light_travel_time` once, rather than looping over
individual coordinates.
