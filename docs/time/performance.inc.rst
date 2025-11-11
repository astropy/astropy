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

Like most of Astropy's classes, |Time| can be array-valued and fully
supports NumPy's `broadcasting <https://numpy.org/doc/stable/user/basics.broadcasting.html>`_ rules.
The best performance will generally be had when making full use of that.
For instance, when one wants to calculate the light travel time for
a large number of sources, rather than loop over them, it is substantially
faster to put all the sources inside a single |SkyCoord| instance and
pass that to :meth:`~astropy.time.Time.light_travel_time`. For a detailed
example, see :ref:`sphx_glr_generated_examples_coordinates_plot_obs-planning.py`.
