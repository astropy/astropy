.. _wcsaxes:

**********************************************
Making plots with coordinate systems (WCSAxes)
**********************************************

WCSAxes is a framework for making plots of Astronomical data in Matplotlib.

Getting started
===============


Using WCSAxes
=============

.. toctree::
   :maxdepth: 1

   getting_started
   ticks_labels_grid
   overlays
   overlaying_coordinate_systems
   slicing_datacubes
   initializing_axes
   controlling_axes
   custom_frames

Notes for developers
====================

This sub-package makes use of image testing with the `pytest-mpl
<https://pypi.python.org/pypi/pytest-mpl/0.6>`_ package. For more information
on writing image tests, see :ref:`image-tests`.

Reference/API
=============

.. automodapi:: astropy.visualization.wcsaxes
   :no-inheritance-diagram:

.. automodapi:: astropy.visualization.wcsaxes.frame
   :no-inheritance-diagram:
