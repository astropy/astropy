.. _astropy-visualization:

********************************************
Data Visualization (`astropy.visualization`)
********************************************

Introduction
============

`astropy.visualization` provides functionality that can be helpful when
visualizing data. This includes a framework for plotting Astronomical images
with coordinates with Matplotlib (previously the standalone **wcsaxes**
package), functionality related to image normalization (including both scaling
and stretching), smart histogram plotting, RGB color image creation from
separate images, and custom plotting styles for Matplotlib.

Using `astropy.visualization`
=============================
.. toctree::
   :maxdepth: 2

   matplotlib_integration.rst
   wcsaxes/index.rst
   normalization.rst
   histogram.rst
   rgb.rst

.. _fits2bitmap:

Scripts
=======

This module includes a command-line script, ``fits2bitmap`` to convert FITS
images to bitmaps, including scaling and stretching of the image. To find out
more about the available options and how to use it, type::

    $ fits2bitmap --help

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. automodapi:: astropy.visualization

.. automodapi:: astropy.visualization.mpl_normalize
