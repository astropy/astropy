.. _astropy-visualization:

********************************************
Data Visualization (`astropy.visualization`)
********************************************

Introduction
============

`astropy.visualization` provides functionality that can be helpful when
visualizing data.
At the moment, functionalities include enhanced histograms,
image normalizing (including both scaling and stretching),
and custom plotting styles for matplotlib.

Using `astropy.visualization`
=============================

.. toctree::
   :maxdepth: 2

   normalization.rst
   histogram.rst


.. _fits2bitmap:

Scripts
=======

This module includes a command-line script, ``fits2bitmap`` to convert FITS
images to bitmaps, including scaling and stretching of the image. To find out
more about the available options and how to use it, type::

    $ fits2bitmap --help

Reference/API
=============

.. automodapi:: astropy.visualization.mpl_style

.. automodapi:: astropy.visualization

.. automodapi:: astropy.visualization.mpl_normalize
