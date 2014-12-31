.. _astropy-visualization:

********************************************
Data Visualization (`astropy.visualization`)
********************************************

Introduction
============

`astropy.visualization` provides functionality that can be helpful when
visualizing data. At the moment, the main functionality is image normalizing
(including both scaling and stretching).

Using `astropy.visualization`
=============================

.. toctree::
   :maxdepth: 2

   normalization.rst


.. _fits2bitmap:

Scripts
=======

This module includes a command-line script, ``fits2bitmap`` to convert FITS
images to bitmaps, including scaling and stretching of the image. To find out
more about the available options and how to use it, tyoe::

    $ fits2bitmap --help

Reference/API
=============

.. automodapi:: astropy.visualization

.. automodapi:: astropy.visualization.mpl_normalize
