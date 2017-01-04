.. _astropy-visualization:

********************************************
Data Visualization (`astropy.visualization`)
********************************************

Introduction
============

`astropy.visualization` provides functionality that can be helpful when
visualizing data. This includes a framework for plotting Astronomical images
with coordinates with Matplotlib (previously the standalone **wcsaxes**
package), functionality related to image normaliation (including both scaling
and stretching), smart histogram plotting, RGB color image creation from
separate images, and custom plotting styles for Matplotlib.

Using `astropy.visualization`
=============================
.. toctree::
   :maxdepth: 2

   wcsaxes/index.rst
   normalization.rst
   histogram.rst
   lupton_rgb.rst

Astropy matplotlib style
=========================

This module contains dictionaries that can be used to set a matplotlib
plotting style.  It is mostly here to allow a consistent plotting style
in tutorials, but can be used to prepare any matplotlib figure.

Using matplotlib version >= 1.5 you can do::

    >>> import matplotlib.pyplot as plt
    >>> from astropy.visualization import astropy_mpl_style
    >>> plt.style.use(astropy_mpl_style)

For older versions of matplotlib the following works::

    >>> import matplotlib as mpl
    >>> from astropy.visualization import astropy_mpl_style
    >>> mpl.rcParams.update(astropy_mpl_style)

This applies the astropy style on top of your existing matplotlib
default parameters. If you want an exactly reproducible plot (again,
this is useful if you are writing teaching material and you want the
plot to come out exactly the same, independent of the users
configuration for example), you should reset the matplotlib settings to
the library defaults *before* applying the astropy style, e.g.::

    >>> import matplotlib as mpl
    >>> from astropy.visualization import astropy_mpl_style
    >>> mpl.rcdefaults()
    >>> mpl.rcParams.update(astropy_mpl_style)

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
