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

   wcsaxes/index.rst
   normalization.rst
   histogram.rst
   rgb.rst

Astropy matplotlib style
=========================

The visualization package contains two dictionaries that can be used to
set the Matplotlib plotting style:

.. data:: astropy_mpl_style

    Improves some settings over the matplotlib default style.

.. data:: astropy_mpl_docs_style

    Matplotlib style used by the Astropy documentation.

To apply the custom style on top of your existing matplotlib style,
perform the following:

Using matplotlib version >= 1.5:

.. NOTE:  skip this doctest because the travis-ci py2.7 test uses
.. matplotlib 1.4.3, which does not support "axes.axisbelow"

.. doctest-skip::

    >>> import matplotlib.pyplot as plt
    >>> from astropy.visualization import astropy_mpl_style
    >>> plt.style.use(astropy_mpl_style)

For older versions of matplotlib:

.. doctest-requires:: matplotlib

    >>> import matplotlib as mpl
    >>> from astropy.visualization import astropy_mpl_style
    >>> mpl.rcParams.update(astropy_mpl_style)

Note that these styles are applied *on top* your existing matplotlib
style.  If you want an exactly reproducible plot (i.e. if you want the
plot to come out exactly the same independent of the user
configuration), you should reset the matplotlib settings to the
defaults *before* applying the astropy style.

Using matplotlib version >= 1.5:

.. NOTE:  skip this doctest because the travis-ci py2.7 test uses
.. matplotlib 1.4.3, which does not support "axes.axisbelow"

.. doctest-skip::

    >>> import matplotlib.pyplot as plt
    >>> from astropy.visualization import astropy_mpl_style
    >>> plt.style.use('default')
    >>> plt.style.use(astropy_mpl_style)

For older versions of matplotlib:

.. doctest-requires:: matplotlib

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
