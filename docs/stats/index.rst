.. _stats:

***************************************
Astrostatistics Tools (`astropy.stats`)
***************************************

Introduction
============

The `astropy.stats` package holds statistical functions or algorithms
used in astronomy.  While the `scipy.stats` and `statsmodels
<http://www.statsmodels.org/stable/index.html>`_ packages contains a
wide range of statistical tools, they are general-purpose packages and
are missing some tools that are particularly useful or specific to
astronomy. This package is intended to provide such functionality,
but *not* to replace `scipy.stats` if its implementation satisfies
astronomers' needs.


Getting Started
===============

A number of different tools are contained in the stats package, and
they can be accessed by importing them::

    >>> from astropy import stats

A full list of the different tools are provided below. Please see the
documentation for their different usages. For example, sigma clipping,
which is a common way to estimate the background of an image, can be
performed with the :func:`~astropy.stats.sigma_clip` function. By
default, the function returns a masked array where outliers are
masked.

Examples
--------

..
  EXAMPLE START
  Sigma Clipping with Astropy Stats sigma_clip Function

To estimate the background of an image::

    >>> data = [1, 5, 6, 8, 100, 5, 3, 2]
    >>> stats.sigma_clip(data, sigma=2, maxiters=5)
    masked_array(data=[1, 5, 6, 8, --, 5, 3, 2],
                 mask=[False, False, False, False,  True, False, False, False],
           fill_value=999999)

..
  EXAMPLE END

..
  EXAMPLE START
  Sigma Clipping with Astropy Stats SigmaClip Class

Alternatively, the :class:`~astropy.stats.SigmaClip` class provides an
object-oriented interface to sigma clipping, which also returns a
masked array by default::

    >>> sigclip = stats.SigmaClip(sigma=2, maxiters=5)
    >>> sigclip(data)
    masked_array(data=[1, 5, 6, 8, --, 5, 3, 2],
                 mask=[False, False, False, False,  True, False, False, False],
           fill_value=999999)

..
  EXAMPLE END

..
  EXAMPLE START
  Calculating Sigma Clipping Statistics

In addition, there are also several convenience functions for making
the calculation of statistics even more convenient. For example,
:func:`~astropy.stats.sigma_clipped_stats` will return the mean,
median, and standard deviation of a sigma-clipped array::

     >>> stats.sigma_clipped_stats(data, sigma=2, maxiters=5)  # doctest: +FLOAT_CMP
     (4.2857142857142856, 5.0, 2.2497165354319457)

There are also tools for calculating :ref:`robust statistics
<stats-robust>`, sampling the data, :ref:`circular statistics
<stats-circular>`, confidence limits, spatial statistics, and adaptive
histograms.

..
  EXAMPLE END

Most tools are fairly self-contained, and include relevant examples in
their docstrings.


Using `astropy.stats`
=====================

More detailed information on using the package is provided on separate pages,
listed below.

.. toctree::
   :maxdepth: 2

   robust.rst
   circ.rst
   ripley.rst
   ../visualization/histogram.rst


Constants
=========

The `astropy.stats` package defines two constants useful for
converting between Gaussian sigma and full width at half maximum
(FWHM):

.. data:: gaussian_sigma_to_fwhm

    Factor with which to multiply Gaussian 1-sigma standard deviation
    to convert it to full width at half maximum (FWHM).

    >>> from astropy.stats import gaussian_sigma_to_fwhm
    >>> gaussian_sigma_to_fwhm  # doctest: +FLOAT_CMP
    2.3548200450309493

.. data:: gaussian_fwhm_to_sigma

    Factor with which to multiply Gaussian full width at half maximum
    (FWHM) to convert it to 1-sigma standard deviation.

    >>> from astropy.stats import gaussian_fwhm_to_sigma
    >>> gaussian_fwhm_to_sigma  # doctest: +FLOAT_CMP
    0.42466090014400953


See Also
========

* :mod:`scipy.stats`
    This SciPy package contains a variety of useful statistical functions
    and classes. The functionality in `astropy.stats` is intended to supplement
    this, *not* replace it.

* `statsmodels <http://www.statsmodels.org/stable/index.html>`_
    The statsmodels package provides functionality for estimating
    different statistical models, tests, and data exploration.

* `astroML <https://www.astroml.org/>`_
    The astroML package is a Python module for machine learning and
    data mining. Some of the tools from this package have been
    migrated here, but there are still a number of tools there that
    are useful for astronomy and statistical analysis.


* :func:`astropy.visualization.hist`
    The :func:`~astropy.stats.histogram` routine and related functionality
    defined here are used within the :func:`astropy.visualization.hist`
    function. For a discussion of these methods for determining histogram
    binnings, see :ref:`astropy-visualization-hist`.

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. automodapi:: astropy.stats
