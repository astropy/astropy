.. _stats:

***************************************
Astrostatistics Tools (`astropy.stats`)
***************************************

Introduction
============

The `astropy.stats` package holds statistical functions or algorithms
used in astronomy.  While the `scipy.stats` and `statsmodels
<https://www.statsmodels.org/stable/index.html>`_ packages contains a
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
performed with the :func:`~astropy.stats.sigma_clip` function.
By default, the function returns a masked array, a type of Numpy array
used for handling missing or invalid entries.  Masked arrays retain the
original data but also store another boolean array of the same shape
where ``True`` indicates that the value is masked. Most Numpy ufuncs
will understand masked arrays and treat them appropriately.
For example, consider the following dataset with a clear outlier::

    >>> import numpy as np
    >>> from astropy.stats import sigma_clip
    >>> x = np.array([1, 0, 0, 1, 99, 0, 0, 1, 0])

The mean is skewed by the outlier::

    >>> x.mean()
    np.float64(11.333333333333334)

Sigma-clipping (3 sigma by default) returns a masked array,
and so functions like ``mean`` will ignore the outlier::

    >>> clipped = sigma_clip(x)
    >>> clipped
    masked_array(data=[1, 0, 0, 1, --, 0, 0, 1, 0],
                 mask=[False, False, False, False,  True, False, False, False,
                       False],
           fill_value=999999)
    >>> clipped.mean()
    np.float64(0.375)

If you need to access the original data directly, you can use the
``data`` property. Combined with the ``mask`` property, you can get the
original outliers, or the values that were not clipped::

    >>> outliers = clipped.data[clipped.mask]
    >>> outliers
    array([99])
    >>> valid = clipped.data[~clipped.mask]
    >>> valid
    array([1, 0, 0, 1, 0, 0, 1, 0])

For more information on masked arrays, including see the
:ref:`numpy.ma module <numpy:maskedarray.generic>`.

Examples
--------

..
  EXAMPLE START
  Sigma Clipping with Astropy Stats sigma_clip Function

To estimate the background of an image::

    >>> data = [1, 5, 6, 8, 100, 5, 3, 2]
    >>> data_clipped = stats.sigma_clip(data, sigma=2, maxiters=5)
    >>> data_clipped
    masked_array(data=[1, 5, 6, 8, --, 5, 3, 2],
                 mask=[False, False, False, False,  True, False, False, False],
           fill_value=999999)
    >>> np.mean(data_clipped)  # doctest: +FLOAT_CMP
    np.float64(4.285714285714286)

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
     (np.float64(4.285714285714286), np.float64(5.0), np.float64(2.249716535431946))

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

Also see :ref:`astropy-visualization-hist`.


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

* `statsmodels <https://www.statsmodels.org/stable/index.html>`_
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

.. toctree::
   :maxdepth: 2

   ref_api
