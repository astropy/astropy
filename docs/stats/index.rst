.. _stats:

***************************************
Astrostatistics Tools (`astropy.stats`)
***************************************

Introduction
============

The `astropy.stats` package holds statistical functions or algorithms used
in astronomy and astropy.

Getting Started
===============

Most tools are fairly self-contained, and include relevant examples in
their docstrings.

Using `astropy.stats`
=====================

More detailed information on using the package is provided on separate pages,
listed below.

.. toctree::
   :maxdepth: 2

   lombscargle.rst

Constants
=========

The `astropy.stats` package defines two constants useful for
converting between Gaussian sigma and full width at half maximum
(FWHM):

.. data:: gaussian_sigma_to_fwhm

    Factor with which to multiply Gaussian 1-sigma standard deviation
    to convert it to full width at half maximum (FWHM).

    >>> from astropy.stats import gaussian_sigma_to_fwhm
    >>> gaussian_sigma_to_fwhm
    2.3548200450309493

.. data:: gaussian_fwhm_to_sigma

    Factor with which to multiply Gaussian full width at half maximum
    (FWHM) to convert it to 1-sigma standard deviation.

    >>> from astropy.stats import gaussian_fwhm_to_sigma
    >>> gaussian_fwhm_to_sigma
    0.42466090014400953

See Also
========

* :mod:`scipy.stats`
    This scipy package contains a variety of useful statistical functions and
    classes.  The functionality in `astropy.stats` is intended to supplement
    this, *not* replace it.

* :func:`astropy.visualization.hist`
    The :func:`~astropy.stats.histogram` routine and related functionality
    defined here are used within the :func:`astropy.visualization.hist`
    function. For a discussion of these methods for determining histogram
    binnings, see :ref:`astropy-visualization-hist`.


Reference/API
=============

.. automodapi:: astropy.stats
