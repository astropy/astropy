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
