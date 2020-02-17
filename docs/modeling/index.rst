.. include:: links.inc

.. _astropy-modeling:

***************************************
Models and Fitting (`astropy.modeling`)
***************************************

Introduction
============

`astropy.modeling` provides a framework for representing models and performing
model evaluation and fitting.  A number of predefined 1-D and 2-D models are
provided and the capability for custom, user defined models is supported.
Different fitting algorithms can be used with any model.  For those fitters
with the capabilities fitting can be done using uncertainties, parameters with
bounds, and priors.

.. note::

    A number of significant changes have been made to the internals that have been
    documented in more detail in :doc:`changes_for_4`. The main change is that
    combining model classes no longer is supported. (Combining model
    instances is still very much supported!)

.. _modeling-using:

Using Modeling
==============

.. toctree::
   :maxdepth: 2

   Models <models.rst>
   Compound Models <compound-models.rst>
   Model Parameters <parameters.rst>
   Fitting <fitting.rst>
   Using Units with Models and Fitting <units.rst>
   Changes in v4.0 <changes_for_4.rst>


.. _getting-started-example:

A Simple Example
================

This simple example illustrates defining a model,
calculating values based on input x values, and using fitting data with a model.

   .. plot::
       :include-source:

       import numpy as np
       import matplotlib.pyplot as plt
       from astropy.modeling import models, fitting

       # define a model for a line
       line_orig = models.Linear1D(slope=1.0, intercept=0.5)

       # generate x, y data non-uniformly spaced in x
       # add noise to y measurements
       npts = 30
       np.random.seed(10)
       x = np.random.uniform(0.0, 10.0, npts)
       y = line_orig(x)
       y += np.random.normal(0.0, 1.5, npts)

       # initialize a linear fitter
       fit = fitting.LinearLSQFitter()

       # initialize a linear model
       line_init = models.Linear1D()

       # fit the data with the fitter
       fitted_line = fit(line_init, x, y)

       # plot the model
       plt.figure()
       plt.plot(x, y, 'ko', label='Data')
       plt.plot(x, fitted_line(x), 'k-', label='Fitted Model')
       plt.xlabel('x')
       plt.ylabel('y')
       plt.legend()

.. _advanced_topics:

Advanced Topics
===============

.. toctree::
   :maxdepth: 2

   Performance Tips <performance.rst>
   Extending Models <new-model.rst>
   Extending Fitters <new-fitter.rst>
   Adding support for units to models <add-units.rst>


Pre-Defined Models
==================

.. To be expanded to include all pre-defined models

Some of the pre-defined models are listed and illustrated.

.. toctree::
   :maxdepth: 2

   1D Models <predef_models1D.rst>
   2D Models <predef_models2D.rst>
   Physical Models <physical_models.rst>
   Polynomial Models <polynomial_models.rst>

Examples
========

.. toctree::
   :maxdepth: 2

   Fitting a line <example-fitting-line>
   example-fitting-constraints
   example-fitting-model-sets

.. TODO list
    fitting with masks
    fitting with priors
    fitting with units
    defining 1d model
    defining 2d model
    fitting 2d model
    defining and using a WCS/gWCS model
    defining and using a Tabular1D model
    statistics functions and how to make your own
    compound models


Reference/API
=============

.. toctree::
   :maxdepth: 1

   reference_api
