.. include:: links.inc

.. _astropy-modeling:

***************************************
Models and Fitting (`astropy.modeling`)
***************************************

Introduction
============

`astropy.modeling` provides a framework for representing models and performing
model evaluation and fitting. The goal of this package is to eventually provide
a rich toolset of models and fitters such that most users will not need to
define new model classes, nor special purpose fitting routines (while making it
reasonably easy to do when necessary). See a
:ref:`quick tutorial <modeling_tutorial>` of the package.

The major components of the package are:

Models
------

Models behave like parametrized functions that can be easily evaluated on
datasets. For example, the `~astropy.modeling.functional_models.Gaussian1D`
model evaluates the Gaussian function
:math:`f(x) = A e^{- \frac{\left(x - x_{0}\right)^{2}}{2 \sigma^{2}}}`.
Single values or arrays can be passed to models to evaluate the function at
these points. Examples of using models can be found in the
:ref:`modeling tutorial <modeling_tutorial>`.

Detailed documentation on models can be found :doc:`here <models>`.
Additional features of models are:

Advanced Model Topics
---------------------

* **Compound models**

  Compound models can be created by combining existing models with
  arithmetic expressions. A quick tutorial of compound models can be found
  :ref:`here <compound-models-intro>`. For more detailed documentation on
  compound models, see :ref:`compound-models`.

* **Model Sets**

  Users can also create a model *set* which is a collection of models of the same
  type, but with different values for their parameters. See
  :ref:`here <model-sets-index>` for a quick tutorial on model sets.

* **Custom models**

  The package is intended to allow users to easily create their own models.
  See documentation :ref:`here<modeling-new-classes>`.

* **Algorithms**

  A description of the algorithms used to evaluate the 1-D polynomial,
  Chebyshev and Legendre polynomial models can be found
  :doc:`here <algorithms>`.

Parameters
----------

This refers to the parameters of a model which are implemented as attributes of
models. For example, the parameters of a ``Gaussian1D`` model are
``amplitude``, ``x_0`` and ``stddev`` which refer to the amplitude, the
position of the center of the peak, and the standard deviation of the function
respectively. Parameters may be single values or arrays. Features of parameters
include:

* **Accessing parameters**

  Parameters can be set when instantiating a model, and can later be easily
  accessed as model attributes as shown :ref:`here <setting-params>`.

* **Constraining parameters**

  Users can constrain parameters to be within a certain range or tie them with
  other parameters using custom functions when fitting models to data. See
  :ref:`here <constraining_parameters>` to view examples of how to modify and
  constrain parameters.

Fitting
-------

Fitters are classes that combine optimization algorithms with statistic
functions to fit user data to models. They can be called as functions by
passing in user data, model data and a model instance. See detailed
documentation on how to fit data to models :doc:`here <fitting>`. Additional
features regarding fitters include:

* **Creating custom fitters**

  Users have the ability to create their own fitters by defining an error
  function and the type of constraints that work with their fitter (if any).
  See documentation :ref:`here<custom-fitter>`.

* **Using custom statistic functions**

  Users can also implement their own statistic functions to be used with
  their custom fitters. See documentation :ref:`here<custom-statistic>`.

Reference/API
=============

 A detailed description of the classes and functions available in the package
 can be found :doc:`here <api>`.

.. note::

    `astropy.modeling` is currently a work-in-progress, and thus it is likely
    there will still be API changes in later versions of Astropy. Backwards
    compatibility support between versions will still be maintained as much as
    possible, but new features and enhancements are coming in future versions.
    If you have specific ideas for how it might be improved, feel free to let
    us know on the `astropy-dev mailing list`_ or at
    http://feedback.astropy.org

`astropy.modeling` documentation
================================

.. toctree::
   :maxdepth: 1

   quick-tutorial
   models
   parameters
   fitting
   compound-models
   new
   bounding-boxes
   algorithms
   api
