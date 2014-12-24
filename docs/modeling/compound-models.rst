.. _compound-models:

Compound Models
===============

.. versionadded:: 1.0

As noted in the :ref:`introduction to the modeling package
<compound-models-intro>`, it is now possible to create new models just by
combining existing models using the arithmetic operators ``+``, ``-``, ``*``,
``/``, and ``**``, as well as by model composition using ``|`` and
concatenation (explained below) with ``&``.


Some terminology
----------------

In discussing the compound model feature, it is useful to be clear about a
few terms where there have been points of confusion:

- The term "model" can refer either to a model *class* or a model *instance*.

  - All models in `astropy.modeling`, whether it represents some
    `function <astropy.modeling.functional_models>`, a
    `rotation <astropy.modeling.rotations>`, etc. is represented in the
    abstract by a model *class*--specifically a subclass of
    `~astropy.modeling.Model`--that encapsulates the routine for evaluating the
    model, a list of its required parameters, and other metadata about the
    model.

  - Per typical object-oriented parlance, a model *instance* is the object
    created when when calling a model class with some arguments--in most cases
    values for the model's parameters.

  A model class, by itself, cannot be used to perform any computation because
  most models, at least, have one or more parameters that must be specified
  before the model can be evaluated on some input data. However, we can still
  get some information about a model class from its representation.  For
  example::

      >>> from astropy.modeling.models import Gaussian1D
      >>> Gaussian1D
      <class 'astropy.modeling.functional_models.Gaussian1D'>
      Name: Gaussian1D
      Inputs: ('x',)
      Outputs: ('y',)
      Fittable parameters: ('amplitude', 'mean', 'stddev')

  We can then create a model *instance* by passing in values for the three
  parameters::

      >>> my_gaussian = Gaussian1D(amplitude=1.0, mean=0, stddev=0.2)
      >>> my_gaussian
      <Gaussian1D(amplitude=1.0, mean=0.0, stddev=0.2)>

  We now have an *instance* of `~astropy.modeling.functional_models.Gaussian1D`
  with all its parameters (and in principle other details like fit constraints)
  filled in so that we can perform calculations with it as though it were a
  function::

      >>> my_gaussian(0.2)  # doctest: +FLOAT_CMP
      0.6065306597126334

  In many cases this document just refers to "models", where the class/instance
  distinction is either irrelevant or clear from context.  But a distinction
  will be made where necessary.

- A *compound model* can be created by combining two or more existing models--
  be they model *instances* or *classes*, and can be models that come with
  Astropy, :doc:`user defined models <new>`, or other compound models--using
  Python expressions consisting of one or more of the supported binary
  operators.

- In some places the term *composite model* is used interchangeably with
  *compound model*.  This can be seen in the cases of the now deprecated
  `~astropy.modeling.SerialCompositeModel` and
  `~astropy.modeling.SummedCompositeModel`.  However, this document uses the
  term *composite model* to refer *only* to the case of a compound model
  created from the functional composition of two or more models using the pipe
  operator ``|`` as explained below.  This distinction is used consistently
  within this document, but it may be helpful to understand the distinction.
