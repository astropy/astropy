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
      >>> my_gaussian  # doctest: +FLOAT_CMP
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


Creating compound models
------------------------

As discussed in the :ref:`introduction to compound models
<compound-models-intro>`, the only way, currently, to create compound models is
to combine existing single models and/or compound models using expressions in
Python with the binary operators ``+``, ``-``, ``*``, ``/``, ``**``, ``|``,
and ``&``, each of which are discussed in the following sections.  The operands
used in these expressions may be model *classes*, or model *instances*.  In
other words, any object for which either ``isinstance(obj, Model)`` or
``issubclass(obj, Model)`` is `True`.

When all models involved in the expression are classes, the result of the
expression is, itself, a class (remember, classes in Python are themselves also
objects just like strings and integers or model instances)::

    >>> TwoGaussians = Gaussian1D + Gaussian1D
    >>> from astropy.modeling import Model
    >>> isinstance(TwoGaussians, Model)
    False
    >>> issubclass(TwoGaussians, Model)
    True

When we inspect the variable ``TwoGaussians`` by printing its representation at
the command prompt we can get some more information about it::

    >>> TwoGaussians
    <class '__main__.CompoundModel...'>
    Name: CompoundModel...
    Inputs: ('x',)
    Outputs: ('y',)
    Fittable parameters: ('amplitude_0', 'mean_0', 'stddev_0', 'amplitude_1', 'mean_1', 'stddev_1')
    Expression: [0] + [1]
    Components: 
        [0]: <class 'astropy.modeling.functional_models.Gaussian1D'>
        Name: Gaussian1D
        Inputs: ('x',)
        Outputs: ('y',)
        Fittable parameters: ('amplitude', 'mean', 'stddev')
    <BLANKLINE>
        [1]: <class 'astropy.modeling.functional_models.Gaussian1D'>
        Name: Gaussian1D
        Inputs: ('x',)
        Outputs: ('y',)
        Fittable parameters: ('amplitude', 'mean', 'stddev')

There are a number of things to point out here:  This model class has six
fittable parameters.  How parameters are handled is discussed further in the
section on :ref:`compound-model-parameters`.  We also see that there is a
listing of the *expression* that was used to create this compound model, which
in this case is summarized as ``[0] + [1]``.  The ``[0]`` and ``[1]`` refer to
the first and second components of the model listed next (in this case both
components are the `~astropy.modeling.functional_models.Gaussian1D` class).

Each component of a compound model is a single, non-compound model.  This is
the case even when performing an existing compound model in a new expression.
The existing compound model is not treated as a single model--instead the
expression represented by that compound model is extended.  An expression
involving two or more compound models results in a new expression that is the
concatenation of all involved models' expressions::

    >>> FourGaussians = TwoGaussians + TwoGaussians
    >>> FourGaussians
    <class '__main__.CompoundModel...'>
    Name: CompoundModel...
    Inputs: ('x',)
    Outputs: ('y',)
    Fittable parameters: ('amplitude_0', 'mean_0', 'stddev_0', ..., 'amplitude_3', 'mean_3', 'stddev_3')
    Expression: [0] + [1] + [2] + [3]
    Components: 
        [0]: <class 'astropy.modeling.functional_models.Gaussian1D'>
        Name: Gaussian1D
        Inputs: ('x',)
        Outputs: ('y',)
        Fittable parameters: ('amplitude', 'mean', 'stddev')
        ...
        [3]: <class 'astropy.modeling.functional_models.Gaussian1D'>
        Name: Gaussian1D
        Inputs: ('x',)
        Outputs: ('y',)
        Fittable parameters: ('amplitude', 'mean', 'stddev')

In a future version it may be possible to "freeze" a compound model, so that
from the user's perspective it is treated as a single model.  However, as this
is the default behavior it is good to be aware of.


Model names
^^^^^^^^^^^

In the last two examples another notable feature of the generated compound
model classes is that the class name, as displayed when printing the class at
the command prompt, is not "TwoGaussians", "FourGaussians", etc.  Instead it is
a generated name consisting of "CompoundModel" followed by an essentially
arbitrary integer that is chosen simply so that every compound model has a
unique default name.  This is a limitation at present, due to the limitation
that it is not generally possible in Python when an object is created by an
expression for it to "know" the name of the variable it will be assigned to, if
any.  It may be possible in the future to work around this in limited cases,
but for now there are a couple workarounds for creating compound model classes
with friendlier names.  The first is to use the
`Model.rename <astropy.modeling.Model.rename>` class method on the result of
the model expression::

    >>> TwoGaussians = (Gaussian1D + Gaussian1D).rename('TwoGaussians')
    >>> TwoGaussians
    <class '__main__.TwoGaussians'>
    Name: TwoGaussians (CompoundModel...)
    ...

This actually takes the generated compound model and creates a light subclass
of it with the desired name.  This does not impose any additional overhead.  An
alternative syntax, which is equivalent to what
`~astropy.modeling.Model.rename` is doing, is to directly use the model
expression as the base class of a new class::

    >>> class TwoGaussians(Gaussian1D + Gaussian1D):
    ...     """A superposition of two Gaussians."""
    ...
    >>> TwoGaussians
    <class '__main__.TwoGaussians'>
    Name: TwoGaussians (CompoundModel...)
    ...

Because the result of the expression ``Gaussian1D + Gaussian1D`` *is* a class,
it can be used directly in the standard class declaration syntax
``class ClassName(Base):`` as the base.  This syntax also has the advantage of
allowing a docstring to be assigned to the new class.  In future versions it
may be possible to customize other aspects of compound model classes in this
way.  Single model classes can also be given custom names by using
`~astropy.modeling.Model.rename`, and model instances can be given names as
well.  This can be used to good effect, for example as shown in the section on
:ref:`compound-model-indexing`.


Compound models with model instances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

So far we have seen how to create compound model *classes* from expressions
involving other model classes.  This is the most "generic" way to create new
models from existing models.  However, many may find it more useful most of the
time, especially when providing an initial guess to a fitter, to create a new
model from a combination of model *instances* with already defined parameter
values.  This can also be done and works mostly the same way::

    >>> both_gaussians = Gaussian1D(1, 0, 0.2) + Gaussian1D(2.5, 0.5, 0.1)
    >>> both_gaussians  # doctest: +FLOAT_CMP
    <CompoundModel...(amplitude_0=1.0, mean_0=0.0, stddev_0=0.2, amplitude_1=2.5, mean_1=0.5, stddev_1=0.1)>

Unlike when a model was created from model classes, this expression does not
directly return a new class; instead it creates a model instance that is ready
to be used for evaluation::

    >>> both_gaussians(0.2)  # doctest: +FLOAT_CMP
    0.6343031510582392

This was found to be much more convenient and natural, in this case, than
returning a class.  It is worth understanding that the way this works under the
hood is to create the compound class, and then immediately instantiate it with
the already known parameter values.  We can see this by checking the type of
``both_gaussians``::

    >>> type(both_gaussians)
    <class '__main__.CompoundModel...'>
    Name: CompoundModel...
    Inputs: ('x',)
    Outputs: ('y',)
    Fittable parameters: ('amplitude_0', 'mean_0', 'stddev_0', 'amplitude_1', 'mean_1', 'stddev_1')
    Expression: [0] + [1]
    Components: 
        [0]: <Gaussian1D(amplitude=1.0, mean=0.0, stddev=0.2)>
    <BLANKLINE>
        [1]: <Gaussian1D(amplitude=2.5, mean=0.5, stddev=0.1)>

It is also possible, and sometimes useful, to make a compound model from a
combinatation of classes *and* instances in the same expression::

    >>> from astropy.modeling.models import Linear1D, Sine1D
    >>> MyModel = Linear1D + Sine1D(amplitude=1, frequency=1)
    >>> MyModel
    <class '__main__.CompoundModel...'>
    Name: CompoundModel...
    Inputs: ('x',)
    Outputs: ('y',)
    Fittable parameters: ('slope_0', 'intercept_0', 'amplitude_1', 'frequency_1')
    Expression: [0] + [1]
    Components: 
        [0]: <class 'astropy.modeling.functional_models.Linear1D'>
        Name: Linear1D
        Inputs: ('x',)
        Outputs: ('y',)
        Fittable parameters: ('slope', 'intercept')
    <BLANKLINE>
        [1]: <Sine1D(amplitude=1.0, frequency=1.0)>

In this case the result is always a class.  However (and this is not
immediately obvious by the representation) the difference is that the
``amplitude`` and ``frequency`` parameters for the
`~astropy.modeling.functional_models.Sine1D` part of the model are
"baked into" the class as default values for those parameters.  So it is
possible to instantiate one of these models by specifying just the ``slope``
and ``intercept`` parameters for the
`~astropy.modeling.functional_models.Linear1D` part of the model::

    >>> my_model = MyModel(1, 0)
    >>> my_model(0.25)  # doctest +FLOAT_CMP
    1.25

This does not prevent the other parameters from being overridden, however::

    >>> my_model = MyModel(slope_0=1, intercept_0=0, frequency_1=2)
    >>> my_model(0.125)  # doctest +FLOAT_CMP
    1.125

In fact, this is currently the only way to use a `polynomial
<astropy.modeling.polynomial>` model in a compound model, because the design of
the polynomial models is currently such that they must be instantiated in order
to specify their polynomial degree.  Because the polynomials are already
designed so that their coefficients all default to zero, this "limitation"
should not have any practical drawbacks.

.. note::

    There is currently a caveat in the example of combining model classes and
    instances, which is that the parameter values of model *instances* are only
    treated as defaults if the expression is written in such a way that all
    model instances are to the right of all model classes.  This limitation
    will be lifted in a later version--in particular, Python 3 offers a lot
    more flexibility with respect to how function arguments are handled.


Operators
---------

Arithmetic operators
^^^^^^^^^^^^^^^^^^^^

Compound models can be created from expressions that include any number of the
arithmetic operators ``+``, ``-``, ``*``, ``/``, and ``**`` which have the same
meanings as they do for other numeric objects in Python.

.. note::

    In the case of division ``/`` always means floating point division--integer
    division and the ``//`` operator are not supported for models).

As demonstrated in previous examples, for models that have a single output
the result of evaluating a model like ``A + B`` is to evaluate ``A`` and
``B`` separately on the given input, and then return the sum of the outputs of
``A`` and ``B``.  This requires that ``A`` and ``B`` take the same number of
inputs and both have a single output.

It is also possible to use arithmetic operators between models with multiple
outputs.  Again, the number of inputs must be the same between the models, as
must be the number of outputs.  In this case the operator is applied to the
operators element-wise, similarly to how arithmetic operators work on two Numpy
arrays.


Model composition
^^^^^^^^^^^^^^^^^

The sixth binary operator that can be used to create compound models is the
composition operator, also known as the "pipe" operator ``|`` (not to be
confused with the boolean "or" operator that this implements for Python numeric
objects).  A model created with the composition operator like ``M = F | G``,
when evaluated, is equivalent to evaluating :math:`g \circ f = g(f(x))`.

.. note::

    The fact that the ``|`` operator has the opposite sense as the functional
    composition operator :math:`\circ` is sometimes a point of confusion.
    This is in part because there is no operator symbol supported in Python
    that corresponds well to this.  The ``|`` operator should instead be read
    like the `pipe operator
    <http://en.wikipedia.org/wiki/Pipeline_%28Unix%29>`_ of UNIX shell syntax:
    It chains together models by piping the output of the left-hand operand to
    the input of the right-hand operand, forming a "pipeline" of models, or
    transformations.

This has different requirements on the inputs/outputs of its operands than do
the arithmetic operators.  For composition all that is required is that the
left-hand model has the same number of outputs as the right-hand model has
inputs.

For simple functional models this is exactly the same as functional
composition, except for the aforementioned caveat about ordering.  For
example:

.. plot::
    :include-source:

    import numpy as np
    from astropy.modeling.models import Redshift, Gaussian1D

    class RedshiftedGaussian(Redshift | Gaussian1D(1, 0.75, 0.1)):
        """Evaluates a Gaussian with optional redshift applied to the input."""

    x = np.linspace(0, 1.2, 100)
    g0 = RedshiftedGaussian(z_0=0)

    plt.plot(x, g0(x), 'g--', lw=2, label='$z=0$')

    for z in (0.2, 0.4, 0.6):
        g = RedshiftedGaussian(z_0=z)
        plt.plot(x, g(x), color=plt.cm.OrRd(z), lw=2,
                 label='$z={0}$'.format(z))

    plt.xlabel('Energy')
    plt.ylabel('Flux')
    plt.legend()


TODO: Example of composite model with rotations.

Model concatenation
^^^^^^^^^^^^^^^^^^^

TODO

.. _compound-model-parameters:

Parameters
----------

TODO

.. _compound-model-indexing:

Indexing and slicing
--------------------

TODO

Advanced mappings
-----------------

TODO
