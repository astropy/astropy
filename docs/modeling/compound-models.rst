.. include:: links.inc

.. _compound-models-intro:

Combining Models
****************

Basics
======

While the Astropy modeling package makes it very easy to define :doc:`new
models <new-model>` either from existing functions, or by writing a
`~astropy.modeling.Model` subclass, an additional way to create new models is
by combining them using arithmetic expressions.  This works with models built
into Astropy, and most user-defined models as well.  For example, it is
possible to create a superposition of two Gaussians like so::

    >>> from astropy.modeling import models
    >>> g1 = models.Gaussian1D(1, 0, 0.2)
    >>> g2 = models.Gaussian1D(2.5, 0.5, 0.1)
    >>> g1_plus_2 = g1 + g2

The resulting object ``g1_plus_2`` is itself a new model.  Evaluating, say,
``g1_plus_2(0.25)`` is the same as evaluating ``g1(0.25) + g2(0.25)``::

    >>> g1_plus_2(0.25)  # doctest: +FLOAT_CMP
    0.5676756958301329
    >>> g1_plus_2(0.25) == g1(0.25) + g2(0.25)
    True

This model can be further combined with other models in new expressions.

These new compound models can also be fitted to data, like most other models
(though this currently requires one of the non-linear fitters):

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling import models, fitting

    # Generate fake data
    np.random.seed(42)
    g1 = models.Gaussian1D(1, 0, 0.2)
    g2 = models.Gaussian1D(2.5, 0.5, 0.1)
    x = np.linspace(-1, 1, 200)
    y = g1(x) + g2(x) + np.random.normal(0., 0.2, x.shape)

    # Now to fit the data create a new superposition with initial
    # guesses for the parameters:
    gg_init = models.Gaussian1D(1, 0, 0.1) + models.Gaussian1D(2, 0.5, 0.1)
    fitter = fitting.SLSQPLSQFitter()
    gg_fit = fitter(gg_init, x, y)

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,5))
    plt.plot(x, y, 'ko')
    plt.plot(x, gg_fit(x))
    plt.xlabel('Position')
    plt.ylabel('Flux')

This works for 1-D models, 2-D models, and combinations thereof, though there
are some complexities involved in correctly matching up the inputs and outputs
of all models used to build a compound model.  You can learn more details in
the :doc:`compound-models` documentation.

Astropy models also support convolution through the function
`~astropy.convolution.convolve_models`, which returns a compound model.

For instance, the convolution of two Gaussian functions is also a Gaussian
function in which the resulting mean (variance) is the sum of the means
(variances) of each Gaussian.

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling import models
    from astropy.convolution import convolve_models

    g1 = models.Gaussian1D(1, -1, 1)
    g2 = models.Gaussian1D(1, 1, 1)
    g3 = convolve_models(g1, g2)

    x = np.linspace(-3, 3, 50)
    plt.plot(x, g1(x), 'k-')
    plt.plot(x, g2(x), 'k-')
    plt.plot(x, g3(x), 'k-')


.. _compound-models:

A comprehensive description
===========================

Some terminology
----------------

It is possible to create new models just by
combining existing models using the arithmetic operators ``+``, ``-``, ``*``,
``/``, and ``**``, or by model composition using ``|`` and
concatenation (explained below) with ``&``, as well as using :func:`~astropy.modeling.fix_inputs`
for :ref:`reducing the number of inputs to a model <model-reduction>`.


In discussing the compound model feature, it is useful to be clear about a
few terms where there have been points of confusion:

- The term "model" can refer either to a model *class* or a model *instance*.

  - All models in `astropy.modeling`, whether it represents some
    `function <astropy.modeling.functional_models>`, a
    `rotation <astropy.modeling.rotations>`, etc., are represented in the
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
      N_inputs: 1
      N_outputs: 1
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

- A *compound model* can be created by combining two or more existing model instances
  which can be models that come with Astropy, :doc:`user defined models <new-model>`, or
  other compound models--using Python expressions consisting of one or more of the
  supported binary operators. The combination of model classes is deprecated and will
  be removed in version 4.0.

- In some places the term *composite model* is used interchangeably with
  *compound model*. However, this document uses the
  term *composite model* to refer *only* to the case of a compound model
  created from the functional composition of two or more models using the pipe
  operator ``|`` as explained below.  This distinction is used consistently
  within this document, but it may be helpful to understand the distinction.


Creating compound models
------------------------

The only way to create compound models is
to combine existing single models and/or compound models using expressions in
Python with the binary operators ``+``, ``-``, ``*``, ``/``, ``**``, ``|``,
and ``&``, each of which is discussed in the following sections.


.. warning:: Creating compound models by combining classes was removed in v4.0.

The result of combining two models is a model instance::

    >>> two_gaussians = Gaussian1D(1.1, 0.1, 0.2) + Gaussian1D(2.5, 0.5, 0.1)
    >>> two_gaussians  # doctest: +FLOAT_CMP
    <CompoundModel...(amplitude_0=1.1, mean_0=0.1, stddev_0=0.2, amplitude_1=2.5, mean_1=0.5, stddev_1=0.1)>

This expression creates a new model instance that is ready to be used for evaluation::

    >>> two_gaussians(0.2)  # doctest: +FLOAT_CMP
    0.9985190841886609

The ``print`` function provides more information about this object::

    >>> print(two_gaussians)
    Model: CompoundModel...
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Expression: [0] + [1]
    Components:
        [0]: <Gaussian1D(amplitude=1.1, mean=0.1, stddev=0.2)>
    <BLANKLINE>
        [1]: <Gaussian1D(amplitude=2.5, mean=0.5, stddev=0.1)>
    Parameters:
        amplitude_0 mean_0 stddev_0 amplitude_1 mean_1 stddev_1
        ----------- ------ -------- ----------- ------ --------
                1.1    0.1      0.2         2.5    0.5      0.1

There are a number of things to point out here:  This model has six
fittable parameters. How parameters are handled is discussed further in the
section on :ref:`compound-model-parameters`.  We also see that there is a
listing of the *expression* that was used to create this compound model, which
in this case is summarized as ``[0] + [1]``.  The ``[0]`` and ``[1]`` refer to
the first and second components of the model listed next (in this case both
components are the `~astropy.modeling.functional_models.Gaussian1D` objects).

Each component of a compound model is a single, non-compound model.  This is
the case even when including an existing compound model in a new expression.
The existing compound model is not treated as a single model--instead the
expression represented by that compound model is extended.  An expression
involving two or more compound models results in a new expression that is the
concatenation of all involved models' expressions::

    >>> four_gaussians = two_gaussians + two_gaussians
    >>> print(four_gaussians)
    Model: CompoundModel...
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Expression: [0] + [1] + [2] + [3]
    Components:
        [0]: <Gaussian1D(amplitude=1.1, mean=0.1, stddev=0.2)>
    <BLANKLINE>
        [1]: <Gaussian1D(amplitude=2.5, mean=0.5, stddev=0.1)>
    <BLANKLINE>
        [2]: <Gaussian1D(amplitude=1.1, mean=0.1, stddev=0.2)>
    <BLANKLINE>
        [3]: <Gaussian1D(amplitude=2.5, mean=0.5, stddev=0.1)>
    Parameters:
        amplitude_0 mean_0 stddev_0 amplitude_1 ... stddev_2 amplitude_3 mean_3 stddev_3
        ----------- ------ -------- ----------- ... -------- ----------- ------ --------
                1.1    0.1      0.2         2.5 ...      0.2         2.5    0.5      0.1


Operators
---------

Arithmetic operators
--------------------

Compound models can be created from expressions that include any
number of the arithmetic operators ``+``, ``-``, ``*``, ``/``, and
``**``, which have the same meanings as they do for other numeric
objects in Python.

.. note::

    In the case of division ``/`` always means floating point division--integer
    division and the ``//`` operator is not supported for models).

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


.. _compound-model-composition:

Model composition
-----------------

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
    <https://en.wikipedia.org/wiki/Pipeline_%28Unix%29>`_ of UNIX shell syntax:
    It chains together models by piping the output of the left-hand operand to
    the input of the right-hand operand, forming a "pipeline" of models, or
    transformations.

This has different requirements on the inputs/outputs of its operands than do
the arithmetic operators.  For composition all that is required is that the
left-hand model has the same number of outputs as the right-hand model has
inputs.

For simple functional models this is exactly the same as functional
composition, except for the aforementioned caveat about ordering.  For
example, to create the following compound model:

.. graphviz::

    digraph {
        in0 [shape="none", label="input 0"];
        out0 [shape="none", label="output 0"];
        redshift0 [shape="box", label="RedshiftScaleFactor"];
        gaussian0 [shape="box", label="Gaussian1D(1, 0.75, 0.1)"];

        in0 -> redshift0;
        redshift0 -> gaussian0;
        gaussian0 -> out0;
    }

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import RedshiftScaleFactor, Gaussian1D

    x = np.linspace(0, 1.2, 100)
    g0 = RedshiftScaleFactor(0) | Gaussian1D(1, 0.75, 0.1)

    plt.figure(figsize=(8, 5))
    plt.plot(x, g0(x), 'g--', label='$z=0$')

    for z in (0.2, 0.4, 0.6):
        g = RedshiftScaleFactor(z) | Gaussian1D(1, 0.75, 0.1)
        plt.plot(x, g(x), color=plt.cm.OrRd(z),
                 label='$z={0}$'.format(z))

    plt.xlabel('Energy')
    plt.ylabel('Flux')
    plt.legend()

If you wish to perform redshifting in the wavelength space instead of energy,
and would also like to conserve flux, here is another way to do it using
model *instances*:

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import RedshiftScaleFactor, Gaussian1D, Scale

    x = np.linspace(1000, 5000, 1000)
    g0 = Gaussian1D(1, 2000, 200)  # No redshift is same as redshift with z=0

    plt.figure(figsize=(8, 5))
    plt.plot(x, g0(x), 'g--', label='$z=0$')

    for z in (0.2, 0.4, 0.6):
        rs = RedshiftScaleFactor(z).inverse  # Redshift in wavelength space
        sc = Scale(1. / (1 + z))  # Rescale the flux to conserve energy
        g = rs | g0 | sc
        plt.plot(x, g(x), color=plt.cm.OrRd(z),
                 label='$z={0}$'.format(z))

    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.legend()

When working with models with multiple inputs and outputs the same idea
applies.  If each input is thought of as a coordinate axis, then this defines a
pipeline of transformations for the coordinates on each axis (though it does
not necessarily guarantee that these transformations are separable).  For
example:

.. graphviz::

    digraph {
        in0 [shape="none", label="input 0"];
        in1 [shape="none", label="input 1"];
        out0 [shape="none", label="output 0"];
        out1 [shape="none", label="output 1"];
        rot0 [shape="box", label="Rotation2D"];
        gaussian0 [shape="box", label="Gaussian2D(1, 0, 0, 0.1, 0.3)"];

        in0 -> rot0;
        in1 -> rot0;
        rot0 -> gaussian0;
        rot0 -> gaussian0;
        gaussian0 -> out0;
        gaussian0 -> out1;
    }

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Rotation2D, Gaussian2D

    x, y = np.mgrid[-1:1:0.01, -1:1:0.01]

    plt.figure(figsize=(8, 2.5))

    for idx, theta in enumerate((0, 45, 90)):
        g = Rotation2D(theta) | Gaussian2D(1, 0, 0, 0.1, 0.3)
        plt.subplot(1, 3, idx + 1)
        plt.imshow(g(x, y), origin='lower')
        plt.xticks([])
        plt.yticks([])
        plt.title('Rotated $ {0}^\circ $'.format(theta))

.. note::

    The above example is a bit contrived in that
    `~astropy.modeling.functional_models.Gaussian2D` already supports an
    optional rotation parameter.  However, this demonstrates how coordinate
    rotation could be added to arbitrary models.

Normally it is not possible to compose, say, a model with two outputs and a
function of only one input::

    >>> from astropy.modeling.models import Rotation2D
    >>> Rotation2D() | Gaussian1D()  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    ModelDefinitionError: Unsupported operands for |: Rotation2D (n_inputs=2, n_outputs=2) and Gaussian1D (n_inputs=1, n_outputs=1); n_outputs for the left-hand model must match n_inputs for the right-hand model.

However, as we will see in the next section,
:ref:`compound-model-concatenation`, provides a means of creating models
that apply transformations to only some of the outputs from a model,
especially when used in concert with :ref:`mappings <compound-model-mappings>`.


.. _compound-model-concatenation:

Model concatenation
-------------------

The concatenation operator ``&``, sometimes also referred to as a "join",
combines two models into a single, fully separable transformation.  That is, it
makes a new model that takes the inputs to the left-hand model, concatenated
with the inputs to the right-hand model, and returns a tuple consisting of the
two models' outputs concatenated together, without mixing in any way.  In other
words, it simply evaluates the two models in parallel--it can be thought of as
something like a tuple of models.

For example, given two coordinate axes, we can scale each coordinate
by a different factor by concatenating two
`~astropy.modeling.functional_models.Scale` models.

.. graphviz::

    digraph {
        in0 [shape="none", label="input 0"];
        in1 [shape="none", label="input 1"];
        out0 [shape="none", label="output 0"];
        out1 [shape="none", label="output 1"];
        scale0 [shape="box", label="Scale(factor=1.2)"];
        scale1 [shape="box", label="Scale(factor=3.4)"];

        in0 -> scale0;
        scale0 -> out0;

        in1 -> scale1;
        scale1 -> out1;
    }

::

    >>> from astropy.modeling.models import Scale
    >>> separate_scales = Scale(factor=1.2) & Scale(factor=3.4)
    >>> separate_scales(1, 2)  # doctest: +FLOAT_CMP
    (1.2, 6.8)

We can also combine concatenation with composition to build chains of
transformations that use both "1D" and "2D" models on two (or more) coordinate
axes:

.. graphviz::

    digraph {
        in0 [shape="none", label="input 0"];
        in1 [shape="none", label="input 1"];
        out0 [shape="none", label="output 0"];
        out1 [shape="none", label="output 1"];
        scale0 [shape="box", label="Scale(factor=1.2)"];
        scale1 [shape="box", label="Scale(factor=3.4)"];
        rot0 [shape="box", label="Rotation2D(90)"];

        in0 -> scale0;
        scale0 -> rot0;

        in1 -> scale1;
        scale1 -> rot0;

        rot0 -> out0;
        rot0 -> out1;
    }

::

    >>> scale_and_rotate = ((Scale(factor=1.2) & Scale(factor=3.4)) |
    ...                     Rotation2D(90))
    >>> scale_and_rotate.n_inputs
    2
    >>> scale_and_rotate.n_outputs
    2
    >>> scale_and_rotate(1, 2)  # doctest: +FLOAT_CMP
    (-6.8, 1.2)

This is of course equivalent to an
`~astropy.modeling.projections.AffineTransformation2D` with the appropriate
transformation matrix::

    >>> from numpy import allclose
    >>> from astropy.modeling.models import AffineTransformation2D
    >>> affine = AffineTransformation2D(matrix=[[0, -3.4], [1.2, 0]])
    >>> # May be small numerical differences due to different implementations
    >>> allclose(scale_and_rotate(1, 2), affine(1, 2))
    True

Other Topics
============

Model names
-----------

In the above two examples another notable feature of the generated compound
model classes is that the class name, as displayed when printing the class at
the command prompt, is not "TwoGaussians", "FourGaussians", etc.  Instead it is
a generated name consisting of "CompoundModel" followed by an essentially
arbitrary integer that is chosen simply so that every compound model has a
unique default name.  This is a limitation at present, due to the limitation
that it is not generally possible in Python when an object is created by an
expression for it to "know" the name of the variable it will be assigned to, if
any.
It is possible to directly assign a name to the compound model instance
by using the `Model.name <astropy.modeling.Model.name>` attribute::

    >>> two_gaussians.name = "TwoGaussians"
    >>> print(two_gaussians)  # doctest: +SKIP
    Model: CompoundModel...
    Name: TwoGaussians
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Expression: [0] + [1]
    Components:
        [0]: <Gaussian1D(amplitude=1.1, mean=0.1, stddev=0.2)>
        <BLANKLINE>
        [1]: <Gaussian1D(amplitude=2.5, mean=0.5, stddev=0.1)>
    Parameters:
        amplitude_0 mean_0 stddev_0 amplitude_1 mean_1 stddev_1
        ----------- ------ -------- ----------- ------ --------
                1.1    0.1      0.2         2.5    0.5      0.1

.. _compound-model-indexing:

Indexing and slicing
--------------------

As seen in some of the previous examples in this document, when creating a
compound model each component of the model is assigned an integer index
starting from zero.  These indices are assigned simply by reading the
expression that defined the model, from left to right, regardless of the order
of operations.  For example::

    >>> from astropy.modeling.models import Const1D
    >>> A = Const1D(1.1, name='A')
    >>> B = Const1D(2.1, name='B')
    >>> C = Const1D(3.1, name='C')
    >>> M = A + B * C
    >>> print(M)
    Model: CompoundModel...
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Expression: [0] + [1] * [2]
    Components:
        [0]: <Const1D(amplitude=1.1, name='A')>
    <BLANKLINE>
        [1]: <Const1D(amplitude=2.1, name='B')>
    <BLANKLINE>
        [2]: <Const1D(amplitude=3.1, name='C')>
    Parameters:
        amplitude_0 amplitude_1 amplitude_2
        ----------- ----------- -----------
                1.1         2.1         3.1


In this example the expression is evaluated ``(B * C) + A``--that is, the
multiplication is evaluated before the addition per usual arithmetic rules.
However, the components of this model are simply read off left to right from
the expression ``A + B * C``, with ``A -> 0``, ``B -> 1``, ``C -> 2``.  If we
had instead defined ``M = C * B + A`` then the indices would be reversed
(though the expression is mathematically equivalent).  This convention is
chosen for simplicity--given the list of components it is not necessary to
jump around when mentally mapping them to the expression.

We can pull out each individual component of the compound model ``M`` by using
indexing notation on it.  Following from the above example, ``M[1]`` should
return the model ``B``::

    >>> M[1]
    <Const1D(amplitude=2.1, name='B')>

We can also take a *slice* of the compound model.  This returns a new compound
model that evaluates the *subexpression* involving the models selected by the
slice.  This follows the same semantics as slicing a `list` or array in Python.
The start point is inclusive and the end point is exclusive.  So a slice like
``M[1:3]`` (or just ``M[1:]``) selects models ``B`` and ``C`` (and all
*operators* between them).  So the resulting model evaluates just the
subexpression ``B * C``::

    >>> print(M[1:])
    Model: CompoundModel
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Expression: [0] * [1]
    Components:
        [0]: <Const1D(amplitude=2.1, name='B')>
    <BLANKLINE>
        [1]: <Const1D(amplitude=3.1, name='C')>
    Parameters:
        amplitude_0 amplitude_1
        ----------- -----------
                2.1         3.1

.. note::

    There is a change in the parameter names of a slice from versions
    prior to 4.0. Previously, the parameter names were identical to that
    of the model being sliced. Now, they are what is expected for a
    compound model of this type apart from the model sliced. That is,
    the sliced model always starts with its own relative index for its
    components, thus the parameter names start with a 0 suffix.

.. note::

    Starting with 4.0, the behavior of slicing is more restrictive than
    previously. For example if::

        m = m1 * m2 + m3

    and one sliced by
    using ``m[1:3]`` previously that would return the model: ``m2 + m3``
    even though there was never any such submodel of m. Starting with 4.0
    a slice must correspond to a submodel (something that corresponds
    to an intermediate result of the computational chain of evaluating
    the compound model). So::

        m1 * m2

    is a submodel (i.e.,``m[:2]``) but
    ``m[1:3]`` is not. Currently this also means that in simpler expressions
    such as::

        m = m1 + m2 + m3 + m4

    where any slice should be valid in
    principle, only slices that include m1 are since it is part of
    all submodules (since the order of evaluation is::

        ((m1 + m2) + m3) + m4

    Anyone creating compound models that wishes submodels to be available
    is advised to use parentheses explicitly  or define intermediate
    models to be used in subsequent expressions so that they can be
    extracted with a slice or simple index depending on the context.
    For example, to make ``m2 + m3`` accessible by slice define ``m`` as::

        m = m1 + (m2 + m3) + m4. In this case ``m[1:3]`` will work.

The new compound model for the subexpression can be evaluated
like any other::

    >>> M[1:](0)  # doctest: +FLOAT_CMP
    6.51

Although the model ``M`` was composed entirely of ``Const1D`` models in this
example, it was useful to give each component a unique name (``A``, ``B``,
``C``) in order to differentiate between them.  This can also be used for
indexing and slicing::

    >>> print(M['B'])
    Model: Const1D
    Name: B
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Parameters:
        amplitude
        ---------
              2.1


In this case ``M['B']`` is equivalent to ``M[1]``.  But by using the name we do
not have to worry about what index that component is in (this becomes
especially useful when combining multiple compound models).  A current
limitation, however, is that each component of a compound model must have a
unique name--if some components have duplicate names then they can only be
accessed by their integer index.

Slicing also works with names.  When using names the start and end points are
*both inclusive*::

    >>> print(M['B':'C'])
    Model: CompoundModel...
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Expression: [0] * [1]
    Components:
        [0]: <Const1D(amplitude=2.1, name='B')>
    <BLANKLINE>
        [1]: <Const1D(amplitude=3.1, name='C')>
    Parameters:
        amplitude_0 amplitude_1
        ----------- -----------
                2.1         3.1

So in this case ``M['B':'C']`` is equivalent to ``M[1:3]``.

.. _compound-model-parameters:

Parameters
----------

A question that frequently comes up when first encountering compound models is
how exactly all the parameters are dealt with.  By now we've seen a few
examples that give some hints, but a more detailed explanation is in order.
This is also one of the biggest areas for possible improvements--the current
behavior is meant to be practical, but is not ideal.  (Some possible
improvements include being able to rename parameters, and providing a means of
narrowing down the number of parameters in a compound model.)

As explained in the general documentation for model :ref:`parameters
<modeling-parameters>`, every model has an attribute called
`~astropy.modeling.Model.param_names` that contains a tuple of all the model's
adjustable parameters.  These names are given in a canonical order that also
corresponds to the order in which the parameters should be specified when
instantiating the model.

The simple scheme used currently for naming parameters in a compound model is
this:  The ``param_names`` from each component model are concatenated with each
other in order from left to right as explained in the section on
:ref:`compound-model-indexing`.  However, each parameter name is appended with
``_<#>``, where ``<#>`` is the index of the component model that parameter
belongs to.  For example::

    >>> Gaussian1D.param_names
    ('amplitude', 'mean', 'stddev')
    >>> (Gaussian1D() + Gaussian1D()).param_names
    ('amplitude_0', 'mean_0', 'stddev_0', 'amplitude_1', 'mean_1', 'stddev_1')

For consistency's sake, this scheme is followed even if not all of the
components have overlapping parameter names::

    >>> from astropy.modeling.models import RedshiftScaleFactor
    >>> (RedshiftScaleFactor() | (Gaussian1D() + Gaussian1D())).param_names
    ('z_0', 'amplitude_1', 'mean_1', 'stddev_1', 'amplitude_2', 'mean_2',
    'stddev_2')

On some level a scheme like this is necessary in order for the compound model
to maintain some consistency with other models with respect to the interface to
its parameters.  However, if one gets lost it is also possible to take
advantage of :ref:`indexing <compound-model-indexing>` to make things easier.
When returning a single component from a compound model the parameters
associated with that component are accessible through their original names, but
are still tied back to the compound model::

    >>> a = Gaussian1D(1, 0, 0.2, name='A')
    >>> b = Gaussian1D(2.5, 0.5, 0.1, name='B')
    >>> m = a + b
    >>> m.amplitude_0
    Parameter('amplitude', value=1.0)

is equivalent to::

    >>> m['A'].amplitude
    Parameter('amplitude', value=1.0)

You can think of these both as different "views" of the same parameter.
Updating one updates the other::

    >>> m.amplitude_0 = 42
    >>> m['A'].amplitude
    Parameter('amplitude', value=42.0)
    >>> m['A'].amplitude = 99
    >>> m.amplitude_0
    Parameter('amplitude', value=99.0)

Note, however, that the original
`~astropy.modeling.functional_models.Gaussian1D` instance ``a`` has been
updated::

    >>> a.amplitude
    Parameter('amplitude', value=99.0)

This is different than the behavior in versions prior to 4.0. Now compound model
parameters share the same Parameter instance as the original model.


.. _compound-model-mappings:

Advanced mappings
-----------------

We have seen in some previous examples how models can be chained together to
form a "pipeline" of transformations by using model :ref:`composition
<compound-model-composition>` and :ref:`concatenation
<compound-model-concatenation>`.  To aid the creation of more complex chains of
transformations (for example for a WCS transformation) a new class of
"`mapping <astropy.modeling.mappings>`" models is provided.

Mapping models do not (currently) take any parameters, nor do they perform any
numeric operation.  They are for use solely with the :ref:`concatenation
<compound-model-concatenation>` (``&``) and :ref:`composition
<compound-model-composition>` (``|``) operators, and can be used to control how
the inputs and outputs of models are ordered, and how outputs from one model
are mapped to inputs of another model in a composition.

Currently there are only two mapping models:
`~astropy.modeling.mappings.Identity`, and (the somewhat generically named)
`~astropy.modeling.mappings.Mapping`.

The `~astropy.modeling.mappings.Identity` mapping simply passes one or more
inputs through, unchanged.  It must be instantiated with an integer specifying
the number of inputs/outputs it accepts.  This can be used to trivially expand
the "dimensionality" of a model in terms of the number of inputs it accepts.
In the section on :ref:`concatenation <compound-model-concatenation>` we saw
an example like::

    >>> m = (Scale(1.2) & Scale(3.4)) | Rotation2D(90)


.. graphviz::

    digraph {
        in0 [shape="none", label="input 0"];
        in1 [shape="none", label="input 1"];
        out0 [shape="none", label="output 0"];
        out1 [shape="none", label="output 1"];
        scale0 [shape="box", label="Scale(factor=1.2)"];
        scale1 [shape="box", label="Scale(factor=3.4)"];
        rot0 [shape="box", label="Rotation2D(90)"];

        in0 -> scale0;
        scale0 -> rot0;

        in1 -> scale1;
        scale1 -> rot0;

        rot0 -> out0;
        rot0 -> out1;
    }

where two coordinate inputs are scaled individually and then rotated into each
other.  However, say we wanted to scale only one of those coordinates.  It
would be fine to simply use ``Scale(1)`` for one them, or any other model that
is effectively a no-op.  But that also adds unnecessary computational overhead,
so we might as well simply specify that that coordinate is not to be scaled or
transformed in any way.  This is a good use case for
`~astropy.modeling.mappings.Identity`:

.. graphviz::

    digraph {
        in0 [shape="none", label="input 0"];
        in1 [shape="none", label="input 1"];
        out0 [shape="none", label="output 0"];
        out1 [shape="none", label="output 1"];
        scale0 [shape="box", label="Scale(factor=1.2)"];
        identity0 [shape="box", label="Identity(1)"];
        rot0 [shape="box", label="Rotation2D(90)"];

        in0 -> scale0;
        scale0 -> rot0;

        in1 -> identity0;
        identity0 -> rot0;

        rot0 -> out0;
        rot0 -> out1;
    }

::

    >>> from astropy.modeling.models import Identity
    >>> m = Scale(1.2) & Identity(1)
    >>> m(1, 2)  # doctest: +FLOAT_CMP
    (1.2, 2.0)


This scales the first input, and passes the second one through unchanged.  We
can use this to build up more complicated steps in a many-axis WCS
transformation.  If for example we had 3 axes and only wanted to scale the
first one:

.. graphviz::

    digraph {
        in0 [shape="none", label="input 0"];
        in1 [shape="none", label="input 1"];
        in2 [shape="none", label="input 2"];
        out0 [shape="none", label="output 0"];
        out1 [shape="none", label="output 1"];
        out2 [shape="none", label="output 2"];
        scale0 [shape="box", label="Scale(1.2)"];
        identity0 [shape="box", label="Identity(2)"];

        in0 -> scale0;
        scale0 -> out0;

        in1 -> identity0;
        in2 -> identity0;
        identity0 -> out1;
        identity0 -> out2;
    }

::

    >>> m = Scale(1.2) & Identity(2)
    >>> m(1, 2, 3)  # doctest: +FLOAT_CMP
    (1.2, 2.0, 3.0)

(Naturally, the last example could also be written out ``Scale(1.2) &
Identity(1) & Identity(1)``.)

The `~astropy.modeling.mappings.Mapping` model is similar in that it does not
modify any of its inputs.  However, it is more general in that it allows inputs
to be duplicated, reordered, or even dropped outright.  It is instantiated with
a single argument: a `tuple`, the number of items of which correspond to the
number of outputs the `~astropy.modeling.mappings.Mapping` should produce.  A
1-tuple means that whatever inputs come in to the
`~astropy.modeling.mappings.Mapping`, only one will be output.  And so on for
2-tuple or higher (though the length of the tuple cannot be greater than the
number of inputs--it will not pull values out of thin air).  The elements of
this mapping are integers corresponding to the indices of the inputs.  For
example, a mapping of ``Mapping((0,))`` is equivalent to ``Identity(1)``--it
simply takes the first (0-th) input and returns it:

.. graphviz::

    digraph G {
        in0 [shape="none", label="input 0"];

        subgraph cluster_A {
            shape=rect;
            color=black;
            label="(0,)";

            a [shape=point, label=""];
        }

        out0 [shape="none", label="output 0"];

        in0 -> a;
        a -> out0;
    }

::

    >>> from astropy.modeling.models import Mapping
    >>> m = Mapping((0,))
    >>> m(1.0)
    1.0

Likewise ``Mapping((0, 1))`` is equivalent to ``Identity(2)``, and so on.
However, `~astropy.modeling.mappings.Mapping` also allows outputs to be
reordered arbitrarily:

.. graphviz::

    digraph G {
        {
            rank=same;
            in0 [shape="none", label="input 0"];
            in1 [shape="none", label="input 1"];
        }

        subgraph cluster_A {
            shape=rect;
            color=black;
            label="(1, 0)";

            {
                rank=same;
                a [shape=point, label=""];
                b [shape=point, label=""];
            }

            {
                rank=same;
                c [shape=point, label=""];
                d [shape=point, label=""];
            }

            a -> c [style=invis];
            a -> d [constraint=false];
            b -> c [constraint=false];
        }

        {
            rank=same;
            out0 [shape="none", label="output 0"];
            out1 [shape="none", label="output 1"];
        }

        in0 -> a;
        in1 -> b;
        c -> out0;
        d -> out1;
    }

::

    >>> m = Mapping((1, 0))
    >>> m(1.0, 2.0)
    (2.0, 1.0)

.. graphviz::

    digraph G {
        {
            rank=same;
            in0 [shape="none", label="input 0"];
            in1 [shape="none", label="input 1"];
            in2 [shape="none", label="input 2"];
        }

        subgraph cluster_A {
            shape=rect;
            color=black;
            label="(1, 0, 2)";

            {
                rank=same;
                a [shape=point, label=""];
                b [shape=point, label=""];
                c [shape=point, label=""];
            }

            {
                rank=same;
                d [shape=point, label=""];
                e [shape=point, label=""];
                f [shape=point, label=""];
            }

            a -> d [style=invis];
            a -> e [constraint=false];
            b -> d [constraint=false];
            c -> f [constraint=false];
        }

        {
            rank=same;
            out0 [shape="none", label="output 0"];
            out1 [shape="none", label="output 1"];
            out2 [shape="none", label="output 2"];
        }

        in0 -> a;
        in1 -> b;
        in2 -> c;
        d -> out0;
        e -> out1;
        f -> out2;
    }

::

    >>> m = Mapping((1, 0, 2))
    >>> m(1.0, 2.0, 3.0)
    (2.0, 1.0, 3.0)

Outputs may also be dropped:

.. graphviz::

    digraph G {
        {
            rank=same;
            in0 [shape="none", label="input 0"];
            in1 [shape="none", label="input 1"];
        }

        subgraph cluster_A {
            shape=rect;
            color=black;
            label="(1,)";

            {
                rank=same;
                a [shape=point, label=""];
                b [shape=point, label=""];
            }

            {
                rank=same;
                c [shape=point, label=""];
            }

            a -> c [style=invis];
            b -> c [constraint=false];
        }

        out0 [shape="none", label="output 0"];

        in0 -> a;
        in1 -> b;
        c -> out0;
    }

::

    >>> m = Mapping((1,))
    >>> m(1.0, 2.0)
    2.0

.. graphviz::

    digraph G {
        {
            rank=same;
            in0 [shape="none", label="input 0"];
            in1 [shape="none", label="input 1"];
            in2 [shape="none", label="input 2"];
        }

        subgraph cluster_A {
            shape=rect;
            color=black;
            label="(0, 2)";

            {
                rank=same;
                a [shape=point, label=""];
                b [shape=point, label=""];
                c [shape=point, label=""];
            }

            {
                rank=same;
                d [shape=point, label=""];
                e [shape=point, label=""];
            }

            a -> d [style=invis];
            a -> d [constraint=false];
            c -> e [constraint=false];
        }

        {
            rank=same;
            out0 [shape="none", label="output 0"];
            out1 [shape="none", label="output 1"];
        }

        in0 -> a;
        in1 -> b;
        in2 -> c;
        d -> out0;
        e -> out1;
    }

::

    >>> m = Mapping((0, 2))
    >>> m(1.0, 2.0, 3.0)
    (1.0, 3.0)

Or duplicated:

.. graphviz::

    digraph G {
        in0 [shape="none", label="input 0"];

        subgraph cluster_A {
            shape=rect;
            color=black;
            label="(0, 0)";

            a [shape=point, label=""];

            {
                rank=same;
                b [shape=point, label=""];
                c [shape=point, label=""];
            }

            a -> b [style=invis];
            a -> b [constraint=false];
            a -> c [constraint=false];
        }

        {
            rank=same;
            out0 [shape="none", label="output 0"];
            out1 [shape="none", label="output 1"];
        }

        in0 -> a;
        b -> out0;
        c -> out1;
    }

::

    >>> m = Mapping((0, 0))
    >>> m(1.0)
    (1.0, 1.0)

.. graphviz::

    digraph G {
        {
            rank=same;
            in0 [shape="none", label="input 0"];
            in1 [shape="none", label="input 1"];
            in2 [shape="none", label="input 2"];
        }

        subgraph cluster_A {
            shape=rect;
            color=black;
            label="(0, 1, 1, 2)";

            {
                rank=same;
                a [shape=point, label=""];
                b [shape=point, label=""];
                c [shape=point, label=""];
            }

            {
                rank=same;
                d [shape=point, label=""];
                e [shape=point, label=""];
                f [shape=point, label=""];
                g [shape=point, label=""];
            }

            a -> d [style=invis];
            a -> d [constraint=false];
            b -> e [constraint=false];
            b -> f [constraint=false];
            c -> g [constraint=false];
        }

        {
            rank=same;
            out0 [shape="none", label="output 0"];
            out1 [shape="none", label="output 1"];
            out2 [shape="none", label="output 2"];
            out3 [shape="none", label="output 3"];
        }

        in0 -> a;
        in1 -> b;
        in2 -> c;
        d -> out0;
        e -> out1;
        f -> out2;
        g -> out3;
    }

::

    >>> m = Mapping((0, 1, 1, 2))
    >>> m(1.0, 2.0, 3.0)
    (1.0, 2.0, 2.0, 3.0)


A complicated example that performs multiple transformations, some separable,
some not, on three coordinate axes might look something like:

.. graphviz::

    digraph G {
        {
            rank=same;
            in0 [shape="none", label="input 0"];
            in1 [shape="none", label="input 1"];
            in2 [shape="none", label="input 2"];
        }

        {
            rank=same;
            poly0 [shape=rect, label="Poly1D(3, c0=1, c3=1)"];
            identity0 [shape=rect, label="Identity(1)"];
            poly1 [shape=rect, label="Poly1D(2, c2=1)"];
        }

        subgraph cluster_A {
            shape=rect;
            color=black;
            label="(0, 2, 1)";

            {
                rank=same;
                a [shape=point, label=""];
                b [shape=point, label=""];
                c [shape=point, label=""];
            }

            {
                rank=same;
                d [shape=point, label=""];
                e [shape=point, label=""];
                f [shape=point, label=""];
            }

            a -> d [style=invis];
            d -> e [style=invis];
            a -> d [constraint=false];
            c -> e [constraint=false];
            b -> f [constraint=false];
        }

        poly2 [shape="rect", label="Poly2D(4, c0_0=1, c1_1=1, c2_2=2)"];
        gaussian0 [shape="rect", label="Gaussian1D(1, 0, 4)"];

        {
            rank=same;
            out0 [shape="none", label="output 0"];
            out1 [shape="none", label="output 1"];
        }

        in0 -> poly0;
        in1 -> identity0;
        in2 -> poly1;
        poly0 -> a;
        identity0 -> b;
        poly1 -> c;
        d -> poly2;
        e -> poly2;
        f -> gaussian0;
        poly2 -> out0;
        gaussian0 -> out1;
    }

::

    >>> from astropy.modeling.models import Polynomial1D as Poly1D
    >>> from astropy.modeling.models import Polynomial2D as Poly2D
    >>> m = ((Poly1D(3, c0=1, c3=1) & Identity(1) & Poly1D(2, c2=1)) |
    ...      Mapping((0, 2, 1)) |
    ...      (Poly2D(4, c0_0=1, c1_1=1, c2_2=2) & Gaussian1D(1, 0, 4)))
    ...
    >>> m(2, 3, 4)  # doctest: +FLOAT_CMP
    (41617.0, 0.7548396019890073)



This expression takes three inputs: :math:`x`, :math:`y`, and :math:`z`.  It
first takes :math:`x \rightarrow x^3 + 1` and :math:`z \rightarrow z^2`.
Then it remaps the axes so that :math:`x` and :math:`z` are passed in to the
`~astropy.modeling.polynomial.Polynomial2D` to evaluate
:math:`2x^2z^2 + xz + 1`, while simultaneously evaluating a Gaussian on
:math:`y`.  The end result is a reduction down to two coordinates.  You can
confirm for yourself that the result is correct.

This opens up the possibility of essentially arbitrarily complex transformation
graphs.  Currently the tools do not exist to make it easy to navigate and
reason about highly complex compound models that use these mappings, but that
is a possible enhancement for future versions.

.. _model-reduction:

Model Reduction
---------------

In order to save much duplication in the construction of complex models, it is
possible to define one complex model that covers all cases where the
variables that distinguish the models are made part of the model's input
variables. The ``fix_inputs`` function allows defining models derived from
the more complex one by setting one or more of the inputs to a constant
value. Examples of this sort of situation arise when working out
the transformations from detector pixel to RA, Dec, and lambda for
spectrographs when the slit locations may be moved (e.g., fiber fed or
commandable slit masks), or different orders may be selected (e.g., Eschelle).
In the case of order, one may have a function of pixel ``x``, ``y``, ``spectral_order``
that map into ``RA``, ``Dec`` and ``lambda``. Without specifying ``spectral_order``, it is
ambiguous what ``RA``, ``Dec`` and ``Lambda`` corresponds to a pixel location. It
is usually possible to define a function of all three inputs. Presuming
this model is ``general_transform`` then ``fix_inputs`` may be used to define
the transform for a specific order as follows:

::
     >>> order1_transform = fix_inputs(general_transform, {'order': 1})  # doctest: +SKIP

creates a new compound model that takes only pixel position and generates
``RA``, ``Dec``, and ``lambda``. The ``fix_inputs`` function can be used to set input
values by position (0 is the first) or by input variable name, and more
than one can be set in the dictionary supplied.

If the input model has a bounding_box, the generated model will have the
bounding for the input coordinate removed.
