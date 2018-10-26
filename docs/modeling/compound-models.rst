.. .. _compound-models:

.. Compound Models
.. ***************

.. .. versionadded:: 1.0

As noted in the :ref:`introduction to the modeling package
<compound-models-intro>`, it is now possible to create new models just by
combining existing models using the arithmetic operators ``+``, ``-``, ``*``,
``/``, and ``**``, as well as by model composition using ``|`` and
concatenation (explained below) with ``&``.


Some terminology
================

..     plt.xlabel('Energy')
..     plt.ylabel('Flux')
..     plt.legend()

.. If you wish to perform redshifting in the wavelength space instead of energy,
.. and would also like to conserve flux, here is another way to do it using
.. model *instances*:

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


.. _compound-model-indexing:

Indexing and slicing
====================

..         in0 -> scale0;
..         scale0 -> rot0;

..         in1 -> identity0;
..         identity0 -> rot0;

..         rot0 -> out0;
..         rot0 -> out1;
..     }

.. ::

..     >>> from astropy.modeling.models import Identity
..     >>> m = Scale(1.2) & Identity(1)
..     >>> m(1, 2)  # doctest: +FLOAT_CMP
..     (1.2, 2.0)


.. This scales the first input, and passes the second one through unchanged.  We
.. can use this to build up more complicated steps in a many-axis WCS
.. transformation.  If for example we had 3 axes and only wanted to scale the
.. first one:

.. .. graphviz::

..     digraph {
..         in0 [shape="none", label="input 0"];
..         in1 [shape="none", label="input 1"];
..         in2 [shape="none", label="input 2"];
..         out0 [shape="none", label="output 0"];
..         out1 [shape="none", label="output 1"];
..         out2 [shape="none", label="output 2"];
..         scale0 [shape="box", label="Scale(1.2)"];
..         identity0 [shape="box", label="Identity(2)"];

..         in0 -> scale0;
..         scale0 -> out0;

..         in1 -> identity0;
..         in2 -> identity0;
..         identity0 -> out1;
..         identity0 -> out2;
..     }

.. ::

..     >>> m = Scale(1.2) & Identity(2)
..     >>> m(1, 2, 3)  # doctest: +FLOAT_CMP
..     (1.2, 2.0, 3.0)

.. (Naturally, the last example could also be written out ``Scale(1.2) &
.. Identity(1) & Identity(1)``.)

.. The `~astropy.modeling.mappings.Mapping` model is similar in that it does not
.. modify any of its inputs.  However, it is more general in that it allows inputs
.. to be duplicated, reordered, or even dropped outright.  It is instantiated with
.. a single argument: a `tuple`, the number of items of which correspond to the
.. number of outputs the `~astropy.modeling.mappings.Mapping` should produce.  A
.. 1-tuple means that whatever inputs come in to the
.. `~astropy.modeling.mappings.Mapping`, only one will be output.  And so on for
.. 2-tuple or higher (though the length of the tuple cannot be greater than the
.. number of inputs--it will not pull values out of thin air).  The elements of
.. this mapping are integers corresponding to the indices of the inputs.  For
.. example, a mapping of ``Mapping((0,))`` is equivalent to ``Identity(1)``--it
.. simply takes the first (0-th) input and returns it:

.. .. graphviz::

..     digraph G {
..         in0 [shape="none", label="input 0"];

..         subgraph cluster_A {
..             shape=rect;
..             color=black;
..             label="(0,)";

..             a [shape=point, label=""];
..         }

..         out0 [shape="none", label="output 0"];

..         in0 -> a;
..         a -> out0;
..     }

.. ::

..     >>> from astropy.modeling.models import Mapping
..     >>> m = Mapping((0,))
..     >>> m(1.0)
..     1.0

.. Likewise ``Mapping((0, 1))`` is equivalent to ``Identity(2)``, and so on.
.. However, `~astropy.modeling.mappings.Mapping` also allows outputs to be
.. reordered arbitrarily:

.. .. graphviz::

..     digraph G {
..         {
..             rank=same;
..             in0 [shape="none", label="input 0"];
..             in1 [shape="none", label="input 1"];
..         }

..         subgraph cluster_A {
..             shape=rect;
..             color=black;
..             label="(1, 0)";

..             {
..                 rank=same;
..                 a [shape=point, label=""];
..                 b [shape=point, label=""];
..             }

..             {
..                 rank=same;
..                 c [shape=point, label=""];
..                 d [shape=point, label=""];
..             }

..             a -> c [style=invis];
..             a -> d [constraint=false];
..             b -> c [constraint=false];
..         }

..         {
..             rank=same;
..             out0 [shape="none", label="output 0"];
..             out1 [shape="none", label="output 1"];
..         }

..         in0 -> a;
..         in1 -> b;
..         c -> out0;
..         d -> out1;
..     }

.. ::

..     >>> m = Mapping((1, 0))
..     >>> m(1.0, 2.0)
..     (2.0, 1.0)

.. .. graphviz::

..     digraph G {
..         {
..             rank=same;
..             in0 [shape="none", label="input 0"];
..             in1 [shape="none", label="input 1"];
..             in2 [shape="none", label="input 2"];
..         }

..         subgraph cluster_A {
..             shape=rect;
..             color=black;
..             label="(1, 0, 2)";

..             {
..                 rank=same;
..                 a [shape=point, label=""];
..                 b [shape=point, label=""];
..                 c [shape=point, label=""];
..             }

..             {
..                 rank=same;
..                 d [shape=point, label=""];
..                 e [shape=point, label=""];
..                 f [shape=point, label=""];
..             }

..             a -> d [style=invis];
..             a -> e [constraint=false];
..             b -> d [constraint=false];
..             c -> f [constraint=false];
..         }

..         {
..             rank=same;
..             out0 [shape="none", label="output 0"];
..             out1 [shape="none", label="output 1"];
..             out2 [shape="none", label="output 2"];
..         }

..         in0 -> a;
..         in1 -> b;
..         in2 -> c;
..         d -> out0;
..         e -> out1;
..         f -> out2;
..     }

.. ::

..     >>> m = Mapping((1, 0, 2))
..     >>> m(1.0, 2.0, 3.0)
..     (2.0, 1.0, 3.0)

.. Outputs may also be dropped:

.. .. graphviz::

..     digraph G {
..         {
..             rank=same;
..             in0 [shape="none", label="input 0"];
..             in1 [shape="none", label="input 1"];
..         }

..         subgraph cluster_A {
..             shape=rect;
..             color=black;
..             label="(1,)";

..             {
..                 rank=same;
..                 a [shape=point, label=""];
..                 b [shape=point, label=""];
..             }

..             {
..                 rank=same;
..                 c [shape=point, label=""];
..             }

..             a -> c [style=invis];
..             b -> c [constraint=false];
..         }

..         out0 [shape="none", label="output 0"];

..         in0 -> a;
..         in1 -> b;
..         c -> out0;
..     }

.. ::

..     >>> m = Mapping((1,))
..     >>> m(1.0, 2.0)
..     2.0

.. .. graphviz::

..     digraph G {
..         {
..             rank=same;
..             in0 [shape="none", label="input 0"];
..             in1 [shape="none", label="input 1"];
..             in2 [shape="none", label="input 2"];
..         }

..         subgraph cluster_A {
..             shape=rect;
..             color=black;
..             label="(0, 2)";

..             {
..                 rank=same;
..                 a [shape=point, label=""];
..                 b [shape=point, label=""];
..                 c [shape=point, label=""];
..             }

..             {
..                 rank=same;
..                 d [shape=point, label=""];
..                 e [shape=point, label=""];
..             }

..             a -> d [style=invis];
..             a -> d [constraint=false];
..             c -> e [constraint=false];
..         }

..         {
..             rank=same;
..             out0 [shape="none", label="output 0"];
..             out1 [shape="none", label="output 1"];
..         }

..         in0 -> a;
..         in1 -> b;
..         in2 -> c;
..         d -> out0;
..         e -> out1;
..     }

.. ::

..     >>> m = Mapping((0, 2))
..     >>> m(1.0, 2.0, 3.0)
..     (1.0, 3.0)

.. Or duplicated:

.. .. graphviz::

..     digraph G {
..         in0 [shape="none", label="input 0"];

..         subgraph cluster_A {
..             shape=rect;
..             color=black;
..             label="(0, 0)";

..             a [shape=point, label=""];

..             {
..                 rank=same;
..                 b [shape=point, label=""];
..                 c [shape=point, label=""];
..             }

..             a -> b [style=invis];
..             a -> b [constraint=false];
..             a -> c [constraint=false];
..         }

..         {
..             rank=same;
..             out0 [shape="none", label="output 0"];
..             out1 [shape="none", label="output 1"];
..         }

..         in0 -> a;
..         b -> out0;
..         c -> out1;
..     }

.. ::

..     >>> m = Mapping((0, 0))
..     >>> m(1.0)
..     (1.0, 1.0)

.. .. graphviz::

..     digraph G {
..         {
..             rank=same;
..             in0 [shape="none", label="input 0"];
..             in1 [shape="none", label="input 1"];
..             in2 [shape="none", label="input 2"];
..         }

..         subgraph cluster_A {
..             shape=rect;
..             color=black;
..             label="(0, 1, 1, 2)";

..             {
..                 rank=same;
..                 a [shape=point, label=""];
..                 b [shape=point, label=""];
..                 c [shape=point, label=""];
..             }

..             {
..                 rank=same;
..                 d [shape=point, label=""];
..                 e [shape=point, label=""];
..                 f [shape=point, label=""];
..                 g [shape=point, label=""];
..             }

..             a -> d [style=invis];
..             a -> d [constraint=false];
..             b -> e [constraint=false];
..             b -> f [constraint=false];
..             c -> g [constraint=false];
..         }

..         {
..             rank=same;
..             out0 [shape="none", label="output 0"];
..             out1 [shape="none", label="output 1"];
..             out2 [shape="none", label="output 2"];
..             out3 [shape="none", label="output 3"];
..         }

..         in0 -> a;
..         in1 -> b;
..         in2 -> c;
..         d -> out0;
..         e -> out1;
..         f -> out2;
..         g -> out3;
..     }

.. ::

..     >>> m = Mapping((0, 1, 1, 2))
..     >>> m(1.0, 2.0, 3.0)
..     (1.0, 2.0, 2.0, 3.0)


.. A complicated example that performs multiple transformations, some separable,
.. some not, on three coordinate axes might look something like:

.. .. graphviz::

..     digraph G {
..         {
..             rank=same;
..             in0 [shape="none", label="input 0"];
..             in1 [shape="none", label="input 1"];
..             in2 [shape="none", label="input 2"];
..         }

..         {
..             rank=same;
..             poly0 [shape=rect, label="Poly1D(3, c0=1, c3=1)"];
..             identity0 [shape=rect, label="Identity(1)"];
..             poly1 [shape=rect, label="Poly1D(2, c2=1)"];
..         }

..         subgraph cluster_A {
..             shape=rect;
..             color=black;
..             label="(0, 2, 1)";

..             {
..                 rank=same;
..                 a [shape=point, label=""];
..                 b [shape=point, label=""];
..                 c [shape=point, label=""];
..             }

..             {
..                 rank=same;
..                 d [shape=point, label=""];
..                 e [shape=point, label=""];
..                 f [shape=point, label=""];
..             }

..             a -> d [style=invis];
..             d -> e [style=invis];
..             a -> d [constraint=false];
..             c -> e [constraint=false];
..             b -> f [constraint=false];
..         }

..         poly2 [shape="rect", label="Poly2D(4, c0_0=1, c1_1=1, c2_2=2)"];
..         gaussian0 [shape="rect", label="Gaussian1D(1, 0, 4)"];

..         {
..             rank=same;
..             out0 [shape="none", label="output 0"];
..             out1 [shape="none", label="output 1"];
..             out2 [shape="none", label="output 2"];
..         }

..         in0 -> poly0;
..         in1 -> identity0;
..         in2 -> poly1;
..         poly0 -> a;
..         identity0 -> b;
..         poly1 -> c;
..         d -> poly2;
..         e -> poly2;
..         f -> gaussian0;
..         poly2 -> out0;
..         poly2 -> out1;
..         gaussian0 -> out2;
..     }

.. ::

..     >>> from astropy.modeling.models import Polynomial1D as Poly1D
..     >>> from astropy.modeling.models import Polynomial2D as Poly2D
..     >>> m = ((Poly1D(3, c0=1, c3=1) & Identity(1) & Poly1D(2, c2=1)) |
..     ...      Mapping((0, 2, 1)) |
..     ...      (Poly2D(4, c0_0=1, c1_1=1, c2_2=2) & Gaussian1D(1, 0, 4)))
..     ...
..     >>> m(2, 3, 4)  # doctest: +FLOAT_CMP
..     (41617.0, 0.7548396019890073)



.. This expression takes three inputs: :math:`x`, :math:`y`, and :math:`z`.  It
.. first takes :math:`x \rightarrow x^3 + 1` and :math:`z \rightarrow z^2`.
.. Then it remaps the axes so that :math:`x` and :math:`z` are passed in to the
.. `~astropy.modeling.polynomial.Polynomial2D` to evaluate
.. :math:`2x^2z^2 + xz + 1`, while simultaneously evaluating a Gaussian on
.. :math:`y`.  The end result is a reduction down to two coordinates.  You can
.. confirm for yourself that the result is correct.

.. This opens up the possibility of essentially arbitrarily complex transformation
.. graphs.  Currently the tools do not exist to make it easy to navigate and
.. reason about highly complex compound models that use these mappings, but that
.. is a possible enhancement for future versions.
