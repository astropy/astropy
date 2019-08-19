.. include:: links.inc

.. _modeling-getting-started-model-sets:

Model sets
==========

In some cases it is necessary to describe many models of the same type but with
different sets of parameter values.  This could be done simply by instantiating
as many instances of a `~astropy.modeling.Model` as are needed.  But that can
be inefficient for a large number of models.  To that end, all model classes in
`astropy.modeling` can also be used to represent a model *set* which is a
collection of models of the same type, but with different values for their
parameters.

To instantiate a model set, use argument ``n_models=N`` where ``N`` is the
number of models in the set when constructing the model.  The value of each
parameter must be a list or array of length ``N``, such that each item in
the array corresponds to one model in the set::

    >>> from astropy.modeling import models
    >>> g = models.Gaussian1D(amplitude=[1, 2], mean=[0, 0],
    ...                       stddev=[0.1, 0.2], n_models=2)
    >>> print(g)
    Model: Gaussian1D
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 2
    Parameters:
        amplitude mean stddev
        --------- ---- ------
              1.0  0.0    0.1
              2.0  0.0    0.2

This is equivalent to two Gaussians with the parameters ``amplitude=1, mean=0,
stddev=0.1`` and ``amplitude=2, mean=0, stddev=0.2`` respectively.  When
printing the model the parameter values are displayed as a table, with each row
corresponding to a single model in the set.

The number of models in a model set can be determined using the `len` builtin::

    >>> len(g)
    2

Single models have a length of 1, and are not considered a model set as such.

When evaluating a model set, by default the input must be the same length as
the number of models, with one input per model::

    >>> g([0, 0.1])  # doctest: +FLOAT_CMP
    array([1.        , 1.76499381])

The result is an array with one result per model in the set.  It is also
possible to broadcast a single value to all models in the set::

    >>> g(0)  # doctest: +FLOAT_CMP
    array([1., 2.])

Model sets are used primarily for fitting, allowing a large number of models of
the same type to be fitted simultaneously (and independently from each other)
to some large set of inputs.  For example, fitting a polynomial to the time
response of each pixel in a data cube.  This can greatly speed up the fitting
process, especially for linear models.
