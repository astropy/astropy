.. _astropy_nddata:

*****************************************
N-dimensional datasets (`astropy.nddata`)
*****************************************

Introduction
============

`astropy.nddata` provides:

+ The `~astropy.nddata.NDDataBase` metaclass to define an astropy-wide
  interface to N-dimensional data sets while allowing flexibility in how those
  datasets are represented internally.
+ The `~astropy.nddata.NDData` class, which provides a basic container for
  N-dimensional datasets.
+ Several mixin classes for adding functionality to `~astropy.nddata.NDData`
  containers.

.. warning::

  `~astropy.nddata` has changed significantly in astropy 1.0. See the section :ref:`nddata_transition` for more information.

Getting started
===============

NDData
------

``NData`` is provided as a basic container to simplify handling of N-dimensional data in astropy and affiliated packages.

An `~astropy.nddata.NDData` object can be instantiated by passing it an
n-dimensional Numpy array::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd = NDData(array)

The underlying Numpy array can be accessed via the ``data`` attribute::

    >>> ndd.data
    array([[[ 0., 0., 0., ...
    ...

Values can be masked using the ``mask`` attribute::

     >>> ndd_masked = NDData(ndd, mask = ndd.data > 0.9)
     INFO: Overwriting NDData's current mask with specified mask [astropy.nddata.nddata]

A mask value of `True` indicates a value that should be ignored, while a mask
value of `False` indicates a valid value.

Similarly, attributes are available to store:

+ generic meta-data, in ``meta``,
+ a unit in ``unit`` and
+ an uncertainty, in ``uncertainty``. Note that the ``uncertainty`` must have
  a string attribute called ``uncertainty_type``.

Mixins for additional functionality
-----------------------------------

Several classes are provided to add functionality to the basic ``NDData`` container. They include:

+ `~astropy.nddata.NDSlicing` to handle slicing of N-dimensional data.
+ `~astropy.nddata.NDArithmetic` to allow arithmetic operations on
  `~astropy.nddata.NDData` objects that include support propagation of
  uncertainties (in limited cases).
+ Placeholder for IO mixin

To use these mixins, create a new class that includes the appropriate mixins as subclasses. For example, to make a class that includes slicing, but not arithmetic or I/O::

    >>> from astropy.nddata import NDData, NDSlicing
    >>> class NDDataSlicable(NDSlicing, NDData): pass

Note that the body of the class need not contain any code; all of the
functionality is provided by the ``NDData`` container and the mixins. The
order of the classes is important because python works from right to left in
determining the order in which methods are resolved.

The class `~astropy.nddata.NDDataArithmetic` is an example of a class which
utilizes mixins *and* adds functionality.

NDDataBase for making new subclasses
------------------------------------

`~astropy.nddata.NDDataBase` is a metaclass provided to support the creation
of objects that support the NDData interface but need the freedom to define
their own ways of storing data, unit, metadata and/or other properties. It
should be used instead of `~astropy.nddata.NDData` as the starting point for
any class for which the `~astropy.nddata.NDData` class is too restrictive.

.. _nddata_transition:

Transition to astropy 1.0
=========================

The most important changes are that:

+ ``NDData`` does not provide a numpy-like interface; to use its data use the
  ``data`` attribute instead.
+ Slicing is no provided in the base `~astropy.nddata.NDData`.
+ Arithmetic is no longer included in the base `~astropy.nddata.NDData` class.

Code that only uses the metadata features of `~astropy.nddata.NDData` should
not need to be modified. Code that uses the arithemtic methods that used to
be included in `~astropy.nddata.NDData` should instead subclass
`~astropy.nddata.NDDataArithmetic`; that class is equivalent to the original
`~astropy.nddata.NDData` class.


Using ``nddata``
================

.. toctree::
   :maxdepth: 2

   nddata.rst
   mixins/index.rst
   subclassing.rst
   decorator.rst

Reference/API
=============

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:

.. automodapi:: astropy.nddata.utils
    :no-inheritance-diagram:
