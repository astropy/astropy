.. _astropy_nddata:

*****************************************
N-dimensional datasets (`astropy.nddata`)
*****************************************

Introduction
============

The `~astropy.nddata` package provides a uniform interface to N-dimensional
datasets in astropy through:

+ The `~astropy.nddata.NDDataBase` metaclass to define an astropy-wide
  interface to N-dimensional data sets while allowing flexibility in
  how those datasets are represented internally.
+ The `~astropy.nddata.NDData` class, which provides a basic container for
  N-dimensional datasets based on `~numpy.ndarray`-like data.
+ Several mixin classes for adding functionality to `~astropy.nddata.NDData`
  containers.
+ Classes for storing uncertainties, like `~astropy.nddata.StdDevUncertainty`
  for normal distributed uncertainties that are given as the sigma or standard
  deviation error. A more detailed description is avaiable in REFERENCE.
+ A decorator, `~astropy.nddata.support_nddata`, for facilitating use of
  `~astropy.nddata` objects  in functions in astropy and affiliated packages.
+ General utility functions (:ref:`nddata_utils`) for array operations.

.. warning::

  `~astropy.nddata` has changed significantly in astropy 1.0. See the section
  :ref:`nddata_transition` for more information.

.. warning::

  `~astropy.nddata` was subject of further changes in astropy xx. See the
  section :ref:`nddata_transition_2` for more information.

Getting started
===============

Of the classes provided by `~astropy.nddata`, the place to start for most
users will be `~astropy.nddata.NDData`, which by default uses a
`~numpy.ndarray`-like object to store the data and provides a fully functional
(if somewhat restricted) implementation. Designers of new classes should
also look at `~astropy.nddata.NDDataBase` before deciding what to subclass
from.

NDData
------

The primary purpose of `~astropy.nddata.NDData` is to act as a *container* for
data, metadata, and other related information like a mask.

An `~astropy.nddata.NDData` object can be instantiated by passing it an
N-dimensional Numpy array::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd = NDData(array)

or something that can be converted to an array::

    >>> ndd2 = NDData([1, 2, 3, 4])

It is also possible to initialize `~astropy.nddata.NDData` with more exotic
objects; see :ref:`nddata_details` for more information.

The underlying Numpy array can be accessed via the ``data`` attribute::

    >>> ndd.data
    array([[[ 0., 0., 0., ...
    ...

Values can be masked using the ``mask`` attribute::

    >>> ndd_masked = NDData(ndd2, mask = ndd2.data > 1.5)
    >>> ndd_masked.mask
    array([False,  True,  True,  True], dtype=bool)

A mask value of `True` indicates a value that should be ignored, while a mask
value of `False` indicates a valid value. Which is adopted from the
`~numpy.ma.MaskedArray` convention.

Similar attributes are available to store:

+ generic meta-data, in ``meta``,
+ a unit for the data values, in ``unit`` and
+ an uncertainty for the data values, in ``uncertainty``. Note that the
  ``uncertainty`` should have a string attribute called ``uncertainty_type``.
+ a world coordinate system in ``wcs`` providing a real-world meaning to the
  gridded coordinates.

Note that a `~astropy.nddata.NDData` object is not sliceable::

    >>> ndd2[1:3]        # doctest: +SKIP
    Traceback (most recent call last):
        ...
    TypeError: 'NDData' object has no attribute '__getitem__'


Mixins for additional functionality
-----------------------------------

Several classes are provided to add functionality to the basic
`~astropy.nddata.NDData` container. They include:

+ `~astropy.nddata.NDSlicingMixin` to handle slicing of N-dimensional data.
+ `~astropy.nddata.NDArithmeticMixin` to allow arithmetic operations on
  `~astropy.nddata.NDData` objects.
+ `~astropy.nddata.NDIOMixin` to use existing astropy functionality for input
  (with the method ``read``) and output (with the method ``write``).

To use these mixins, create a new class that includes the appropriate mixins
as subclasses. For example, to make a class that includes slicing, but not
arithmetic or I/O::

    >>> from astropy.nddata import NDData, NDSlicingMixin
    >>> class NDDataSliceable(NDSlicingMixin, NDData): pass

Note that the body of the class need not contain any code; all of the
functionality is provided by the ``NDData`` container and the mixins. The
order of the classes is important because python works from right to left in
determining the order in which methods are resolved.

``NDDataSliceable`` is initialized the same way that `~astropy.nddata.NDData`
is::

    >>> ndd_sliceable = NDDataSliceable([1, 2, 3, 4])

but can be sliced::

    >>> ndd_sliceable[1:3]
    NDDataSliceable([2, 3])

The class `~astropy.nddata.NDDataArray` is an example of a class which
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

The nddata package underwent substantial revision as a result of `APE 7`_;
please see that APE for an extensive discussion of the motivation and the
changes.

The most important changes are that:

+ ``NDData`` does not provide a numpy-like interface; to use its data use the
  ``data`` attribute instead.
+ Slicing is no provided in the base `~astropy.nddata.NDData`.
+ Arithmetic is no longer included in the base `~astropy.nddata.NDData` class.

Code that only uses the metadata features of `~astropy.nddata.NDData` should
not need to be modified.

Code that uses the arithemtic methods that used to be included in
`~astropy.nddata.NDData` and relied on it to behave like a numpy array should
instead subclass `~astropy.nddata.NDDataArray`; that class is equivalent to
the original `~astropy.nddata.NDData` class.

.. _nddata_transition_2:

Transition to astropy xx
========================

The `~astropy.nddata` package was subject to another revision mainly to fix
existing bugs but also to implemented new functionality while trying to keep
the previous API or only to extend it.

These changes included:

+ `~astropy.nddata.NDDataBase` now enforces no restrictions. Previous
  implementations forced
  the ``uncertainty`` to match certain requirements, these were dropped.
+ `~astropy.nddata.NDData` underwent several  bugfixes and documentation
  updates. But several
  small new functionalities are now included like creating an instance now
  allows for a ``copy`` parameter (default is ``False``) and is less
  restrictive on the ``uncertainty`` (only a Warning is issued if the
  uncertainty does not meet the requirements instead of an Exception).
+ `~astropy.nddata.NDSlicingMixin` is now more restrictive while slicing
  ``meta`` attributes
  like ``uncertainty``, ``wcs`` and ``mask``. But the Mixin was refactored in
  a way that allows subclasses to customize slicing attributes in a different
  way. Visit :ref:`nddata/subclassing` for more information.
+ `~astropy.nddata.NDArithmeticMixin` now allows to customize what attributes
  are computed
  during arithmetic operations. As with `~astropy.nddata.NDSlicingMixin` these
  are now
  arranged in a way that subclasses can easily manipulate the way one
  attribute is handled. Visit :ref:`nddata/subclassing` for more information.
+ `~astropy.nddata.NDUncertainty` underwent a major revision and is probably
  not backwards
  compatible while `~astropy.nddata.StdDevUncertainty` has kept most of it's
  API. For more
  information look at REFERENCE.
+ `~astropy.nddata.StdDevUncertainty` now handles ``units`` correctly and
  allows to compute
  the propagation with correlation. But only if the correlation is an input
  parameter. Evaluating the correlation itself is still not possible.


Using ``nddata``
================

.. toctree::
   :maxdepth: 2

   nddata.rst
   decorator.rst
   mixins/index.rst
   subclassing.rst
   utils.rst

Reference/API
=============

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:

.. automodapi:: astropy.nddata.utils
    :no-inheritance-diagram:

.. _APE 7: https://github.com/astropy/astropy-APEs/blob/master/APE7.rst
