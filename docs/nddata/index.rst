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
  gridded N-dimensional datasets.
+ Several mixin classes for adding functionality to `~astropy.nddata.NDData`
  containers.
+ A decorator, `~astropy.nddata.support_nddata`, for facilitating use of
  `~astropy.nddata` objects  in functions in astropy and affiliated packages.
+ General utility functions (:ref:`nddata_utils`) for array operations.

.. warning::

  `~astropy.nddata` has changed significantly in astropy 1.0. See the section
  :ref:`nddata_transition` for more information. Some more changes were done in
  astropy x.xx, see more details in :ref:`nddata_transition2`.

Getting started
===============

Of the classes provided by `~astropy.nddata`, the place to start for most
users will be `~astropy.nddata.NDData`, which by default uses a numpy array to
store the data. Designers of new classes should also look at
`~astropy.nddata.NDDataBase` before deciding what to subclass from.

NDData
------

The primary purpose of `~astropy.nddata.NDData` is to act as a *container* for
data, metadata, and other related information like a mask.

An `~astropy.nddata.NDData` object can be instantiated by passing it an
n-dimensional Numpy array::

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

     >>> ndd_masked = NDData(ndd, mask = ndd.data > 0.9)

A mask value of `True` indicates a value that should be ignored, while a mask
value of `False` indicates a valid value.


Similar attributes are available to store:

+ generic meta-data, in ``meta``,
+ a unit for the data values, in ``unit`` and
+ an uncertainty for the data values, in ``uncertainty``. Note that the
  ``uncertainty`` must have a attribute called ``uncertainty_type`` otherwise
  it will be converted to an `~astropy.nddata.UnknownUncertainty`.

Note that a `~astropy.nddata.NDData` object is not sliceable::

    >>> ndd2[1:3]        # doctest: +SKIP
    Traceback (most recent call last):
        ...
    TypeError: 'NDData' object has no attribute '__getitem__'



Mixins for additional functionality
-----------------------------------

Several classes are provided to add functionality to the basic ``NDData``
container. They include:

+ `~astropy.nddata.NDSlicingMixin` to handle slicing of N-dimensional data.
+ `~astropy.nddata.NDArithmeticMixin` to allow arithmetic operations on
  `~astropy.nddata.NDData` objects that include support propagation of
  uncertainties (in limited cases).
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

``NDDataSliceable`` is initialized the same way that `~astropy.nddata.NDData` is::

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

.. _nddata_transition2:

Transition to astropy x.xx
==========================

The implementation of `APE 7`_ where `~astropy.nddata.NDData` underwent a major
revision had some bugs and minor shortcomings;
So the interface was changed in some small ways and some bugs were fixed:

`~astropy.nddata.NDDataBase`:

+ `NDDataBase.uncertainty.getter`: is now abstract.
+ `NDDataBase.uncertainty.setter`: is moved to `~astropy.nddata.NDData`.

`~astropy.nddata.NDData`:

+ ``NDData`` now implements the ``uncertainty`` getter and setter from the
  former `~astropy.nddata.NDDataBase` with some minor modifications.

+ ``NDData.uncertainty.setter``:

  + wraps an ``uncertainty`` inside an `~astropy.nddata.UnknownUncertainty` if
    no ``uncertainty_type`` attribute is present in the uncertainty.
  + the uncertainty_type of the uncertainty must not be a string
    anymore. But it is still recommended.
  + did not set the ``parent_nddata`` of the ``uncertainty`` if the uncertainty
    was a subclass of `~astropy.nddata.NDUncertainty`. This is now fixed.

+ ``NDData.__init__``:

  + got an additional optional parameter ``copy`` which is
    by default False. If it is True all other parameters are copied before saving
    them as attributes. If False the parameters are only copied if there is no
    way of saving them as reference.
  + the ``data`` parameter allows now also ``Masked_Quantity`` objects
  + if the ``data`` is another instance or subclass of
    `~astropy.nddata.NDData` all implicit properties are checked because
    subclasses might implement another set of restrictions.
  + if an implicit and explicit attribute is provided, the
    explicit one is kept (without conversion) and a warning is issued. For
    example if ``data`` is a `~astropy.units.Quantity` and also a ``unit`` is
    given. The unit of the instance will be the ``unit`` parameter.
  + allowed a ``data`` parameter with a ``mask`` without checking that it had a
    ``data`` attribute.

Not backwards-compatible changes:

+ Subclasses of `~astropy.nddata.NDDataBase` need to implement an
  ``uncertainty`` getter and setter. Subclasses of `~astropy.nddata.NDData`
  are not affected.


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
