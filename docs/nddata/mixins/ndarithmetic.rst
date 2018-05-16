.. _nddata_arithmetic:

NDData Arithmetic
*****************

Introduction
============

`~astropy.nddata.NDDataRef` implements the following arithmetic operations:

- addition: :meth:`~astropy.nddata.NDArithmeticMixin.add`
- subtraction: :meth:`~astropy.nddata.NDArithmeticMixin.subtract`
- multiplication: :meth:`~astropy.nddata.NDArithmeticMixin.multiply`
- division: :meth:`~astropy.nddata.NDArithmeticMixin.divide`

Using basic arithmetic methods
==============================

Using the standard arithmetic methods requires that the first operand
is an `~astropy.nddata.NDDataRef` instance

    >>> from astropy.nddata import NDDataRef
    >>> import numpy as np
    >>> ndd1 = NDDataRef([1, 2, 3, 4])

while the requirement for the second operand is simply: It must be convertible
to the first operand. It can be a number::

    >>> ndd1.add(3)
    NDDataRef([4, 5, 6, 7])

or a `list`::

    >>> ndd1.subtract([1,1,1,1])
    NDDataRef([0, 1, 2, 3])

a `numpy.ndarray`::

    >>> ndd1.multiply(np.arange(4, 8))
    NDDataRef([ 4, 10, 18, 28])
    >>> ndd1.divide(np.arange(1,13).reshape(3,4))  # a 3 x 4 numpy array  # doctest: +FLOAT_CMP
    NDDataRef([[1.        , 1.        , 1.        , 1.        ],
               [0.2       , 0.33333333, 0.42857143, 0.5       ],
               [0.11111111, 0.2       , 0.27272727, 0.33333333]])

here broadcasting takes care of the different dimensions. Also several other
classes are possible.

Using arithmetic classmethods
=============================

Here both operands don't need to be `~astropy.nddata.NDDataRef`-like::

    >>> NDDataRef.add(1, 3)
    NDDataRef(4)

or to wrap the result of an arithmetic operation between two Quantities::

    >>> import astropy.units as u
    >>> ndd = NDDataRef.multiply([1,2] * u.m, [10, 20] * u.cm)
    >>> ndd  # doctest: +FLOAT_CMP
    NDDataRef([10., 40.])
    >>> ndd.unit
    Unit("cm m")

or taking the inverse of a `~astropy.nddata.NDDataRef` object::

    >>> NDDataRef.divide(1, ndd1)  # doctest: +FLOAT_CMP
    NDDataRef([1.        , 0.5       , 0.33333333, 0.25      ])


Possible operands
-----------------

The possible types of input for operands are:

+ scalars of any type
+ lists containing numbers (or nested lists)
+ numpy arrays
+ numpy masked arrays
+ astropy quantities
+ other nddata classes or subclasses

Advanced options
================

The normal python operators ``+``, ``-``, ... are not implemented because
the methods provide several options how to proceed with the additional
attributes.

data, unit
----------

For ``data`` and ``unit`` there are no parameters. Every arithmetic
operation lets the `astropy.units.Quantity`-framework evaluate the result
or fail and abort the operation.

Adding two NDData objects with the same unit works::

    >>> ndd1 = NDDataRef([1,2,3,4,5], unit='m')
    >>> ndd2 = NDDataRef([100,150,200,50,500], unit='m')

    >>> ndd = ndd1.add(ndd2)
    >>> ndd.data  # doctest: +FLOAT_CMP
    array([101., 152., 203.,  54., 505.])
    >>> ndd.unit
    Unit("m")

Adding two NDData objects with compatible units also works::

    >>> ndd1 = NDDataRef(ndd1, unit='pc')
    INFO: overwriting NDData's current unit with specified unit. [astropy.nddata.nddata]
    >>> ndd2 = NDDataRef(ndd2, unit='lyr')
    INFO: overwriting NDData's current unit with specified unit. [astropy.nddata.nddata]

    >>> ndd = ndd1.subtract(ndd2)
    >>> ndd.data  # doctest: +FLOAT_CMP
    array([ -29.66013938,  -43.99020907,  -58.32027876,  -11.33006969,
           -148.30069689])
    >>> ndd.unit
    Unit("pc")

this will keep by default the unit of the first operand. However units will
not be decomposed during division::

    >>> ndd = ndd2.divide(ndd1)
    >>> ndd.data  # doctest: +FLOAT_CMP
    array([100.        ,  75.        ,  66.66666667,  12.5       , 100.        ])
    >>> ndd.unit
    Unit("lyr / pc")

mask
----

The ``handle_mask`` parameter for the arithmetic operations implements what the
resulting mask will be. There are several options.

- ``None``, the result will have no ``mask``::

      >>> ndd1 = NDDataRef(1, mask=True)
      >>> ndd2 = NDDataRef(1, mask=False)
      >>> ndd1.add(ndd2, handle_mask=None).mask is None
      True

- ``"first_found"`` or ``"ff"``, the result will have the mask of the first
  operand or if that is None the mask of the second operand::

      >>> ndd1 = NDDataRef(1, mask=True)
      >>> ndd2 = NDDataRef(1, mask=False)
      >>> ndd1.add(ndd2, handle_mask="first_found").mask
      True
      >>> ndd3 = NDDataRef(1)
      >>> ndd3.add(ndd2, handle_mask="first_found").mask
      False

- a function (or an arbitary callable) that takes at least two arguments.
  For example `numpy.logical_or` is the default::

      >>> ndd1 = NDDataRef(1, mask=np.array([True, False, True, False]))
      >>> ndd2 = NDDataRef(1, mask=np.array([True, False, False, True]))
      >>> ndd1.add(ndd2).mask
      array([ True, False,  True,  True]...)

  This defaults to ``"first_found"`` in case only one ``mask`` is not None::

      >>> ndd1 = NDDataRef(1)
      >>> ndd2 = NDDataRef(1, mask=np.array([True, False, False, True]))
      >>> ndd1.add(ndd2).mask
      array([ True, False, False,  True]...)

  Custom functions are also possible::

      >>> def take_alternating_values(mask1, mask2, start=0):
      ...     result = np.zeros(mask1.shape, dtype=np.bool)
      ...     result[start::2] = mask1[start::2]
      ...     result[start+1::2] = mask2[start+1::2]
      ...     return result

  This function is non-sense but let's see how it performs::

      >>> ndd1 = NDDataRef(1, mask=np.array([True, False, True, False]))
      >>> ndd2 = NDDataRef(1, mask=np.array([True, False, False, True]))
      >>> ndd1.add(ndd2, handle_mask=take_alternating_values).mask
      array([ True, False,  True,  True]...)

  Additional parameters can be given by prefixing them with ``mask_``
  (which will be stripped before passing it to the function)::

      >>> ndd1.add(ndd2, handle_mask=take_alternating_values, mask_start=1).mask
      array([False, False, False, False]...)
      >>> ndd1.add(ndd2, handle_mask=take_alternating_values, mask_start=2).mask
      array([False, False,  True,  True]...)

meta
----

The ``handle_meta`` parameter for the arithmetic operations implements what the
resulting meta will be. The options are the same as for the ``mask``:

- If ``None`` the resulting ``meta`` will be an empty `collections.OrderedDict`.

      >>> ndd1 = NDDataRef(1, meta={'object': 'sun'})
      >>> ndd2 = NDDataRef(1, meta={'object': 'moon'})
      >>> ndd1.add(ndd2, handle_meta=None).meta
      OrderedDict()

  For ``meta`` this is the default so you don't need to pass it in this case::

      >>> ndd1.add(ndd2).meta
      OrderedDict()

- If ``"first_found"`` or ``"ff"`` the resulting meta will be the meta of the
  first operand or if that contains no keys the meta of the second operand is
  taken.

      >>> ndd1 = NDDataRef(1, meta={'object': 'sun'})
      >>> ndd2 = NDDataRef(1, meta={'object': 'moon'})
      >>> ndd1.add(ndd2, handle_meta='ff').meta
      {'object': 'sun'}

- If it's a ``callable`` it must take at least two arguments. Both ``meta``
  attributes will be passed to this function (even if one or both of them are
  empty) and the callable evaluates the result's meta. For example just a
  function that merges these two::

      >>> # It's expected with arithmetics that the result is not a reference,
      >>> # so we need to copy
      >>> from copy import deepcopy

      >>> def combine_meta(meta1, meta2):
      ...     if not meta1:
      ...         return deepcopy(meta2)
      ...     elif not meta2:
      ...         return deepcopy(meta1)
      ...     else:
      ...         meta_final = deepcopy(meta1)
      ...         meta_final.update(meta2)
      ...         return meta_final

      >>> ndd1 = NDDataRef(1, meta={'time': 'today'})
      >>> ndd2 = NDDataRef(1, meta={'object': 'moon'})
      >>> ndd1.subtract(ndd2, handle_meta=combine_meta).meta # doctest: +SKIP
      {'object': 'moon', 'time': 'today'}

  Here again additional arguments for the function can be passed in using
  the prefix ``meta_`` (which will be stripped away before passing it to this)
  function. See the description for the mask-attribute for further details.

wcs
^^^

The ``compare_wcs`` argument will determine what the result's ``wcs`` will be
or if the operation should be forbidden. The possible values are identical to
``mask`` and ``meta``:

- If ``None`` the resulting ``wcs`` will be an empty ``None``.

      >>> ndd1 = NDDataRef(1, wcs=0)
      >>> ndd2 = NDDataRef(1, wcs=1)
      >>> ndd1.add(ndd2, compare_wcs=None).wcs is None
      True

- If ``"first_found"`` or ``"ff"`` the resulting wcs will be the wcs of the
  first operand or if that is None the meta of the second operand is
  taken.

      >>> ndd1 = NDDataRef(1, wcs=1)
      >>> ndd2 = NDDataRef(1, wcs=0)
      >>> ndd1.add(ndd2, compare_wcs='ff').wcs
      1

- If it's a ``callable`` it must take at least two arguments. Both ``wcs``
  attributes will be passed to this function (even if one or both of them are
  None) and the callable should return ``True`` if these wcs are identical
  (enough) to allow the arithmetic operation or ``False`` if the arithmetic
  operation should be aborted with a ``ValueError``. If ``True`` the ``wcs``
  are identical and the first one is used for the result::

      >>> def compare_wcs_scalar(wcs1, wcs2, allowed_deviation=0.1):
      ...     if wcs1 is None and wcs2 is None:
      ...         return True  # both have no WCS so they are identical
      ...     if wcs1 is None or wcs2 is None:
      ...         return False  # one has WCS, the other doesn't not possible
      ...     else:
      ...         return abs(wcs1 - wcs2) < allowed_deviation

      >>> ndd1 = NDDataRef(1, wcs=1)
      >>> ndd2 = NDDataRef(1, wcs=1)
      >>> ndd1.subtract(ndd2, compare_wcs=compare_wcs_scalar).wcs
      1

  Additional arguments can be passed in prefixing them with ``wcs_`` (this
  prefix will be stripped away before passing it to the function)::

      >>> ndd1 = NDDataRef(1, wcs=1)
      >>> ndd2 = NDDataRef(1, wcs=2)
      >>> ndd1.subtract(ndd2, compare_wcs=compare_wcs_scalar, wcs_allowed_deviation=2).wcs
      1

  If one is using `~astropy.wcs.WCS` objects a very handy function to use might
  be::

      >>> def wcs_compare(wcs1, wcs2, *args, **kwargs):
      ...     return wcs1.wcs.compare(wcs2.wcs, *args, **kwargs)

  see :meth:`astropy.wcs.Wcsprm.compare` for the arguments this comparison
  allows.

uncertainty
-----------

The ``propagate_uncertainties`` argument can be used to turn the propagation
of uncertainties on or off.

- If ``None`` the result will have no uncertainty::

      >>> from astropy.nddata import StdDevUncertainty
      >>> ndd1 = NDDataRef(1, uncertainty=StdDevUncertainty(0))
      >>> ndd2 = NDDataRef(1, uncertainty=StdDevUncertainty(1))
      >>> ndd1.add(ndd2, propagate_uncertainties=None).uncertainty is None
      True

- If ``False`` the result will have the first found uncertainty.

  .. note::
      Setting ``propagate_uncertainties=False`` is not generally not
      recommended.

- If ``True`` both uncertainties must be ``NDUncertainty`` subclasses that
  implement propagation. This is possible for
  `~astropy.nddata.StdDevUncertainty`::

      >>> ndd1 = NDDataRef(1, uncertainty=StdDevUncertainty([10]))
      >>> ndd2 = NDDataRef(1, uncertainty=StdDevUncertainty([10]))
      >>> ndd1.add(ndd2, propagate_uncertainties=True).uncertainty  # doctest: +FLOAT_CMP
      StdDevUncertainty([14.14213562])

uncertainty with correlation
----------------------------

If ``propagate_uncertainties`` is ``True`` you can give also an argument
for ``uncertainty_correlation``. `~astropy.nddata.StdDevUncertainty` cannot
keep track of it's correlations by itself but it can evaluate the correct
resulting uncertainty if the correct ``correlation`` is given.

The default (``0``) represents uncorrelated while ``1`` means correlated and
``-1`` anti-correlated. If given a `numpy.ndarray` it should represent the
element-wise correlation coefficient.

For example without correlation subtracting a `~astropy.nddata.NDDataRef`
instance from itself results in a non-zero uncertainty::

    >>> ndd1 = NDDataRef(1, uncertainty=StdDevUncertainty([10]))
    >>> ndd1.subtract(ndd1, propagate_uncertainties=True).uncertainty  # doctest: +FLOAT_CMP
    StdDevUncertainty([14.14213562])

Given a correlation of ``1`` because they clearly correlate gives the
correct uncertainty of ``0``::

    >>> ndd1 = NDDataRef(1, uncertainty=StdDevUncertainty([10]))
    >>> ndd1.subtract(ndd1, propagate_uncertainties=True,
    ...               uncertainty_correlation=1).uncertainty  # doctest: +FLOAT_CMP
    StdDevUncertainty([0.])

which would be consistent with the equivalent operation ``ndd1 * 0``::

    >>> ndd1.multiply(0, propagate_uncertainties=True).uncertainty
    StdDevUncertainty([0])

.. warning::
    The user needs to calculate or know the appropriate value or array manually
    and pass it to ``uncertainty_correlation``. The implementation follows
    general first order error propagation formulas, see for example:
    `Wikipedia <https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas>`_.

You can also give element-wise correlations::

    >>> ndd1 = NDDataRef([1,1,1,1], uncertainty=StdDevUncertainty([1,1,1,1]))
    >>> ndd2 = NDDataRef([2,2,2,2], uncertainty=StdDevUncertainty([2,2,2,2]))
    >>> ndd1.add(ndd2,uncertainty_correlation=np.array([1,0.5,0,-1])).uncertainty  # doctest: +FLOAT_CMP
    StdDevUncertainty([3.        , 2.64575131, 2.23606798, 1.        ])

The correlation ``np.array([1, 0.5, 0, -1])`` would indicate that the first
element is fully correlated, the second element partially correlates while
element 3 is uncorrelated and 4 is anti-correlated.

uncertainty with unit
---------------------

`~astropy.nddata.StdDevUncertainty` implements correct error propagation even
if the unit of the data differs from the unit of the uncertainty::

    >>> ndd1 = NDDataRef([10], unit='m', uncertainty=StdDevUncertainty([10], unit='cm'))
    >>> ndd2 = NDDataRef([20], unit='m', uncertainty=StdDevUncertainty([10]))
    >>> ndd1.subtract(ndd2, propagate_uncertainties=True).uncertainty  # doctest: +FLOAT_CMP
    StdDevUncertainty([10.00049999])

but it needs to be convertible to the unit for the data.
