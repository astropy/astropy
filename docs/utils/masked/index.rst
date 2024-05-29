.. _utils-masked:

**************************************
Masked Values (`astropy.utils.masked`)
**************************************

Often, data sets are incomplete or corrupted and it would be handy to be able
to mask certain values.  Astropy provides a |Masked| class to help represent
such data sets.

.. warning:: |Masked| is experimental! While we hope basic usage will remain
   similar, we are not yet sure whether it will not be necessary to change it
   to make things work throughout Astropy. This also means that comments and
   suggestions for improvements are especially welcome!

.. note:: |Masked| is similar to Numpy's :class:`~numpy.ma.MaskedArray`,
   but it supports subclasses much better and also has some important
   :ref:`differences in behaviour <utils-masked-vs-numpy-maskedarray>`.
   As a result, the behaviour of functions inside `numpy.ma` is poorly
   defined, and one should instead use regular ``numpy`` functions, which
   are overridden to work properly with masks (with non-obvious
   choices documented in `astropy.utils.masked.function_helpers`; please
   report numpy functions that do not work properly with |Masked| values!).

Usage
=====

Astropy |Masked| instances behave like `~numpy.ndarray` or subclasses such as
|Quantity| but with a mask associated, which is propagated in operations such
as addition, etc.::

  >>> import numpy as np
  >>> from astropy import units as u
  >>> from astropy.utils.masked import Masked
  >>> ma = Masked([1., 2., 3.], mask=[False, False, True])
  >>> ma
  MaskedNDArray([1., 2., ——])
  >>> mq = ma * u.m
  >>> mq + 25 * u.cm
  <MaskedQuantity [1.25, 2.25,  ———] m>

You can get the values without the mask using
`~astropy.utils.masked.Masked.unmasked`, or, if you need to control what
should be substituted for any masked values, with
:meth:`~astropy.utils.masked.Masked.filled`::

  >>> mq.unmasked
  <Quantity [1., 2., 3.] m>
  >>> mq.filled(fill_value=-75*u.cm)
  <Quantity [ 1.  ,  2.  , -0.75] m>

For reductions such as sums, the mask propagates as if the sum was
done directly::

  >>> ma = Masked([[0., 1.], [2., 3.]], mask=[[False, True], [False, False]])
  >>> ma.sum(axis=-1)
  MaskedNDArray([——, 5.])
  >>> ma.sum()
  MaskedNDArray(——)

You might wonder why masked elements are propagated, instead of just being
skipped (as is done in `~numpy.ma.MaskedArray`; see :ref:`below
<utils-masked-vs-numpy-maskedarray>`).  The rationale is that this leaves a
sum which is generally not useful unless one knows the number of masked
elements.  In contrast, for sample properties such as the mean, for which the
number of elements are counted, it seems natural to simply omit the masked
elements from the calculation::

  >> ma.mean(-1)
  MaskedNDArray([0.0, 2.5])


.. _utils-masked-vs-numpy-maskedarray:

Differences from `numpy.ma.MaskedArray`
=======================================

|Masked| differs from `~numpy.ma.MaskedArray` in a number of ways.  In usage,
a major difference is that most operations act on the masked values, i.e., no
effort is made to preserve values.  For instance, compare::

  >>> np_ma = np.ma.MaskedArray([1., 2., 3.], mask=[False, True, False])
  >>> (np_ma + 1).data
  array([2., 2., 4.])
  >>> (Masked(np_ma) + 1).unmasked
  array([2., 3., 4.])

The main reason for this decision is that for some masked subclasses, like
masked |Quantity|, keeping the original value makes no sense (e.g., consider
dividing a length by a time: if the unit of a masked quantity is changing, why
should its value not change?).  But it also helps to keep the implementation
considerably simpler, as the |Masked| class now primarily has to deal with
propagating the mask rather than deciding what to do with values.

A second difference is that for reductions, the mask propagates as it would
have if the operations were done on the individual elements::

  >>> np_ma.prod()
  np.float64(3.0)
  >>> np_ma[0] * np_ma[1] * np_ma[2]
  masked
  >>> Masked(np_ma).prod()
  MaskedNDArray(——)

The rationale for this becomes clear again by thinking about subclasses like a
masked |Quantity|.  For instance, consider an array ``s`` of lengths with
shape ``(N, 3)``, in which the last axis represents width, height, and depth.
With this, you could compute corresponding volumes by taking the product of
the values in the last axis, ``s.prod(axis=-1)``. But if masked elements were
skipped, the physical dimension of entries in the result would depend how many
elements were masked, which is something |Quantity| could not represent (and
would be rather surprising!).  As noted above, however, masked elements are
skipped for operations for which this is well defined, such as for getting the
mean and other sample properties such as the variance and standard deviation.

A third difference is more conceptual.  For `~numpy.ma.MaskedArray`, the
instance that is created is a masked version of the unmasked instance, i.e.,
`~numpy.ma.MaskedArray` remembers that is has wrapped a subclass like
|Quantity|, but does not share any of its methods.  Hence, even though the
resulting class looks reasonable at first glance, it does not work as expected::

  >>> q = [1., 2.] * u.m
  >>> np_mq = np.ma.MaskedArray(q, mask=[False, True])
  >>> np_mq
  masked_Quantity(data=[1.0, --],
                  mask=[False,  True],
            fill_value=1e+20)
  >>> np_mq.unit
  Traceback (most recent call last):
  ...
  AttributeError: 'MaskedArray' object has no attribute 'unit'...
  >>> np_mq / u.s
  <Quantity [1., 2.] 1 / s>

In contrast, |Masked| is always wrapped around the data properly, i.e., a
``MaskedQuantity`` is a quantity which has masked values, but with a unit that
is never masked.  Indeed, one can see this from the class hierarchy::

  >>> mq.__class__.__mro__
  (<class 'astropy.utils.masked.core.MaskedQuantity'>,
   <class 'astropy.units.quantity.Quantity'>,
   <class 'astropy.utils.masked.core.MaskedNDArray'>,
   <class 'astropy.utils.masked.core.Masked'>,
   <class 'astropy.utils.shapes.NDArrayShapeMethods'>,
   <class 'numpy.ndarray'>,
   <class 'object'>)

This choice has made the implementation much simpler: |Masked| only has to
worry about how to deal with masked values, while |Quantity| can worry just
about unit propagation, etc.  Indeed, an experiment showed that applying
|Masked| to `~astropy.table.Column` (which is a subclass of `~numpy.ndarray`),
the result is a new ``MaskedColumn`` that "just works", with no need for the
overrides and special-casing that were needed to make `~numpy.ma.MaskedArray`
work with `~astropy.table.Column`.  (Because the behaviour does change
somewhat, however, we chose not to replace the existing implementation.)

In some respects, rather than think of |Masked| as similar to
`~numpy.ma.MaskedArray`, it may be more useful to think of |Masked| as similar
to marking bad elements in arrays with NaN (not-a-number).  Like those NaN,
the mask just propagates, except that for some operations like taking the mean
the equivalence of `~numpy.nanmean` is used.

Reference/API
=============

.. automodapi:: astropy.utils.masked

.. automodapi:: astropy.utils.masked.function_helpers
