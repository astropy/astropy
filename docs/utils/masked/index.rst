.. _utils-masked:

*****************************************************
Masked Values and Quantities (`astropy.utils.masked`)
*****************************************************

Often, data sets are incomplete or corrupted and it would be handy to be able
to mask certain values.  Astropy provides a |Masked| class to help represent
such data sets.

.. note:: |Masked| is similar to Numpy's :class:`~numpy.ma.MaskedArray`,
   but it supports subclasses much better and also has some important
   :ref:`differences in behaviour <utils-masked-vs-numpy-maskedarray>`.

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

You can mask or unmask individual elements by setting them to
`~numpy.ma.masked` or `~numpy.ma.nomask`::

  >>> mq.mask
  array([False, False,  True])
  >>> mq[:] = np.ma.nomask
  >>> mq[2] = np.ma.masked
  >>> mq.mask
  array([False, False,  True])
  >>> mq
  <MaskedQuantity [1., 2., ——] m>

These same procedures also work for higher-level classes like |Time| and
|SkyCoord|, which use |Masked| under the hood.

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

Numpy functions work as expected on |Masked| instances, with non-obvious
choices documented in `astropy.utils.masked.function_helpers` (please report
numpy functions that do not work properly with |Masked| values!). For example,
:func:`~astropy.utils.masked.function_helpers.nansum` does not propagate
masked elements, but instead replaces them with zero, and returns an unmasked
instance::

  >> np.nansum(ma, axis=-1)
  array([0., 5.])

.. _utils-masked-vs-numpy-maskedarray:

Differences from `~numpy.ma.MaskedArray`
========================================

|Masked| differs from `~numpy.ma.MaskedArray` in a number of ways, which we
detail below.  Overall, it may be helpful to think of |Masked| not as a
replacement of `~numpy.ma.MaskedArray`, but just as a way of marking bad
elements, as one might do without needing a different class by setting them to
NaN (not-a-number).  Like those NaN, the mask just propagates, except that for
some operations like taking the mean the equivalent of `~numpy.nanmean` is
used.

Values under masked are operated on
-----------------------------------

A difference in usage is that most operations act on the masked values,
i.e., no effort is made to preserve values.  For instance, compare::

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

Masked values are not skipped in reductions
-------------------------------------------

In reductions, the mask propagates as it would have if the operations were
done on the individual elements::

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

.. _utils-masked-mask-setting:

Setting the mask attribute replaces it
--------------------------------------

If one sets the mask attribute of a `~numpy.ma.MaskedArray`, it will
attempt to change the mask inplace::

  >>> np_ma = np.ma.MaskedArray([1., 2., 3.], mask=[False, True, False])
  >>> np_ma_mask_ref = np_ma.mask
  >>> np_ma.mask = False
  >>> np_ma_mask_ref
  array([False, False, False])

In contrast, if one sets the mask on a |Masked| class, it just sets it::

  >>> ma = Masked([1., 2., 3.], mask=[False, True, False])
  >>> ma_mask_ref = ma.mask
  >>> ma.mask = False
  >>> ma.mask
  array([False, False, False])
  >>> ma_mask_ref
  array([False,  True, False])

This has a consequence for setting the mask on a slice: for
`~numpy.ma.MaskedArray` it propagates back, but for |Masked| it does not::

  >>> np_ma = np.ma.MaskedArray([1., 2., 3.], mask=[False, True, False])
  >>> np_ma_view = np_ma[2:3]
  >>> np_ma_view.mask = True
  >>> np_ma_view
  masked_array(data=[--],
               mask=[ True],
         fill_value=1e+20,
              dtype=float64)
  >>> np_ma
  masked_array(data=[1.0, --, --],
               mask=[False,  True,  True],
         fill_value=1e+20)

  >>> ma = Masked([1., 2., 3.], mask=[False, True, False])
  >>> ma_view = ma[2:3]
  >>> ma_view.mask = True
  >>> ma_view
  MaskedNDArray([——])
  >>> ma
  MaskedNDArray([1., ——, 3.])

In order for the mask to be set in-place, one should do it explicitly::

  >>> ma[1:2].mask[...] = True
  >>> ma.mask
  array([False,  True, False])

The reason for not attempting to propagate is partially just that assignment
should be just that, assignment. But also that it is tricky to get right.
Indeed, also for `~numpy.ma.MaskedArray` it does not always work::

  >>> np_ma[0].mask = True
  Traceback (most recent call last):
  ...
  AttributeError: 'numpy.float64' object has no attribute 'mask'...

.. note:: We recommend not dealing with the mask directly but setting
  the instance to `~numpy.ma.masked` or `~numpy.ma.nomask`, as described
  above. This is also the only way to mask values for the
  higher-level classes such as |Time| and |SkyCoord|.

Numpy functions work as expected
--------------------------------

For `~numpy.ma.MaskedArray`, a number of regular numpy functions do not work
properly, and instead one has to use variants from the ``np.ma`` namespace.
For |Masked|, numpy functions do work as expected (but those under the
``np.ma`` namespace typically do not).

Masked subclasses behave like the subclass
------------------------------------------

A more conceptual difference is that for `~numpy.ma.MaskedArray`, the
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

In contrast, |Masked| is always wrapped around the data proper, i.e., a
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

Reference/API
=============

.. automodapi:: astropy.utils.masked

.. automodapi:: astropy.utils.masked.function_helpers
