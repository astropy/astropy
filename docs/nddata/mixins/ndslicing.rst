.. _nddata_slicing:

Slicing and Indexing NDData
***************************

Introduction
============

This page only deals with peculiarities applying to
`~astropy.nddata.NDData`-like classes. For a tutorial about slicing/indexing see the
`python documentation <https://docs.python.org/3/tutorial/introduction.html#lists>`_
and `numpy documentation <https://docs.scipy.org/doc/numpy/reference/arrays.indexing.html>`_.

.. warning::
    `~astropy.nddata.NDData` and `~astropy.nddata.NDDataRef` enforce almost no
    restrictions on the properties so it might happen that some **valid but
    unusual** combination of properties always results in an IndexError or
    incorrect results. In this case see :ref:`nddata_subclassing` on how to
    customize slicing for a particular property.


Slicing NDDataRef
=================

Unlike `~astropy.nddata.NDData` the class `~astropy.nddata.NDDataRef`
implements slicing or indexing. The result will be wrapped inside the same
class as the sliced object.

Getting one element::

    >>> import numpy as np
    >>> from astropy.nddata import NDDataRef

    >>> data = np.array([1, 2, 3, 4])
    >>> ndd = NDDataRef(data)
    >>> ndd[1]
    NDDataRef(2)

Getting a sliced portion of the original::

    >>> ndd[1:3]  # Get element 1 (inclusive) to 3 (exclusive)
    NDDataRef([2, 3])

This will return a reference (and as such **not a copy**) of the original
properties so changing a slice will affect the original::

    >>> ndd_sliced = ndd[1:3]
    >>> ndd_sliced.data[0] = 5
    >>> ndd_sliced
    NDDataRef([5, 3])
    >>> ndd
    NDDataRef([1, 5, 3, 4])

except you indexed only one element (for example ``ndd_sliced = ndd[1]``). Then
the element is a scalar and changes will not propagate to the original.

Slicing NDDataRef including attributes
======================================

In case a ``wcs``, ``mask`` or ``uncertainty`` is present this attribute will
be sliced too::

    >>> from astropy.nddata import StdDevUncertainty
    >>> data = np.array([1, 2, 3, 4])
    >>> mask = data > 2
    >>> uncertainty = StdDevUncertainty(np.sqrt(data))
    >>> wcs = np.ones(4)
    >>> ndd = NDDataRef(data, mask=mask, uncertainty=uncertainty, wcs=wcs)
    >>> ndd_sliced = ndd[1:3]

    >>> ndd_sliced.data
    array([2, 3])

    >>> ndd_sliced.mask
    array([False,  True]...)

    >>> ndd_sliced.uncertainty  # doctest: +FLOAT_CMP
    StdDevUncertainty([1.41421356, 1.73205081])

    >>> ndd_sliced.wcs  # doctest: +FLOAT_CMP
    array([1., 1.])

but ``unit`` and ``meta`` will be unaffected.

If any of the attributes is set but doesn't implement slicing an info will be
printed and the property will be kept as is::

    >>> data = np.array([1, 2, 3, 4])
    >>> mask = False
    >>> uncertainty = StdDevUncertainty(0)
    >>> wcs = {'a': 5}
    >>> ndd = NDDataRef(data, mask=mask, uncertainty=uncertainty, wcs=wcs)
    >>> ndd_sliced = ndd[1:3]
    INFO: uncertainty cannot be sliced. [astropy.nddata.mixins.ndslicing]
    INFO: mask cannot be sliced. [astropy.nddata.mixins.ndslicing]
    INFO: wcs cannot be sliced. [astropy.nddata.mixins.ndslicing]

    >>> ndd_sliced.mask
    False

Example: Remove masked data
===========================

.. warning::
    If you are using a `~astropy.wcs.WCS` object as ``wcs`` this will **NOT**
    be possible. But you could work around it, i.e. set it to ``None`` before
    slicing.

By convention the ``mask`` attribute indicates if a point is valid or invalid.
So we are able to get all valid data points by slicing with the mask::

    >>> data = np.array([[1,2,3],[4,5,6],[7,8,9]])
    >>> mask = np.array([[0,1,0],[1,1,1],[0,0,1]], dtype=bool)
    >>> uncertainty = StdDevUncertainty(np.sqrt(data))
    >>> ndd = NDDataRef(data, mask=mask, uncertainty=uncertainty)
    >>> # don't forget that ~ or you'll get the invalid points
    >>> ndd_sliced = ndd[~ndd.mask]
    >>> ndd_sliced
    NDDataRef([1, 3, 7, 8])

    >>> ndd_sliced.mask
    array([False, False, False, False]...)

    >>> ndd_sliced.uncertainty  # doctest: +FLOAT_CMP
    StdDevUncertainty([1.        , 1.73205081, 2.64575131, 2.82842712])

or all invalid points::

    >>> ndd_sliced = ndd[ndd.mask] # without the ~ now!
    >>> ndd_sliced
    NDDataRef([2, 4, 5, 6, 9])

    >>> ndd_sliced.mask
    array([ True,  True,  True,  True,  True]...)

    >>> ndd_sliced.uncertainty  # doctest: +FLOAT_CMP
    StdDevUncertainty([1.41421356, 2.        , 2.23606798, 2.44948974, 3.        ])

.. note::
    The result of this kind of indexing (boolean indexing) will always be
    one-dimensional!
