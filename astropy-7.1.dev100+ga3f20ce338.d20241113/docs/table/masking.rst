.. _masking_and_missing_values:

Masking and Missing Values
**************************

The `astropy.table` package provides support for masking and missing values in
a table by using the ``numpy.ma`` `masked array
<https://numpy.org/doc/stable/reference/maskedarray.html>`_ package to define
masked columns and by supporting :ref:`mixin_columns` that provide masking.
This allows handling tables with missing or invalid entries in much the same
manner as for standard (unmasked) tables. It is useful to be familiar with the
`masked array documentation
<https://numpy.org/doc/stable/reference/maskedarray.generic.html>`_
when using masked tables within `astropy.table`.

In a nutshell, the concept is to define a boolean mask that mirrors
the structure of a column data array. Wherever a mask value is
`True`, the corresponding entry is considered to be missing or invalid.
Operations involving column or row access and slicing are unchanged.
The key difference is that arithmetic or reduction operations involving
columns or column slices follow the rules for `operations
on masked arrays
<https://numpy.org/doc/stable/reference/maskedarray.generic.html#operations-on-masked-arrays>`_.

.. Note::

   Reduction operations like :func:`numpy.sum` or :func:`numpy.mean` follow the
   convention of ignoring masked (invalid) values. This differs from
   the behavior of the floating point ``NaN``, for which the sum of an
   array including one or more ``NaN's`` will result in ``NaN``.

   For more information see NumPy Enhancement Proposals `24
   <https://numpy.org/neps/nep-0024-missing-data-2.html>`_, `25
   <https://numpy.org/neps/nep-0025-missing-data-3.html>`_, and `26
   <https://numpy.org/neps/nep-0026-missing-data-summary.html>`_.

Table Creation
==============

A masked table can be created in several ways:

**Create a table with one or more columns as a MaskedColumn object**

  >>> from astropy.table import Table, Column, MaskedColumn
  >>> a = MaskedColumn([1, 2], name='a', mask=[False, True], dtype='i4')
  >>> b = Column([3, 4], name='b', dtype='i8')
  >>> Table([a, b])
  <Table length=2>
    a     b
  int32 int64
  ----- -----
      1     3
     --     4

The |MaskedColumn| is the masked analog of the |Column| class and provides the
interface for creating and manipulating a column of masked data. The
|MaskedColumn| class inherits from :class:`numpy.ma.MaskedArray`, in contrast
to |Column| which inherits from |ndarray|. This distinction is the main reason
there are different classes for these two cases.

Notice that masked entries in the table output are shown as ``--``.

**Create a table with one or more columns as a NumPy MaskedArray**

  >>> import numpy as np
  >>> a = np.ma.array([1, 2])
  >>> b = [3, 4]
  >>> t = Table([a, b], names=('a', 'b'))

**Create a table from list data containing numpy.ma.masked**

You can use the `numpy.ma.masked` constant to indicate masked or invalid data::

  >>> a = [1.0, np.ma.masked]
  >>> b = [np.ma.masked, 'val']
  >>> Table([a, b], names=('a', 'b'))
  <Table length=2>
    a     b
  float64 str3
  ------- ----
      1.0   --
      --  val

Initializing from lists with embedded `numpy.ma.masked` elements is
considerably slower than using :func:`numpy.ma.array` or |MaskedColumn|
directly, so if performance is a concern you should use the latter methods if
possible.

**Add a MaskedColumn object to an existing table**

  >>> t = Table([[1, 2]], names=['a'])
  >>> b = MaskedColumn([3, 4], mask=[True, False])
  >>> t['b'] = b

**Add a new row to an existing table and specify a mask argument**

  >>> a = Column([1, 2], name='a')
  >>> b = Column([3, 4], name='b')
  >>> t = Table([a, b])
  >>> t.add_row([3, 6], mask=[True, False])

**Create a new table object and specify masked=True**

If ``masked=True`` is provided when creating the table then every column will
be created as a |MaskedColumn|, and new columns will always be added as a
|MaskedColumn|.

  >>> Table([(1, 2), (3, 4)], names=('a', 'b'), masked=True, dtype=('i4', 'i8'))
  <Table masked=True length=2>
    a     b
  int32 int64
  ----- -----
      1     3
      2     4

**Convert an existing table to a masked table**

  >>> t = Table([[1, 2], ['x', 'y']])  # standard (unmasked) table
  >>> t = Table(t, masked=True, copy=False)  # convert to masked table

This operation will convert every |Column| to |MaskedColumn| and ensure that any
subsequently added columns are masked.

Table Access
============

Nearly all of the standard methods for accessing and modifying data
columns, rows, and individual elements also apply to masked tables.

There is a difference however regarding the |Row| objects that are obtained by
indexing a single row of a table. For standard tables, two such rows can be
compared for equality, but for masked tables this comparison will produce an
exception::

  >>> t[0] == t[1]
  Traceback (most recent call last):
  ...
  ValueError: Unable to compare rows for masked table due to numpy.ma bug

Masking and Filling
===================

Both the |Table| and |MaskedColumn| classes provide attributes and methods to
support manipulating tables with missing or invalid data.

Mask
----

.. EXAMPLE START: Manipulating Tables with Missing Data using Masks

The mask for a column can be viewed and modified via the ``mask`` attribute::

  >>> t = Table([(1, 2), (3, 4)], names=('a', 'b'), masked=True)
  >>> t['a'].mask = [False, True]  # Modify column mask (boolean array)
  >>> t['b'].mask = [True, False]  # Modify column mask (boolean array)
  >>> print(t)
   a   b
  --- ---
    1  --
   --   4

Masked entries are shown as ``--`` when the table is printed. You can
view the mask directly, either at the column or table level::

  >>> t['a'].mask
  array([False,  True]...)

  >>> t.mask
  <Table length=2>
    a     b
   bool  bool
  ----- -----
  False  True
   True False

To get the indices of masked elements, use an expression like::

  >>> t['a'].mask.nonzero()[0]  # doctest: +SKIP
  array([1])

.. EXAMPLE END

Filling
-------

.. EXAMPLE START: Manipulating Tables with Missing Data by Filling Masked Values

The entries which are masked (i.e., missing or invalid) can be replaced with
specified fill values. Filling a |MaskedColumn| produces a |Column|. Each
column in a masked table has a ``fill_value`` attribute that specifies the
default fill value for that column. To perform the actual replacement operation
the :meth:`~astropy.table.Table.filled` method is called. This takes an
optional argument which can override the default column ``fill_value``
attribute.
::

  >>> t['a'].fill_value = -99
  >>> t['b'].fill_value = 33

  >>> print(t.filled())
   a   b
  --- ---
    1  33
  -99   4

  >>> print(t['a'].filled())
   a
  ---
    1
  -99

  >>> print(t['a'].filled(999))
   a
  ---
    1
  999

  >>> print(t.filled(1000))
   a    b
  ---- ----
     1 1000
  1000    4

.. EXAMPLE END
