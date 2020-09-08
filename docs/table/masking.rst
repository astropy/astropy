.. include:: references.txt

.. _masking_and_missing_values:

Masking and missing values
**************************

The `astropy.table` package provides support for masking and missing
values in a table by using the ``numpy.ma`` masked array package to
define masked columns and by supporting :ref:`mixin_columns` that provide masking.
This allows handling tables with missing or invalid entries in much
the same manner as for standard (unmasked) tables.  It
is useful to be familiar with the `masked array
<https://numpy.org/doc/stable/reference/maskedarray.generic.html>`_
documentation when using masked tables within `astropy.table`.

In a nutshell, the concept is to define a boolean mask that mirrors
the structure of a column data array.  Wherever a mask value is
`True`, the corresponding entry is considered to be missing or invalid.
Operations involving column or row access and slicing are unchanged.
The key difference is that arithmetic or reduction operations involving
columns or column slices follow the rules for `operations
on masked arrays
<https://numpy.org/doc/stable/reference/maskedarray.generic.html#operations-on-masked-arrays>`_.

.. Important:: Changes in astropy 4.0

   In astropy 4.0 the behavior of masked tables was changed in a way that
   could impact program functionality.  See :ref:`table-masked-4.0` for details.

.. Note::

   Reduction operations like `numpy.sum` or `numpy.mean` follow the
   convention of ignoring masked (invalid) values.  This differs from
   the behavior of the floating point ``NaN``, for which the sum of an
   array including one or more ``NaN's`` will result in ``NaN``.

   See `<https://numpy.org/neps/>`_ for information on NumPy
   Enhancement Proposals 24, 25, and 26.

Table creation
===============

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

The |MaskedColumn| is the masked analog of the |Column| class and
provides the interface for creating and manipulating a column of
masked data.  The |MaskedColumn| class inherits from
`numpy.ma.MaskedArray`, in contrast to |Column| which inherits from
`numpy.ndarray`.  This distinction is the main reason there are
different classes for these two cases.

Notice that masked entries in the table output are shown as ``--``.

**Create a table with one or more columns as a numpy MaskedArray**

  >>> import numpy as np
  >>> a = np.ma.array([1, 2])
  >>> b = [3, 4]
  >>> t = Table([a, b], names=('a', 'b'))

**Create a table from list data containing `numpy.ma.masked`**

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

Initializing from lists with embedded `numpy.ma.masked` elements is considerably
slower than using `numpy.ma.array` or |MaskedColumn| directly, so if performance
is a concern you should use the latter methods if possible.

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
a |MaskedColumn|.

  >>> Table([(1, 2), (3, 4)], names=('a', 'b'), masked=True, dtype=('i4', 'i8'))
  <Table masked=True length=2>
    a     b
  int32 int64
  ----- -----
      1     3
      2     4

Notice the table attributes ``mask`` and ``fill_value`` that are
available for a masked table.

**Convert an existing table to a masked table**

  >>> t = Table([[1, 2], ['x', 'y']])  # standard (unmasked) table
  >>> t = Table(t, masked=True, copy=False)  # convert to masked table

This operation will convert every |Column| to |MaskedColumn| and ensure that any
subsequently added columns are masked.

Table access
============

Nearly all the of standard methods for accessing and modifying data
columns, rows, and individual elements also apply to masked tables.

There are two minor differences for the |Row| object that is obtained by
indexing a single row of a table:

- For standard tables, two such rows can be compared for equality, but
  in masked tables this comparison will produce an exception.

Both of these differences are due to issues in the underlying
`numpy.ma.MaskedArray` implementation.

Masking and filling
====================

Both the |Table| and |MaskedColumn| classes provide
attributes and methods to support manipulating tables with missing or
invalid data.

Mask
----

The mask for a column can be viewed and modified via the ``mask`` attribute::

  >>> t = Table([(1, 2), (3, 4)], names=('a', 'b'), masked=True)
  >>> t['a'].mask = [False, True]  # Modify column mask (boolean array)
  >>> t['b'].mask = [True, False]  # Modify column mask (boolean array)
  >>> print(t)
   a   b
  --- ---
    1  --
   --   4

Masked entries are shown as ``--`` when the table is printed.  You can
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

To get the indices of masked elements use an expression like::

  >>> t['a'].mask.nonzero()[0]  # doctest: +SKIP
  array([1])


Filling
-------

The entries which are masked (i.e. missing or invalid) can be replaced
with specified fill values.  In this case the |MaskedColumn| or masked
|Table| will be converted to a standard |Column| or table. Each column
in a masked table has a ``fill_value`` attribute that specifies the
default fill value for that column.  To perform the actual replacement
operation the ``filled()`` method is called.  This takes an optional
argument which can override the default column ``fill_value``
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

.. _table-masked-4.0:

Masking change in astropy 4.0
=============================

In astropy 4.0 a change was introduced in the behavior of |Table| that impacts the
handling of masked columns.

Prior to 4.0, in order to include one or more |MaskedColumn| columns in a table, it was
required that *every* column be masked, even those with no missing or masked data.  This
was holdover from the original implementation of |Table| that used a numpy structured
array as the underlying container for the column data.  Since astropy 1.0 the |Table|
object is simply an ordered dictionary of columns (:ref:`table_implementation_details`)
and there is no requirement that column types be homogenous.

Starting with 4.0, a |Table| can contain both |Column| and |MaskedColumn| columns, and
by default the column type is determined solely by the data for each column.

The details of this change are discussed in the sections below.

.. Note::

   For most applications, even those with masked column data, we now recommend using
   the default |Table| behavior which allows heterogenous column types.  This implies
   creating tables *without* specifying the ``masked`` keyword argument.

Meaning of the ``masked`` table attribute
-----------------------------------------

The |Table| object has a ``masked`` attribute which determines the table behavior when
adding a new column:

- ``masked=True`` : non-mixin columns or data are always converted to |MaskedColumn|, and
  mixin columns have a ``mask`` attribute added if necessary.
- ``masked=False`` : each column is added based on the type or contents of the data.

The behavior associated with the ``masked`` attribute has *not changed* in version 4.0.
What has changed is that from 4.0 onward a table with ``masked=False`` may contain
|MaskedColumn| columns.

It is important to recognize that the ``masked`` attribute for a table does not imply
whether any of the column data are actually masked.  A table can have ``masked=True`` but
not have any masked elements in any table column.  Starting with version 4.0 there are two
table properties which give more useful information about masking:

- ``has_masked_columns`` : table has at least one |MaskedColumn| column.  This does *not*
  check if any data values are actually masked.
- ``has_masked_values`` : table has one or more column data values which are masked.
  This may be relatively slow for large tables as it requires checking the mask
  values of each column.

Starting with version 4.0 the term "masked table" should be reserved for the narrow and
less-common case of a table created with ``masked=True``.  In most cases there should be
no need worry about "masked" or "unmasked" at the table level, but instead focus on the
individual columns.

Auto-upgrade to masked
----------------------

Prior to version 4.0, adding a |MaskedColumn| or a new row with masked elements to a table
with ``masked=False`` would set ``masked=True`` and automatically "upgrade" other
columns to be masked.  In many cases this upgrade of the other columns was unnecessary and
an annoyance.

Starting with 4.0, new columns are added using the column type which is appropriate for
the data.  For instance, if a numpy masked array is added, then that will turn into a
|MaskedColumn|, but no other columns will be affected and the ``masked`` attribute will
remain as ``False``.

A commonly-encountered implication of this change is that tables read with
`~astropy.table.Table.read` will *always* have ``masked=False``, and only columns with
masked values will be |MaskedColumn|.  Prior to 4.0 if the input table had any masked
values then the returned table would have ``masked=True`` and all |MaskedColumn| columns.
An example is in the next section.

Recovering the pre-4.0 behavior
-------------------------------

For code that requires every existing or newly-added column to be masked, it is now
required to explicitly specify ``masked=True`` when creating the table.  Previously the
table would be auto-upgraded to use |MaskedColumn| for all columns as soon as the first
masked column was added.  If the table already exists (e.g. after using
`~astropy.table.Table.read` to read a data file), then one needs to make a new table:

.. doctest-skip::

  >> dat = Table.read('data.fits')
  >> dat = Table(dat, masked=True, copy=False)  # Convert to masked table
  >> dat['new_column'] = [1, 2, 3, 4, 5]  # Will be added as a MaskedColumn

For most applications this should not be necessary, and the preferred idiom is the more
explicit version below:

.. doctest-skip::

  >> dat = Table.read('data.fits')
  >> dat['new_column'] = np.ma.MaskedArray([1, 2, 3, 4, 5])
