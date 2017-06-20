.. _modify_table:

.. include:: references.txt

Modifying a table
*****************

The data values within a |Table| object can be modified in much the same manner
as for `numpy` structured arrays by accessing columns or rows of data and
assigning values appropriately.  A key enhancement provided by the |Table| class
is the ability to easily modify the structure of the table: one can add or
remove columns, and add new rows of data.

Quick overview
==============

The code below shows the basics of modifying a table and its data.


**Make a table**
::

  >>> from astropy.table import Table
  >>> import numpy as np
  >>> arr = np.arange(15).reshape(5, 3)
  >>> t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})

**Modify data values**
::

  >>> t['a'][:] = [1, -2, 3, -4, 5]  # Set all column values
  >>> t['a'][2] = 30                 # Set row 2 of column 'a'
  >>> t[1] = (8, 9, 10)              # Set all row values
  >>> t[1]['b'] = -9                 # Set column 'b' of row 1
  >>> t[0:3]['c'] = 100              # Set column 'c' of rows 0, 1, 2

Note that ``table[row][column]`` assignments will not work with
`numpy` "fancy" ``row`` indexing (in that case ``table[row]`` would be
a *copy* instead of a *view*).  "Fancy" `numpy` indices include a
`list`, `numpy.ndarray`, or `tuple` of `numpy.ndarray` (e.g. the
return from `numpy.where`)::

  >>> t[[1, 2]]['a'] = [3., 5.]             # doesn't change table t
  >>> t[np.array([1, 2])]['a'] = [3., 5.]   # doesn't change table t
  >>> t[np.where(t['a'] > 3)]['a'] = 3.     # doesn't change table t

Instead use ``table[column][row]`` order::

  >>> t['a'][[1, 2]] = [3., 5.]
  >>> t['a'][np.array([1, 2])] = [3., 5.]
  >>> t['a'][np.where(t['a'] > 3)] = 3.

You can also modify data columns with ``unit`` set in a way that follows
the conventions of `~astropy.units.Quantity` by using the
:attr:`~astropy.table.Column.quantity` property::

  >>> from astropy import units as u
  >>> tu = Table([[1, 2.5]], names=('a',))
  >>> tu['a'].unit = u.m
  >>> tu['a'].quantity[:] = [1, 2] * u.km
  >>> tu['a']
  <Column name='a' dtype='float64' unit='m' length=2>
  1000.0
  2000.0

**Add a column or columns**

A single column can be added to a table using syntax like adding a dict value.
The value on the right hand side can be a list or array
of the correct size, or a scalar value that will be broadcast::

  >>> t['d1'] = np.arange(5)
  >>> t['d2'] = [1, 2, 3, 4, 5]
  >>> t['d3'] = 6  # all 5 rows set to 6

For more explicit control the :meth:`~astropy.table.Table.add_column` and
:meth:`~astropy.table.Table.add_columns` methods can be used to add one or multiple
columns to a table.  In both cases the new columns must be specified as |Column| or
|MaskedColumn| objects::

  >>> from astropy.table import Column
  >>> aa = Column(np.arange(5), name='aa')
  >>> t.add_column(aa, index=0)  # Insert before the first table column
  >>> bb = Column(np.arange(5))
  >>> t.add_column(bb, name='bb')  # Append unnamed column to the table with 'bb' as name

  # Make a new table with the same number of rows and add columns to original table
  >>> t2 = Table(np.arange(25).reshape(5, 5), names=('e', 'f', 'g', 'h', 'i'))
  >>> t.add_columns(t2.columns.values())

Finally, columns can also be added from
:class:`~astropy.units.Quantity` objects, which automatically sets the
``.unit`` attribute on the column:

  >>> from astropy import units as u
  >>> t['d'] = np.arange(1., 6.) * u.m
  >>> t['d']
  <Column name='d' dtype='float64' unit='m' length=5>
  1.0
  2.0
  3.0
  4.0
  5.0

**Remove columns**
::

  >>> t.remove_column('f')
  >>> t.remove_columns(['aa', 'd1', 'd2', 'd3', 'e'])
  >>> del t['g']
  >>> del t['h', 'i']
  >>> t.keep_columns(['a', 'b'])

**Replace a column**

One can entirely replace an existing column with a new column by setting the
column to any object that could be used to initialize a table column (e.g.  a
list or numpy array).  For example, one could change the data type of the ``a``
column from ``int`` to ``float`` using::

  >>> t['a'] = t['a'].astype(float)

If the right hand side value is not column-like, then an in-place update
using broadcasting will be done, e.g.::

  >>> t['a'] = 1  # Internally does t['a'][:] = 1

.. Note ::

   Prior to astropy version 1.3, assignment as shown above performed
   an in-place update of the existing column values and it was not possible
   to change the data type in this way.  Prior to 1.3 it was necessary
   to use the :meth:`~astropy.table.Table.replace_column` method in this case.
   See the section `API change in replacing columns`_ for additional information.

**Rename columns**
::

  >>> t.rename_column('a', 'a_new')
  >>> t['b'].name = 'b_new'

**Add a row of data**
::

  >>> t.add_row([-8, -9])

**Remove rows**
::

  >>> t.remove_row(0)
  >>> t.remove_rows(slice(4, 5))
  >>> t.remove_rows([1, 2])

**Sort by one more more columns**
::

  >>> t.sort('b_new')
  >>> t.sort(['a_new', 'b_new'])

**Reverse table rows**
::

  >>> t.reverse()

**Modify meta-data**
::

  >>> t.meta['key'] = 'value'

**Select or reorder columns**

A new table with a subset or reordered list of columns can be
created as shown in the following example::

  >>> t = Table(arr, names=('a', 'b', 'c'))
  >>> t_acb = t['a', 'c', 'b']

Another way to do the same thing is to provide a list or tuple
as the item as shown below::

  >>> new_order = ['a', 'c', 'b']  # List or tuple
  >>> t_acb = t[new_order]

Caveats
=======

Modifying the table data and properties is fairly straightforward.  One thing
to keep in mind is that adding a row *may* require a new copy in memory of the
table data.  This depends on the detailed layout of Python objects in memory
and cannot be reliably controlled.  In some cases it may be possible to build a
table row by row in less than O(N**2) time but you cannot count on it.

Another subtlety to keep in mind are cases where the return value of an
operation results in a new table in memory versus a view of the existing
table data.  As an example, imagine trying to set two table elements
using column selection with ``t['a', 'c']`` in combination with row index selection::

  >>> t = Table([[1, 2], [3, 4], [5, 6]], names=('a', 'b', 'c'))
  >>> t['a', 'c'][1] = (100, 100)
  >>> print(t)
   a   b   c
  --- --- ---
    1   3   5
    2   4   6

This might be surprising because the data values did not change and there
was no error.  In fact what happened is that ``t['a', 'c']`` created a
new temporary table in memory as a *copy* of the original and then updated
row 1 of the copy.  The original ``t`` table was unaffected and the new
temporary table disappeared once the statement was complete.  The takeaway
is to pay attention to how certain operations are performed one step at
a time.

.. _table-replace-1_3:

API change in replacing columns
===============================

Astropy version 1.3 introduces an API change in the way that the following
behaves::

  >>> t = Table([[1, 2, 3]], names=['a'])
  >>> t['a'] = [10.5, 20.5, 30.5]

Prior to 1.3 this always did an in-place replacement of the data values so that
the ``t['a']`` column object reference was maintained.  However, since the
original data type was integer in this case, the replaced values would silently
be converted to integer by truncation.

Starting with astropy 1.3 the operation shown above does a *complete
replacement* of the column object.  In this case it makes a new column
object with float values by internally calling
``t.replace_column('a', [10.5, 20.5, 30.5])``.  In general this behavior
is more consistent with Python and Pandas behavior, but there is potential
for somewhat subtle bugs in code that was written that expects the pre-1.3
in-place behavior.

**Examples**
::

  >>> t = Table([[1, 2, 3]], names=['a'])
  >>> t['a'].description = 'My data column'

  # Sliced column gets replaced
  >>> t2 = t[:2]  # Make a slice

  # In astropy 1.3 the following emits a warning about replacing a slice.
  >>> t2['a'] = [10, 20]  # doctest: +SKIP

  >>> list(t['a'])  # Outputs [10, 20, 3] prior to astropy 1.3.
  [1, 2, 3]

  # Column reference count changes
  >>> ta = t['a']  # Make a reference to the original column
  >>> t['a'] = [10, 20, 30]
  >>> t['a'] is ta  # Outputs True prior to astropy 1.3
  False

  # Column attributes change
  >>> print(t['a'].description)  # Outputs 'My data column' prior to astropy 1.3
  None

**Replicating pre-1.3 behavior**

If the pre-1.3 in-place behavior is required in code, it is straightforward
to achieve this.  Simply replace::

  t[colname] = value

with::

  t[colname][:] = value

As a *temporary* measure, or if the problematic code is in a package
that cannot be modified, one can modify the ``table.replace_inplace``
configuration variable::

  from astropy import table
  table.conf.replace_inplace = True

This will entirely revert to the pre-1.3 behavior.  This configuration option
will be deprecated and then subsequently removed in future releases, so it is
meant only as a stop-gap to provide time to appropriately update all code.

**Finding the source of problems**

In order to find potential problems related to the API change, the
configuration option ``table.conf.replace_warnings`` controls a set of warnings
that are emitted under certain circumstances when a table column is replaced.
This option must be set to a list that includes zero or more of the
following string values:

``always`` :
  Print a warning every time a column gets replaced via the
  setitem syntax (i.e. ``t['a'] = new_col``).

``slice`` :
  Print a warning when a column that appears to be a slice of
  a parent column is replaced.

``refcount`` :
  Print a warning when the Python reference count for the
  column changes.  This indicates that a stale object exists that might
  be used elsewhere in the code and give unexpected results.

``attributes`` :
  Print a warning if any of the standard column attributes changed.

The default value for the ``table.conf.replace_warnings`` option is
``['slice']``.
