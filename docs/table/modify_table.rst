.. include:: references.txt

.. _modify_table:

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

For more explicit control, the :meth:`~astropy.table.Table.add_column` and
:meth:`~astropy.table.Table.add_columns` methods can be used to add one or
multiple columns to a table. In both cases the new column(s) can be specified as
a list, an array (including |Column| or |MaskedColumn|), or a scalar::

  >>> from astropy.table import Column
  >>> t.add_column(np.arange(5), name='aa', index=0)  # Insert before first table column
  >>> t.add_column(1.0, name='bb')  # Add column of all 1.0 to end of table
  >>> c = Column(np.arange(5), name='e')
  >>> t.add_column(c, index=0)  # Add Column using the existing column name 'e'
  >>> t.add_columns([[1, 2, 3, 4, 5], ['v', 'w', 'x', 'y', 'z']], names=['h', 'i'])

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

  >>> t.remove_column('d1')
  >>> t.remove_columns(['aa', 'd2', 'e'])
  >>> del t['d3']
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

In-place versus replace column update
=====================================

Consider this code snippet::

  >>> t = Table([[1, 2, 3]], names=['a'])
  >>> t['a'] = [10.5, 20.5, 30.5]

There are a couple of ways this could be handled.  It could update the existing array
values in-place (truncating to integer), or it could replace the entire column with a new
column based on the supplied data values.

The answer for astropy (since version 1.3) is that the operation shown above does a *complete
replacement* of the column object.  In this case it makes a new column
object with float values by internally calling
``t.replace_column('a', [10.5, 20.5, 30.5])``.  In general this behavior
is more consistent with Python and Pandas behavior.

**Forcing in-place update**

It is straightforward to force an in-place update of a column as follows::

  t[colname][:] = value

**Finding the source of problems**

In order to find potential problems related to the replacing columns, there is a
configuration option ``table.conf.replace_warnings``.  This controls a set of warnings
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
``[]`` (no warnings).
