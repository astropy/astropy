.. _modify_table:

.. include:: references.txt

Modifying a table
-----------------

The data values within a |Table| object can be modified in much the same manner
as for `numpy` structured arrays by accessing columns or rows of data and
assigning values appropriately.  A key enhancement provided by the |Table| class
is the ability to easily modify the structure of the table: one can add or
remove columns, and add new rows of data.

Quick overview
^^^^^^^^^^^^^^

The code below shows the basics of modifying a table and its data.


**Make a table**
::

  >>> from astropy.table import Table, Column
  >>> import numpy as np
  >>> arr = np.arange(15).reshape(5, 3)
  >>> t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})

**Modify data values**
::

  >>> t['a'] = [1, -2, 3, -4, 5]  # Set all column values
  >>> t['a'][2] = 30              # Set row 2 of column 'a'
  >>> t[1] = (8, 9, 10)           # Set all row values
  >>> t[1]['b'] = -9              # Set column 'b' of row 1
  >>> t[0:3]['c'] = 100           # Set column 'c' of rows 0, 1, 2

**Add a column or columns**

A single column can be added to a table using syntax like adding a dict value.
The value on the right hand side can be a list or array
of the correct size, or a scalar value that will be broadcast::

  >>> t['d1'] = np.arange(5)
  >>> t['d2'] = [1, 2, 3, 4, 5]
  >>> t['d3'] = 6  # all 5 rows set to 6

For more explicit control the :func:`~astropy.table.table.add_column` and
:func:`~astropy.table.table.add_columns` functions can be used to add one or multiple
columns to a table.  In both cases the new columns must be specified as |Column| or
|MaskedColumn| objects with the ``name`` defined::

  >>> aa = Column(np.arange(5), name='aa')
  >>> t.add_column(aa, index=0)  # Insert before the first table column

  # Make a new table with the same number of rows and add columns to original table
  >>> t2 = Table(np.arange(25).reshape(5, 5), names=('e', 'f', 'g', 'h', 'i'))
  >>> t.add_columns(t2.columns.values())

**Remove columns**
::

  >>> t.remove_column('f')
  >>> t.remove_columns(['aa', 'd1', 'd2', 'd3', 'e'])
  >>> del t['g']
  >>> del t['h', 'i']
  >>> t.keep_columns(['a', 'b'])

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

**Reorder columns**

.. note::

  In this example, it is important that `neworder` is a tuple, and not a
  list, slice, or `~numpy.ndarray`.

::

  >>> t_acb = t['a','c','b']
  >>> neworder = ('a','c','b')
  >>> t_acb = t[neworder]


Caveats
^^^^^^^

Modifying the table data and properties is fairly straightforward.  There are
only a few things to keep in mind:

- The data type for a column cannot be changed in place.  In order to do this
  you must make a copy of the table with the column type changed appropriately.
- Adding or removing a column will generate a new copy
  in memory of all the data.  If the table is very large this may be slow.
- Adding a row *may* require a new copy in memory of the table data.  This
  depends on the detailed layout of Python objects in memory and cannot be
  reliably controlled.  In some cases it may be possible to build a table
  row by row in less than O(N**2) time but you cannot count on it.

Another subtlety to keep in mind are cases where the return value of an
operation results in a new table in memory versus a view of the existing
table data.  As an example, imagine trying to set two table elements
using column selection with ``t['a', 'c']`` in combination with row index selection::

  >>> t = Table([[1, 2], [3, 4], [5, 6]], names=('a', 'b', 'c'))
  >>> t['a', 'c'][1] = (100, 100)
  >>> print t
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

