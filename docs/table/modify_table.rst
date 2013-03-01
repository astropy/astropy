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

The :func:`~astropy.table.table.add_column` and :func:`~astropy.table.table.add_columns`
functions can be used to add one or multiple columns to a table.  In both cases the new
columns must be specified as |Column| or |MaskedColumn| objects.
::

  >>> c = Column(np.arange(5), name='d')
  >>> t.add_column(c)

  # Make a new table with the same number of rows and add columns to original table
  >>> t2 = Table(np.arange(25).reshape(5, 5), names=('e', 'f', 'g', 'h', 'i'))
  >>> t.add_columns(t2.columns.values())

**Remove columns**
::

  >>> t.remove_column('f')
  >>> t.remove_columns(['d', 'e'])
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
