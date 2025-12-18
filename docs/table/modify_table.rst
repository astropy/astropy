.. _modify_table:

Modifying a Table
*****************

The data values within a |Table| object can be modified in much the same manner
as for ``numpy`` `structured arrays
<https://numpy.org/doc/stable/user/basics.rec.html>`_ by accessing columns or
rows of data and assigning values appropriately. A key enhancement provided by
the |Table| class is the ability to modify the structure of the table: you can
add or remove columns, and add new rows of data.

Quick Overview
==============

The code below shows the basics of modifying a table and its data.

Examples
--------

.. EXAMPLE START: Making a Table and Modifying Data

**Make a table**
::

  >>> from astropy.table import Table
  >>> import numpy as np
  >>> arr = np.arange(15).reshape(5, 3)
  >>> t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})

**Modify data values**
::

  >>> t['a'][:] = [1, -2, 3, -4, 5]  # Set all values of column 'a'
  >>> t['a'][2] = 30                 # Set row 2 of column 'a'
  >>> t[1] = (8, 9, 10)              # Set all values of row 1
  >>> t[1]['b'] = -9                 # Set column 'b' of row 1
  >>> t[0:3]['c'] = 100              # Set column 'c' of rows 0, 1, 2

Note that ``table[row][column]`` assignments will not work with ``numpy``
"fancy" ``row`` indexing (in that case ``table[row]`` would be a *copy* instead
of a *view*). "Fancy" ``numpy`` indices include a :class:`list`, |ndarray|, or
:class:`tuple` of |ndarray| (e.g., the return from :func:`numpy.where`)::

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

.. note::

  The best way to combine the functionality of the |Table| and |Quantity|
  classes is to use a |QTable|. See :ref:`quantity_and_qtable` for more
  information.

.. EXAMPLE END

**Add a column or columns**

.. EXAMPLE START: Adding Columns to Tables

A single column can be added to a table using syntax like adding a key-value
pair to a :class:`dict`. The value on the right hand side can be a
:class:`list` or |ndarray| of the correct size, or a scalar value that will be
`broadcast <https://numpy.org/doc/stable/user/basics.broadcasting.html>`_::

  >>> t['d1'] = np.arange(5)
  >>> t['d2'] = [1, 2, 3, 4, 5]
  >>> t['d3'] = 6  # all 5 rows set to 6

For more explicit control, the :meth:`~astropy.table.Table.add_column` and
:meth:`~astropy.table.Table.add_columns` methods can be used to add one or
multiple columns to a table. In both cases the new column(s) can be specified as
a :class:`list`, |ndarray|, |Column|, |MaskedColumn|, or a scalar::

  >>> from astropy.table import Column
  >>> t.add_column(np.arange(5), name='aa', index=0)  # Insert before first table column
  >>> t.add_column(1.0, name='bb')  # Add column of all 1.0 to end of table
  >>> c = Column(np.arange(5), name='e')
  >>> t.add_column(c, index=0)  # Add Column using the existing column name 'e'
  >>> t.add_columns([[1, 2, 3, 4, 5], ['v', 'w', 'x', 'y', 'z']], names=['h', 'i'])

Finally, columns can also be added from |Quantity| objects, which automatically
sets the ``unit`` attribute on the column (but you might find it more
convenient to add a |Quantity| to a |QTable| instead, see
:ref:`quantity_and_qtable` for details)::

  >>> from astropy import units as u
  >>> t['d'] = np.arange(1., 6.) * u.m
  >>> t['d']
  <Column name='d' dtype='float64' unit='m' length=5>
  1.0
  2.0
  3.0
  4.0
  5.0

.. EXAMPLE END

**Remove columns**

.. EXAMPLE START: Removing Columns from Tables

To remove a column from a table::

  >>> t.remove_column('d1')
  >>> t.remove_columns(['aa', 'd2', 'e'])
  >>> del t['d3']
  >>> del t['h', 'i']
  >>> t.keep_columns(['a', 'b'])

.. EXAMPLE END

**Replace a column**

.. EXAMPLE START: Replacing Columns in Tables

You can entirely replace an existing column with a new column by setting the
column to any object that could be used to initialize a table column (e.g.,  a
:class:`list` or |ndarray|). For example, you could change the data type of the
``a`` column from ``int`` to ``float`` using::

  >>> t['a'] = t['a'].astype(float)

If the right-hand side value is not column-like, then an in-place update using
`broadcasting <https://numpy.org/doc/stable/user/basics.broadcasting.html>`_
will be done, for example::

  >>> t['a'] = 1  # Internally does t['a'][:] = 1

.. EXAMPLE END

**Perform a dictionary-style update**

It is possible to perform a dictionary-style update, which adds new columns to
the table and replaces existing ones::

  >>> t1 = Table({'name': ['foo', 'bar'], 'val': [0., 0.]}, meta={'n': 2})
  >>> t2 = Table({'val': [1., 2.], 'val2': [10., 10.]}, meta={'id': 0})
  >>> t1 |= t2
  >>> t1
  <Table length=2>
  name   val     val2
  str3 float64 float64
  ---- ------- -------
   foo     1.0    10.0
   bar     2.0    10.0

When using ``|=``, the other object does not need to be a |Table|, it can be
anything that can be used for :ref:`construct_table` with a compatible number
of rows::

  >>> t1 = Table({'name': ['foo', 'bar'], 'val': [0., 0.]}, meta={'n': 2})
  >>> d = dict({'val': [1., 2.], 'val2': [10., 10.]})
  >>> t1 |= d
  >>> t1
  <Table length=2>
  name   val     val2
  str3 float64 float64
  ---- ------- -------
   foo     1.0    10.0
   bar     2.0    10.0

It is also possible to use the ``|`` operator to merge multiple |Table| instances
into a new table::

  >>> from astropy.table import QTable
  >>> t1 = Table({'name': ['foo', 'bar'], 'val': [0., 0.]}, meta={'n': 2})
  >>> t2 = QTable({'val': [1., 2.], 'val2': [10., 10.]}, meta={'id': 0})
  >>> t3 = t1 | t2  # Create a new table as result of update
  >>> t3
  <Table length=2>
  name   val     val2
  str3 float64 float64
  ---- ------- -------
   foo     1.0    10.0
   bar     2.0    10.0

``|`` and ``|=`` also take care of silently :ref:`merging_metadata`::

  >>> t3.meta
  {'n': 2, 'id': 0}

The columns in the updated |Table| are going to be copies of the originals. If
you need them to be references you can use the
:meth:`~astropy.table.Table.update` method with ``copy=False``, see :ref:`copy_versus_reference`
for details.

**Ensure the existence of a column**

|Table| has a :meth:`~astropy.table.Table.setdefault` method, which is
analogous to :meth:`dict.setdefault`.
It adds a column with a given name to the table if such a column is not in the
table already.
The default value passed to the method will be validated and, if necessary,
converted.
Either way the (possibly just inserted) column in the table is returned::

  >>> t0 = Table({"a": ["Ham", "Spam"]})
  >>> t0
  <Table length=2>
   a
  str4
  ----
   Ham
  Spam
  >>> t0.setdefault("a", ["Breakfast"])  # Existing column
  <Column name='a' dtype='str4' length=2>
   Ham
  Spam
  >>> t0.setdefault("approved", False)  # New column
  <Column name='approved' dtype='bool' length=2>
  False
  False
  >>> t0
  <Table length=2>
   a   approved
  str4   bool
  ---- --------
   Ham    False
  Spam    False

**Rename columns**

.. EXAMPLE START: Renaming Columns in Tables

To rename a column::

  >>> t.rename_column('a', 'a_new')
  >>> t['b'].name = 'b_new'

To rename multiple columns at once::

  >>> t.rename_columns(['a_new', 'b_new'], ['a', 'b'])

.. EXAMPLE END

**Add a row of data**

.. EXAMPLE START: Adding a Row of Data to a Table

To add a row::

  >>> t.add_row([-8, -9])

.. EXAMPLE END

**Remove rows**

.. EXAMPLE START: Removing Rows of Data from Tables

To remove a row::

  >>> t.remove_row(0)
  >>> t.remove_rows(slice(4, 5))
  >>> t.remove_rows([1, 2])

.. EXAMPLE END

**Sort by one or more columns**

.. EXAMPLE START: Sorting Columns in Tables

To sort columns::

  >>> t.sort('b')
  >>> t.sort(['a', 'b'])

.. EXAMPLE END

**Reverse table rows**

.. EXAMPLE START: Reversing Table Rows

To reverse the order of table rows::

  >>> t.reverse()

.. EXAMPLE END

**Modify metadata**

.. EXAMPLE START: Modifying Metadata in Tables

To modify metadata::

  >>> t.meta['key'] = 'value'

.. EXAMPLE END

**Select or reorder columns**

.. EXAMPLE START: Selecting or Reordering Columns in Tables

A new table with a subset or reordered list of columns can be
created as shown in the following example::

  >>> t = Table(arr, names=('a', 'b', 'c'))
  >>> t_acb = t['a', 'c', 'b']

Another way to do the same thing is to provide a list or tuple
as the item, as shown below::

  >>> new_order = ['a', 'c', 'b']  # List or tuple
  >>> t_acb = t[new_order]

.. EXAMPLE END

Caveats
=======

Modifying the table data and properties is fairly clear-cut, but one thing
to keep in mind is that adding a row *may* require a new copy in memory of the
table data. This depends on the detailed layout of Python objects in memory
and cannot be reliably controlled. In some cases it may be possible to build a
table row by row in less than O(N**2) time but you cannot count on it.

Another subtlety to keep in mind is that in some cases the return value of an
operation results in a new table in memory while in other cases it results in a
view of the existing table data. As an example, imagine trying to set two table
elements using column selection with ``t['a', 'c']`` in combination with row
index selection::

  >>> t = Table([[1, 2], [3, 4], [5, 6]], names=('a', 'b', 'c'))
  >>> t['a', 'c'][1] = (100, 100)
  >>> print(t)
   a   b   c
  --- --- ---
    1   3   5
    2   4   6

This might be surprising because the data values did not change and there
was no error. In fact, what happened is that ``t['a', 'c']`` created a
new temporary table in memory as a *copy* of the original and then updated the
first row of the copy. The original ``t`` table was unaffected and the new
temporary table disappeared once the statement was complete. The takeaway
is to pay attention to how certain operations are performed one step at
a time.

.. _table-replace-1_3:

In-Place Versus Replace Column Update
=====================================

Consider this code snippet::

  >>> t = Table([[1, 2, 3]], names=['a'])
  >>> t['a'] = [10.5, 20.5, 30.5]

There are a couple of ways this could be handled. It could update the existing
array values in-place (truncating to integer), or it could replace the entire
column with a new column based on the supplied data values.

The answer for ``astropy`` is that the operation shown above does a *complete
replacement* of the column object. In this case it makes a new column object
with float values by internally calling ``t.replace_column('a', [10.5, 20.5,
30.5])``. In general this behavior is more consistent with Python and `pandas
<https://pandas.pydata.org>`_ behavior.

**Forcing in-place update**

It is possible to force an in-place update of a column as follows::

  t[colname][:] = value

**Finding the source of problems**

In order to find potential problems related to replacing columns, there is the
option `astropy.table.conf.replace_warnings
<astropy.table.Conf.replace_warnings>` in the :ref:`astropy_config`. This
controls a set of warnings that are emitted under certain circumstances when a
table column is replaced. This option must be set to a list that includes zero
or more of the following string values:

``always`` :
  Print a warning every time a column gets replaced via the
  ``__setitem__()`` syntax (i.e., ``t['a'] = new_col``).

``slice`` :
  Print a warning when a column that appears to be a :class:`slice` of
  a parent column is replaced.

``refcount`` :
  Print a warning when the Python reference count for the
  column changes. This indicates that a stale object exists that might
  be used elsewhere in the code and give unexpected results.

``attributes`` :
  Print a warning if any of the standard column attributes changed.

The default value for the ``table.conf.replace_warnings`` option is
``[]`` (no warnings).
