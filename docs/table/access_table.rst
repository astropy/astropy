.. _access_table:

.. include:: references.txt

=============================
Accessing a table
=============================

Accessing the table properties and data is straightforward and is generally consistent with
the basic interface for `numpy` structured arrays.

Quick overview
==============

For the impatient, the code below shows the basics of accessing table data.
Where relevant there is a comment about what sort of object.  Except where
noted, the table access returns objects that can be modified in order to
update table data or properties.
In cases where is returned and how
the data contained in that object relate to the original table data
(i.e. whether it is a copy or reference, see :ref:`copy_versus_reference`).
  
**Make table**
::
  
  from astropy.table import Table
  import numpy as np
  
  arr = np.arange(15).reshape(5, 3)
  t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})
  
**Table properties**
::
  
  t.columns   # Dict of table columns (read-only)
  t.colnames  # List of column names (read-only) *DEV NOTE: not strictly true right now)*
  t.meta      # Dict of meta-data
  len(t)      # Number of table rows
  
**Access table data**
::
  
  t['a']       # Column 'a'
  t['a'][1]    # Row 1 of column 'a'
  t[1]         # Row obj for with row 1 values
  t[1]['a']    # Column 'a' of row 1
  t[2:5]       # Table object with rows 2:5
  t[[1, 3, 4]]  # Table object with rows 1, 3, 4 (copy)
  t[np.array([1, 3, 4])]  # Table object with rows 1, 3, 4 (copy)
  t['a', 'c']  # Table with cols 'a', 'c' (copy)
  dat = np.array(t)  # Copy table data to numpy structured array object

Details
================

For all the following examples it is assumed that the table has been created as below::

  >>> from astropy.table import Table, Column
  >>> import numpy as np
  
  >>> arr = np.arange(15).reshape(5, 3)
  >>> t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})
  >>> t['a'].format = "%6.3f"  # print as a float with 3 digits after decimal point
  >>> t['a'].units = 'm sec^-1'
  >>> t['a'].description = 'unladen swallow velocity'

**Accessing table properties**

The code below shows accessing the table columns as a |TableColumns| object,
getting the column names, table meta-data, and number of table rows.  The table
meta-data is simply an ordered dictionary (OrderedDict_) by default.
::

  >>> t.columns
  <TableColumns names=('a','b','c')>
  
  >>> t.colnames
  ['a', 'b', 'c']
  
  >>> t.meta  # Dict of meta-data
  {'keywords': {'key1': 'val1'}}
  
  >>> len(t)
  5


**Accessing table data**

As expected one can access a table column by name and get an element from that
column with a numerical index::

  >>> t['a']  # Column 'a'
  <Column name='a' units='m sec^-1' format='%6.3f' description='unladen swallow velocity'>
  array([ 0,  3,  6,  9, 12])

  >>> t['a'][1]  # Row 1 of column 'a'
  3
  
When a table column is printed, either with ``print`` or via the ``str()``
built-in function, it is formatted according to the ``format`` attribute (see
:ref:`table_format_string`)::

  >>> print t['a'].description, t['a']
  unladen swallow velocity  0.000,  3.000,  6.000,  9.000, 12.000

Likewise a table row and a column from that row can be selected::

  >>> t[1]  # Row object corresponding to row 1
  <Row 1 of table
   values=(3, 4, 5)
   dtype=[('a', '<i8'), ('b', '<i8'), ('c', '<i8')]>
  
  >>> t[1]['a']  # Column 'a' of row 1
  3
  
A |Row| object has the same columns and meta-data as its parent table::

  >>> t[1].columns
  <TableColumns names=('a','b','c')>
  
  >>> t[1].colnames
  ['a', 'b', 'c']

Slicing a table returns a new table object which references to the original
data within the slice region (See :ref:`copy_versus_reference`).  The table
meta-data and column definitions are copied.
::

  >>> t[2:5]  # Table object with rows 2:5 (reference)
  <Table rows=3 names=('a','b','c')>
  array([(6, 7, 8), (9, 10, 11), (12, 13, 14)], 
        dtype=[('a', '<i8'), ('b', '<i8'), ('c', '<i8')])
  
It is possible to select table rows with an array of indexes or by providing
specifying multiple column names.  This returns a copy of the original table
for the selected rows.  ::

  >>> t[[1, 3, 4]]  # Table object with rows 1, 3, 4 (copy)
  <Table rows=3 names=('a','b','c')>
  array([(3, 4, 5), (9, 10, 11), (12, 13, 14)], 
        dtype=[('a', '<i8'), ('b', '<i8'), ('c', '<i8')])
  
  >>> t[np.array([1, 3, 4])]  # Table object with rows 1, 3, 4 (copy)
  <Table rows=3 names=('a','b','c')>
  array([(3, 4, 5), (9, 10, 11), (12, 13, 14)], 
        dtype=[('a', '<i8'), ('b', '<i8'), ('c', '<i8')])
  
  >>> t['a', 'c']  # Table with cols 'a', 'c' (copy)
  <Table rows=5 names=('a','c')>
  array([(0, 2), (3, 5), (6, 8), (9, 11), (12, 14)], 
        dtype=[('a', '<i8'), ('c', '<i8')])

Finally, one can access the underlying table data as a native `numpy`
structured array by creating a copy or reference with ``np.array``::

  >>> data = np.array(t)  # copy of data in t as a structured array
  >>> data = np.array(t, copy=False)  # reference to data in t
