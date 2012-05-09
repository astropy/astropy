.. _astropy_table:

.. include:: references.txt


Data Tables (`astropy.table`)
===============================================

Introduction
------------

`astropy.table` provides functionality for storing and manipulating
heterogenous tables of data in a way that is familiar to `numpy` users.  A few
notable features of this package are:

* Initialize a table from a wide variety of input data structures and types.
* Modify a table by adding or removing columns, changing column names,
  or adding new rows of data.
* Include table and column metadata as flexible data structures.
* Specify a description, units and output formatting for columns.
* Create a new table by selecting rows or columns from a table.
* Full support for multidimensional columns.
* Create a table by referencing (not copying) an existing `numpy` table.

Currently `astropy.table` is used when reading an ASCII table using
`astropy.io.ascii`.  Future releases of AstroPy are expected to use the |Table|
class for other subpackages such as `astropy.io.vo` and `astropy.io.fits`.

Getting Started
----------------

The basic workflow for creating a table, accessing table elements,
and modifying the table is shown below.  These examples show a very simple
case, while the full `astropy.table` documentation is available from the 
`Using Tables`_ section.

First create a simple table with three columns of data named ``a``, ``b``,
and ``c``.  These columns have integer, float, and string values respectively::

  >>> from astropy.table import Table, Column
  >>> a = [1, 4, 5]
  >>> b = [2.0, 5.0, 8.2]
  >>> c = ['x', 'y', 'z']
  >>> t = Table([a, b, c], names=('a', 'b', 'c'), meta={'name': 'first table'})
  >>> t
  <Table rows=3 names=('a','b','c')>
  array([(1, 2.0, 'x'), (4, 5.0, 'y'), (5, 8.2, 'z')], 
        dtype=[('a', '<i8'), ('b', '<f8'), ('c', '|S1')])

Now examine some high-level information about the table::

  >>> t.colnames
  ['a', 'b', 'c']
  >>> len(t)
  3
  >>> t.meta
  {'name': 'first table'}

Access the data by column or row using familiar `numpy` structured array syntax::
  
  >>> t['a']       # Column 'a'
  <Column name='a' units=None format=None description=None>
  array([1, 4, 5])

  >>> t['a'][1]    # Row 1 of column 'a'
  4

  >>> t[1]         # Row obj for with row 1 values
  <Row 1 of table
   values=(4, 5.0, 'y')
   dtype=[('a', '<i8'), ('b', '<f8'), ('c', '|S1')]>

  >>> t[1]['a']    # Column 'a' of row 1
  4

One can retreive a subset of a table by rows (using a slice) or
columns (using column names), where the subset is returned as a new table::

  >>> t[0:2]       # Table object with rows 0 and 1
  <Table rows=2 names=('a','b','c')>
  array([(1, 2.0, 'x'), (4, 5.0, 'y')], 
        dtype=[('a', '<i8'), ('b', '<f8'), ('c', '|S1')])

  >>> t['a', 'c']  # Table with cols 'a', 'c'
  <Table rows=3 names=('a','c')>
  array([(1, 'x'), (4, 'y'), (5, 'z')], 
        dtype=[('a', '<i8'), ('c', '|S1')])

Modifying table values in place is flexible and works as one would expect::

  >>> t['a'] = [-1, -2, -3]       # Set all column values
  >>> t['a'][2] = 30              # Set row 2 of column 'a'
  >>> t[1] = (8, 9.0, "W")        # Set all row values
  >>> t[1]['b'] = -9              # Set column 'b' of row 1
  >>> t[0:2]['b'] = 100.0         # Set column 'c' of rows 0 and 1
  >>> t
  <Table rows=3 names=('a','b','c')>
  array([(-1, 100.0, 'n'), (8, 100.0, 'W'), (30, 8.2, 'z')], 
        dtype=[('a', '<i8'), ('b', '<f8'), ('c', '|S1')])

Add, remove, and rename columns with the following::

  >>> t.add_column(Column('d', [1, 2, 3])))
  >>> t.remove_column('c')
  >>> t.rename_column('a', 'A')
  >>> t.colnames
  ['A', 'b', 'd']
  
Lastly, adding a new row of data to the table is as follows::

  >>> t.add_row([-8, -9, 10])
  >>> len(t)
  4

Using Tables
------------

The details of using `astropy.table` are provided in the following sections:

.. toctree::
   :maxdepth: 2

   construct_table.rst
   access_table.rst
   modify_table.rst
   masking.rst
   table_api.rst
