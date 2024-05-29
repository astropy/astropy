.. |add_index| replace:: :func:`~astropy.table.Table.add_index`
.. |index_mode| replace:: :func:`~astropy.table.Table.index_mode`

.. _table-indexing:

Table Indexing
**************

Once a |Table| has been created, it is possible to create indices on one or
more columns of the table. An index internally sorts the rows of a table based
on the index column(s), allowing for element retrieval by column value and
improved performance for certain table operations.

Creating an Index
=================

.. EXAMPLE START: Creating Indexes on Table Columns

To create an index on a table, use the |add_index| method::

   >>> from astropy.table import Table
   >>> t = Table([(2, 3, 2, 1), (8, 7, 6, 5)], names=('a', 'b'))
   >>> t.add_index('a')

The optional argument ``unique`` may be specified to create an index with
uniquely valued elements.

To create a composite index on multiple columns, pass a list of columns
instead::

   >>> t.add_index(['a', 'b'])

In particular, the first index created using the
|add_index| method is considered the default index or the "primary key." To
retrieve an index from a table, use the `~astropy.table.Table.indices`
property::

   >>> t.indices['a']
   <SlicedIndex original=True index=<Index columns=('a',) data=<SortedArray length=4>
    a  rows
   --- ----
     1    3
     2    0
     2    2
     3    1>>
   >>> t.indices['a', 'b']
   <SlicedIndex original=True index=<Index columns=('a', 'b') data=<SortedArray length=4>
    a   b  rows
   --- --- ----
     1   5    3
     2   6    2
     2   8    0
     3   7    1>>

.. EXAMPLE END

Row Retrieval using Indices
===========================

.. EXAMPLE START: Retrieving Table Rows using Indices

Row retrieval can be accomplished using two table properties:
`~astropy.table.Table.loc` and `~astropy.table.Table.iloc`. The
`~astropy.table.Table.loc` property can be indexed either by column value,
range of column values (*including* the bounds), or a :class:`list` or
|ndarray| of column values::

   >>> t = Table([(1, 2, 3, 4), (10, 1, 9, 9)], names=('a', 'b'), dtype=['i8', 'i8'])
   >>> t.add_index('a')
   >>> t.loc[2]  # the row(s) where a == 2
   <Row index=1>
     a     b
   int64 int64
   ----- -----
       2     1
   >>> t.loc[[1, 4]]  # the row(s) where a in [1, 4]
   <Table length=2>
     a     b
   int64 int64
   ----- -----
       1    10
       4     9
   >>> t.loc[1:3]  # the row(s) where a in [1, 2, 3]
   <Table length=3>
     a     b
   int64 int64
   ----- -----
       1    10
       2     1
       3     9
   >>> t.loc[:]
   <Table length=4>
     a     b
   int64 int64
   ----- -----
       1    10
       2     1
       3     9
       4     9

Note that by default, `~astropy.table.Table.loc` uses the primary index, which
here is column ``'a'``. To use a different index, pass the indexed column name
before the retrieval data::

   >>> t.add_index('b')
   >>> t.loc['b', 8:10]
   <Table length=3>
     a     b
   int64 int64
   ----- -----
       3     9
       4     9
       1    10

The property `~astropy.table.Table.iloc` works similarly, except that the
retrieval information must be either an integer or a :class:`slice`, and
relates to the sorted order of the index rather than column values. For
example::

   >>> t.iloc[0] # smallest row by value 'a'
   <Row index=0>
     a     b
   int64 int64
   ----- -----
       1    10
   >>> t.iloc['b', 1:] # all but smallest value of 'b'
   <Table length=3>
     a     b
   int64 int64
   ----- -----
       3     9
       4     9
       1    10

.. EXAMPLE END

Effects on Performance
======================

Table operations change somewhat when indices are present, and there are a
number of factors to consider when deciding whether the use of indices will
improve performance. In general, indexing offers the following advantages:

* Table grouping and sorting based on indexed column(s) both become faster.
* Retrieving values by index is faster than custom searching.

There are certain caveats, however:

* Creating an index requires time and memory.
* Table modifications become slower due to automatic index updates.
* Slicing a table becomes slower due to index relabeling.

See `here
<https://nbviewer.jupyter.org/github/mdmueller/astropy-notebooks/blob/master/table/indexing-profiling.ipynb>`_
for an IPython notebook profiling various aspects of table indexing.

Index Modes
===========

The |index_mode| method allows for some flexibility in the behavior of table
indexing by allowing the user to enter a specific indexing mode via a context
manager. There are currently three indexing modes: ``'freeze'``,
``'copy_on_getitem'``, and ``'discard_on_copy'``.

.. EXAMPLE START: Table Indexing with the "freeze" Index Mode

The ``'freeze'`` mode prevents automatic index updates whenever a column of the
index is modified, and all indices refresh themselves after the context ends::

  >>> with t.index_mode('freeze'):
  ...    t['a'][0] = 0
  ...    print(t.indices['a']) # unmodified
  <SlicedIndex original=True index=<Index columns=('a',) data=<SortedArray length=4>
   a  rows
  --- ----
    1    0
    2    1
    3    2
    4    3>>
  >>> print(t.indices['a']) # modified
  <SlicedIndex original=True index=<Index columns=('a',) data=<SortedArray length=4>
   a  rows
  --- ----
    0    0
    2    1
    3    2
    4    3>>

.. EXAMPLE END

.. EXAMPLE START: Table Indexing with the "copy_on_getitem" Index Mode

The ``'copy_on_getitem'`` mode forces columns to copy and relabel their indices
upon slicing. In the absence of this mode, table slices will preserve
indices while column slices will not::

  >>> ca = t['a'][[1, 3]]
  >>> ca.info.indices
  []
  >>> with t.index_mode('copy_on_getitem'):
  ...     ca = t['a'][[1, 3]]
  ...     print(ca.info.indices)
  [<SlicedIndex original=True index=<Index columns=('a',) data=<SortedArray length=2>
   a  rows
  --- ----
    2    0
    4    1>>]

.. EXAMPLE END

.. EXAMPLE START: Table Indexing with the "discard_on_copy" Index Mode

The ``'discard_on_copy'`` mode prevents indices from being copied whenever a
column or table is copied::

  >>> t2 = Table(t)
  >>> t2.indices['a']
  <SlicedIndex original=True index=<Index columns=('a',) data=<SortedArray length=4>
   a  rows
  --- ----
    0    0
    2    1
    3    2
    4    3>>
  >>> with t.index_mode('discard_on_copy'):
  ...    t2 = Table(t)
  ...    print(t2.indices)
  []

.. EXAMPLE END

Updating Rows using Indices
===========================

.. EXAMPLE START: Updating Table Rows using Indices

Row updates can be accomplished by assigning the table property
`~astropy.table.Table.loc` a complete row or a list of rows::

   >>> t = Table([('w', 'x', 'y', 'z'), (10, 1, 9, 9)], names=('a', 'b'), dtype=['str', 'i8'])
   >>> t.add_index('a')
   >>> t.loc['x']
   <Row index=1>
    a     b
   str1 int64
   ---- -----
      x     1
   >>> t.loc['x'] = ['a', 12]
   >>> t
   <Table length=4>
    a     b
   str1 int64
   ---- -----
      w    10
      a    12
      y     9
      z     9
   >>> t.loc[['w', 'y']]
   <Table length=2>
    a     b
   str1 int64
   ---- -----
      w    10
      y     9
   >>> t.loc[['w', 'z']] = [['b', 23], ['c', 56]]
   >>> t
   <Table length=4>
    a     b
   str1 int64
   ---- -----
      b    23
      a    12
      y     9
      c    56

.. EXAMPLE END

Retrieving the Location of Rows using Indices
=============================================

.. EXAMPLE START: Retrieving the Location of Table Rows using Indices

Retrieval of the location of rows can be accomplished using a table property:
`~astropy.table.Table.loc_indices`. The `~astropy.table.Table.loc_indices`
property can be indexed either by column value, range of column values
(*including* the bounds), or a :class:`list` or |ndarray| of column values::

   >>> t = Table([('w', 'x', 'y', 'z'), (10, 1, 9, 9)], names=('a', 'b'), dtype=['str', 'i8'])
   >>> t.add_index('a')
   >>> t.loc_indices['x']
   np.int64(1)

.. EXAMPLE END

Engines
=======

When creating an index via |add_index|, the keyword argument ``engine`` may be
specified to use a particular indexing engine. The available engines are:

* `~astropy.table.SortedArray`, a sorted array engine using an underlying
  sorted |Table|.
* `~astropy.table.SCEngine`, a sorted list engine using the `Sorted Containers
  <https://pypi.org/project/sortedcontainers/>`_ package.
* `~astropy.table.BST`, a Python-based binary search tree engine (not recommended).

The SCEngine depends on the ``sortedcontainers`` dependency. The most important takeaway is that
`~astropy.table.SortedArray` (the default engine) is usually best, although
`~astropy.table.SCEngine` may be more appropriate for an index created on an
empty column since adding new values is quicker.

The `~astropy.table.BST` engine demonstrates a simple pure Python implementation
of a search tree engine, but the performance is poor for larger tables. This
is available in the code largely as an implementation reference.
