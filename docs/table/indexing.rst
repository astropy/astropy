.. include:: references.txt
.. |add_index| replace:: :func:`~astropy.table.Table.add_index`
.. |index_mode| replace:: :func:`~astropy.table.Table.index_mode`

Table indexing
--------------

Once a |Table| has been created, it is possible to create indexes on one or
more columns of the table. An index internally sorts the rows of a table based
on the index column(s), allowing for element retrieval by column value and
improved performance for certain table operations.

.. Warning::

   The table indexing engine is new and is not yet considered stable.
   It is recommended to avoid using this engine in production code for now.

Creating an index
^^^^^^^^^^^^^^^^^

To create an index on a table, use the |add_index| method::

   >>> from astropy.table import Table
   >>> t = Table([(2, 3, 2, 1), (8, 7, 6, 5)], names=('a', 'b'))
   >>> t.add_index('a')

The optional argument "unique" may be specified to create an index with
uniquely valued elements.

To create a composite index on multiple columns, pass a list of columns
instead::

   >>> t.add_index(['a', 'b'])

In particular, the first index created using the
|add_index| method is considered the default index or the "primary key". To
retrieve an index from a table, use the `~astropy.table.Table.indices` property::

   >>> t.indices['a']
    a  rows
   --- ----
     1    3
     2    0
     2    2
     3    1
   >>> t.indices['a', 'b']
    a   b  rows
   --- --- ----
     1   5    3
     2   6    2
     2   8    0
     3   7    1



Row retrieval using indices
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Row retrieval can be accomplished using two table properties: `~astropy.table.Table.loc` and
`~astropy.table.Table.iloc`. The `~astropy.table.Table.loc` property can be indexed either by column value, range of
column values (*including* the bounds), or a list or ndarray of column values::

   >>> t = Table([(1, 2, 3, 4), (10, 1, 9, 9)], names=('a', 'b'))
   >>> t.add_index('a')
   >>> t.loc[2]
   <Row index=1>
     a     b
   int64 int64
   ----- -----
       2     1
   >>> t.loc[[1, 4]]
   <Table length=2>
     a     b  
   int64 int64
   ----- -----
       1    10
       4     9
   >>> t.loc[1:3]
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


Note that by default, `~astropy.table.Table.loc` uses the primary index, which here is column
'a'. To use a different index, pass the indexed column name before the
retrieval data::

   >>> t.add_index('b')
   >>> t.loc['b', 8:10]
   <Table length=3>
     a     b  
   int64 int64
   ----- -----
       3     9
       4     9
       1    10

The property `~astropy.table.Table.iloc` works similarly, except that the retrieval information must
be either an int or a slice, and relates to the sorted order of the index
rather than column values. For example::

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

Effects on performance
^^^^^^^^^^^^^^^^^^^^^^
Table operations change somewhat when indices are present, and there are a
number of factors to consider when deciding whether the use of indices will
improve performance. In general, indexing offers the following advantages:

* Table grouping and sorting based on indexed column(s) become faster
* Retrieving values by index is faster than custom searching

There are certain caveats, however:

* Creating an index requires time and memory
* Table modifications become slower due to automatic index updates
* Slicing a table becomes slower due to index relabeling

See `here <http://nbviewer.ipython.org/github/mdmueller/astropy-notebooks/blob/master/table/indexing-profiling.ipynb>`_ for an IPython notebook profiling various aspects of table indexing.

Index modes
^^^^^^^^^^^
The |index_mode| method allows for some flexibility in the behavior of table
indexing by allowing the user to enter a specific indexing mode via a context manager. There are
currently three indexing modes: *freeze*, *copy_on_getitem*, and
*discard_on_copy*. The *freeze* mode prevents automatic index updates whenever
a column of the index is modified, and all indices refresh themselves after the
context ends::

  >>> with t.index_mode('freeze'):
  ...    t['a'][0] = 0
  ...    print(t.indices['a']) # unmodified
   a  rows
  --- ----
    1    0
    2    1
    3    2
    4    3
  >>> print(t.indices['a']) # modified
   a  rows
  --- ----
    0    0
    2    1
    3    2
    4    3

The *copy_on_getitem* mode forces columns to copy and relabel their indices upon
slicing. In the absence of this mode, table slices will preserve
indices while column slices will not::

  >>> t['a'][[1, 3]].info.indices
  []
  >>> with t.index_mode('copy_on_getitem'):
  ...    print(t['a'][[1, 3]].info.indices)
  [ a  rows
  --- ----
    2    0
    4    1]

The *discard_on_copy* mode prevents indices from being copied whenever a column
or table is copied::

  >>> t2 = Table(t)
  >>> t2.indices['a']
   a  rows
  --- ----
    0    0
    2    1
    3    2
    4    3
  >>> t2.indices['b']
   b  rows
  --- ----
    1    1
    9    2
    9    3
   10    0
  >>> with t.index_mode('discard_on_copy'):
  ...    t2 = Table(t)
  ...    print(t2.indices)
  []

Engines
^^^^^^^
When creating an index via |add_index|, the keyword argument "engine" may be
specified to use a particular indexing engine. The available engines are

* `~astropy.table.SortedArray`, a sorted array engine using an underlying
  sorted Table
* `~astropy.table.FastRBT`, a C-based red-black tree engine
* `~astropy.table.FastBST`, a C-based binary search tree engine
* `~astropy.table.BST`, a Python-based binary search tree engine

Note that FastRBT and FastBST depend on the bintrees dependency; without this
dependency, both classes default to `~astropy.table.BST`. For a comparison of
engine performance, see `this IPython notebook
<http://nbviewer.ipython.org/github/mdmueller/astropy-notebooks/blob/master/table/indexing-profiling.ipynb>`_. Probably
the most important takeaway is that `~astropy.table.SortedArray` (the default
engine) is usually best, although `~astropy.table.FastRBT` may be more
appropriate for an index created on an empty column since adding new values is quicker.
