.. _construct_table:

Constructing a Table
********************

There is great deal of flexibility in the way that a table can be initially
constructed. Details on the inputs to the |Table| and |QTable|
constructors are in the `Initialization Details`_ section. However, the
best way to understand how to make a table is by example.

Examples
========

Setup
-----

For the following examples you need to import the |QTable|, |Table|, and
|Column| classes along with the :ref:`astropy-units` package and the ``numpy``
package::

  >>> from astropy.table import QTable, Table, Column
  >>> from astropy import units as u
  >>> import numpy as np

Creating from Scratch
---------------------

.. EXAMPLE START: Creating an Astropy Table from Scratch

A |Table| can be created without any initial input data or even without any
initial columns. This is useful for building tables dynamically if the initial
size, columns, or data are not known.

.. Note::
   Adding rows requires making a new copy of the entire
   table each time, so in the case of large tables this may be slow.
   On the other hand, adding columns is fast.

::

  >>> t = Table()
  >>> t['a'] = [1, 4]
  >>> t['b'] = [2.0, 5.0]
  >>> t['c'] = ['x', 'y']

  >>> t = Table(names=('a', 'b', 'c'), dtype=('f4', 'i4', 'S2'))
  >>> t.add_row((1, 2.0, 'x'))
  >>> t.add_row((4, 5.0, 'y'))

  >>> t = Table(dtype=[('a', 'f4'), ('b', 'i4'), ('c', 'S2')])

If your data columns have physical units associated with them then we
recommend using the |QTable| class. This will allow the column to be
stored in the table as a native |Quantity| and bring the full power of
:ref:`astropy-units` to the table. See :ref:`quantity_and_qtable` for details.
::

  >>> t = QTable()
  >>> t['a'] = [1, 4]
  >>> t['b'] = [2.0, 5.0] * u.cm / u.s
  >>> t['c'] = ['x', 'y']
  >>> type(t['b'])
  <class 'astropy.units.quantity.Quantity'>

.. EXAMPLE END

List of Columns
---------------

.. EXAMPLE START: Creating an Astropy Table from a List of Columns

A typical case is where you have a number of data columns with the same length
defined in different variables. These might be Python lists or ``numpy`` arrays
or a mix of the two. These can be used to create a |Table| by putting the column
data variables into a Python list. In this case the column names are not
defined by the input data, so they must either be set using the ``names``
keyword or they will be automatically generated as ``col<N>``.

::

  >>> a = np.array([1, 4], dtype=np.int32)
  >>> b = [2.0, 5.0]
  >>> c = ['x', 'y']
  >>> t = Table([a, b, c], names=('a', 'b', 'c'))
  >>> t
  <Table length=2>
    a      b     c
  int32 float64 str1
  ----- ------- ----
      1     2.0    x
      4     5.0    y

.. EXAMPLE END

**Make a new table using columns from the first table**

Once you have a |Table|, then you can make a new table by selecting columns
and putting them into a Python list (e.g., ``[ t['c'], t['a'] ]``)::

  >>> Table([t['c'], t['a']])
  <Table length=2>
   c     a
  str1 int32
  ---- -----
     x     1
     y     4

**Make a new table using expressions involving columns**

The |Column| object is derived from the standard |ndarray| and can be used
directly in arithmetic expressions. This allows for a compact way of making a
new table with modified column values::

  >>> Table([t['a']**2, t['b'] + 10])
  <Table length=2>
    a      b
  int32 float64
  ----- -------
      1    12.0
     16    15.0


**Different types of column data**

The list input method for |Table| is very flexible since you can use a mix
of different data types to initialize a table::

  >>> a = (1., 4.)
  >>> b = np.array([[2, 3], [5, 6]], dtype=np.int64)  # vector column
  >>> c = Column(['x', 'y'], name='axis')
  >>> d = u.Quantity([([1., 2., 3.], [.1, .2, .3]),
  ...                 ([4., 5., 6.], [.4, .5, .6])], 'm,m/s')
  >>> QTable([a, b, c, d])
  <QTable length=2>
    col0    col1   axis          col3 [f0, f1]
                                    (m, m / s)
  float64 int64[2] str1     (float64[3], float64[3])
  ------- -------- ---- -------------------------------
      1.0   2 .. 3    x ([1., 2., 3.], [0.1, 0.2, 0.3])
      4.0   5 .. 6    y ([4., 5., 6.], [0.4, 0.5, 0.6])

Notice that in the third column the existing column name ``'axis'`` is used.

Dict of Columns
---------------

.. EXAMPLE START: Creating an Astropy Table from a Dictionary of Columns

A :class:`dict` of column data can be used to initialize a |Table|::

  >>> arr = {'a': np.array([1, 4], dtype=np.int32),
  ...        'b': [2.0, 5.0],
  ...        'c': ['x', 'y']}
  >>>
  >>> Table(arr)
  <Table length=2>
    a      b     c
  int32 float64 str1
  ----- ------- ----
      1     2.0    x
      4     5.0    y

.. EXAMPLE END

**Specify the column order and optionally the data types**
::

  >>> Table(arr, names=('a', 'c', 'b'), dtype=('f8', 'U2', 'i4'))
  <Table length=2>
     a     c     b
  float64 str2 int32
  ------- ---- -----
      1.0    x     2
      4.0    y     5

**Different types of column data**

The input column data can be any data type that can initialize a |Column|
object::

  >>> arr = {'a': (1., 4.),
  ...        'b': np.array([[2, 3], [5, 6]], dtype=np.int64),
  ...        'c': Column(['x', 'y'], name='axis')}
  >>> Table(arr, names=('a', 'b', 'c'))
  <Table length=2>
     a       b      c
  float64 int64[2] str1
  ------- -------- ----
      1.0   2 .. 3    x
      4.0   5 .. 6    y

Notice that the key ``'c'`` takes precedence over the existing column name
``'axis'`` in the third column. Also see that the ``'b'`` column is a vector
column where each row element is itself a two-element array.

**Renaming columns is not possible**
::

  >>> Table(arr, names=('a_new', 'b_new', 'c_new'))
  Traceback (most recent call last):
    ...
  KeyError: 'a_new'

Row Data
--------

Row-oriented data can be used to create a table using the ``rows``
keyword argument.

**List or tuple of data records**

If you have row-oriented input data such as a list of records, you
need to use the ``rows`` keyword to create a table::

  >>> data_rows = [(1, 2.0, 'x'),
  ...              (4, 5.0, 'y'),
  ...              (5, 8.2, 'z')]
  >>> t = Table(rows=data_rows, names=('a', 'b', 'c'))
  >>> print(t)
   a   b   c
  --- --- ---
    1 2.0   x
    4 5.0   y
    5 8.2   z

**List of dict objects**

You can also initialize a table with row values. This is constructed as a
list of :class:`dict` objects. The keys determine the column names::

  >>> data = [{'a': 5, 'b': 10},
  ...         {'a': 15, 'b': 20}]
  >>> t = Table(rows=data)
  >>> print(t)
   a   b
  --- ---
    5  10
   15  20

If there are missing keys in one or more rows then the corresponding values
will be marked as missing (masked)::

  >>> t = Table(rows=[{'a': 5, 'b': 10}, {'a': 15, 'c': 50}])
  >>> print(t)
   a   b   c
  --- --- ---
    5  10  --
   15  --  50

You can specify the column order with the ``names`` argument::

  >>> data = [{'a': 5, 'b': 10},
  ...         {'a': 15, 'b': 20}]
  >>> t = Table(rows=data, names=('b', 'a'))
  >>> print(t)
   b   a
  --- ---
   10   5
   20  15

If ``names`` are not provided then column ordering will be determined by the
first :class:`dict` if it contains values for all the columns, or by sorting
the column names alphabetically if it doesn't::

  >>> data = [{'b': 10, 'c': 7, 'a': 5},
  ...         {'a': 15, 'c': 35, 'b': 20}]
  >>> t = Table(rows=data)
  >>> print(t)
   b   c   a
  --- --- ---
   10   7   5
   20  35  15
  >>> data = [{'b': 10, 'c': 7, },
  ...         {'a': 15, 'c': 35, 'b': 20}]
  >>> t = Table(rows=data)
  >>> print(t)
   a   b   c
  --- --- ---
   --  10   7
   15  20  35

**Single row**

You can also make a new table from a single row of an existing table::

  >>> a = [1, 4]
  >>> b = [2.0, 5.0]
  >>> t = Table([a, b], names=('a', 'b'))
  >>> t2 = Table(rows=t[1])

Remember that a |Row| has effectively a zero length compared to the
newly created |Table| which has a length of one. This is similar to
the difference between a scalar ``1`` (length 0) and an array such as
``np.array([1])`` with length 1.

.. Note::

   In the case of input data as a list of dicts or a single |Table| row, you
   can supply the data as the ``data`` argument since these forms
   are always unambiguous. For example, ``Table([{'a': 1}, {'a': 2}])`` is
   accepted. However, a list of records must always be provided using the
   ``rows`` keyword, otherwise it will be interpreted as a list of columns.

NumPy Structured Array
----------------------

The `structured array <https://numpy.org/doc/stable/user/basics.rec.html>`_ is
the standard mechanism in ``numpy`` for storing heterogeneous table data. Most
scientific I/O packages that read table files (e.g., `astropy.io.fits`,
`astropy.io.votable`, and `asciitable
<https://cxc.harvard.edu/contrib/asciitable/>`_) will return the table in an
object that is based on the structured array. A structured array can be
created using::

  >>> arr = np.array([(1, 2.0, 'x'),
  ...                 (4, 5.0, 'y')],
  ...                dtype=[('a', 'i4'), ('b', 'f8'), ('c', 'U2')])

From ``arr`` it is possible to create the corresponding |Table| object::

  >>> Table(arr)
  <Table length=2>
    a      b     c
  int32 float64 str2
  ----- ------- ----
      1     2.0    x
      4     5.0    y

Note that in the above example and most of the following examples we are
creating a table and immediately asking the interactive Python interpreter to
print the table to see what we made. In real code you might do something like::

  >>> table = Table(arr)
  >>> print(table)
   a   b   c
  --- --- ---
    1 2.0   x
    4 5.0   y

**New column names**

The column names can be changed from the original values by providing the
``names`` argument::

  >>> Table(arr, names=('a_new', 'b_new', 'c_new'))
  <Table length=2>
  a_new  b_new  c_new
  int32 float64  str2
  ----- ------- -----
      1     2.0     x
      4     5.0     y

**New data types**

The data type for each column can likewise be changed with ``dtype``::

  >>> Table(arr, dtype=('f4', 'i4', 'U4'))
  <Table length=2>
     a      b    c
  float32 int32 str4
  ------- ----- ----
      1.0     2    x
      4.0     5    y

  >>> Table(arr, names=('a_new', 'b_new', 'c_new'), dtype=('f4', 'i4', 'U4'))
  <Table length=2>
   a_new  b_new c_new
  float32 int32  str4
  ------- ----- -----
      1.0     2     x
      4.0     5     y

NumPy Homogeneous Array
-----------------------

A ``numpy`` 1D array is treated as a single row table where each element of the
array corresponds to a column::

  >>> Table(np.array([1, 2, 3]), names=['a', 'b', 'c'], dtype=('i8', 'i8', 'i8'))
  <Table length=1>
    a     b     c
  int64 int64 int64
  ----- ----- -----
      1     2     3

A ``numpy`` 2D array (where all elements have the same type) can also be
converted into a |Table|. In this case the column names are not specified by
the data and must either be provided by the user or will be automatically
generated as ``col<N>`` where ``<N>`` is the column number.

**Basic example with automatic column names**
::

  >>> arr = np.array([[1, 2, 3],
  ...                 [4, 5, 6]], dtype=np.int32)
  >>> Table(arr)
  <Table length=2>
   col0  col1  col2
  int32 int32 int32
  ----- ----- -----
      1     2     3
      4     5     6

**Column names and types specified**
::

  >>> Table(arr, names=('a_new', 'b_new', 'c_new'), dtype=('f4', 'i4', 'U4'))
  <Table length=2>
   a_new  b_new c_new
  float32 int32  str4
  ------- ----- -----
      1.0     2     3
      4.0     5     6

**Referencing the original data**

It is possible to reference the original data as long as the data types are not
changed::

  >>> t = Table(arr, copy=False)

See the `Copy versus Reference`_ section for more information.

**Python arrays versus NumPy arrays as input**

There is a slightly subtle issue that is important to understand about the way
that |Table| objects are created. Any data input that looks like a Python
:class:`list` (including a :class:`tuple`) is considered to be a list of
columns. In contrast, a homogeneous |ndarray| input is interpreted as a list of
rows::

  >>> arr = [[1, 2, 3],
  ...        [4, 5, 6]]
  >>> np_arr = np.array(arr)

  >>> print(Table(arr))    # Two columns, three rows
  col0 col1
  ---- ----
     1    4
     2    5
     3    6

  >>> print(Table(np_arr))  # Three columns, two rows
  col0 col1 col2
  ---- ---- ----
     1    2    3
     4    5    6

This dichotomy is needed to support flexible list input while retaining the
natural interpretation of 2D ``numpy`` arrays where the first index corresponds
to data "rows" and the second index corresponds to data "columns."

From an Existing Table
----------------------

.. EXAMPLE START: Creating an Astropy Table from an Existing Table

A new table can be created by selecting a subset of columns in an existing
table::

  >>> t = Table(names=('a', 'b', 'c'))
  >>> t['c', 'b', 'a']  # Makes a copy of the data
  <Table length=0>
     c       b       a
  float64 float64 float64
  ------- ------- -------

An alternate way is to use the ``columns`` attribute (explained in the
`TableColumns`_ section) to initialize a new table. This lets you choose
columns by their numerical index or name and supports slicing syntax::

  >>> Table(t.columns[0:2])
  <Table length=0>
     a       b
  float64 float64
  ------- -------

  >>> Table([t.columns[0], t.columns['c']])
  <Table length=0>
     a       c
  float64 float64
  ------- -------

To create a copy of an existing table that is empty (has no rows)::

 >>> t = Table([[1.0, 2.3], [2.1, 3]], names=['x', 'y'])
 >>> t
 <Table length=2>
    x       y
 float64 float64
 ------- -------
     1.0     2.1
     2.3     3.0

 >>> tcopy = t[:0].copy()
 >>> tcopy
 <Table length=0>
    x       y
 float64 float64
 ------- -------

.. EXAMPLE END

Empty Array of a Known Size
---------------------------

.. EXAMPLE START: Creating an Astropy Table from an Empty Array

If you do know the size that your table will be, but do not know the values in
advance, you can create a zeroed |ndarray| and build the |Table| from it::

  >>> N = 3
  >>> dtype = [('a', 'i4'), ('b', 'f8'), ('c', 'bool')]
  >>> t = Table(data=np.zeros(N, dtype=dtype))
  >>> t
  <Table length=3>
    a      b      c
  int32 float64  bool
  ----- ------- -----
      0     0.0 False
      0     0.0 False
      0     0.0 False

For example, you can then fill in this table row by row with values extracted
from another table, or generated on the fly::

  >>> for i in range(len(t)):
  ...     t[i] = (i, 2.5*i, i % 2)
  >>> t
  <Table length=3>
    a      b      c
  int32 float64  bool
  ----- ------- -----
      0     0.0 False
      1     2.5  True
      2     5.0 False

.. EXAMPLE END

SkyCoord
--------

A |SkyCoord| object can be converted to a |QTable| using its
:meth:`~astropy.coordinates.SkyCoord.to_table` method. For details and examples
see :ref:`skycoord-table-conversion`.

Pandas DataFrame
----------------

The section on :ref:`pandas` gives details on how to initialize a |Table| using
a :class:`pandas.DataFrame` via the :func:`~astropy.table.Table.from_pandas`
class method. This provides a convenient way to take advantage of the many I/O
and table manipulation methods in `pandas <https://pandas.pydata.org/>`_.

Comment Lines
-------------

.. EXAMPLE START: Adding Comment Lines in an ASCII File

Comment lines in an ASCII file can be added via the ``'comments'`` key in the
table's metadata. The following will insert two comment lines in the output
ASCII file unless ``comment=False`` is explicitly set in ``write()``::

  >>> import sys
  >>> from astropy.table import Table
  >>> t = Table(names=('a', 'b', 'c'), dtype=('f4', 'i4', 'S2'))
  >>> t.add_row((1, 2.0, 'x'))
  >>> t.meta['comments'] = ['Here is my explanatory text. This is awesome.',
  ...                       'Second comment line.']
  >>> t.write(sys.stdout, format='ascii')
  # Here is my explanatory text. This is awesome.
  # Second comment line.
  a b c
  1.0 2 x

.. EXAMPLE END

Initialization Details
======================

A table object is created by initializing a |Table| class
object with the following arguments, all of which are optional:

``data`` : |ndarray|, :class:`dict`, :class:`list`, |Table|, or table-like object, optional
    Data to initialize table.
``masked`` : :class:`bool`, optional
    Specify whether the table is masked.
``names`` : :class:`list`, optional
    Specify column names.
``dtype`` : :class:`list`, optional
    Specify column data types.
``meta`` : :class:`dict`, optional
    Metadata associated with the table.
``copy`` : :class:`bool`, optional
    Copy the input data. If the input is a |Table| the ``meta`` is always
    copied regardless of the ``copy`` parameter.
    Default is `True`.
``rows`` : |ndarray|, :class:`list` of lists, optional
    Row-oriented data for table instead of ``data`` argument.
``copy_indices`` : :class:`bool`, optional
    Copy any indices in the input data. Default is `True`.
``units`` : :class:`list`, :class:`dict`, optional
    List or dict of units to apply to columns.
``descriptions`` : :class:`list`, :class:`dict`, optional
    List or dict of descriptions to apply to columns.
``**kwargs`` : :class:`dict`, optional
    Additional keyword args when converting table-like object.

The following subsections provide further detail on the values and options for
each of the keyword arguments that can be used to create a new |Table| object.

data
----

The |Table| object can be initialized with several different forms
for the ``data`` argument.

**NumPy ndarray (structured array)**
    The base column names are the field names of the ``data`` structured
    array. The ``names`` list (optional) can be used to select
    particular fields and/or reorder the base names. The ``dtype`` list
    (optional) must match the length of ``names`` and is used to
    override the existing ``data`` types.

**NumPy ndarray (homogeneous)**
    If the ``data`` is a one-dimensional |ndarray| then it is treated as a
    single row table where each element of the array corresponds to a column.

    If the ``data`` is an at least two-dimensional |ndarray|, then the first
    (left-most) index corresponds to row number (table length) and the
    second index corresponds to column number (table width). Higher
    dimensions get absorbed in the shape of each table cell.

    If provided, the ``names`` list must match the "width" of the ``data``
    argument. The default for ``names`` is to auto-generate column names
    in the form ``col<N>``. If provided, the ``dtype`` list overrides the
    base column types and must match the length of ``names``.

**dict-like**
    The keys of the ``data`` object define the base column names. The
    corresponding values can be |Column| objects, ``numpy`` arrays, or list-
    like objects. The ``names`` list (optional) can be used to select
    particular fields and/or reorder the base names. The ``dtype`` list
    (optional) must match the length of ``names`` and is used to override
    the existing or default data types.

**list-like**
    Each item in the ``data`` list provides a column of data values and
    can be a |Column| object, |ndarray|, or list-like object. The
    ``names`` list defines the name of each column. The names will be
    auto-generated if not provided (either with the ``names`` argument or
    by |Column| objects). If provided, the ``names`` argument must match the
    number of items in the ``data`` list. The optional ``dtype`` list
    will override the existing or default data types and must match
    ``names`` in length.

**list-of-dicts**
    Similar to Python's built-in :class:`csv.DictReader`, each item in the
    ``data`` list provides a row of data values and must be a :class:`dict`.
    The key values in each :class:`dict` define the column names. The ``names``
    argument may be supplied to specify column ordering. If ``names`` are not
    provided then column ordering will be determined by the first :class:`dict`
    if it contains values for all the columns, or by sorting the column names
    alphabetically if it does not. The ``dtype`` list may be specified, and
    must correspond to the order of output columns.

**Table-like object**
    If another table-like object has a ``__astropy_table__()`` method then
    that object can be used to directly create a |Table|. See the
    `table-like objects`_ section for details.

**None**
    Initialize a zero-length table. If ``names`` and optionally ``dtype``
    are provided, then the corresponding columns are created.

names
-----

The ``names`` argument provides a way to specify the table column names or
override the existing ones. By default, the column names are either taken from
existing names (for |ndarray| or |Table| input) or auto-generated as
``col<N>``. If ``names`` is provided, then it must be a list with the same
length as the number of columns. Any list elements with value `None` fall back
to the default name.

In the case where ``data`` is provided as a :class:`dict` of columns, the
``names`` argument can be supplied to specify the order of columns. The
``names`` list must then contain each of the keys in the ``data``
:class:`dict`.

dtype
-----

The ``dtype`` argument provides a way to specify the table column data types or
override the existing types. By default, the types are either taken from
existing types (for |ndarray| or |Table| input) or auto-generated by the
:func:`numpy.array` routine. If ``dtype`` is provided then it must be a list
with the same length as the number of columns. The values must be valid
:class:`numpy.dtype` initializers or `None`. Any list elements with value
`None` fall back to the default type.

meta
----

The ``meta`` argument is an object that contains metadata associated with the
table. It is recommended that this object be a :class:`dict` or
:class:`~collections.OrderedDict`, but the only firm requirement is that it can
be copied with the standard library :func:`copy.deepcopy` routine. By
default, ``meta`` is an empty :class:`~collections.OrderedDict`.

copy
----

In the case where ``data`` is either an |ndarray| object, a :class:`dict`, or
an existing |Table|, it is possible to use a reference to the existing data by
setting ``copy=False``. This has the advantage of reducing memory use and being
faster. However, you should take care because any modifications to the new
|Table| data will also be seen in the original input data. See the `Copy versus
Reference`_ section for more information.

rows
----

This argument allows for providing data as a sequence of rows, in contrast
to the ``data`` keyword, which generally assumes data are a sequence of columns.
The `Row data`_ section provides details.

copy_indices
------------

If you are initializing a |Table| from another |Table| that makes use of
:ref:`table-indexing`, then this option allows copying that table *without*
copying the indices by setting ``copy_indices=False``. By default, the indices
are copied.

units
-----

This allows for setting the unit for one or more columns at the time of
creating the table. The input can be either a list of unit values corresponding
to each of the columns in the table (using `None` or ``''`` for no unit), or a
:class:`dict` that provides the unit for specified column names. For example::

  >>> dat = [[1, 2], ['hello', 'world']]
  >>> qt = QTable(dat, names=['a', 'b'], units=(u.m, None))
  >>> qt = QTable(dat, names=['a', 'b'], units={'a': u.m})

See :ref:`quantity_and_qtable` for why we used a |QTable| here instead of a
|Table|.

descriptions
------------

This allows for setting the description for one or more columns at the time of
creating the table. The input can be either a list of description values
corresponding to each of the columns in the table (using `None` for no
description), or a :class:`dict` that provides the description for specified
column names. This works in the same way as the ``units`` example above.

.. _copy_versus_reference:

Copy versus Reference
=====================

Normally when a new |Table| object is created, the input data are *copied*.
This ensures that if the new table elements are modified then the original data
will not be affected. However, when creating a table from an existing |Table|,
a |ndarray| object (structured or homogeneous) or a :class:`dict`, it is
possible to disable copying so that a memory reference to the original data is
used instead. This has the advantage of being faster and using less memory.
However, caution must be exercised because the new table data and original data
will be linked, as shown below::

  >>> arr = np.array([(1, 2.0, 'x'),
  ...                 (4, 5.0, 'y')],
  ...                dtype=[('a', 'i8'), ('b', 'f8'), ('c', 'S2')])
  >>> print(arr['a'])  # column "a" of the input array
  [1 4]
  >>> t = Table(arr, copy=False)
  >>> t['a'][1] = 99
  >>> print(arr['a'])  # arr['a'] got changed when we modified t['a']
  [ 1 99]

Note that when referencing the data it is not possible to change the data types
since that operation requires making a copy of the data. In this case an error
occurs::

  >>> t = Table(arr, copy=False, dtype=('f4', 'i4', 'S4'))
  Traceback (most recent call last):
    ...
  ValueError: Cannot specify dtype when copy=False

Another caveat to using referenced data is that if you add a new row to the
table, the reference to the original data array is lost and the table will now
instead hold a copy of the original values (in addition to the new row).

Column and TableColumns Classes
===============================

There are two classes, |Column| and |TableColumns|, that are useful when
constructing new tables.

Column
------

A |Column| object can be created as follows, where in all cases the column
``name`` should be provided as a keyword argument and you can optionally provide
these values:

``data`` : :class:`list`, |ndarray| or `None`
    Column data values.
``dtype`` : :class:`numpy.dtype` compatible value
    Data type for column.
``description`` : :class:`str`
    Full description of column.
``unit`` : :class:`str`
    Physical unit.
``format`` : :class:`str` or function
    `Format specifier`_ for outputting column values.
``meta`` : :class:`dict`
    Metadata associated with the column.

Initialization Options
^^^^^^^^^^^^^^^^^^^^^^

The column data values, shape, and data type are specified in one of two ways:

**Provide data but not length or shape**

  Examples::

    col = Column([1, 2], name='a')  # shape=(2,)
    col = Column([[1, 2], [3, 4]], name='a')  # shape=(2, 2)
    col = Column([1, 2], name='a', dtype=float)
    col = Column(np.array([1, 2]), name='a')
    col = Column(['hello', 'world'], name='a')

  The ``dtype`` argument can be any value which is an acceptable fixed-size
  data type initializer for a :class:`numpy.dtype`. See the reference for
  `data type objects
  <https://numpy.org/doc/stable/reference/arrays.dtypes.html>`_. Examples
  include:

  - Python non-string type (:class:`float`, :class:`int`, :class:`bool`).
  - ``numpy`` non-string type (e.g., ``np.float32``, ``np.int64``).
  - ``numpy.dtype`` array-protocol type strings (e.g., ``'i4'``, ``'f8'``, ``'U15'``).

  If no ``dtype`` value is provided, then the type is inferred using
  :func:`numpy.array`. When ``data`` is provided then the ``shape``
  and ``length`` arguments are ignored.

**Provide length and optionally shape, but not data**

  Examples::

    col = Column(name='a', length=5)
    col = Column(name='a', dtype=int, length=10, shape=(3,4))

  The default ``dtype`` is ``np.float64``. The ``shape`` argument is the array
  shape of a single cell in the column. The default ``shape`` is ``()`` which means
  a single value in each element.

.. note::

   After setting the type for a column, that type cannot be changed.
   If data values of a different type are assigned to the column then they
   will be cast to the existing column type.

.. _table_format_string:

Format Specifier
^^^^^^^^^^^^^^^^

The format specifier controls the output of column values when a table or column
is printed or written to an ASCII table. In the simplest case, it is a string
that can be passed to Python's built-in :func:`format` function. For more
complicated formatting, one can also give "old style" or "new style"
format strings, or even a function:

**Plain format specification**

This type of string specifies directly how the value should be formatted
using a `format specification mini-language
<https://docs.python.org/3/library/string.html#formatspec>`_ that is
quite similar to C.

   ``".4f"`` will give four digits after the decimal in float format, or

   ``"6d"`` will give integers in six-character fields.

**Old style format string**

This corresponds to syntax like ``"%.4f" % value`` as documented in
`printf-style String Formatting
<https://docs.python.org/3/library/stdtypes.html#printf-style-string-formatting>`_.

   ``"%.4f"`` to print four digits after the decimal in float format, or

   ``"%6d"`` to print an integer in a six-character wide field.

**New style format string**

This corresponds to syntax like ``"{:.4f}".format(value)`` as documented in
`format string syntax
<https://docs.python.org/3/library/string.html#format-string-syntax>`_.

   ``"{:.4f}"`` to print four digits after the decimal in float format, or

   ``"{:6d}"`` to print an integer in a six-character wide field.

Note that in either format string case any Python string that formats exactly
one value is valid, so ``{:.4f} angstroms`` or ``Value: %12.2f`` would both
work.

**Function**

.. EXAMPLE START: Initialization Options for Column Objects

The greatest flexibility can be achieved by setting a formatting function. This
function must accept a single argument (the value) and return a string. In the
following example this is used to make a LaTeX ready output::

    >>> t = Table([[1,2],[1.234e9,2.34e-12]], names = ('a','b'))
    >>> def latex_exp(value):
    ...     val = f'{value:8.2}'
    ...     mant, exp = val.split('e')
    ...     # remove leading zeros
    ...     exp = exp[0] + exp[1:].lstrip('0')
    ...     return f'$ {mant} \\times 10^{{ {exp} }}$'
    >>> t['b'].format = latex_exp
    >>> t['a'].format = '.4f'
    >>> import sys
    >>> t.write(sys.stdout, format='latex')
    \begin{table}
    \begin{tabular}{cc}
    a & b \\
    1.0000 & $  1.2 \times 10^{ +9 }$ \\
    2.0000 & $  2.3 \times 10^{ -12 }$ \\
    \end{tabular}
    \end{table}

.. EXAMPLE END

TableColumns
------------

Each |Table| object has an attribute ``columns`` which is an ordered dictionary
that stores all of the |Column| objects in the table (see also the `Column`_
section). Technically, the ``columns`` attribute is a |TableColumns| object,
which is an enhanced ordered dictionary that provides easier ways to select
multiple columns. There are a few key points to remember:

- A |Table| can be initialized from a |TableColumns| object (``copy`` is always
  `True`).
- Selecting multiple columns from a |TableColumns| object returns another
  |TableColumns| object.
- Selecting one column from a |TableColumns| object returns a |Column|.

There are a few different ways to select columns from a |TableColumns| object:

**Select columns by name**
::

  >>> t = Table(names=('a', 'b', 'c', 'd'))

  >>> t.columns['d', 'c', 'b']
  <TableColumns names=('d','c','b')>

**Select columns by index slicing**
::

  >>> t.columns[0:2]  # Select first two columns
  <TableColumns names=('a','b')>

  >>> t.columns[::-1]  # Reverse column order
  <TableColumns names=('d','c','b','a')>

**Select single columns by index or name**
::

  >>> t.columns[1]  # Choose a column by index
  <Column name='b' dtype='float64' length=0>

  >>> t.columns['b']  # Choose a column by name
  <Column name='b' dtype='float64' length=0>

.. _subclassing_table:

Subclassing Table
=================

For some applications it can be useful to subclass the |Table| class in order
to introduce specialized behavior. Here we address two particular use cases
for subclassing: adding custom table attributes and changing the behavior of
internal class objects.

.. _table-custom-attributes:

Adding Custom Table Attributes
------------------------------

One simple customization that can be useful is adding new attributes to
the table object.  There is nothing preventing setting an attribute on an
existing table object, for example ``t.foo = 'hello'``.  However, this attribute
would be ephemeral because it will be lost if the table is sliced, copied, or
pickled. Instead, you can add persistent attributes as shown in this example::

  from astropy.table import Table, TableAttribute

  class MyTable(Table):
      foo = TableAttribute()
      bar = TableAttribute(default=[])
      baz = TableAttribute(default=1)

  t = MyTable([[1, 2]], foo='foo')
  t.bar.append(2.0)
  t.baz = 'baz'

Some key points:

- A custom attribute can be set when the table is created or using
  the usual syntax for setting an object attribute.
- A custom attribute always has a default value, either explicitly set
  in the class definition or `None`.
- The attribute values are stored in the table ``meta`` dictionary. This is
  the mechanism by which they are persistent through copy, slice, and
  serialization such as pickling or writing to an :ref:`ecsv_format` file.

Changing Behavior of Internal Class Objects
-------------------------------------------

It is also possible to change the behavior of the internal class objects which
are contained or created by a |Table|. This includes rows, columns, formatting,
and the columns container. In order to do this the subclass needs to declare
what class to use (if it is different from the built-in version). This is done
by specifying one or more of the class attributes ``Row``, ``Column``,
``MaskedColumn``, ``TableColumns``, or ``TableFormatter``.

The following trivial example overrides all of these with do-nothing
subclasses, but in practice you would override only the necessary
subcomponents::

  >>> from astropy.table import Table, Row, Column, MaskedColumn, TableColumns, TableFormatter

  >>> class MyRow(Row): pass
  >>> class MyColumn(Column): pass
  >>> class MyMaskedColumn(MaskedColumn): pass
  >>> class MyTableColumns(TableColumns): pass
  >>> class MyTableFormatter(TableFormatter): pass

  >>> class MyTable(Table):
  ...     """
  ...     Custom subclass of astropy.table.Table
  ...     """
  ...     Row = MyRow  # Use MyRow to create a row object
  ...     Column = MyColumn  # Column
  ...     MaskedColumn = MyMaskedColumn  # Masked Column
  ...     TableColumns = MyTableColumns  # Ordered dict holding Column objects
  ...     TableFormatter = MyTableFormatter  # Controls table output


Example
^^^^^^^

.. EXAMPLE START: Subclassing the Table Class

As a more practical example, suppose you have a table of data with a certain
set of fixed columns, but you also want to carry an arbitrary dictionary of
parameters for each row and then access those values using the same item access
syntax as if they were columns. It is assumed here that the extra parameters
are contained in a ``numpy`` object-dtype column named ``params``::

  >>> from astropy.table import Table, Row
  >>> class ParamsRow(Row):
  ...    """
  ...    Row class that allows access to an arbitrary dict of parameters
  ...    stored as a dict object in the ``params`` column.
  ...    """
  ...    def __getitem__(self, item):
  ...        if item not in self.colnames:
  ...            return super().__getitem__('params')[item]
  ...        else:
  ...            return super().__getitem__(item)
  ...
  ...    def keys(self):
  ...        out = [name for name in self.colnames if name != 'params']
  ...        params = [key.lower() for key in sorted(self['params'])]
  ...        return out + params
  ...
  ...    def values(self):
  ...        return [self[key] for key in self.keys()]

Now we put this into action with a trivial |Table| subclass::

  >>> class ParamsTable(Table):
  ...     Row = ParamsRow

First make a table and add a couple of rows::

  >>> t = ParamsTable(names=['a', 'b', 'params'], dtype=['i', 'f', 'O'])
  >>> t.add_row((1, 2.0, {'x': 1.5, 'y': 2.5}))
  >>> t.add_row((2, 3.0, {'z': 'hello', 'id': 123123}))
  >>> print(t)
   a   b             params
  --- --- ----------------------------
    1 2.0         {'x': 1.5, 'y': 2.5}
    2 3.0 {'z': 'hello', 'id': 123123}

Now see what we have from our specialized ``ParamsRow`` object::

  >>> t[0]['y']
  2.5
  >>> t[1]['id']
  123123
  >>> t[1].keys()
  ['a', 'b', 'id', 'z']
  >>> t[1].values()
  [2, 3.0, 123123, 'hello']

To make this example really useful, you might want to override
``Table.__getitem__()`` in order to allow table-level access to the parameter
fields. This might look something like::

  class ParamsTable(table.Table):
      Row = ParamsRow

      def __getitem__(self, item):
          if isinstance(item, str):
              if item in self.colnames:
                  return self.columns[item]
              else:
                  # If item is not a column name then create a new MaskedArray
                  # corresponding to self['params'][item] for each row.  This
                  # might not exist in some rows so mark as masked (missing) in
                  # those cases.
                  mask = np.zeros(len(self), dtype=np.bool_)
                  item = item.upper()
                  values = [params.get(item) for params in self['params']]
                  for ii, value in enumerate(values):
                      if value is None:
                          mask[ii] = True
                          values[ii] = ''
                  return self.MaskedColumn(name=item, data=values, mask=mask)

          # ... and then the rest of the original __getitem__ ...

.. EXAMPLE END

Columns and Quantities
======================

.. EXAMPLE START: Handling Astropy Column and Quantity Objects within Tables

``astropy`` `~astropy.units.Quantity` objects can be handled within tables in
two complementary ways. The first method stores the `~astropy.units.Quantity`
object natively within the table via the "mixin" column protocol. See the
sections on :ref:`mixin_columns` and :ref:`quantity_and_qtable` for details,
but in brief, the key difference is using the `~astropy.table.QTable` class to
indicate that a `~astropy.units.Quantity` should be stored natively within the
table::

  >>> from astropy.table import QTable
  >>> from astropy import units as u
  >>> t = QTable()
  >>> t['velocity'] = [3, 4] * u.m / u.s
  >>> type(t['velocity'])
  <class 'astropy.units.quantity.Quantity'>

For new code that is quantity-aware we recommend using `~astropy.table.QTable`,
but this may not be possible in all situations (particularly when interfacing
with legacy code that does not handle quantities) and there are
:ref:`details_and_caveats` that apply. In this case, use the
`~astropy.table.Table` class, which will convert a `~astropy.units.Quantity` to
a `~astropy.table.Column` object with a ``unit`` attribute::

  >>> from astropy.table import Table
  >>> t = Table()
  >>> t['velocity'] = [3, 4] * u.m / u.s
  >>> type(t['velocity'])
  <class 'astropy.table.column.Column'>
  >>> t['velocity'].unit
  Unit("m / s")

To learn more about using standard `~astropy.table.Column` objects with defined
units, see the :ref:`columns_with_units` section.

.. EXAMPLE END

.. _Table-like Objects:

Table-like Objects
==================

In order to improve interoperability between different table classes, an
``astropy`` |Table| object can be created directly from any other table-like
object that provides an ``__astropy_table__()`` method. In this case the
``__astropy_table__()`` method will be called as follows::

  >>> data = SomeOtherTableClass({'a': [1, 2], 'b': [3, 4]})  # doctest: +SKIP
  >>> t = QTable(data, copy=False, strict_copy=True)  # doctest: +SKIP

Internally the following call will be made to ask the ``data`` object
to return a representation of itself as an ``astropy`` |Table|, respecting
the ``copy`` preference of the original call to ``QTable()``::

  data.__astropy_table__(cls, copy, **kwargs)

Here ``cls`` is the |Table| class or subclass that is being instantiated
(|QTable| in this example), ``copy`` indicates whether a copy of the values in
``data`` should be provided, and ``**kwargs`` are any extra keyword arguments
which are not valid |Table| ``_init_()`` keyword arguments. In the example
above, ``strict_copy=True`` would end up in ``**kwargs`` and get passed to
``__astropy_table__()``.

If ``copy`` is `True` then the ``__astropy_table__()`` method must ensure that
a copy of the original data is returned. If ``copy`` is `False` then a
reference to the table data should be returned if possible. If it is not
possible (e.g., the original data are in a Python list or must be otherwise
transformed in memory) then ``__astropy_table__()`` method is free to either
return a copy or else raise an exception. This choice depends on the preference
of the implementation. The implementation might choose to allow an additional
keyword argument (e.g., ``strict_copy`` which gets passed via ``**kwargs``) to
control the behavior in this case.

As a concise example, imagine a dict-based table class. (Note that |Table|
already can be initialized from a dict-like object, so this is a bit contrived
but does illustrate the principles involved.) Please pay attention to the
method signature::

  def __astropy_table__(self, cls, copy, **kwargs):

Your class implementation of this must use the ``**kwargs`` technique for
catching keyword arguments at the end. This is to ensure future compatibility
in case additional keywords are added to the internal ``table =
data.__astropy_table__(cls, copy)`` call. Including ``**kwargs`` will prevent
breakage in this case. ::

  class DictTable(dict):
      """
      Trivial "table" class that just uses a dict to hold columns.
      This does not actually implement anything useful that makes
      this a table.

      The non-standard ``strict_copy=False`` keyword arg here will be passed
      via the **kwargs of Table __init__().
      """

      def __astropy_table__(self, cls, copy, strict_copy=False, **kwargs):
          """
          Return an astropy Table of type ``cls``.

          Parameters
          ----------
          cls : type
               Astropy ``Table`` class or subclass.
          copy : bool
               Copy input data (True) or return a reference (False).
          strict_copy : bool, optional
               Raise an exception if copy is False but reference is not
               possible.
          **kwargs : dict, optional
               Additional keyword args (ignored currently).
          """
          if kwargs:
              warnings.warn(f'unexpected keyword args {kwargs}')

          cols = list(self.values())
          names = list(self.keys())

          # If returning a reference to existing data (copy=False) and
          # strict_copy=True, make sure that each column is a numpy ndarray.
          # If a column is a Python list or tuple then it must be copied for
          # representation in an astropy Table.

          if not copy and strict_copy:
              for name, col in zip(names, cols):
                  if not isinstance(col, np.ndarray):
                      raise ValueError(f'cannot have copy=False because column {name} is '
                                       'not an ndarray')

          return cls(cols, names=names, copy=copy)
