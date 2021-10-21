.. _access_table:

Accessing a Table
*****************

Accessing table properties and data is generally consistent with the basic
interface for ``numpy`` `structured arrays
<https://numpy.org/doc/stable/user/basics.rec.html>`_.

Basics
======

For a quick overview, the code below shows the basics of accessing table data.
Where relevant, there is a comment about what sort of object is returned.
Except where noted, table access returns objects that can be modified in order
to update the original table data or properties. See also the section on
:ref:`copy_versus_reference` to learn more about this topic.

**Make a table**
::

  from astropy.table import Table
  import numpy as np

  arr = np.arange(15).reshape(5, 3)
  t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})

**Table properties**
::

  t.columns   # Dict of table columns (access by column name, index, or slice)
  t.colnames  # List of column names
  t.meta      # Dict of meta-data
  len(t)      # Number of table rows

**Access table data**
::

  t['a']       # Column 'a'
  t['a'][1]    # Row 1 of column 'a'
  t[1]         # Row 1
  t[1]['a']    # Column 'a' of row 1
  t[2:5]       # Table object with rows 2:5
  t[[1, 3, 4]]  # Table object with rows 1, 3, 4 (copy)
  t[np.array([1, 3, 4])]  # Table object with rows 1, 3, 4 (copy)
  t[[]]        # Same table definition but with no rows of data
  t['a', 'c']  # Table with cols 'a', 'c' (copy)
  dat = np.array(t)  # Copy table data to numpy structured array object
  t['a'].quantity  # an astropy.units.Quantity for Column 'a'
  t['a'].to('km')  # an astropy.units.Quantity for Column 'a' in units of kilometers
  t.columns[1]  # Column 1 (which is the 'b' column)
  t.columns[0:2]  # New table with columns 0 and 1

.. Note::
   Although they appear nearly equivalent, there is a factor of two performance
   difference between ``t[1]['a']`` (slower, because an intermediate |Row|
   object gets created) versus ``t['a'][1]`` (faster). Always use the latter
   when possible.

**Print table or column**
::

  print(t)     # Print formatted version of table to the screen
  t.pprint()   # Same as above
  t.pprint(show_unit=True)  # Show column unit
  t.pprint(show_name=False)  # Do not show column names
  t.pprint_all() # Print full table no matter how long / wide it is (same as t.pprint(max_lines=-1, max_width=-1))

  t.more()  # Interactively scroll through table like Unix "more"

  print(t['a'])    # Formatted column values
  t['a'].pprint()  # Same as above, with same options as Table.pprint()
  t['a'].more()    # Interactively scroll through column
  t['a', 'c'].pprint()  # Print columns 'a' and 'c' of table

  lines = t.pformat()  # Formatted table as a list of lines (same options as pprint)
  lines = t['a'].pformat()  # Formatted column values as a list


Details
=======

For all of the following examples it is assumed that the table has been created
as follows::

  >>> from astropy.table import Table, Column
  >>> import numpy as np
  >>> import astropy.units as u

  >>> arr = np.arange(15, dtype=np.int32).reshape(5, 3)
  >>> t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})
  >>> t['a'].format = "{:.3f}"  # print with 3 digits after decimal point
  >>> t['a'].unit = 'm sec^-1'
  >>> t['a'].description = 'unladen swallow velocity'
  >>> print(t)
       a      b   c
    m sec^-1
    -------- --- ---
       0.000   1   2
       3.000   4   5
       6.000   7   8
       9.000  10  11
      12.000  13  14

.. Note::

   In the example above the ``format``, ``unit``, and ``description``
   attributes of the |Column| were set directly. For :ref:`mixin_columns` like
   |Quantity| you must set via the ``info`` attribute, for example,
   ``t['a'].info.format = "{:.3f}"``. You can use the ``info`` attribute with
   |Column| objects as well, so the general solution that works with any table
   column is to set via the ``info`` attribute. See :ref:`mixin_attributes` for
   more information.

.. _table-summary-information:

Summary Information
-------------------

You can get summary information about the table as follows::

  >>> t.info
  <Table length=5>
  name dtype   unit   format       description
  ---- ----- -------- ------ ------------------------
     a int32 m sec^-1 {:.3f} unladen swallow velocity
     b int32
     c int32

If called as a function then you can supply an ``option`` that specifies
the type of information to return. The built-in ``option`` choices are
``'attributes'`` (column attributes, which is the default) or ``'stats'``
(basic column statistics). The ``option`` argument can also be a list
of available options::

  >>> t.info('stats')  # doctest: +FLOAT_CMP
  <Table length=5>
  name mean   std   min max
  ---- ---- ------- --- ---
     a    6 4.24264   0  12
     b    7 4.24264   1  13
     c    8 4.24264   2  14

  >>> t.info(['attributes', 'stats'])  # doctest: +FLOAT_CMP
  <Table length=5>
  name dtype   unit   format       description        mean   std   min max
  ---- ----- -------- ------ ------------------------ ---- ------- --- ---
     a int32 m sec^-1 {:.3f} unladen swallow velocity    6 4.24264   0  12
     b int32                                             7 4.24264   1  13
     c int32                                             8 4.24264   2  14

Columns also have an ``info`` property that has the same behavior and
arguments, but provides information about a single column::

  >>> t['a'].info
  name = a
  dtype = int32
  unit = m sec^-1
  format = {:.3f}
  description = unladen swallow velocity
  class = Column
  n_bad = 0
  length = 5

  >>> t['a'].info('stats')  # doctest: +FLOAT_CMP
  name = a
  mean = 6
  std = 4.24264
  min = 0
  max = 12
  n_bad = 0
  length = 5


Accessing Properties
--------------------

The code below shows accessing the table columns as a |TableColumns| object,
getting the column names, table metadata, and number of table rows. The table
metadata is an `~collections.OrderedDict` by default.
::

  >>> t.columns
  <TableColumns names=('a','b','c')>

  >>> t.colnames
  ['a', 'b', 'c']

  >>> t.meta  # Dict of meta-data
  {'keywords': {'key1': 'val1'}}

  >>> len(t)
  5


Accessing Data
--------------

As expected you can access a table column by name and get an element from that
column with a numerical index::

  >>> t['a']  # Column 'a'
  <Column name='a' dtype='int32' unit='m sec^-1' format='{:.3f}' description='unladen swallow velocity' length=5>
   0.000
   3.000
   6.000
   9.000
  12.000


  >>> t['a'][1]  # Row 1 of column 'a'
  3

When a table column is printed, it is formatted according to the ``format``
attribute (see :ref:`table_format_string`). Note the difference between the
column representation above and how it appears via ``print()`` or ``str()``::

  >>> print(t['a'])
     a
  m sec^-1
  --------
     0.000
     3.000
     6.000
     9.000
    12.000

Likewise a table row and a column from that row can be selected::

  >>> t[1]  # Row object corresponding to row 1
  <Row index=1>
     a       b     c
  m sec^-1
   int32   int32 int32
  -------- ----- -----
     3.000     4     5

  >>> t[1]['a']  # Column 'a' of row 1
  3

A |Row| object has the same columns and metadata as its parent table::

  >>> t[1].columns
  <TableColumns names=('a','b','c')>

  >>> t[1].meta
  {'keywords': {'key1': 'val1'}}

Slicing a table returns a new table object with references to the original
data within the slice region (See :ref:`copy_versus_reference`). The table
metadata and column definitions are copied.
::

  >>> t[2:5]  # Table object with rows 2:5 (reference)
  <Table length=3>
     a       b     c
  m sec^-1
   int32   int32 int32
  -------- ----- -----
     6.000     7     8
     9.000    10    11
    12.000    13    14

It is possible to select table rows with an array of indexes or by specifying
multiple column names. This returns a copy of the original table for the
selected rows or columns.  ::

  >>> print(t[[1, 3, 4]])  # Table object with rows 1, 3, 4 (copy)
       a      b   c
    m sec^-1
    -------- --- ---
       3.000   4   5
       9.000  10  11
      12.000  13  14


  >>> print(t[np.array([1, 3, 4])])  # Table object with rows 1, 3, 4 (copy)
       a      b   c
    m sec^-1
    -------- --- ---
       3.000   4   5
       9.000  10  11
      12.000  13  14


  >>> print(t['a', 'c'])  # or t[['a', 'c']] or t[('a', 'c')]
  ...                     # Table with cols 'a', 'c' (copy)
       a      c
    m sec^-1
    -------- ---
       0.000   2
       3.000   5
       6.000   8
       9.000  11
      12.000  14

We can select rows from a table using conditionals to create boolean masks. A
table indexed with a boolean array will only return rows where the mask array
element is `True`. Different conditionals can be combined using the bitwise
operators.  ::

  >>> mask = (t['a'] > 4) & (t['b'] > 8)  # Table rows where column a > 4
  >>> print(t[mask])                      # and b > 8
  ...
       a      b   c
    m sec^-1
    -------- --- ---
       9.000  10  11
      12.000  13  14

Finally, you can access the underlying table data as a native ``numpy``
structured array by creating a copy or reference with :func:`numpy.array`::

  >>> data = np.array(t)  # copy of data in t as a structured array
  >>> data = np.array(t, copy=False)  # reference to data in t


Table Equality
--------------

We can check table data equality using two different methods:

- The ``==`` comparison operator. This returns a `True` or `False` for
  each row if the *entire row* matches. This is the same as the behavior of
  ``numpy`` structured arrays.
- Table :meth:`~astropy.table.Table.values_equal` to compare table values
  element-wise. This returns a boolean `True` or `False` for each table
  *element*, so you get a `~astropy.table.Table` of values.

Examples
^^^^^^^^

.. EXAMPLE START: Checking Table Equality

To check table equality::

  >>> t1 = Table(rows=[[1, 2, 3],
  ...                  [4, 5, 6],
  ...                  [7, 7, 9]], names=['a', 'b', 'c'])
  >>> t2 = Table(rows=[[1, 2, -1],
  ...                  [4, -1, 6],
  ...                  [7, 7, 9]], names=['a', 'b', 'c'])

  >>> t1 == t2
  array([False, False,  True])

  >>> t1.values_equal(t2)  # Compare to another table
  <Table length=3>
   a     b     c
  bool  bool  bool
  ---- ----- -----
  True  True False
  True False  True
  True  True  True

  >>> t1.values_equal([2, 4, 7])  # Compare to an array column-wise
  <Table length=3>
    a     b     c
   bool  bool  bool
  ----- ----- -----
  False  True False
   True False False
   True  True False

  >>> t1.values_equal(7)  # Compare to a scalar column-wise
  <Table length=3>
    a     b     c
   bool  bool  bool
  ----- ----- -----
  False False False
  False False False
   True  True False

.. EXAMPLE END

Formatted Printing
------------------

The values in a table or column can be printed or retrieved as a formatted
table using one of several methods:

- `print()` function.
- `Table.more() <astropy.table.Table.more>` or `Column.more()
  <astropy.table.Column.more>` methods to interactively scroll through
  table values.
- `Table.pprint() <astropy.table.Table.pprint>` or `Column.pprint()
  <astropy.table.Column.pprint>` methods to print a formatted version of
  the table to the screen.
- `Table.pformat() <astropy.table.Table.pformat>` or `Column.pformat()
  <astropy.table.Column.pformat>` methods to return the formatted table
  or column as a list of fixed-width strings. This could be used as a quick way
  to save a table.

These methods use :ref:`table_format_string`
if available and strive to make the output readable.
By default, table and column printing will
not print the table larger than the available interactive screen size. If the
screen size cannot be determined (in a non-interactive environment or on
Windows) then a default size of 25 rows by 80 columns is used. If a table is
too large, then rows and/or columns are cut from the middle so it fits.

Example
^^^^^^^

.. EXAMPLE START: Printing Formatted Tables

To print a formatted table::

  >>> arr = np.arange(3000).reshape(100, 30)  # 100 rows x 30 columns array
  >>> t = Table(arr)
  >>> print(t)
  col0 col1 col2 col3 col4 col5 col6 ... col23 col24 col25 col26 col27 col28 col29
  ---- ---- ---- ---- ---- ---- ---- ... ----- ----- ----- ----- ----- ----- -----
     0    1    2    3    4    5    6 ...    23    24    25    26    27    28    29
    30   31   32   33   34   35   36 ...    53    54    55    56    57    58    59
    60   61   62   63   64   65   66 ...    83    84    85    86    87    88    89
    90   91   92   93   94   95   96 ...   113   114   115   116   117   118   119
   120  121  122  123  124  125  126 ...   143   144   145   146   147   148   149
   150  151  152  153  154  155  156 ...   173   174   175   176   177   178   179
   180  181  182  183  184  185  186 ...   203   204   205   206   207   208   209
   210  211  212  213  214  215  216 ...   233   234   235   236   237   238   239
   240  241  242  243  244  245  246 ...   263   264   265   266   267   268   269
   270  271  272  273  274  275  276 ...   293   294   295   296   297   298   299
   ...  ...  ...  ...  ...  ...  ... ...   ...   ...   ...   ...   ...   ...   ...
  2670 2671 2672 2673 2674 2675 2676 ...  2693  2694  2695  2696  2697  2698  2699
  2700 2701 2702 2703 2704 2705 2706 ...  2723  2724  2725  2726  2727  2728  2729
  2730 2731 2732 2733 2734 2735 2736 ...  2753  2754  2755  2756  2757  2758  2759
  2760 2761 2762 2763 2764 2765 2766 ...  2783  2784  2785  2786  2787  2788  2789
  2790 2791 2792 2793 2794 2795 2796 ...  2813  2814  2815  2816  2817  2818  2819
  2820 2821 2822 2823 2824 2825 2826 ...  2843  2844  2845  2846  2847  2848  2849
  2850 2851 2852 2853 2854 2855 2856 ...  2873  2874  2875  2876  2877  2878  2879
  2880 2881 2882 2883 2884 2885 2886 ...  2903  2904  2905  2906  2907  2908  2909
  2910 2911 2912 2913 2914 2915 2916 ...  2933  2934  2935  2936  2937  2938  2939
  2940 2941 2942 2943 2944 2945 2946 ...  2963  2964  2965  2966  2967  2968  2969
  2970 2971 2972 2973 2974 2975 2976 ...  2993  2994  2995  2996  2997  2998  2999
  Length = 100 rows

.. EXAMPLE END

more() method
^^^^^^^^^^^^^

In order to browse all rows of a table or column use the `Table.more()
<astropy.table.Table.more>` or `Column.more() <astropy.table.Column.more>`
methods. These let you interactively scroll through the rows much like the Unix
``more`` command. Once part of the table or column is displayed the supported
navigation keys are:

|  **f, space** : forward one page
|  **b** : back one page
|  **r** : refresh same page
|  **n** : next row
|  **p** : previous row
|  **<** : go to beginning
|  **>** : go to end
|  **q** : quit browsing
|  **h** : print this help

pprint() method
^^^^^^^^^^^^^^^

In order to fully control the print output use the `Table.pprint()
<astropy.table.Table.pprint>` or `Column.pprint()
<astropy.table.Column.pprint>` methods. These have keyword arguments
``max_lines``, ``max_width``, ``show_name``, and ``show_unit``, with meanings
as shown below::

  >>> arr = np.arange(3000, dtype=float).reshape(100, 30)
  >>> t = Table(arr)
  >>> t['col0'].format = '%e'
  >>> t['col0'].unit = 'km**2'
  >>> t['col29'].unit = 'kg sec m**-2'

  >>> t.pprint(max_lines=8, max_width=40)
      col0     ...    col29
      km2      ... kg sec m**-2
  ------------ ... ------------
  0.000000e+00 ...         29.0
           ... ...          ...
  2.940000e+03 ...       2969.0
  2.970000e+03 ...       2999.0
  Length = 100 rows

  >>> t.pprint(max_lines=8, max_width=40, show_unit=False)
      col0     ... col29
  ------------ ... ------
  0.000000e+00 ...   29.0
           ... ...    ...
  2.940000e+03 ... 2969.0
  2.970000e+03 ... 2999.0
  Length = 100 rows

  >>> t.pprint(max_lines=8, max_width=40, show_name=False)
      km2      ... kg sec m**-2
  ------------ ... ------------
  0.000000e+00 ...         29.0
  3.000000e+01 ...         59.0
           ... ...          ...
  2.940000e+03 ...       2969.0
  2.970000e+03 ...       2999.0
  Length = 100 rows

In order to force printing all values regardless of the output length or width
use :meth:`~astropy.table.Table.pprint_all`, which is equivalent to setting
``max_lines`` and ``max_width`` to ``-1`` in :meth:`~astropy.table.Table.pprint`.
:meth:`~astropy.table.Table.pprint_all` takes the same arguments as :meth:`~astropy.table.Table.pprint`.
For the wide table in this example you see six lines of wrapped output like the
following::

  >>> t.pprint_all(max_lines=8)  # doctest: +SKIP
      col0         col1     col2   col3   col4   col5   col6   col7   col8   col9  col10  col11  col12  col13  col14  col15  col16  col17  col18  col19  col20  col21  col22  col23  col24  col25  col26  col27  col28     col29
      km2                                                                                                                                                                                                               kg sec m**-2
  ------------ ----------- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------------
  0.000000e+00    1.000000    2.0    3.0    4.0    5.0    6.0    7.0    8.0    9.0   10.0   11.0   12.0   13.0   14.0   15.0   16.0   17.0   18.0   19.0   20.0   21.0   22.0   23.0   24.0   25.0   26.0   27.0   28.0         29.0
           ...         ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...          ...
  2.940000e+03 2941.000000 2942.0 2943.0 2944.0 2945.0 2946.0 2947.0 2948.0 2949.0 2950.0 2951.0 2952.0 2953.0 2954.0 2955.0 2956.0 2957.0 2958.0 2959.0 2960.0 2961.0 2962.0 2963.0 2964.0 2965.0 2966.0 2967.0 2968.0       2969.0
  2.970000e+03 2971.000000 2972.0 2973.0 2974.0 2975.0 2976.0 2977.0 2978.0 2979.0 2980.0 2981.0 2982.0 2983.0 2984.0 2985.0 2986.0 2987.0 2988.0 2989.0 2990.0 2991.0 2992.0 2993.0 2994.0 2995.0 2996.0 2997.0 2998.0       2999.0
  Length = 100 rows

For columns, the syntax and behavior of :func:`~astropy.table.Column.pprint` is
the same except that there is no ``max_width`` keyword argument::

  >>> t['col3'].pprint(max_lines=8)
   col3
  ------
     3.0
    33.0
     ...
  2943.0
  2973.0
  Length = 100 rows

Column alignment
^^^^^^^^^^^^^^^^

Individual columns have the ability to be aligned in a number of different
ways for an enhanced viewing experience::

  >>> t1 = Table()
  >>> t1['long column name 1'] = [1, 2, 3]
  >>> t1['long column name 2'] = [4, 5, 6]
  >>> t1['long column name 3'] = [7, 8, 9]
  >>> t1['long column name 4'] = [700000, 800000, 900000]
  >>> t1['long column name 2'].info.format = '<'
  >>> t1['long column name 3'].info.format = '0='
  >>> t1['long column name 4'].info.format = '^'
  >>> t1.pprint()
   long column name 1 long column name 2 long column name 3 long column name 4
  ------------------ ------------------ ------------------ ------------------
                   1 4                  000000000000000007       700000
                   2 5                  000000000000000008       800000
                   3 6                  000000000000000009       900000

Conveniently, alignment can be handled another way — by passing a list to the
keyword argument ``align``::

  >>> t1 = Table()
  >>> t1['column1'] = [1, 2, 3]
  >>> t1['column2'] = [2, 4, 6]
  >>> t1.pprint(align=['<', '0='])
  column1 column2
  ------- -------
  1       0000002
  2       0000004
  3       0000006

It is also possible to set the alignment of all columns with a single
string value::

  >>> t1.pprint(align='^')
  column1 column2
  ------- -------
     1       2
     2       4
     3       6

The fill character for justification can be set as a prefix to the
alignment character (see `Format Specification Mini-Language
<https://docs.python.org/3/library/string.html#format-specification-mini-language>`_
for additional explanation). This can be done both in the ``align`` argument
and in the column ``format`` attribute. Note the interesting interaction below::

  >>> t1 = Table([[1.0, 2.0], [1, 2]], names=['column1', 'column2'])

  >>> t1['column1'].format = '#^.2f'
  >>> t1.pprint()
  column1 column2
  ------- -------
  ##1.00#       1
  ##2.00#       2

Now if we set a global align, it seems like our original column format
got lost::

  >>> t1.pprint(align='!<')
  column1 column2
  ------- -------
  1.00!!! 1!!!!!!
  2.00!!! 2!!!!!!

The way to avoid this is to explicitly specify the alignment strings
for every column and use `None` where the column format should be
used::

  >>> t1.pprint(align=[None, '!<'])
  column1 column2
  ------- -------
  ##1.00# 1!!!!!!
  ##2.00# 2!!!!!!

pformat() method
^^^^^^^^^^^^^^^^

In order to get the formatted output for manipulation or writing to a file use
the `Table.pformat() <astropy.table.Table.pformat>` or `Column.pformat()
<astropy.table.Column.pformat>` methods. These behave just as for
:meth:`~astropy.table.Table.pprint` but return a list corresponding to each
formatted line in the :meth:`~astropy.table.Table.pprint` output. The
:meth:`~astropy.table.Table.pformat_all` method can be used to return a list
for all lines in the |Table|.

  >>> lines = t['col3'].pformat(max_lines=8)

Hiding columns
^^^^^^^^^^^^^^

The |Table| class has functionality to selectively show or hide certain columns
within the table when using any of the print methods. This can be useful for
columns that are very wide or else "uninteresting" for various reasons. The
specification of which columns are outputted is associated with the table itself
so that it persists through slicing, copying, and serialization (e.g. saving to
:ref:`ecsv_format`). One use case is for specialized table subclasses that
contain auxiliary columns that are not typically useful to the user.

The specification of which columns to include when printing is handled through
two complementary |Table| attributes:

- `~astropy.table.Table.pprint_include_names`: column names to include, where
  the default value of `None` implies including all columns.
- `~astropy.table.Table.pprint_exclude_names`: column names to exclude, where
  the default value of `None` implies excluding no columns.

Typically you should use just one of the two attributes at a time. However,
both can be set at once and the set of columns that actually gets printed
is conceptually expressed in this pseudo-code::

  include_names = (set(table.pprint_include_names() or table.colnames)
                   - set(table.pprint_exclude_names() or ())

Examples
""""""""
Let's start with defining a simple table with one row and six columns::

  >>> from astropy.table.table_helpers import simple_table
  >>> t = simple_table(size=1, cols=6)
  >>> print(t)
  a   b   c   d   e   f
  --- --- --- --- --- ---
  1 1.0   c   4 4.0   f

Now you can get the value of the ``pprint_include_names`` attribute by calling
it as a function, and then include some names for printing::

  >>> print(t.pprint_include_names())
  None
  >>> t.pprint_include_names = ('a', 'c', 'e')
  >>> print(t.pprint_include_names())
  ('a', 'c', 'e')
  >>> print(t)
   a   c   e
  --- --- ---
    1   c 4.0

Now you can instead exclude some columns from printing. Note that for both
include and exclude, you can add column names that do not exist in the table.
This allows pre-defining the attributes before the table has been fully
constructed.
::

  >>> t.pprint_include_names = None  # Revert to printing all columns
  >>> t.pprint_exclude_names = ('a', 'c', 'e', 'does-not-exist')
  >>> print(t)
   b   d   f
  --- --- ---
  1.0   4   f

Next you can ``add`` or ``remove`` names from the attribute::

  >>> t = simple_table(size=1, cols=6)  # Start with a fresh table
  >>> t.pprint_exclude_names.add('b')  # Single name
  >>> t.pprint_exclude_names.add(['d', 'f'])  # List or tuple of names
  >>> t.pprint_exclude_names.remove('f')  # Single name or list/tuple of names
  >>> t.pprint_exclude_names()
  ('b', 'd')

Finally, you can temporarily set the attributes within a `context manager
<https://docs.python.org/3/reference/datamodel.html#context-managers>`_. For
example::

  >>> t = simple_table(size=1, cols=6)
  >>> t.pprint_include_names = ('a', 'b')
  >>> print(t)
   a   b
  --- ---
    1 1.0

  >>> # Show all (for pprint_include_names the value of None => all columns)
  >>> with t.pprint_include_names.set(None):
  ...     print(t)
   a   b   c   d   e   f
  --- --- --- --- --- ---
    1 1.0   c   4 4.0   f

The specification of names for these attributes can include Unix-style globs
like ``*`` and ``?``. See `fnmatch` for details (and in particular how to
escape those characters if needed). For example::

  >>> t = Table()
  >>> t.pprint_exclude_names = ['boring*']
  >>> t['a'] = [1]
  >>> t['b'] = ['b']
  >>> t['boring_ra'] = [122.0]
  >>> t['boring_dec'] = [89.9]
  >>> print(t)
   a   b
  --- ---
    1   b

Multidimensional columns
^^^^^^^^^^^^^^^^^^^^^^^^

If a column has more than one dimension then each element of the column is
itself an array. In the example below there are three rows, each of which is a
``2 x 2`` array. The formatted output for such a column shows only the first
and last value of each row element and indicates the array dimensions in the
column name header::

  >>> t = Table()
  >>> arr = [ np.array([[ 1,  2],
  ...                   [10, 20]]),
  ...         np.array([[ 3,  4],
  ...                   [30, 40]]),
  ...         np.array([[ 5,  6],
  ...                   [50, 60]]) ]
  >>> t['a'] = arr
  >>> t['a'].shape
  (3, 2, 2)
  >>> t.pprint()
  a [2,2]
  -------
  1 .. 20
  3 .. 40
  5 .. 60

In order to see all of the data values for a multidimensional column use the
column representation. This uses the standard ``numpy`` mechanism for printing
any array::

  >>> t['a'].data
  array([[[ 1,  2],
          [10, 20]],
         [[ 3,  4],
          [30, 40]],
         [[ 5,  6],
          [50, 60]]])

.. _columns_with_units:

Columns with Units
^^^^^^^^^^^^^^^^^^

.. note::

  |Table| and |QTable| instances handle entries with units differently. The
  following describes |Table|. :ref:`quantity_and_qtable` explains how a
  |QTable| differs from a |Table|.

A |Column| object with units within a standard |Table| has certain
quantity-related conveniences available. To begin with, it can be converted
explicitly to a |Quantity| object via the
:attr:`~astropy.table.Column.quantity` property and the
:meth:`~astropy.table.Column.to` method::

  >>> data = [[1., 2., 3.], [40000., 50000., 60000.]]
  >>> t = Table(data, names=('a', 'b'))
  >>> t['a'].unit = u.m
  >>> t['b'].unit = 'km/s'
  >>> t['a'].quantity  # doctest: +FLOAT_CMP
  <Quantity [1., 2., 3.] m>
  >>> t['b'].to(u.kpc/u.Myr)  # doctest: +FLOAT_CMP
  <Quantity [40.9084866 , 51.13560825, 61.3627299 ] kpc / Myr>

Note that the :attr:`~astropy.table.Column.quantity` property is actually
a *view* of the data in the column, not a copy. Hence, you can set the
values of a column in a way that respects units by making in-place
changes to the :attr:`~astropy.table.Column.quantity` property::

  >>> t['b']
  <Column name='b' dtype='float64' unit='km / s' length=3>
  40000.0
  50000.0
  60000.0

  >>> t['b'].quantity[0] = 45000000*u.m/u.s
  >>> t['b']
  <Column name='b' dtype='float64' unit='km / s' length=3>
  45000.0
  50000.0
  60000.0

Even without explicit conversion, columns with units can be treated like a
|Quantity| in *some* arithmetic expressions (see the warning below for caveats
to this)::

  >>> t['a'] + .005*u.km  # doctest: +FLOAT_CMP
  <Quantity [6., 7., 8.] m>
  >>> from astropy.constants import c
  >>> (t['b'] / c).decompose()  # doctest: +FLOAT_CMP
  <Quantity [0.15010384, 0.16678205, 0.20013846]>

.. warning::

  |Table| columns do *not* always behave the same as |Quantity|. |Table|
  columns act more like regular ``numpy`` arrays unless either explicitly
  converted to a |Quantity| or combined with a |Quantity| using an arithmetic
  operator. For example, the following does not work in the way you would
  expect::

    >>> data = [[30, 90]]
    >>> t = Table(data, names=('angle',))
    >>> t['angle'].unit = 'deg'
    >>> np.sin(t['angle'])  # doctest: +FLOAT_CMP
    <Column name='angle' dtype='float64' unit='deg' length=2>
    -0.988031624093
     0.893996663601

  This is wrong both in that it says the result is in degrees, *and*
  `~numpy.sin` treated the values as radians rather than degrees. If at all in
  doubt that you will get the right result, the safest choice is to either use
  |QTable| or to explicitly convert to |Quantity|::

    >>> np.sin(t['angle'].quantity)  # doctest: +FLOAT_CMP
    <Quantity [0.5, 1. ]>

.. _bytestring-columns-python-3:

Bytestring Columns
^^^^^^^^^^^^^^^^^^

Using bytestring columns (``numpy`` ``'S'`` dtype) is possible
with ``astropy`` tables since they can be compared with the natural
Python string (``str``) type. See `The bytes/str dichotomy in Python 3
<https://eli.thegreenplace.net/2012/01/30/the-bytesstr-dichotomy-in-python-3>`_
for a very brief overview of the difference.

The standard method of representing strings in ``numpy`` is via the
unicode ``'U'`` dtype. The problem is that this requires 4 bytes per
character, and if you have a very large number of strings this could
fill memory and impact performance. A very common use case is that these
strings are actually ASCII and can be represented with 1 byte per character.
In ``astropy`` it is possible to work directly and conveniently with
bytestring data in |Table| and |Column| operations.

Note that the bytestring issue is a particular problem when dealing with HDF5
files, where character data are read as bytestrings (``'S'`` dtype) when using
the :ref:`table_io`. Since HDF5 files are frequently used to store very large
datasets, the memory bloat associated with conversion to ``'U'`` dtype is
unacceptable.


Examples
""""""""

.. EXAMPLE START: Bytestring Data in Astropy Tables

The examples below illustrate dealing with bytestring data in ``astropy``::

    >>> t = Table([['abc', 'def']], names=['a'], dtype=['S'])

    >>> t['a'] == 'abc'  # Gives expected answer
    array([ True, False])

    >>> t['a'] == b'abc'  # Still gives expected answer
    array([ True, False])

    >>> t['a'][0] == 'abc'  # Expected answer
    True

    >>> t['a'][0] == b'abc'  # Cannot compare to bytestring
    False

    >>> t['a'][0] = 'bä'
    >>> t
    <Table length=2>
      a
    bytes3
    ------
        bä
       def

    >>> t['a'] == 'bä'
    array([ True, False])

.. doctest-skip::

    >>> # Round trip unicode strings through HDF5
    >>> t.write('test.hdf5', format='hdf5', path='data', overwrite=True)
    >>> t2 = Table.read('test.hdf5', format='hdf5', path='data')
    >>> t2
    <Table length=2>
     col0
    bytes3
    ------
        bä
       def

.. EXAMPLE END
