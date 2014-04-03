.. _access_table:

.. include:: references.txt

Accessing a table
-----------------

Accessing the table properties and data is straightforward and is generally consistent with
the basic interface for `numpy` structured arrays.

Quick overview
^^^^^^^^^^^^^^

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

  t.columns   # Dict of table columns
  t.colnames  # List of column names
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

**Print table or column**
::

  print t      # Print formatted version of table to the screen
  t.pprint()   # Same as above
  t.pprint(show_unit=True)  # Show column unit
  t.pprint(show_name=False)  # Do not show column names
  t.pprint(max_lines=-1, max_width=-1)  # Print full table no matter how long / wide it is

  t.more()  # Interactively scroll through table like Unix "more"

  print t['a'] # Formatted column values
  t['a'].pprint()  # Same as above, with same options as Table.pprint()
  t['a'].more()  # Interactively scroll through column

  lines = t.pformat()  # Formatted table as a list of lines (same options as pprint)
  lines = t['a'].pformat()  # Formatted column values as a list


Details
^^^^^^^

For all the following examples it is assumed that the table has been created as below::

  >>> from astropy.table import Table, Column
  >>> import numpy as np

  >>> arr = np.arange(15).reshape(5, 3)
  >>> t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})
  >>> t['a'].format = "%6.3f"  # print as a float with 3 digits after decimal point
  >>> t['a'].unit = 'm sec^-1'
  >>> t['a'].description = 'unladen swallow velocity'
  >>> print t
       a      b   c
    m sec^-1
    -------- --- ---
       0.000   1   2
       3.000   4   5
       6.000   7   8
       9.000  10  11
      12.000  13  14

Accessing properties
""""""""""""""""""""

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


Accessing data
""""""""""""""

As expected one can access a table column by name and get an element from that
column with a numerical index::

  >>> t['a']  # Column 'a'
  <Column name='a' unit='m sec^-1' format='%6.3f' description='unladen swallow velocity'>
  array([ 0,  3,  6,  9, 12])

  >>> t['a'][1]  # Row 1 of column 'a'
  3

When a table column is printed, either with ``print`` or via the ``str()``
built-in function, it is formatted according to the ``format`` attribute (see
:ref:`table_format_string`)::

  >>> print(t['a'])
      a
    ------
     0.000
     3.000
     6.000
     9.000
    12.000

Likewise a table row and a column from that row can be selected::

  >>> t[1]  # Row object corresponding to row 1
  <Row 1 of table
   values=(3, 4, 5)
   dtype=[('a', '<i4'), ('b', '<i8'), ('c', '<i8')]>

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

It is possible to select table rows with an array of indexes or by specifying
multiple column names.  This returns a copy of the original table for the
selected rows or columns.  ::

  >>> print t[[1, 3, 4]]  # Table object with rows 1, 3, 4 (copy)
       a      b   c
    m sec^-1
    -------- --- ---
       3.000   4   5
       9.000  10  11
      12.000  13  14


  >>> print t[np.array([1, 3, 4])]  # Table object with rows 1, 3, 4 (copy)
       a      b   c
    m sec^-1
    -------- --- ---
       3.000   4   5
       9.000  10  11
      12.000  13  14


  >>> print t['a', 'c']  # or t[['a', 'c']] or t[('a', 'c')]
  ...                    # Table with cols 'a', 'c' (copy)
       a      c
    m sec^-1
    -------- ---
       0.000   2
       3.000   5
       6.000   8
       9.000  11
      12.000  14

Finally, one can access the underlying table data as a native `numpy`
structured array by creating a copy or reference with ``np.array``::

  >>> data = np.array(t)  # copy of data in t as a structured array
  >>> data = np.array(t, copy=False)  # reference to data in t


Formatted printing
""""""""""""""""""

The values in a table or column can be printed or retrieved as a formatted
table using one of several methods:

- `print` statement (Python 2) or `print()` function (Python 3).
- Table :meth:`~astropy.table.Table.more` or Column
  :meth:`~astropy.table.table.Column.more` methods to interactively scroll
  through table values.
- Table :meth:`~astropy.table.Table.pprint` or Column
  :func:`~astropy.table.table.Column.pprint` methods to print a formatted version of
  the table to the screen.
- Table :meth:`~astropy.table.Table.pformat` or Column
  :func:`~astropy.table.table.Column.pformat` methods to return the formatted table
  or column as a list of fixed-width strings.  This could be used as a quick
  way to save a table.

These methods use :ref:`table_format_string`
if available and strive to make the output readable.
By default, table and column printing will
not print the table larger than the available interactive screen size.  If the
screen size cannot be determined (in a non-interactive environment or on
Windows) then a default size of 25 rows by 80 columns is used.  If a table is
too large then rows and/or columns are cut from the middle so it fits.  For example::

  >>> arr = np.arange(3000).reshape(100, 30)  # 100 rows x 30 columns array
  >>> t = Table(arr)
  >>> print t
  col0 col1 col2 col3 col4 col5 col6 ... col24 col25 col26 col27 col28 col29
  ---- ---- ---- ---- ---- ---- ---- ... ----- ----- ----- ----- ----- -----
     0    1    2    3    4    5    6 ...    24    25    26    27    28    29
    30   31   32   33   34   35   36 ...    54    55    56    57    58    59
    60   61   62   63   64   65   66 ...    84    85    86    87    88    89
    90   91   92   93   94   95   96 ...   114   115   116   117   118   119
   120  121  122  123  124  125  126 ...   144   145   146   147   148   149
   150  151  152  153  154  155  156 ...   174   175   176   177   178   179
   180  181  182  183  184  185  186 ...   204   205   206   207   208   209
   210  211  212  213  214  215  216 ...   234   235   236   237   238   239
   240  241  242  243  244  245  246 ...   264   265   266   267   268   269
   ...  ...  ...  ...  ...  ...  ... ...   ...   ...   ...   ...   ...   ...
  2760 2761 2762 2763 2764 2765 2766 ...  2784  2785  2786  2787  2788  2789
  2790 2791 2792 2793 2794 2795 2796 ...  2814  2815  2816  2817  2818  2819
  2820 2821 2822 2823 2824 2825 2826 ...  2844  2845  2846  2847  2848  2849
  2850 2851 2852 2853 2854 2855 2856 ...  2874  2875  2876  2877  2878  2879
  2880 2881 2882 2883 2884 2885 2886 ...  2904  2905  2906  2907  2908  2909
  2910 2911 2912 2913 2914 2915 2916 ...  2934  2935  2936  2937  2938  2939
  2940 2941 2942 2943 2944 2945 2946 ...  2964  2965  2966  2967  2968  2969
  2970 2971 2972 2973 2974 2975 2976 ...  2994  2995  2996  2997  2998  2999

more() method
'''''''''''''

In order to browse all rows of a table or column use the Table
:meth:`~astropy.table.Table.more` or Column :func:`~astropy.table.table.Column.more`
methods.  These let you interactively scroll through the rows much like the
linux ``more`` command.  Once part of the table or column is displayed the
supported navigation keys are:

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
'''''''''''''''

In order to fully control the print output use the Table
:meth:`~astropy.table.Table.pprint` or Column
:func:`~astropy.table.table.Column.pprint` methods.  These have keyword
arguments ``max_lines``, ``max_width``, ``show_name``, ``show_unit`` with
meaning as shown below::

  >>> arr = np.arange(3000, dtype=float).reshape(100, 30)
  >>> t = Table(arr)
  >>> t['col0'].format = '%e'
  >>> t['col1'].format = '%.6f'
  >>> t['col0'].unit = 'km**2'
  >>> t['col29'].unit = 'kg sec m**-2'

  >>> t.pprint(max_lines=8, max_width=40)
        col0     ...    col29
        km2      ... kg sec m**-2
    ------------ ... ------------
    0.000000e+00 ...         29.0
    3.000000e+01 ...         59.0
             ... ...          ...
    2.940000e+03 ...       2969.0
    2.970000e+03 ...       2999.0


  >>> t.pprint(max_lines=8, max_width=40, show_unit=True)
      col0     ...    col29
      km2      ... kg sec m**-2
  ------------ ... ------------
  0.000000e+00 ...         29.0
  3.000000e+01 ...         59.0
           ... ...          ...
  2.940000e+03 ...       2969.0
  2.970000e+03 ...       2999.0

  >>> t.pprint(max_lines=8, max_width=40, show_name=False)
        km2      ... kg sec m**-2
    ------------ ... ------------
    0.000000e+00 ...         29.0
    3.000000e+01 ...         59.0
    6.000000e+01 ...         89.0
             ... ...          ...
    2.940000e+03 ...       2969.0
    2.970000e+03 ...       2999.0

In order to force printing all values regardless of the output length or width
set ``max_lines`` or ``max_width`` to ``-1``, respectively.  For the wide
table in this example one sees 6 lines of wrapped output like the following::

  >>> t.pprint(max_lines=6, max_width=-1)  # doctest: +SKIP
      col0         col1     col2   col3   col4   col5   col6   col7   col8   col
  9  col10  col11  col12  col13  col14  col15  col16  col17  col18  col19  col20
    col21  col22  col23  col24  col25  col26  col27  col28  col29
  ------------ ----------- ------ ------ ------ ------ ------ ------ ------ ----
  -- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ -----
  - ------ ------ ------ ------ ------ ------ ------ ------ ------
  0.000000e+00    1.000000    2.0    3.0    4.0    5.0    6.0    7.0    8.0    9
  .0   10.0   11.0   12.0   13.0   14.0   15.0   16.0   17.0   18.0   19.0   20.
  0   21.0   22.0   23.0   24.0   25.0   26.0   27.0   28.0   29.0
  3.000000e+01   31.000000   32.0   33.0   34.0   35.0   36.0   37.0   38.0   39
  .0   40.0   41.0   42.0   43.0   44.0   45.0   46.0   47.0   48.0   49.0   50.
  0   51.0   52.0   53.0   54.0   55.0   56.0   57.0   58.0   59.0
           ...         ...    ...    ...    ...    ...    ...    ...    ...    .
  ..    ...    ...    ...    ...    ...    ...    ...    ...    ...    ...    ..
  .    ...    ...    ...    ...    ...    ...    ...    ...    ...
  2.970000e+03 2971.000000 2972.0 2973.0 2974.0 2975.0 2976.0 2977.0 2978.0 2979
  .0 2980.0 2981.0 2982.0 2983.0 2984.0 2985.0 2986.0 2987.0 2988.0 2989.0 2990.
  0 2991.0 2992.0 2993.0 2994.0 2995.0 2996.0 2997.0 2998.0 2999.0

For columns the syntax and behavior of
:func:`~astropy.table.table.Column.pprint` is the same except that there is no
``max_width`` keyword argument::

  >>> t['col3'].pprint(max_lines=8)
   col3
  ------
     3.0
    33.0
    63.0
     ...
  2943.0
  2973.0

pformat() method
''''''''''''''''

In order to get the formatted output for manipulation or writing to a file use
the Table :meth:`~astropy.table.Table.pformat` or Column
:func:`~astropy.table.table.Column.pformat` methods.  These behave just as for
:meth:`~astropy.table.Table.pprint` but return a list corresponding to each formatted line in the
:meth:`~astropy.table.Table.pprint` output.

  >>> lines = t['col3'].pformat(max_lines=8)
  >>> lines
  [' col3 ', '------', '   3.0', '  33.0', '  63.0', '   ...', '2943.0', '2973.0']

Multidimensional columns
''''''''''''''''''''''''

If a column has more than one dimension then each element of the column is
itself an array.  In the example below there are 3 rows, each of which is a
``2 x 2`` array.  The formatted output for such a column shows only the first
and last value of each row element and indicates the array dimensions in the
column name header::

  >>> from astropy.table import Table, Column
  >>> import numpy as np
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

In order to see all the data values for a multidimensional column use the
column representation.  This uses the standard `numpy` mechanism for printing
any array::

  >>> t['a']
  <Column name='a' unit=None format=None description=None>
  array([[[ 1,  2],
          [10, 20]],
         [[ 3,  4],
          [30, 40]],
         [[ 5,  6],
          [50, 60]]])
