.. |join| replace:: :func:`~astropy.table.join`

.. _table_operations:

Table Operations
****************

In this section we describe high-level operations that can be used to generate
a new table from one or more input tables. This includes:

=======================

.. list-table::
   :header-rows: 1
   :widths: 28 52 20

   * - Documentation
     - Description
     - Function
   * - `Grouped operations`_
     - Group tables and columns by keys
     - :func:`~astropy.table.Table.group_by`
   * - `Binning`_
     - Binning tables
     - :func:`~astropy.table.Table.group_by`
   * - `Stack vertically`_
     - Concatenate input tables along rows
     - :func:`~astropy.table.vstack`
   * - `Stack horizontally`_
     - Concatenate input tables along columns
     - :func:`~astropy.table.hstack`
   * - `Join`_
     - Database-style join of two tables
     - |join|
   * - `Unique rows`_
     - Unique table rows by keys
     - :func:`~astropy.table.unique`
   * - `Set difference`_
     - Set difference of two tables
     - :func:`~astropy.table.setdiff`
   * - `Table diff`_
     - Generic difference of two simple tables
     - :func:`~astropy.utils.diff.report_diff_values`


.. _grouped-operations:

Grouped Operations
------------------

.. EXAMPLE START: Grouped Operations in Tables

Sometimes in a table or table column there are natural groups within the dataset
for which it makes sense to compute some derived values. A minimal example is a
list of objects with photometry from various observing runs::

  >>> from astropy.table import Table
  >>> obs = Table.read("""name    obs_date    mag_b  mag_v
  ...                     M31     2012-01-02  17.0   17.5
  ...                     M31     2012-01-02  17.1   17.4
  ...                     M101    2012-01-02  15.1   13.5
  ...                     M82     2012-02-14  16.2   14.5
  ...                     M31     2012-02-14  16.9   17.3
  ...                     M82     2012-02-14  15.2   15.5
  ...                     M101    2012-02-14  15.0   13.6
  ...                     M82     2012-03-26  15.7   16.5
  ...                     M101    2012-03-26  15.1   13.5
  ...                     M101    2012-03-26  14.8   14.3
  ...                     """, format='ascii')
  >>> # Make sure magnitudes are printed with one digit after the decimal point
  >>> obs['mag_b'].info.format = '{:.1f}'
  >>> obs['mag_v'].info.format = '{:.1f}'

.. EXAMPLE END

Table Groups
^^^^^^^^^^^^

Now suppose we want the mean magnitudes for each object. We first group the data
by the ``name`` column with the :func:`~astropy.table.Table.group_by` method.
This returns a new table sorted by ``name`` which has a ``groups`` property
specifying the unique values of ``name`` and the corresponding table rows::

  >>> obs_by_name = obs.group_by('name')
  >>> print(obs_by_name)  # doctest: +SKIP
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
  M101 2012-01-02  15.1  13.5  << First group (index=0, key='M101')
  M101 2012-02-14  15.0  13.6
  M101 2012-03-26  15.1  13.5
  M101 2012-03-26  14.8  14.3
   M31 2012-01-02  17.0  17.5  << Second group (index=4, key='M31')
   M31 2012-01-02  17.1  17.4
   M31 2012-02-14  16.9  17.3
   M82 2012-02-14  16.2  14.5  << Third group (index=7, key='M83')
   M82 2012-02-14  15.2  15.5
   M82 2012-03-26  15.7  16.5
                               << End of groups (index=10)
  >>> print(obs_by_name.groups.keys)
  name
  ----
  M101
   M31
   M82
  >>> print(obs_by_name.groups.indices)
  [ 0  4  7 10]

The ``groups`` property is the portal to all grouped operations with tables and
columns. It defines how the table is grouped via an array of the unique row key
values and the indices of the group boundaries for those key values. The groups
here correspond to the row slices ``0:4``, ``4:7``, and ``7:10`` in the
``obs_by_name`` table.

The initial argument (``keys``) for the :func:`~astropy.table.Table.group_by`
function can take a number of input data types:

- Single string value with a table column name (as shown above)
- List of string values with table column names
- Another |Table| or |Column| with same length as table
- ``numpy`` structured array with same length as table
- ``numpy`` homogeneous array with same length as table

In all cases the corresponding row elements are considered as a :class:`tuple`
of values which form a key value that is used to sort the original table and
generate the required groups.

As an example, to get the average magnitudes for each object on each observing
night, we would first group the table on both ``name`` and ``obs_date`` as
follows::

  >>> print(obs.group_by(['name', 'obs_date']).groups.keys)
  name  obs_date
  ---- ----------
  M101 2012-01-02
  M101 2012-02-14
  M101 2012-03-26
   M31 2012-01-02
   M31 2012-02-14
   M82 2012-02-14
   M82 2012-03-26


Manipulating Groups
^^^^^^^^^^^^^^^^^^^

.. EXAMPLE START: Manipulating Groups in Tables

Once you have applied grouping to a table then you can access the individual
groups or subsets of groups. In all cases this returns a new grouped table.
For instance, to get the subtable which corresponds to the second group
(index=1) do::

  >>> print(obs_by_name.groups[1])
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
   M31 2012-01-02  17.0  17.5
   M31 2012-01-02  17.1  17.4
   M31 2012-02-14  16.9  17.3

To get the first and second groups together use a :class:`slice`::

  >>> groups01 = obs_by_name.groups[0:2]
  >>> print(groups01)
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
  M101 2012-01-02  15.1  13.5
  M101 2012-02-14  15.0  13.6
  M101 2012-03-26  15.1  13.5
  M101 2012-03-26  14.8  14.3
   M31 2012-01-02  17.0  17.5
   M31 2012-01-02  17.1  17.4
   M31 2012-02-14  16.9  17.3
  >>> print(groups01.groups.keys)
  name
  ----
  M101
   M31

You can also supply a ``numpy`` array of indices or a boolean mask to select
particular groups, for example::

  >>> mask = obs_by_name.groups.keys['name'] == 'M101'
  >>> print(obs_by_name.groups[mask])
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
  M101 2012-01-02  15.1  13.5
  M101 2012-02-14  15.0  13.6
  M101 2012-03-26  15.1  13.5
  M101 2012-03-26  14.8  14.3

You can iterate over the group subtables and corresponding keys with::

  >>> for key, group in zip(obs_by_name.groups.keys, obs_by_name.groups):
  ...     print(f'****** {key["name"]} *******')
  ...     print(group)
  ...     print('')
  ...
  ****** M101 *******
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
  M101 2012-01-02  15.1  13.5
  M101 2012-02-14  15.0  13.6
  M101 2012-03-26  15.1  13.5
  M101 2012-03-26  14.8  14.3
  ****** M31 *******
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
   M31 2012-01-02  17.0  17.5
   M31 2012-01-02  17.1  17.4
   M31 2012-02-14  16.9  17.3
  ****** M82 *******
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
   M82 2012-02-14  16.2  14.5
   M82 2012-02-14  15.2  15.5
   M82 2012-03-26  15.7  16.5

.. EXAMPLE END

Column Groups
^^^^^^^^^^^^^

Like |Table| objects, |Column| objects can also be grouped for subsequent
manipulation with grouped operations. This can apply both to columns within a
|Table| or bare |Column| objects.

As for |Table|, the grouping is generated with the
:func:`~astropy.table.Table.group_by` method. The difference here is that
there is no option of providing one or more column names since that
does not make sense for a |Column|.

Examples
~~~~~~~~

.. EXAMPLE START: Grouping Column Objects in Tables

To generate grouping in columns::

  >>> from astropy.table import Column
  >>> import numpy as np
  >>> c = Column([1, 2, 3, 4, 5, 6], name='a')
  >>> key_vals = np.array(['foo', 'bar', 'foo', 'foo', 'qux', 'qux'])
  >>> cg = c.group_by(key_vals)

  >>> for key, group in zip(cg.groups.keys, cg.groups):
  ...     print(f'****** {key} *******')
  ...     print(group)
  ...     print('')
  ...
  ****** bar *******
   a
  ---
    2
  ****** foo *******
   a
  ---
    1
    3
    4
  ****** qux *******
   a
  ---
    5
    6

.. EXAMPLE END

Aggregation
^^^^^^^^^^^

Aggregation is the process of applying a specified reduction function to the
values within each group for each non-key column. This function must accept a
|ndarray| as the first argument and return a single scalar value. Common
function examples are :func:`numpy.sum`, :func:`numpy.mean`, and
:func:`numpy.std`.

For the example grouped table ``obs_by_name`` from above, we compute the group
means with the :meth:`~astropy.table.groups.TableGroups.aggregate` method::

  >>> obs_mean = obs_by_name.groups.aggregate(np.mean)  # doctest: +SHOW_WARNINGS
  AstropyUserWarning: Cannot aggregate column 'obs_date' with type '<U10'
  >>> print(obs_mean)
  name mag_b mag_v
  ---- ----- -----
  M101  15.0  13.7
   M31  17.0  17.4
   M82  15.7  15.5

It seems the magnitude values were successfully averaged, but what about the
:class:`~astropy.utils.exceptions.AstropyUserWarning`? Since the ``obs_date``
column is a string-type array, the :func:`numpy.mean` function failed and
raised an exception.  Any time this happens
:meth:`~astropy.table.groups.TableGroups.aggregate` will issue a warning and
then drop that column from the output result. Note that the ``name`` column is
one of the ``keys`` used to determine the grouping so it is automatically
ignored from aggregation.

.. EXAMPLE START: Performing Aggregation on Grouped Tables

From a grouped table it is possible to select one or more columns on which
to perform the aggregation::

  >>> print(obs_by_name['mag_b'].groups.aggregate(np.mean))
  mag_b
  -----
   15.0
   17.0
   15.7

The order of the columns can be specified too::

  >>> print(obs_by_name['name', 'mag_v', 'mag_b'].groups.aggregate(np.mean))
  name mag_v mag_b
  ---- ----- -----
  M101  13.7  15.0
   M31  17.4  17.0
   M82  15.5  15.7


A single column of data can be aggregated as well::

  >>> c = Column([1, 2, 3, 4, 5, 6], name='a')
  >>> key_vals = np.array(['foo', 'bar', 'foo', 'foo', 'qux', 'qux'])
  >>> cg = c.group_by(key_vals)
  >>> cg_sums = cg.groups.aggregate(np.sum)
  >>> for key, cg_sum in zip(cg.groups.keys, cg_sums):
  ...     print(f'Sum for {key} = {cg_sum}')
  ...
  Sum for bar = 2
  Sum for foo = 8
  Sum for qux = 11

.. EXAMPLE END

If the specified function has a :meth:`numpy.ufunc.reduceat` method, this will
be called instead. This can improve the performance by a factor of 10 to 100
(or more) for large unmasked tables or columns with many relatively small
groups.  It also allows for the use of certain ``numpy`` functions which
normally take more than one input array but also work as reduction functions,
like `numpy.add`.  The ``numpy`` functions which should take advantage of using
:meth:`numpy.ufunc.reduceat` include:

- `numpy.add`
- `numpy.arctan2`
- `numpy.bitwise_and`
- `numpy.bitwise_or`
- `numpy.bitwise_xor`
- `numpy.copysign`
- `numpy.divide`
- `numpy.equal`
- `numpy.floor_divide`
- `numpy.fmax`
- `numpy.fmin`
- `numpy.fmod`
- `numpy.greater_equal`
- `numpy.greater`
- `numpy.hypot`
- `numpy.left_shift`
- `numpy.less_equal`
- `numpy.less`
- `numpy.logaddexp2`
- `numpy.logaddexp`
- `numpy.logical_and`
- `numpy.logical_or`
- `numpy.logical_xor`
- `numpy.maximum`
- `numpy.minimum`
- `numpy.mod`
- `numpy.multiply`
- `numpy.not_equal`
- `numpy.power`
- `numpy.remainder`
- `numpy.right_shift`
- `numpy.subtract`
- `numpy.true_divide`

In special cases, :func:`numpy.sum` and :func:`numpy.mean` are substituted with
their respective ``reduceat`` methods.

Filtering
^^^^^^^^^

Table groups can be filtered by means of the
:meth:`~astropy.table.groups.TableGroups.filter` method. This is done by
supplying a function which is called for each group. The function
which is passed to this method must accept two arguments:

- ``table`` : |Table| object
- ``key_colnames`` : list of columns in ``table`` used as keys for grouping

It must then return either `True` or `False`.

Example
~~~~~~~

.. EXAMPLE START: Filtering Table Groups

The following will select all table groups with only positive values in the non-
key columns::

  >>> def all_positive(table, key_colnames):
  ...     colnames = [name for name in table.colnames if name not in key_colnames]
  ...     for colname in colnames:
  ...         if np.any(table[colname] <= 0):
  ...             return False
  ...     return True

An example of using this function is::

  >>> t = Table.read(""" a   b    c
  ...                   -2  7.0   2
  ...                   -2  5.0   1
  ...                    1  3.0  -5
  ...                    1 -2.0  -6
  ...                    1  1.0   7
  ...                    0  4.0   4
  ...                    3  3.0   5
  ...                    3 -2.0   6
  ...                    3  1.0   7""", format='ascii')
  >>> tg = t.group_by('a')
  >>> t_positive = tg.groups.filter(all_positive)
  >>> for group in t_positive.groups:
  ...     print(group)
  ...     print('')
  ...
   a   b   c
  --- --- ---
   -2 7.0   2
   -2 5.0   1
  <BLANKLINE>
   a   b   c
  --- --- ---
    0 4.0   4

As can be seen only the groups with ``a == -2`` and ``a == 0`` have all
positive values in the non-key columns, so those are the ones that are selected.

Likewise a grouped column can be filtered with the
:meth:`~astropy.table.groups.ColumnGroups.filter`, method but in this case the
filtering function takes only a single argument which is the column group. It
still must return either `True` or `False`. For example::

  def all_positive(column):
      return np.all(column > 0)

.. EXAMPLE END

.. _table_binning:

Binning
-------

A common tool in analysis is to bin a table based on some reference value.
Examples:

- Photometry of a binary star in several bands taken over a
  span of time which should be binned by orbital phase.
- Reducing the sampling density for a table by combining
  100 rows at a time.
- Unevenly sampled historical data which should binned to
  four points per year.

All of these examples of binning a table can be accomplished using
`grouped operations`_. The examples in that section are focused on the
case of discrete key values such as the name of a source. In this
section we show a concise yet powerful way of applying grouped operations to
accomplish binning on key values such as time, phase, or row number.

The common theme in all of these cases is to convert the key value array into
a new float- or int-valued array whose values are identical for rows in the same
output bin.

Example
^^^^^^^

.. EXAMPLE START: Binning a Table using Grouped Operations

As an example, we generate a fake light curve::

  >>> year = np.linspace(2000.0, 2010.0, 200)  # 200 observations over 10 years
  >>> period = 1.811
  >>> y0 = 2005.2
  >>> mag = 14.0 + 1.2 * np.sin(2 * np.pi * (year - y0) / period)
  >>> phase = ((year - y0) / period) % 1.0
  >>> dat = Table([year, phase, mag], names=['year', 'phase', 'mag'])

Now we make an array that will be used for binning the data by 0.25 year
intervals::

  >>> year_bin = np.trunc(year / 0.25)

This has the property that all samples in each 0.25 year bin have the same
value of ``year_bin``. Think of ``year_bin`` as the bin number for ``year``.
Then do the binning by grouping and immediately aggregating with
:func:`numpy.mean`.

  >>> dat_grouped = dat.group_by(year_bin)
  >>> dat_binned = dat_grouped.groups.aggregate(np.mean)

We can plot the results with ``plt.plot(dat_binned['year'], dat_binned['mag'],
'.')``. Alternately, we could bin into 10 phase bins::

  >>> phase_bin = np.trunc(phase / 0.1)
  >>> dat_grouped = dat.group_by(phase_bin)
  >>> dat_binned = dat_grouped.groups.aggregate(np.mean)

This time, try plotting with ``plt.plot(dat_binned['phase'],
dat_binned['mag'])``.

.. EXAMPLE END

.. _stack-vertically:

Stack Vertically
----------------

The |Table| class supports stacking tables vertically with the
:func:`~astropy.table.vstack` function. This process is also commonly known as
concatenating or appending tables in the row direction. It corresponds roughly
to the :func:`numpy.vstack` function.

Examples
^^^^^^^^

.. EXAMPLE START: Stacking (or Concatenating) Tables Vertically

Suppose we have two tables of observations with several column names in
common::

  >>> from astropy.table import Table, vstack
  >>> obs1 = Table.read("""name    obs_date    mag_b  logLx
  ...                      M31     2012-01-02  17.0   42.5
  ...                      M82     2012-10-29  16.2   43.5
  ...                      M101    2012-10-31  15.1   44.5""", format='ascii')

  >>> obs2 = Table.read("""name    obs_date    logLx
  ...                      NGC3516 2011-11-11  42.1
  ...                      M31     1999-01-05  43.1
  ...                      M82     2012-10-30  45.0""", format='ascii')

Now we can stack these two tables::

  >>> print(vstack([obs1, obs2]))
    name   obs_date  mag_b logLx
  ------- ---------- ----- -----
      M31 2012-01-02  17.0  42.5
      M82 2012-10-29  16.2  43.5
     M101 2012-10-31  15.1  44.5
  NGC3516 2011-11-11    --  42.1
      M31 1999-01-05    --  43.1
      M82 2012-10-30    --  45.0

Notice that the ``obs2`` table is missing the ``mag_b`` column, so in the
stacked output table those values are marked as missing. This is the default
behavior and corresponds to ``join_type='outer'``. There are two other allowed
values for the ``join_type`` argument, ``'inner'`` and ``'exact'``::

  >>> print(vstack([obs1, obs2], join_type='inner'))
    name   obs_date  logLx
  ------- ---------- -----
      M31 2012-01-02  42.5
      M82 2012-10-29  43.5
     M101 2012-10-31  44.5
  NGC3516 2011-11-11  42.1
      M31 1999-01-05  43.1
      M82 2012-10-30  45.0

  >>> print(vstack([obs1, obs2], join_type='exact'))  # doctest: +IGNORE_EXCEPTION_DETAIL
  Traceback (most recent call last):
    ...
  TableMergeError: Inconsistent columns in input arrays (use 'inner'
  or 'outer' join_type to allow non-matching columns)

In the case of ``join_type='inner'``, only the common columns (the intersection)
are present in the output table. When ``join_type='exact'`` is specified, then
:func:`~astropy.table.vstack` requires that all of the input tables have
exactly the same column names.

More than two tables can be stacked by supplying a longer list of tables::

  >>> obs3 = Table.read("""name    obs_date    mag_b  logLx
  ...                      M45     2012-02-03  15.0   40.5""", format='ascii')
  >>> print(vstack([obs1, obs2, obs3]))
    name   obs_date  mag_b logLx
  ------- ---------- ----- -----
      M31 2012-01-02  17.0  42.5
      M82 2012-10-29  16.2  43.5
     M101 2012-10-31  15.1  44.5
  NGC3516 2011-11-11    --  42.1
      M31 1999-01-05    --  43.1
      M82 2012-10-30    --  45.0
      M45 2012-02-03  15.0  40.5

See also the sections on `Merging metadata`_ and `Merging column attributes`_
for details on how these characteristics of the input tables are merged in the
single output table. Note also that you can use a single table |Row| instead of
a full table as one of the inputs.

.. EXAMPLE END

.. _stack-horizontally:

Stack Horizontally
------------------

The |Table| class supports stacking tables horizontally (in the column-wise
direction) with the :func:`~astropy.table.hstack` function. It corresponds
roughly to the :func:`numpy.hstack` function.

Examples
^^^^^^^^

.. EXAMPLE START: Stacking (or Concatenating) Tables Horizontally

Suppose we have the following two tables::

  >>> from astropy.table import Table, hstack
  >>> t1 = Table.read("""a   b    c
  ...                    1   foo  1.4
  ...                    2   bar  2.1
  ...                    3   baz  2.8""", format='ascii')
  >>> t2 = Table.read("""d     e
  ...                    ham   eggs
  ...                    spam  toast""", format='ascii')

Now we can stack these two tables horizontally::

  >>> print(hstack([t1, t2]))
   a   b   c   d     e
  --- --- --- ---- -----
    1 foo 1.4  ham  eggs
    2 bar 2.1 spam toast
    3 baz 2.8   --    --

As with :func:`~astropy.table.vstack`, there is an optional ``join_type``
argument that can take values ``'inner'``, ``'exact'``, and ``'outer'``. The
default is ``'outer'``, which effectively takes the union of available rows and
masks out any missing values. This is illustrated in the example above. The
other options give the intersection of rows, where ``'exact'`` requires that
all tables have exactly the same number of rows::

  >>> print(hstack([t1, t2], join_type='inner'))
   a   b   c   d     e
  --- --- --- ---- -----
    1 foo 1.4  ham  eggs
    2 bar 2.1 spam toast

  >>> print(hstack([t1, t2], join_type='exact'))  # doctest: +IGNORE_EXCEPTION_DETAIL
  Traceback (most recent call last):
    ...
  TableMergeError: Inconsistent number of rows in input arrays (use 'inner' or
  'outer' join_type to allow non-matching rows)

More than two tables can be stacked by supplying a longer list of tables. The
example below also illustrates the behavior when there is a conflict in the
input column names (see the section on `Column renaming`_ for details)::

  >>> t3 = Table.read("""a    b
  ...                    M45  2012-02-03""", format='ascii')
  >>> print(hstack([t1, t2, t3]))
  a_1 b_1  c   d     e   a_3    b_3
  --- --- --- ---- ----- --- ----------
    1 foo 1.4  ham  eggs M45 2012-02-03
    2 bar 2.1 spam toast  --         --
    3 baz 2.8   --    --  --         --

The metadata from the input tables is merged by the process described in the
`Merging metadata`_ section. Note also that you can use a single table |Row|
instead of a full table as one of the inputs.

.. EXAMPLE END

.. _stack-depthwise:

Stack Depth-Wise
----------------

The |Table| class supports stacking columns within tables depth-wise using the
:func:`~astropy.table.dstack` function. It corresponds roughly to running the
:func:`numpy.dstack` function on the individual columns matched by name.

Examples
^^^^^^^^

.. EXAMPLE START: Stacking (or Concatenating) Tables Depth-Wise

Suppose we have tables of data for sources giving information on the enclosed
source counts for different PSF fractions::

  >>> from astropy.table import Table, dstack
  >>> src1 = Table.read("""psf_frac  counts
  ...                      0.10        45
  ...                      0.50        90
  ...                      0.90       120
  ...                      """, format='ascii')

  >>> src2 = Table.read("""psf_frac  counts
  ...                      0.10       200
  ...                      0.50       300
  ...                      0.90       350
  ...                      """, format='ascii')

Now we can stack these two tables depth-wise to get a single table with the
characteristics of both sources::

  >>> srcs = dstack([src1, src2])
  >>> print(srcs)
   psf_frac    counts
  float64[2]  int64[2]
  ---------- ----------
  0.1 .. 0.1  45 .. 200
  0.5 .. 0.5  90 .. 300
  0.9 .. 0.9 120 .. 350

In this case the counts for the first source are accessible as
``srcs['counts'][:, 0]``, and likewise the second source counts are
``srcs['counts'][:, 1]``.

For this function the length of all input tables must be the same. This
function can accept ``join_type`` and ``metadata_conflicts`` just like the
:func:`~astropy.table.vstack` function. The ``join_type`` argument controls how
to handle mismatches in the columns of the input table.

See also the sections on `Merging metadata`_ and `Merging column attributes`_
for details on how these characteristics of the input tables are merged in the
single output table. Note also that you can use a single table |Row| instead of
a full table as one of the inputs.

.. EXAMPLE END

.. _table-join:

Join
----

The |Table| class supports the `database join
<https://en.wikipedia.org/wiki/Join_(SQL)>`_ operation. This provides a flexible
and powerful way to combine tables based on the values in one or more key
columns.

Examples
^^^^^^^^

.. EXAMPLE START: Combining Tables using the Database Join Operation

Suppose we have two tables of observations, the first with B and V magnitudes
and the second with X-ray luminosities of an overlapping (but not identical)
sample::

  >>> from astropy.table import Table, join
  >>> optical = Table.read("""name    obs_date    mag_b  mag_v
  ...                         M31     2012-01-02  17.0   16.0
  ...                         M82     2012-10-29  16.2   15.2
  ...                         M101    2012-10-31  15.1   15.5""", format='ascii')
  >>> xray = Table.read("""   name    obs_date    logLx
  ...                         NGC3516 2011-11-11  42.1
  ...                         M31     1999-01-05  43.1
  ...                         M82     2012-10-29  45.0""", format='ascii')

The |join| method allows you to merge these two tables into a single table based
on matching values in the "key columns". By default, the key columns are the set
of columns that are common to both tables. In this case the key columns are
``name`` and ``obs_date``. We can find all of the observations of the same
object on the same date as follows::

  >>> opt_xray = join(optical, xray)
  >>> print(opt_xray)
  name  obs_date  mag_b mag_v logLx
  ---- ---------- ----- ----- -----
   M82 2012-10-29  16.2  15.2  45.0

We can perform the match by ``name`` only by providing the ``keys`` argument,
which can be either a single column name or a list of column names::

  >>> print(join(optical, xray, keys='name'))
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
   M31 2012-01-02  17.0  16.0 1999-01-05  43.1
   M82 2012-10-29  16.2  15.2 2012-10-29  45.0

This output table has all of the observations that have both optical and X-ray
data for an object (M31 and M82). Notice that since the ``obs_date`` column
occurs in both tables, it has been split into two columns, ``obs_date_1`` and
``obs_date_2``. The values are taken from the "left" (``optical``) and "right"
(``xray``) tables, respectively.

.. EXAMPLE END

Different Join Options
^^^^^^^^^^^^^^^^^^^^^^

The table joins so far are known as "inner" joins and represent the strict
intersection of the two tables on the key columns.

.. EXAMPLE START: Table Join Options

If you want to make a new table which has *every* row from the left table and
includes matching values from the right table when available, this is known as a
left join::

  >>> print(join(optical, xray, join_type='left'))
  name  obs_date  mag_b mag_v logLx
  ---- ---------- ----- ----- -----
  M101 2012-10-31  15.1  15.5    --
   M31 2012-01-02  17.0  16.0    --
   M82 2012-10-29  16.2  15.2  45.0

Two of the observations do not have X-ray data, as indicated by the ``--`` in
the table. You might be surprised that there is no X-ray data for M31 in the
output. Remember that the default matching key includes both ``name`` and
``obs_date``. Specifying the key as only the ``name`` column gives::

  >>> print(join(optical, xray, join_type='left', keys='name'))
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
  M101 2012-10-31  15.1  15.5         --    --
   M31 2012-01-02  17.0  16.0 1999-01-05  43.1
   M82 2012-10-29  16.2  15.2 2012-10-29  45.0

Likewise you can construct a new table with every row of the right table and
matching left values (when available) using ``join_type='right'``.

To make a table with the union of rows from both tables do an "outer" join::

  >>> print(join(optical, xray, join_type='outer'))
    name   obs_date  mag_b mag_v logLx
  ------- ---------- ----- ----- -----
     M101 2012-10-31  15.1  15.5    --
      M31 1999-01-05    --    --  43.1
      M31 2012-01-02  17.0  16.0    --
      M82 2012-10-29  16.2  15.2  45.0
  NGC3516 2011-11-11    --    --  42.1

In all the above cases the output join table will be sorted by the key
column(s) and in general will not preserve the row order of the input tables.

Finally, you can do a "Cartesian" join, which is the Cartesian product of all
available rows. In this case there are no key columns (and supplying the
``keys`` argument is an error)::

  >>> print(join(optical, xray, join_type='cartesian'))
  name_1 obs_date_1 mag_b mag_v  name_2 obs_date_2 logLx
  ------ ---------- ----- ----- ------- ---------- -----
     M31 2012-01-02  17.0  16.0 NGC3516 2011-11-11  42.1
     M31 2012-01-02  17.0  16.0     M31 1999-01-05  43.1
     M31 2012-01-02  17.0  16.0     M82 2012-10-29  45.0
     M82 2012-10-29  16.2  15.2 NGC3516 2011-11-11  42.1
     M82 2012-10-29  16.2  15.2     M31 1999-01-05  43.1
     M82 2012-10-29  16.2  15.2     M82 2012-10-29  45.0
    M101 2012-10-31  15.1  15.5 NGC3516 2011-11-11  42.1
    M101 2012-10-31  15.1  15.5     M31 1999-01-05  43.1
    M101 2012-10-31  15.1  15.5     M82 2012-10-29  45.0

.. EXAMPLE END

Non-Identical Key Column Names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. EXAMPLE START: Joining Tables with Unique Key Column Names

To use the |join| function with non-identical key column names, use the
``keys_left`` and ``keys_right`` arguments. In the following example one table
has a ``'name'`` column while the other has an ``'obj_id'`` column::

  >>> optical = Table.read("""name    obs_date    mag_b  mag_v
  ...                         M31     2012-01-02  17.0   16.0
  ...                         M82     2012-10-29  16.2   15.2
  ...                         M101    2012-10-31  15.1   15.5""", format='ascii')
  >>> xray_1 = Table.read("""obj_id    obs_date    logLx
  ...                        NGC3516 2011-11-11  42.1
  ...                        M31     1999-01-05  43.1
  ...                        M82     2012-10-29  45.0""", format='ascii')

In order to perform a match based on the names of the objects, do the
following::

  >>> print(join(optical, xray_1, keys_left='name', keys_right='obj_id'))
  name obs_date_1 mag_b mag_v obj_id obs_date_2 logLx
  ---- ---------- ----- ----- ------ ---------- -----
   M31 2012-01-02  17.0  16.0    M31 1999-01-05  43.1
   M82 2012-10-29  16.2  15.2    M82 2012-10-29  45.0

The ``keys_left`` and ``keys_right`` arguments can also take a list of column
names or even a list of column-like objects. The latter case allows specifying
the matching key column values independent of the tables being joined.

.. EXAMPLE END

Identical Key Values
^^^^^^^^^^^^^^^^^^^^

.. EXAMPLE START: Joining Tables with Identical Key Values

The |Table| join operation works even if there are multiple rows with identical
key values. For example, the following tables have multiple rows for the column
``'key'``::

  >>> from astropy.table import Table, join
  >>> left = Table([[0, 1, 1, 2], ['L1', 'L2', 'L3', 'L4']], names=('key', 'L'))
  >>> right = Table([[1, 1, 2, 4], ['R1', 'R2', 'R3', 'R4']], names=('key', 'R'))
  >>> print(left)
  key  L
  --- ---
    0  L1
    1  L2
    1  L3
    2  L4
  >>> print(right)
  key  R
  --- ---
    1  R1
    1  R2
    2  R3
    4  R4

Doing an outer join on these tables shows that what is really happening is a
`Cartesian product <https://en.wikipedia.org/wiki/Cartesian_product>`_. For
each matching key, every combination of the left and right tables is
represented. When there is no match in either the left or right table, the
corresponding column values are designated as missing::

  >>> print(join(left, right, join_type='outer'))
  key  L   R
  --- --- ---
    0  L1  --
    1  L2  R1
    1  L2  R2
    1  L3  R1
    1  L3  R2
    2  L4  R3
    4  --  R4

An inner join is the same but only returns rows where there is a key match in
both the left and right tables::

  >>> print(join(left, right, join_type='inner'))
  key  L   R
  --- --- ---
    1  L2  R1
    1  L2  R2
    1  L3  R1
    1  L3  R2
    2  L4  R3

Conflicts in the input table names are handled by the process described in the
section on `Column renaming`_. See also the sections on `Merging metadata`_ and
`Merging column attributes`_ for details on how these characteristics of the
input tables are merged in the single output table.

.. EXAMPLE END

Merging Details
---------------

When combining two or more tables there is the need to merge certain
characteristics in the inputs and potentially resolve conflicts. This
section describes the process.

Column Renaming
^^^^^^^^^^^^^^^

In cases where the input tables have conflicting column names, there
is a mechanism to generate unique output column names. There are two
keyword arguments that control the renaming behavior:

``table_names``
    List of strings that provide names for the tables being joined.
    By default this is ``['1', '2', ...]``, where the numbers correspond to
    the input tables.

``uniq_col_name``
    String format specifier with a default value of ``'{col_name}_{table_name}'``.

This is best understood by example using the ``optical`` and ``xray`` tables
in the |join| example defined previously::

  >>> print(join(optical, xray, keys='name',
  ...            table_names=['OPTICAL', 'XRAY'],
  ...            uniq_col_name='{table_name}_{col_name}'))
  name OPTICAL_obs_date mag_b mag_v XRAY_obs_date logLx
  ---- ---------------- ----- ----- ------------- -----
   M31       2012-01-02  17.0  16.0    1999-01-05  43.1
   M82       2012-10-29  16.2  15.2    2012-10-29  45.0

.. _merging_metadata:

Merging Metadata
^^^^^^^^^^^^^^^^

|Table| objects can have associated metadata:

- ``Table.meta``: table-level metadata as an ordered dictionary
- ``Column.meta``: per-column metadata as an ordered dictionary

The table operations described here handle the task of merging the metadata in
the input tables into a single output structure. Because the metadata can be
arbitrarily complex there is no unique way to do the merge. The current
implementation uses a recursive algorithm with four rules:

- :class:`dict` elements are merged by keys.
- Conflicting :class:`list` or :class:`tuple` elements are concatenated.
- Conflicting :class:`dict` elements are merged by recursively calling the
  merge function.
- Conflicting elements that are not :class:`list`, :class:`tuple`, or
  :class:`dict` will follow the following rules:

    - If both metadata values are identical, the output is set to this value.
    - If one of the conflicting metadata values is `None`, the other value is
      picked.
    - If both metadata values are different and neither is `None`, the one for
      the last table in the list is picked.

By default, a warning is emitted in the last case (both metadata values are not
`None`). The warning can be silenced or made into an exception using the
``metadata_conflicts`` argument to :func:`~astropy.table.hstack`,
:func:`~astropy.table.vstack`, or
:func:`~astropy.table.join`. The ``metadata_conflicts`` option can be set to:

- ``'silent'`` – no warning is emitted, the value for the last table is silently
  picked.
- ``'warn'`` – a warning is emitted, the value for the last table is picked.
- ``'error'`` – an exception is raised.

The default strategies for merging metadata can be augmented or customized by
defining subclasses of the `~astropy.utils.metadata.MergeStrategy` base class.
In most cases you will also use
:func:`~astropy.utils.metadata.enable_merge_strategies` for enabling the custom
strategies. The linked documentation strings provide details.

Merging Column Attributes
^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the table and column ``meta`` attributes, the column attributes
``unit``, ``format``, and ``description`` are merged by going through the input
tables in order and taking the last value which is defined (i.e., is not
`None`).

Example
~~~~~~~

.. EXAMPLE START: Merging Column Attributes in a Table

To merge column attributes ``unit``, ``format``, or ``description``::

  >>> from astropy.table import Column, Table, vstack
  >>> col1 = Column([1], name='a')
  >>> col2 = Column([2], name='a', unit='cm')
  >>> col3 = Column([3], name='a', unit='m')
  >>> t1 = Table([col1])
  >>> t2 = Table([col2])
  >>> t3 = Table([col3])
  >>> out = vstack([t1, t2, t3])  # doctest: +SHOW_WARNINGS
  MergeConflictWarning: In merged column 'a' the 'unit' attribute does
  not match (cm != m).  Using m for merged output
  >>> out['a'].unit
  Unit("m")

The rules for merging are the same as for `Merging metadata`_, and the
``metadata_conflicts`` option also controls the merging of column attributes.

.. EXAMPLE END

.. _astropy-table-join-functions:

Joining Coordinates and Custom Join Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Source catalogs that have |SkyCoord| coordinate columns can be joined using
cross-matching of the coordinates with a specified distance threshold. This is
a special case of a more general problem of "fuzzy" matching of key column
values, where instead of an exact match we require only an approximate match.
This is supported using the ``join_funcs`` argument.

.. warning::

   The coordinate and distance table joins discussed in this section are most
   applicable in the case where the relevant entries in at least one of the
   tables are all separated from one another by more than twice the join
   distance. If this is not satisfied then the join results may be unexpected.

   This is a consequence of the algorithm which effectively finds clusters of
   nearby points (an "equivalence class") and assigns a unique cluster
   identifier to each entry in both tables. This assumes the join matching
   function is a transitive relation where ``join_func(A, B)`` and
   ``join_func(B, C)`` implies ``join_func(A, C)``. With multiple matches on
   both left and right sides it is possible for the cluster of points having a
   single cluster identifier to expand in size beyond the distance threshold.

   Users should be especially aware of this issue if additional join keys
   are provided beyond the ``join_funcs``. The code does not do a "pre-join"
   on the other keys, so the possibility of having overlaps within the distance
   in both tables is higher.

Example
~~~~~~~

.. EXAMPLE START: Joining a Table on Coordinates

To join two tables on a |SkyCoord| key column we use the ``join_funcs`` keyword
to supply a :class:`dict` of functions that specify how to match a particular
key column by name. In the example below we are joining on the ``sc`` column,
so we provide the following argument::

  join_funcs={'sc': join_skycoord(0.2 * u.deg)}

This tells |join| to match the ``sc`` key column using the join function
:func:`~astropy.table.join_skycoord` with a matching distance threshold of 0.2
deg. Under the hood this calls
:meth:`~astropy.coordinates.SkyCoord.search_around_sky` or
:meth:`~astropy.coordinates.SkyCoord.search_around_3d` to do the
cross-matching. The default is to use
:meth:`~astropy.coordinates.SkyCoord.search_around_sky` (angle) matching, but
:meth:`~astropy.coordinates.SkyCoord.search_around_3d` (length or
dimensionless) is also available. This is specified using the ``distance_func``
argument of :func:`~astropy.table.join_skycoord`, which can also be a function
that matches the input and output API of
:meth:`~astropy.coordinates.SkyCoord.search_around_sky`.

Now we show the whole process:

..  doctest-requires:: scipy

  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u
  >>> from astropy.table import Table, join, join_skycoord

..  doctest-requires:: scipy

  >>> sc1 = SkyCoord([0, 1, 1.1, 2], [0, 0, 0, 0], unit='deg')
  >>> sc2 = SkyCoord([1.05, 0.5, 2.1], [0, 0, 0], unit='deg')

..  doctest-requires:: scipy

  >>> t1 = Table([sc1, [0, 1, 2, 3]], names=['sc', 'idx'])
  >>> t2 = Table([sc2, [0, 1, 2]], names=['sc', 'idx'])

..  doctest-requires:: scipy

  >>> t12 = join(t1, t2, keys='sc', join_funcs={'sc': join_skycoord(0.2 * u.deg)})
  >>> print(t12)
  sc_id   sc_1  idx_1   sc_2   idx_2
        deg,deg       deg,deg
  ----- ------- ----- -------- -----
      1 1.0,0.0     1 1.05,0.0     0
      1 1.1,0.0     2 1.05,0.0     0
      2 2.0,0.0     3  2.1,0.0     2

The joined table has matched the sources within 0.2 deg and created a new
column ``sc_id`` with a unique identifier for each source.

.. EXAMPLE END

You might be wondering what is happening in the join function defined above,
especially if you are interested in defining your own such function. This could
be done in order to allow fuzzy word matching of tables, for example joining
tables of people by name where the names do not always match exactly.

The first thing to note here is that the :func:`~astropy.table.join_skycoord`
function actually returns a function itself. This allows specifying a variable
match distance via a function enclosure. The requirement of the join function
is that it accepts two arguments corresponding to the two key columns, and
returns a tuple of ``(ids1, ids2)``. These identifiers correspond to the
identification of each column entry with a unique matched source.

..  doctest-requires:: scipy

    >>> join_func = join_skycoord(0.2 * u.deg)
    >>> join_func(sc1, sc2)  # Associate each coordinate with unique source ID
    (array([3, 1, 1, 2]), array([1, 4, 2]))

If you would like to write your own fuzzy matching function, we suggest starting
from the source code for :func:`~astropy.table.join_skycoord` or
:func:`~astropy.table.join_distance`.

Join on Distance
~~~~~~~~~~~~~~~~

The example above focused on joining on a |SkyCoord|, but you can also join on
a generic distance between column values using the
:func:`~astropy.table.join_distance` join function. This can apply to 1D or 2D
(vector) columns. This will look very similar to the coordinates example, but
here there is a bit more flexibility. The matching is done using
:class:`scipy.spatial.cKDTree` and
:meth:`scipy.spatial.cKDTree.query_ball_tree`, and the behavior of these can be
controlled via the ``kdtree_args`` and ``query_args`` arguments, respectively.

.. _unique-rows:

Unique Rows
-----------

Sometimes it makes sense to use only rows with unique key columns or even
fully unique rows from a table. This can be done using the above described
:meth:`~astropy.table.Table.group_by` method and ``groups`` attribute, or with
the :func:`~astropy.table.unique` convenience function. The
:func:`~astropy.table.unique` function returns a sorted table containing the
first row for each unique ``keys`` column value. If no ``keys`` is provided, it
returns a sorted table containing all of the fully unique rows.

Example
^^^^^^^

.. EXAMPLE START: Grouping Unique Rows in Tables

An example of a situation where you might want to use rows with unique key
columns is a list of objects with photometry from various observing
runs. Using ``'name'`` as the only ``keys``, it returns with the first
occurrence of each of the three targets::

  >>> from astropy import table
  >>> obs = table.Table.read("""name    obs_date    mag_b  mag_v
  ...                           M31     2012-01-02  17.0   17.5
  ...                           M82     2012-02-14  16.2   14.5
  ...                           M101    2012-01-02  15.1   13.5
  ...                           M31     2012-01-02  17.1   17.4
  ...                           M101    2012-01-02  15.1   13.5
  ...                           M82     2012-02-14  16.2   14.5
  ...                           M31     2012-02-14  16.9   17.3
  ...                           M82     2012-02-14  15.2   15.5
  ...                           M101    2012-02-14  15.0   13.6
  ...                           M82     2012-03-26  15.7   16.5
  ...                           M101    2012-03-26  15.1   13.5
  ...                           M101    2012-03-26  14.8   14.3
  ...                           """, format='ascii')
  >>> unique_by_name = table.unique(obs, keys='name')
  >>> print(unique_by_name)
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
  M101 2012-01-02  15.1  13.5
   M31 2012-01-02  17.0  17.5
   M82 2012-02-14  16.2  14.5

Using multiple columns as ``keys``::

  >>> unique_by_name_date = table.unique(obs, keys=['name', 'obs_date'])
  >>> print(unique_by_name_date)
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
  M101 2012-01-02  15.1  13.5
  M101 2012-02-14  15.0  13.6
  M101 2012-03-26  15.1  13.5
   M31 2012-01-02  17.0  17.5
   M31 2012-02-14  16.9  17.3
   M82 2012-02-14  16.2  14.5
   M82 2012-03-26  15.7  16.5

.. EXAMPLE END

.. _set-difference:

Set Difference
--------------

A set difference will tell you the elements that are contained in the first set
but not in the other. This concept can be applied to rows of a table by using
the :func:`~astropy.table.setdiff` function. You provide the function with two
input tables and it will return all rows in the first table which do not occur
in the second table.

The optional ``keys`` parameter specifies the names of columns that are used to
match table rows. This can be a subset of the full list of columns, but both
the first and second tables must contain all columns specified by ``keys``.
If not provided, then ``keys`` defaults to all column names in the first table.

If no different rows are found, the :func:`~astropy.table.setdiff` function
will return an empty table.

Example
^^^^^^^

.. EXAMPLE START: Using Set Difference in Tables

The example below illustrates finding the set difference of two observation
lists using a common subset of the columns in two tables.::

  >>> from astropy.table import Table, setdiff
  >>> cat_1 = Table.read("""name    obs_date    mag_b  mag_v
  ...                       M31     2012-01-02  17.0   16.0
  ...                       M82     2012-10-29  16.2   15.2
  ...                       M101    2012-10-31  15.1   15.5""", format='ascii')
  >>> cat_2 = Table.read("""   name    obs_date    logLx
  ...                          NGC3516 2011-11-11  42.1
  ...                          M31     2012-01-02  43.1
  ...                          M82     2012-10-29  45.0""", format='ascii')
  >>> sdiff = setdiff(cat_1, cat_2, keys=['name', 'obs_date'])
  >>> print(sdiff)
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
  M101 2012-10-31  15.1  15.5

In this example there is a column in the first table that is not
present in the second table, so the ``keys`` parameter must be used to specify
the desired column names.

.. EXAMPLE END

.. _table-diff:

Table Diff
----------

To compare two tables, you can use
:func:`~astropy.utils.diff.report_diff_values`, which would produce a report
identical to :ref:`FITS diff <io-fits-differs>`.

Example
^^^^^^^

.. EXAMPLE START: Using Table Diff to Compare Tables

The example below illustrates finding the difference between two tables::

  >>> from astropy.table import Table
  >>> from astropy.utils.diff import report_diff_values
  >>> import sys
  >>> cat_1 = Table.read("""name    obs_date    mag_b  mag_v
  ...                       M31     2012-01-02  17.0   16.0
  ...                       M82     2012-10-29  16.2   15.2
  ...                       M101    2012-10-31  15.1   15.5""", format='ascii')
  >>> cat_2 = Table.read("""name    obs_date    mag_b  mag_v
  ...                       M31     2012-01-02  17.0   16.5
  ...                       M82     2012-10-29  16.2   15.2
  ...                       M101    2012-10-30  15.1   15.5
  ...                       NEW     2018-05-08   nan    9.0""", format='ascii')
  >>> identical = report_diff_values(cat_1, cat_2, fileobj=sys.stdout)
       name  obs_date  mag_b mag_v
       ---- ---------- ----- -----
    a>  M31 2012-01-02  17.0  16.0
     ?                           ^
    b>  M31 2012-01-02  17.0  16.5
     ?                           ^
        M82 2012-10-29  16.2  15.2
    a> M101 2012-10-31  15.1  15.5
     ?               ^
    b> M101 2012-10-30  15.1  15.5
     ?               ^
    b>  NEW 2018-05-08   nan   9.0
  >>> identical
  False

.. EXAMPLE END
