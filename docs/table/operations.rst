.. include:: references.txt
.. |join| replace:: :func:`~astropy.table.join`

.. _table_operations:

Table operations
-----------------

In this section we describe higher-level operations that can be used to generate a new
table from one or more input tables.  This includes:

=======================

.. list-table::
   :header-rows: 1
   :widths: 28 52 20

   * - Documentation
     - Description
     - Function
   * - `Grouped operations`_
     - Group tables and columns by keys
     - `~astropy.table.Table.group_by`
   * - `Stack vertically`_
     - Concatenate input tables along rows
     - `~astropy.table.vstack`
   * - `Stack horizontally`_
     - Concatenate input tables along columns
     - `~astropy.table.hstack`
   * - `Join`_
     - Database-style join of two tables
     - `~astropy.table.join`
   * - `Unique rows`_
     - Unique table rows by keys
     - `~astropy.table.unique`


.. _grouped-operations:

Grouped operations
^^^^^^^^^^^^^^^^^^

Sometimes in a table or table column there are natural groups within the dataset for which
it makes sense to compute some derived values.  A simple example is a list of objects with
photometry from various observing runs::

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

Table groups
~~~~~~~~~~~~~~

Now suppose we want the mean magnitudes for each object.  We first group the data by the
``name`` column with the :func:`~astropy.table.Table.group_by` method.  This returns
a new table sorted by ``name`` which has a ``groups`` property specifying the unique
values of ``name`` and the corresponding table rows::

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

The ``groups`` property is the portal to all grouped operations with tables and columns.
It defines how the table is grouped via an array of the unique row key values and the
indices of the group boundaries for those key values.  The groups here correspond to the
row slices ``0:4``, ``4:7``, and ``7:10`` in the ``obs_by_name`` table.

The initial argument (``keys``) for the `~astropy.table.Table.group_by` function
can take a number of input data types:

- Single string value with a table column name (as shown above)
- List of string values with table column names
- Another |Table| or |Column| with same length as table
- Numpy structured array with same length as table
- Numpy homogeneous array with same length as table

In all cases the corresponding row elements are considered as a tuple of values which
form a key value that is used to sort the original table and generate
the required groups.

As an example, to get the average magnitudes for each object on each observing
night, we would first group the table on both ``name`` and ``obs_date`` as follows::

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


Manipulating groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have applied grouping to a table then you can easily access the individual
groups or subsets of groups.  In all cases this returns a new grouped table.
For instance to get the sub-table which corresponds to the second group (index=1)
do::

  >>> print(obs_by_name.groups[1])
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
   M31 2012-01-02  17.0  17.5
   M31 2012-01-02  17.1  17.4
   M31 2012-02-14  16.9  17.3

To get the first and second groups together use a slice::

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

You can also supply a numpy array of indices or a boolean mask to select particular
groups, e.g.::

  >>> mask = obs_by_name.groups.keys['name'] == 'M101'
  >>> print(obs_by_name.groups[mask])
  name  obs_date  mag_b mag_v
  ---- ---------- ----- -----
  M101 2012-01-02  15.1  13.5
  M101 2012-02-14  15.0  13.6
  M101 2012-03-26  15.1  13.5
  M101 2012-03-26  14.8  14.3

One can iterate over the group sub-tables and corresponding keys with::

  >>> from itertools import izip
  >>> for key, group in izip(obs_by_name.groups.keys, obs_by_name.groups):
  ...     print('****** {0} *******'.format(key['name']))
  ...     print(group)
  ...     print
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

Column Groups
~~~~~~~~~~~~~~

Like |Table| objects, |Column| objects can also be grouped for subsequent
manipulation with grouped operations.  This can apply both to columns within a
|Table| or bare |Column| objects.

As for |Table|, the grouping is generated with the
`~astropy.table.Table.group_by` method.  The difference here is that
there is no option of providing one or more column names since that
doesn't make sense for a |Column|.

Examples::

  >>> from astropy.table import Column
  >>> import numpy as np
  >>> c = Column([1, 2, 3, 4, 5, 6], name='a')
  >>> key_vals = np.array(['foo', 'bar', 'foo', 'foo', 'qux', 'qux'])
  >>> cg = c.group_by(key_vals)

  >>> for key, group in izip(cg.groups.keys, cg.groups):
  ...     print('****** {0} *******'.format(key))
  ...     print(group)
  ...     print
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


Aggregation
~~~~~~~~~~~~~~

Aggregation is the process of applying a
specified reduction function to the values within each group for each
non-key column.  This function must accept a numpy array as the first
argument and return a single scalar value.  Common function examples are
`numpy.sum`, `numpy.mean`, and `numpy.std`.

For the example grouped table ``obs_by_name`` from above we compute the group means with
the `~astropy.table.groups.TableGroups.aggregate` method::

  >>> obs_mean = obs_by_name.groups.aggregate(np.mean)  # doctest: +SKIP
  WARNING: Cannot aggregate column 'obs_date' [astropy.table.groups]
  >>> print(obs_mean)  # doctest: +SKIP
  name mag_b mag_v
  ---- ----- ------
  M101  15.0 13.725
   M31  17.0   17.4
   M82  15.7   15.5

It seems the magnitude values were successfully averaged, but what
about the WARNING?  Since the ``obs_date`` column is a string-type
array, the `numpy.mean` function failed and raised an exception.
Any time this happens then `~astropy.table.groups.TableGroups.aggregate`
will issue a warning and then
drop that column from the output result.  Note that the ``name``
column is one of the ``keys`` used to determine the grouping so
it is automatically ignored from aggregation.

From a grouped table it is possible to select one or more columns on which
to perform the aggregation::

  >>> print(obs_by_name['mag_b'].groups.aggregate(np.mean))
  mag_b
  -----
   15.0
   17.0
   15.7

  >>> print(obs_by_name['name', 'mag_v', 'mag_b'].groups.aggregate(np.mean))
  name mag_v  mag_b
  ---- ------ -----
  M101 13.725  15.0
   M31   17.4  17.0
   M82   15.5  15.7

A single column of data can be aggregated as well::

  >>> c = Column([1, 2, 3, 4, 5, 6], name='a')
  >>> key_vals = np.array(['foo', 'bar', 'foo', 'foo', 'qux', 'qux'])
  >>> cg = c.group_by(key_vals)
  >>> cg_sums = cg.groups.aggregate(np.sum)
  >>> for key, cg_sum in izip(cg.groups.keys, cg_sums):
  ...     print('Sum for {0} = {1}'.format(key, cg_sum))
  ...
  Sum for bar = 2
  Sum for foo = 8
  Sum for qux = 11

If the specified function has a `numpy.ufunc.reduceat` method, this will be called instead.
This can improve the performance by a factor of 10 to 100 (or more) for large unmasked
tables or columns with many relatively small groups.  It also allows for the use of
certain numpy functions which normally take more than one input array but also work as
reduction functions, like `numpy.add`.  The numpy functions which should take advantage of
using `numpy.ufunc.reduceat` include:

`numpy.add`, `numpy.arctan2`, `numpy.bitwise_and`, `numpy.bitwise_or`, `numpy.bitwise_xor`,
`numpy.copysign`, `numpy.divide`, `numpy.equal`, `numpy.floor_divide`, `numpy.fmax`,
`numpy.fmin`, `numpy.fmod`, `numpy.greater_equal`, `numpy.greater`, `numpy.hypot`,
`numpy.left_shift`, `numpy.less_equal`, `numpy.less`, `numpy.logaddexp2`,
`numpy.logaddexp`, `numpy.logical_and`, `numpy.logical_or`, `numpy.logical_xor`,
`numpy.maximum`, `numpy.minimum`, `numpy.mod`, `numpy.multiply`, `numpy.not_equal`,
`numpy.power`, `numpy.remainder`, `numpy.right_shift`, `numpy.subtract` and `numpy.true_divide`.

As special cases `numpy.sum` and `numpy.mean` are substituted with their
respective reduceat methods.


Filtering
~~~~~~~~~~

Table groups can be filtered by means of the
`~astropy.table.groups.TableGroups.filter` method.  This is done by
supplying a function which is called for each group.  The function
which is passed to this method must accept two arguments:

- ``table`` : |Table| object
- ``key_colnames`` : list of columns in ``table`` used as keys for grouping

It must then return either `True` or `False`.  As an example, the following
will select all table groups with only positive values in the non-key columns::

  >>> def all_positive(table, key_colnames):
  ...     colnames = [name for name in table.colnames if name not in key_colnames]
  ...     for colname in colnames:
  ...         if np.any(table[colname] < 0):
  ...             return False
  ...     return True

An example of using this function is::

  >>> t = Table.read(""" a   b    c
  ...                   -2  7.0   0
  ...                   -2  5.0   1
  ...                    1  3.0  -5
  ...                    1 -2.0  -6
  ...                    1  1.0   7
  ...                    0  0.0   4
  ...                    3  3.0   5
  ...                    3 -2.0   6
  ...                    3  1.0   7""", format='ascii')
  >>> tg = t.group_by('a')
  >>> t_positive = tg.groups.filter(all_positive)
  >>> for group in t_positive.groups:
  ...     print(group)
  ...     print
  ...
   a   b   c
  --- --- ---
   -2 7.0   0
   -2 5.0   1
  <BLANKLINE>
   a   b   c
  --- --- ---
    0 0.0   4

As can be seen only the groups with ``a == -2`` and ``a == 0`` have all positive values
in the non-key columns, so those are the ones that are selected.

Likewise a grouped column can be filtered with the
`~astropy.table.groups.ColumnGroups.filter`, method but in this case the filtering
function takes only a single argument which is the column group.  It still must return
either `True` or `False`.  For example::

  def all_positive(column):
      if np.any(column < 0):
          return False
      return True

.. _stack-vertically:

Stack vertically
^^^^^^^^^^^^^^^^^^^^

The |Table| class supports stacking tables vertically with the
`~astropy.table.vstack` function.  This process is also commonly known as
concatenating or appending tables in the row direction.  It corresponds roughly
to the `numpy.vstack` function.

For example, suppose one has two tables of observations with several
column names in common::

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

Notice that the ``obs2`` table is missing the ``mag_b`` column, so in the stacked output
table those values are marked as missing.  This is the default behavior and corresponds to
``join_type='outer'``.  There are two other allowed values for the ``join_type`` argument,
``'inner'`` and ``'exact'``::

  >>> print(vstack([obs1, obs2], join_type='inner'))
    name   obs_date  logLx
  ------- ---------- -----
      M31 2012-01-02  42.5
      M82 2012-10-29  43.5
     M101 2012-10-31  44.5
  NGC3516 2011-11-11  42.1
      M31 1999-01-05  43.1
      M82 2012-10-30  45.0

  >>> print(vstack([obs1, obs2], join_type='exact'))
  Traceback (most recent call last):
    ...
  TableMergeError: Inconsistent columns in input arrays (use 'inner'
  or 'outer' join_type to allow non-matching columns)

In the case of ``join_type='inner'``, only the common columns (the intersection) are
present in the output table.  When ``join_type='exact'`` is specified then
`~astropy.table.vstack` requires that all the input tables
have exactly the same column names.

More than two tables can be stacked by supplying a list of table objects::

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

See also the sections on `Merging metadata`_ and `Merging column
attributes`_ for details on how these characteristics of the input tables are merged in
the single output table.  Note also that you can use a single table row instead of a
full table as one of the inputs.

.. _stack-horizontally:

Stack horizontally
^^^^^^^^^^^^^^^^^^^^^

The |Table| class supports stacking tables horizontally (in the column-wise direction) with the
`~astropy.table.hstack` function.    It corresponds roughly
to the `numpy.hstack` function.

For example, suppose one has the following two tables::

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

As with `~astropy.table.vstack`, there is an optional ``join_type`` argument
that can take values ``'inner'``, ``'exact'``, and ``'outer'``.  The default is
``'outer'``, which effectively takes the union of available rows and masks out any missing
values.  This is illustrated in the example above.  The other options give the
intersection of rows, where ``'exact'`` requires that all tables have exactly the same
number of rows::

  >>> print(hstack([t1, t2], join_type='inner'))
   a   b   c   d     e
  --- --- --- ---- -----
    1 foo 1.4  ham  eggs
    2 bar 2.1 spam toast

  >>> print(hstack([t1, t2], join_type='exact'))
  Traceback (most recent call last):
    ...
  TableMergeError: Inconsistent number of rows in input arrays (use 'inner' or
  'outer' join_type to allow non-matching rows)

More than two tables can be stacked by supplying a list of table objects.  The example
below also illustrates the behavior when there is a conflict in the input column names
(see the section on `Column renaming`_ for details)::

  >>> t3 = Table.read("""a    b
  ...                    M45  2012-02-03""", format='ascii')
  >>> print(hstack([t1, t2, t3]))
  a_1 b_1  c   d     e   a_3    b_3
  --- --- --- ---- ----- --- ----------
    1 foo 1.4  ham  eggs M45 2012-02-03
    2 bar 2.1 spam toast  --         --
    3 baz 2.8   --    --  --         --


The metadata from the input tables is merged by the process described in the `Merging
metadata`_ section.  Note also that you can use a single table row instead of a
full table as one of the inputs.

.. _table-join:

Join
^^^^^^^^^^^^^^

The |Table| class supports the `database join <http://en.wikipedia.org/wiki/Join_(SQL)>`_
operation.  This provides a flexible and powerful way to combine tables based on the
values in one or more key columns.

For example, suppose one has two tables of observations, the first with B and V magnitudes
and the second with X-ray luminosities of an overlapping (but not identical) sample::

  >>> from astropy.table import Table, join
  >>> optical = Table.read("""name    obs_date    mag_b  mag_v
  ...                         M31     2012-01-02  17.0   16.0
  ...                         M82     2012-10-29  16.2   15.2
  ...                         M101    2012-10-31  15.1   15.5""", format='ascii')
  >>> xray = Table.read("""   name    obs_date    logLx
  ...                         NGC3516 2011-11-11  42.1
  ...                         M31     1999-01-05  43.1
  ...                         M82     2012-10-29  45.0""", format='ascii')

The |join| method allows one to merge these two tables into a single table based on
matching values in the "key columns".  By default the key columns are the set of columns
that are common to both tables.  In this case the key columns are ``name`` and
``obs_date``.  We can find all the observations of the same object on the same date as
follows::

  >>> opt_xray = join(optical, xray)
  >>> print(opt_xray)
  name  obs_date  mag_b mag_v logLx
  ---- ---------- ----- ----- -----
   M82 2012-10-29  16.2  15.2  45.0

We can perform the match only by ``name`` by providing the ``keys`` argument, which can be
either a single column name or a list of column names::

  >>> print(join(optical, xray, keys='name'))
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
   M31 2012-01-02  17.0  16.0 1999-01-05  43.1
   M82 2012-10-29  16.2  15.2 2012-10-29  45.0

This output table has all observations that have both optical and X-ray data for an object
(M31 and M82).  Notice that since the ``obs_date`` column occurs in both tables it has
been split into two columns, ``obs_date_1`` and ``obs_date_2``.  The values are taken from
the "left" (``optical``) and "right" (``xray``) tables, respectively.


Different join options
~~~~~~~~~~~~~~~~~~~~~~

The table joins so far are known as "inner" joins and represent the strict intersection of
the two tables on the key columns.

If one wants to make a new table which has *every* row from the left table and includes
matching values from the right table when available, this is known as a left join::

  >>> print(join(optical, xray, join_type='left'))
  name  obs_date  mag_b mag_v logLx
  ---- ---------- ----- ----- -----
  M101 2012-10-31  15.1  15.5    --
   M31 2012-01-02  17.0  16.0    --
   M82 2012-10-29  16.2  15.2  45.0

Two of the observations do not have X-ray data, as indicated by the "--" in the table.
When there are any missing values the output will be a masked table.  You might be
surprised that there is no X-ray data for M31 in the output.  Remember that the default
matching key includes both ``name`` and ``obs_date``.  Specifying the key as only the
``name`` column gives::

  >>> print(join(optical, xray, join_type='left', keys='name'))
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
  M101 2012-10-31  15.1  15.5         --    --
   M31 2012-01-02  17.0  16.0 1999-01-05  43.1
   M82 2012-10-29  16.2  15.2 2012-10-29  45.0

Likewise one can construct a new table with every row of the right table and matching left
values (when available) using ``join_type='right'``.

Finally, to make a table with the union of rows from both tables do an "outer" join::

  >>> print(join(optical, xray, join_type='outer'))
    name   obs_date  mag_b mag_v logLx
  ------- ---------- ----- ----- -----
     M101 2012-10-31  15.1  15.5    --
      M31 1999-01-05    --    --  43.1
      M31 2012-01-02  17.0  16.0    --
      M82 2012-10-29  16.2  15.2  45.0
  NGC3516 2011-11-11    --    --  42.1


Non-identical key column names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The |join| function requires the key column names to be identical in the
two tables. However, in the following one table has a ``'name'`` column
while the other has an ``'obj_id'`` column::

  >>> optical = Table.read("""name    obs_date    mag_b  mag_v
  ...                         M31     2012-01-02  17.0   16.0
  ...                         M82     2012-10-29  16.2   15.2
  ...                         M101    2012-10-31  15.1   15.5""", format='ascii')
  >>> xray_1 = Table.read("""   obj_id    obs_date    logLx
  ...                           NGC3516 2011-11-11  42.1
  ...                           M31     1999-01-05  43.1
  ...                           M82     2012-10-29  45.0""", format='ascii')

In order to perform a match based on the names of the objects, one has to
temporarily rename one of the columns mentioned above, right before creating
the new table::

  >>> xray_1.rename_column('obj_id', 'name')
  >>> opt_xray_1 = join(optical, xray_1, keys='name')
  >>> xray_1.rename_column('name', 'obj_id')
  >>> print(opt_xray_1)
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
  M31 2012-01-02  17.0  16.0 1999-01-05  43.1
  M82 2012-10-29  16.2  15.2 2012-10-29  45.0

The original ``xray_1`` table remains unchanged after the operation::

  >>> print(xray_1)
  obj_id  obs_date  logLx
  ------- ---------- -----
  NGC3516 2011-11-11  42.1
      M31 1999-01-05  43.1
      M82 2012-10-29  45.0


Identical key values
~~~~~~~~~~~~~~~~~~~~

The |Table| join operation works even if there are multiple rows with identical key
values.  For example the following tables have multiple rows for the key column ``x``::

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

Doing an outer join on these tables shows that what is really happening is a `Cartesian
product <http://en.wikipedia.org/wiki/Cartesian_product>`_.  For each matching key, every
combination of the left and right tables is represented.  When there is no match in either
the left or right table, the corresponding column values are designated as missing.

.. doctest-skip:: win32

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

.. note::

   The output table is sorted on the key columns, but when there are rows with identical
   keys the output order in the non-key columns is not guaranteed to be identical across
   installations.  In the example above the order within the four rows with ``key == 1``
   can vary.

An inner join is the same but only returns rows where there is a key match in both the
left and right tables:

.. doctest-skip:: win32

  >>> print(join(left, right, join_type='inner'))
  key  L   R
  --- --- ---
    1  L2  R1
    1  L2  R2
    1  L3  R1
    1  L3  R2
    2  L4  R3

Conflicts in the input table names are handled by the process described in the section on
`Column renaming`_.  See also the sections on `Merging metadata`_ and `Merging column
attributes`_ for details on how these characteristics of the input tables are merged in
the single output table.

Merging details
^^^^^^^^^^^^^^^^^^^^

When combining two or more tables there is the need to merge certain
characteristics in the inputs and potentially resolve conflicts.  This
section describes the process.

Column renaming
~~~~~~~~~~~~~~~~~


In cases where the input tables have conflicting column names, there
is a mechanism to generate unique output column names.  There are two
keyword arguments that control the renaming behavior:

``table_names``
    Two-element list of strings that provide a name for the tables being joined.
    By default this is ``['1', '2', ...]``, where the numbers correspond to
    the input tables.

``uniq_col_name``
    String format specifier with a default value of ``'{col_name}_{table_name}'``.

This is most easily understood by example using the ``optical`` and ``xray`` tables
in the |join| example defined previously::

  >>> print(join(optical, xray, keys='name',
  ...            table_names=['OPTICAL', 'XRAY'],
  ...            uniq_col_name='{table_name}_{col_name}'))
  name OPTICAL_obs_date mag_b mag_v XRAY_obs_date logLx
  ---- ---------------- ----- ----- ------------- -----
   M31       2012-01-02  17.0  16.0    1999-01-05  43.1
   M82       2012-10-29  16.2  15.2    2012-10-29  45.0


Merging metadata
~~~~~~~~~~~~~~~~~~~

|Table| objects can have associated metadata:

- ``Table.meta``: table-level metadata as an ordered dictionary
- ``Column.meta``: per-column metadata as an ordered dictionary

The table operations described here handle the task of merging the metadata in the input
tables into a single output structure.  Because the metadata can be arbitrarily complex
there is no unique way to do the merge.  The current implementation uses a simple
recursive algorithm with four rules:

- `dict` elements are merged by keys
- Conflicting `list` or `tuple` elements are concatenated
- Conflicting `dict` elements are merged by recursively calling the merge function
- Conflicting elements that are not both `list`, `tuple`, or `dict` will follow the following rules:
    - If both metadata values are identical, the output is set to this value
    - If one of the conflicting metadata values is `None`, the other value is picked
    - If both metadata values are different and neither is `None`, the one for the last table in the list is picked

By default, a warning is emitted in the last case (both metadata values are not
`None`). The warning can be silenced or made into an exception using the
``metadata_conflicts`` argument to :func:`~astropy.table.hstack`,
:func:`~astropy.table.vstack`, or
:func:`~astropy.table.join`. The ``metadata_conflicts`` option can be set to:

- ``'silent'`` - no warning is emitted, the value for the last table is silently picked
- ``'warn'`` - a warning is emitted, the value for the last table is picked
- ``'error'`` - an exception is raised

Merging column attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the table and column ``meta`` attributes, the column attributes ``unit``,
``format``, and ``description`` are merged by going through the input tables in
order and taking the first value which is defined (i.e. is not None).  For example::

  >>> from astropy.table import Column, Table, vstack
  >>> col1 = Column([1], name='a')
  >>> col2 = Column([2], name='a', unit='cm')
  >>> col3 = Column([3], name='a', unit='m')
  >>> t1 = Table([col1])
  >>> t2 = Table([col2])
  >>> t3 = Table([col3])
  >>> out = vstack([t1, t2, t3])  # doctest: +SKIP
  WARNING: MergeConflictWarning: In merged column 'a' the 'unit' attribute does
  not match (cm != m).  Using m for merged output [astropy.table.operations]
  >>> out['a'].unit  # doctest: +SKIP
  Unit("m")

The rules for merging are as for `Merging metadata`_, and the
``metadata_conflicts`` option also controls the merging of column attributes.


.. _unique-rows:

Unique rows
^^^^^^^^^^^

Sometimes it makes sense to use only rows with unique key columns or even
fully unique rows from a table. This can be done using the above described
:func:`~astropy.table.Table.group_by` method and ``groups`` attribute, or
with the `~astropy.table.unique` convenience method. The
`~astropy.table.unique` method returns with a sorted table containing the
first row for each unique ``keys`` column value. If no ``keys`` is provided
it returns with a sorted table containing all the fully unique rows.

A simple example is a list of objects with photometry from various observing
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
