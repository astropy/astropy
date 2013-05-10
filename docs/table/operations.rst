.. include:: references.txt
.. |join| replace:: :func:`~astropy.table.operations.join`

.. _table_operations:

Table operations
-----------------

In this section we describe higher-level operations that can be used to generate a new
table from two or more input tables.  This includes:

=======================

.. list-table::
   :header-rows: 1
   :widths: 28 52 20

   * - Documentation
     - Description
     - Function
   * - `Stack vertically`_
     - Concatenate input tables along rows
     - `~astropy.table.operations.vstack`
   * - `Stack horizontally`_
     - Concatenate input tables along columns
     - `~astropy.table.operations.hstack`
   * - `Join`_
     - Database-style join of two tables
     - `~astropy.table.operations.join`


Stack vertically
^^^^^^^^^^^^^^^^^^^^

The |Table| class supports stacking tables vertically with the
`~astropy.table.operations.vstack` function.  This process is also commonly known as
concatenating or appending tables in the row direction.  It corresponds roughly
to the `numpy.vstack` function.

For example, suppose one has two tables of observations with several
column names in common::

  >>> from astropy.table import Table, vstack
  >>> obs1 = Table.read("""name    obs_date    mag_b  logLx
                           M31     2012-01-02  17.0   42.5
                           M82     2012-10-29  16.2   43.5
                           M101    2012-10-31  15.1   44.5""", format='ascii')

  >>> obs2 = Table.read("""name    obs_date    logLx
                           NGC3516 2011-11-11  42.1
                           M31     1999-01-05  43.1
                           M82     2012-10-30  45.0""", format='ascii')

Now we can stack these two tables::

  >>> print vstack([obs1, obs2])
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

  >>> print vstack([obs1, obs2], join_type='inner')
    name   obs_date  logLx
  ------- ---------- -----
      M31 2012-01-02  42.5
      M82 2012-10-29  43.5
     M101 2012-10-31  44.5
  NGC3516 2011-11-11  42.1
      M31 1999-01-05  43.1
      M82 2012-10-30  45.0

  >>> print vstack([obs1, obs2], join_type='exact')
  ...
  TableMergeError: Inconsistent columns in input arrays (use 'inner' or
  'outer' join_type to allow non-matching columns)

In the case of ``join_type='inner'`, only the common columns (the intersection) are
present in the output table.  When ``join_type='exact'`` is specified then
`~astropy.table.operations.vstack` requires that all the input tables
have exactly the same column names.

More than two tables can be stacked by supplying a list of table objects::

  >>> obs3 = Table.read("""name    obs_date    mag_b  logLx
                           M45     2012-02-03  15.0   40.5""", format='ascii')
  >>> print vstack([obs1, obs2, obs3])
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
the single output table.

Stack horizontally
^^^^^^^^^^^^^^^^^^^^^

The |Table| class supports stacking tables horizontally (in the column-wise direction) with the
`~astropy.table.operations.hstack` function.    It corresponds roughly
to the `numpy.hstack` function.

For example, suppose one has the following two tables::

  >>> from astropy.table import Table, hstack
  >>> t1 = Table.read("""a   b    c
                         1   foo  1.4
                         2   bar  2.1
                         3   baz  2.8""", format='ascii')

  >>> t2 = Table.read("""d     e
                         ham   eggs
                         spam  toast""", format='ascii')

Now we can stack these two tables horizontally::

  >>> print hstack([t1, t2])
   a   b   c   d     e
  --- --- --- ---- -----
    1 foo 1.4  ham  eggs
    2 bar 2.1 spam toast
    3 baz 2.8   --    --

As with `~astropy.table.operations.vstack`, there is an optional ``join_type`` argument
that can take values `'inner'``, ``'exact'``, and ``'outer'``.  The default is
``'outer'``, which effectively takes the union of available rows and masks out any missing
values.  This is illustrated in the example above.  The other options give the
intersection of rows, where ``'exact'`` requires that all tables have exactly the same
number of rows::

  >>> print hstack([t1, t2], join_type='inner')
   a   b   c   d     e
  --- --- --- ---- -----
    1 foo 1.4  ham  eggs
    2 bar 2.1 spam toast

  >>> print hstack([t1, t2], join_type='exact')
  ...
  TableMergeError: Inconsistent number of rows in input arrays (use 'inner' or
  'outer' join_type to allow non-matching columns)

More than two tables can be stacked by supplying a list of table objects.  The example
below also illustrates the behavior when there is a conflict in the input column names
(see the section on `Column renaming`_ for details)::

  >>> t3 = Table.read("""a    b
                         M45  2012-02-03""", format='ascii')
  >>> print hstack([t1, t2, t3])
  a_1 b_1  c   d     e   a_3    b_3
  --- --- --- ---- ----- --- ----------
    1 foo 1.4  ham  eggs M45 2012-02-03
    2 bar 2.1 spam toast  --         --
    3 baz 2.8   --    --  --         --


The metadata from the input tables is merged by the process described in the `Merging
metadata`_ section.

Join
^^^^^^^^^^^^^^

The |Table| class supports the `database join <http://en.wikipedia.org/wiki/Join_(SQL)>`_
operation.  This provides a flexible and powerful way to combine tables based on the
values in one or more key columns.

For example, suppose one has two tables of observations, the first with B and V magnitudes
and the second with X-ray luminosities of an overlapping (but not identical) sample::

  >>> from astropy.table import Table, join
  >>> optical = Table.read("""name    obs_date    mag_b  mag_v
                              M31     2012-01-02  17.0   16.0
                              M82     2012-10-29  16.2   15.2
                              M101    2012-10-31  15.1   15.5""", format='ascii')

  >>> xray = Table.read("""   name    obs_date    logLx
                              NGC3516 2011-11-11  42.1
                              M31     1999-01-05  43.1
                              M82     2012-10-29  45.0""", format='ascii')

The |join| method allows one to merge these two tables into a single table based on
matching values in the "key columns".  By default the key columns are the set of columns
that are common to both tables.  In this case the key columns are ``name`` and
``obs_date``.  We can find all the observations of the same object on the same date as
follows::

  >>> opt_xray = join(optical, xray)
  >>> print opt_xray
  name  obs_date  mag_b mag_v logLx
  ---- ---------- ----- ----- -----
   M82 2012-10-29  16.2  15.2  45.0

We can perform the match only by ``name`` by providing the ``keys`` argument, which can be
either a single column name or a list of column names::

  >>> print join(optical, xray, keys='name')
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
   M31 2012-01-02  17.0  16.0 1999-01-05  43.1
   M82 2012-10-29  16.2  15.2 2012-10-29  45.0

This output table has all observations that have both optical and X-ray data for an object
(M31 and M82).  Notice that since the ``obs_date`` column occurs in both tables it has
been split into two columns, ``obs_date_1`` and ``obs_date_2``.  The values are taken from
the "left" (``optical``) and "right" (``xray``) tables, respectively.

The table joins so far are known as "inner" joins and represent the strict intersection of
the two tables on the key columns.

If one wants to make a new table which has *every* row from the left table and includes
matching values from the right table when available, this is known as a left join::

  >>> print join(optical, xray, join_type='left')
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

  >>> print join(optical, xray, join_type='left', keys='name')
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
  M101 2012-10-31  15.1  15.5         --    --
   M31 2012-01-02  17.0  16.0 2012-01-05  43.1
   M82 2012-10-29  16.2  15.2 2012-10-29  45.0

Likewise one can construct a new table with every row of the right table and matching left
values (when available) using ``join_type='right'``.

Finally, to make a table with the union of rows from both tables do an "outer" join::

  >>> print join(optical, xray, join_type='outer')
    name   obs_date  mag_b mag_v logLx
  ------- ---------- ----- ----- -----
     M101 2012-10-31  15.1  15.5    --
      M31 2012-01-02  17.0  16.0    --
      M31 2012-01-05    --    --  43.1
      M82 2012-10-29  16.2  15.2  45.0
  NGC3516 2011-11-11    --    --  42.1


Identical keys
~~~~~~~~~~~~~~

The |Table| join operation works even if there are multiple rows with identical key
values.  For example the following tables have multiple rows for the key column ``x``::

  >>> from astropy.table import Table, join
  >>> left = Table([[0, 1, 1, 2], ['L1', 'L2', 'L3', 'L4']], names=('key', 'L'))
  >>> right = Table([[1, 1, 2, 4], ['R1', 'R2', 'R3', 'R4']], names=('key', 'R'))
  >>> print left
  key  L
  --- ---
    0  L1
    1  L2
    1  L3
    2  L4
  >>> print right
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

  >>> print join(left, right, join_type='outer')
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
left and right tables::

  >>> print join(left, right, join_type='inner')
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

  >>> print join(optical, xray, keys='name',
                 table_names=['OPTICAL', 'XRAY'],
                 uniq_col_name='{table_name}_{col_name}')
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

- Dict elements are merged by keys
- Conflicting list or tuple elements are concatenated
- Conflicting dict elements are merged by recursively calling the merge function
- Conflicting elements that are not both list, tuple, or dict results in an exception

Merging column attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the table and column ``meta`` attributes, the column attributes ``units``,
``format``, and ``description`` are merged by going through the input tables in
order and taking the first value which is defined (i.e. is not None).  For example::

  >>> from astropy.table import Column, Table, vstack
  >>> col1 = Column([1], name='a')
  >>> col2 = Column([2], name='a', units='cm')
  >>> col3 = Column([3], name='a', units='m')
  >>> t1 = Table([col1])
  >>> t2 = Table([col2])
  >>> t3 = Table([col3])
  >>> out = vstack([t1, t2, t3])
  WARNING: MergeConflictWarning: In merged column 'a' the 'units' attribute does
  not match (cm != m).  Using cm for merged output [astropy.table.table]
  >>> out['a'].units
  Unit("cm")

In this case there was a conflict so a warning is shown.
