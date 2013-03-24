.. _table_utilities:

.. include:: references.txt

Table operations
-----------------

In this section we describe higher-level operations that can be used to generate a new table from two or more input tables.

Join
^^^^^^^^^^^^^^

The |Table| class supports the `database join <wikipedia ref>`_ operation.  This
provides a flexible and powerful way to combine tables based on the 
values in one or more key columns.

For example, suppose two tables of observations, the first with B and V magnitudes and the second with X-ray luminosities of an overlapping (but not identical) sample::

  >>> from astropy.io import ascii
  >>> optical = ascii.read("""name    obs_date    mag_b  mag_v
                              M31     2012-01-02  17.0   16.0
                              M82     2012-10-29  16.2   15.2
                              M101    2012-10-31  15.1   15.5""")

  >>> xray = ascii.read("""   name    obs_date    logLx
                              NGC3516 2011-11-11  42.1
                              M31     1999-01-05  43.1
                              M82     2012-10-29  45.0""")

The '~astropy.table.table.join` method allows one to merge these two tables into a single table based on matching values in the "key columns".  By default the key columns are the set of columns that are common to both tables.  In this case the key columns are``name`` and ``obs_date``.  We can find all the observations of the same object on the same date as follows::

  >>> print optical.join(xray)
  name  obs_date  mag_b mag_v logLx
  ---- ---------- ----- ----- -----
   M82 2012-10-29  16.2  15.2  45.0

If the observations do not need to be simultaneous then we can choose to match
only by ``name`` by providing the ``keys`` argument::

  >>> print optical.join(xray, keys='name')
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
   M31 2012-01-02  17.0  16.0 1999-01-05  43.1
   M82 2012-10-29  16.2  15.2 2012-10-29  45.0

This output table has all observations that have both optical and X-ray data for an object (M31 and M82).  Notice that since the ``obs_date`` column occurs in both tables it has been split into two columns, ``obs_date_1`` and ``obs_date_2``.  The values are taken from the "left" (``optical``) and "right" (``xray``) tables, respectively.

The table joins so far are known as "inner" joins and represent the strict intersection of the two tables on the key columns.

If one wants to make a new table which has *every* row from the left table and includes matching values from the right table when available, this is known as a left join::

  >>> print optical.join(xray, join_type='left')
  name  obs_date  mag_b mag_v logLx
  ---- ---------- ----- ----- -----
  M101 2012-10-31  15.1  15.5    --
   M31 2012-01-02  17.0  16.0    --
   M82 2012-10-29  16.2  15.2  45.0

Two of the observations do not have X-ray data, as indicated by the "--" in the table.  When there are any missing values the output will be a masked table.  You might be surprised that there is no X-ray data for M31 in the output.  Remember that the default matching key includes both ``name`` and ``obs_date``.  Specifying the key as only the ``name`` column gives::

  >>> print optical.join(xray, join_type='left', keys='name')
  name obs_date_1 mag_b mag_v obs_date_2 logLx
  ---- ---------- ----- ----- ---------- -----
  M101 2012-10-31  15.1  15.5         --    --
   M31 2012-01-02  17.0  16.0 2012-01-05  43.1
   M82 2012-10-29  16.2  15.2 2012-10-29  45.0

Likewise one can construct a new table with every row of the right table and matching left values (when available) using ``join_type='right'``.

Finally, to make a table with the union of rows from both tables do an "outer" join::

  >>> print optical.join(xray, join_type='outer')
    name   obs_date  mag_b mag_v logLx
  ------- ---------- ----- ----- -----
     M101 2012-10-31  15.1  15.5    --
      M31 2012-01-02  17.0  16.0    --
      M31 2012-01-05    --    --  43.1
      M82 2012-10-29  16.2  15.2  45.0
  NGC3516 2011-11-11    --    --  42.1


Identical keys
~~~~~~~~~~~~~~

The |Table| join operation works even if there are multiple rows with identical key values.  For example the following tables have multiple rows for the key column ``x``::

  >>> from astropy.table import Table
  >>> left = Table([[0, 1, 1, 2], ['L1', 'L2', 'L3', 'L4']], names=('x', 'y'))
  >>> right = Table([[1, 1, 2, 4], ['R1', 'R2', 'R3', 'R4']], names=('x', 'z'))
  >>> print left
   x   y 
  --- ---
    0  L1
    1  L2
    1  L3
    2  L4
  >>> print right
   x   z 
  --- ---
    1  R1
    1  R2
    2  R3
    4  R4

Doing an outer join on these tables shows that what is really happening is a "cartesian" join.  For each matching key, every combination of the left and right tables is represented.  When there is no match in either the left or right table, the corresponding column values are designated as missing.

  >>> print left.join(right, join_type='outer')
   x   y   z 
  --- --- ---
    0  L1  --
    1  L2  R1
    1  L2  R2
    1  L3  R1
    1  L3  R2
    2  L4  R3
    4  --  R4

An inner join is the same but only returns rows where there is a key match in both the left and right tables::

  >>> print left.join(right, join_type='inner')
   x   y   z 
  --- --- ---
    1  L2  R1
    1  L2  R2
    1  L3  R1
    1  L3  R2
    2  L4  R3
