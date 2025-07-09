.. doctest-skip-all

.. _polars:

Interfacing with the Polars Package
***********************************

The `polars <https://pola.rs/>`__ package is a high-performance, multi-threaded
DataFrame library for table-like data, which complements the
:class:`~astropy.table.Table` class in ``astropy``.

To enable interoperability, the |Table| class includes two methods:
:meth:`~astropy.table.Table.to_polars` and :meth:`~astropy.table.Table.from_polars`,
which convert between the |Table| and ``polars.DataFrame`` classes.

Basic Example
-------------

.. EXAMPLE START: Interfacing Tables with the Polars Package

To demonstrate, we can create a minimal table::

    >>> from astropy.table import Table
    >>> t = Table()
    >>> t['a'] = [1, 2, 3, 4]
    >>> t['b'] = ['a', 'b', 'c', 'd']

Which we can then convert to a ``polars.DataFrame``::

    >>> df = t.to_polars()
    >>> df
    shape: (4, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ str │
    ├─────┼─────┤
    │ 1   ┆ a   │
    │ 2   ┆ b   │
    │ 3   ┆ c   │
    │ 4   ┆ d   │
    └─────┴─────┘
    >>> type(df)
    <class 'polars.dataframe.frame.DataFrame'>

It is also possible to create a table from a ``polars.DataFrame``::

    >>> t2 = Table.from_polars(df)
    >>> t2
    <Table length=4>
      a      b
    int64 string8
    ----- -------
        1       a
        2       b
        3       c
        4       d

.. EXAMPLE END

Details
-------

Conversions to and from ``polars`` are subject to a few limitations:

* The ``polars.DataFrame`` structure does not support multidimensional
  columns. |Table| objects with such columns cannot be converted to
  ``polars.DataFrame``.

* Masked tables are partially supported. ``polars`` represents missing values
  using a sentinel (`None`) similar to ``pandas``, but the original value
  under the mask is lost.

* Tables with :ref:`mixin_columns` such as `~astropy.time.Time`,
  `~astropy.coordinates.SkyCoord`, and |Quantity| can be converted, but with
  *loss of information or fidelity*. These objects are broken down into one or
  more primitive columns (e.g., float or string), and units or metadata are not
  preserved.

The following example demonstrates a more complex table with masked and mixin columns:

.. EXAMPLE START: Interfacing Tables with the Polars Package (Complex Example)

First, we create a table with masked and mixin columns::

    >>> import numpy as np
    >>> from astropy.table import MaskedColumn, QTable
    >>> from astropy.time import Time
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> t = QTable()
    >>> t['a'] = MaskedColumn([1, 2, 3], mask=[False, True, False])
    >>> t['b'] = MaskedColumn([1.0, 2.0, 3.0], mask=[False, False, True])
    >>> t['c'] = MaskedColumn(["a", "b", "c"], mask=[True, False, False])
    >>> t['tm'] = Time(["2021-01-01", "2021-01-02", "2021-01-03"])
    >>> t['sc'] = SkyCoord(ra=[1, 2, 3] * u.deg, dec=[4, 5, 6] * u.deg)
    >>> t['q'] = [1, 2, 3] * u.m

    >>> t
    <QTable length=3>
      a      b     c              tm              sc       q
                                               deg,deg     m
    int64 float64 str1           Time          SkyCoord float64
    ----- ------- ---- ----------------------- -------- -------
        1     1.0   -- 2021-01-01 00:00:00.000  1.0,4.0     1.0
       --     2.0    b 2021-01-02 00:00:00.000  2.0,5.0     2.0
        3      --    c 2021-01-03 00:00:00.000  3.0,6.0     3.0

Now we convert this table to a ``polars.DataFrame``::

    >>> df = t.to_polars()
    >>> df
    shape: (3, 7)
    ┌──────┬──────┬──────┬─────────────────────┬───────┬────────┬─────┐
    │ a    ┆ b    ┆ c    ┆ tm                  ┆ sc.ra ┆ sc.dec ┆ q   │
    │ ---  ┆ ---  ┆ ---  ┆ ---                 ┆ ---   ┆ ---    ┆ --- │
    │ i64  ┆ f64  ┆ str  ┆ str                 ┆ f64   ┆ f64    ┆ f64 │
    ├──────┼──────┼──────┼─────────────────────┼───────┼────────┼─────┤
    │ 1    ┆ 1.0  ┆ null ┆ 2021-01-01T00:00:00Z ┆ 1.0   ┆ 4.0    ┆ 1.0 │
    │ null ┆ 2.0  ┆ b    ┆ 2021-01-02T00:00:00Z ┆ 2.0   ┆ 5.0    ┆ 2.0 │
    │ 3    ┆ null ┆ c    ┆ 2021-01-03T00:00:00Z ┆ 3.0   ┆ 6.0    ┆ 3.0 │
    └──────┴──────┴──────┴─────────────────────┴───────┴────────┴─────┘

    >>> df.schema
    {'a': pl.Int64, 'b': pl.Float64, 'c': pl.Utf8, 'tm': pl.Utf8,
     'sc.ra': pl.Float64, 'sc.dec': pl.Float64, 'q': pl.Float64}

Key observations:

- Masked values are represented as `null` in ``polars``.
- The `~astropy.time.Time` column is converted to an ISO string column.
- The `~astropy.coordinates.SkyCoord` column is split into separate `ra` and `dec` float columns, losing angular units and SkyCoord semantics.
- The `~astropy.units.Quantity` column becomes a plain float column without unit tracking.

Now convert back to a table::

    >>> t_df = QTable.from_polars(df)
    >>> t_df
    <QTable length=3>
      a      b     c              tm            sc.ra   sc.dec    q
    int64 float64 str1           str           float64 float64 float64
    ----- ------- ---- --------------------- -------- -------- -------
        1     1.0       2021-01-01T00:00:00Z     1.0      4.0     1.0
              2.0    b  2021-01-02T00:00:00Z     2.0      5.0     2.0
        3            c  2021-01-03T00:00:00Z     3.0      6.0     3.0

The `~astropy.time.Time` column is not automatically parsed and is returned as a string.
Other mixin columns are not restored and require manual reconstruction if needed.

Finally, observe that masked values are not preserved through this round-trip::

    >>> for nm in 'a', 'b', 'c':
    ...     print(t_df[nm].data)
    [1 -- 3]
    [1.0 2.0 nan]
    ['' 'b' 'c']

.. EXAMPLE END
