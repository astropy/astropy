.. doctest-skip-all

.. _pandas:

Interfacing with the Pandas Package
***********************************

The `pandas <https://pandas.pydata.org/>`__ package is a package for high
performance data analysis of table-like structures that is complementary to the
:class:`~astropy.table.Table` class in ``astropy``.

In order to exchange data between the :class:`~astropy.table.Table` class and
the :class:`pandas.DataFrame` class (the main data structure in ``pandas``),
the |Table| class includes two methods, :meth:`~astropy.table.Table.to_pandas`
and :meth:`~astropy.table.Table.from_pandas`.

Basic Example
-------------

.. EXAMPLE START: Interfacing Tables with the Pandas Package

To demonstrate, we can create a minimal table::

    >>> from astropy.table import Table
    >>> t = Table()
    >>> t['a'] = [1, 2, 3, 4]
    >>> t['b'] = ['a', 'b', 'c', 'd']

Which we can then convert to a :class:`~pandas.DataFrame`::

    >>> df = t.to_pandas()
    >>> df
       a  b
    0  1  a
    1  2  b
    2  3  c
    3  4  d
    >>> type(df)
    <class 'pandas.core.frame.DataFrame'>

It is also possible to create a table from a :class:`~pandas.DataFrame`::

    >>> t2 = Table.from_pandas(df)
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
The conversions to and from ``pandas`` are subject to the following caveats:

* The :class:`~pandas.DataFrame` structure does not support multidimensional
  columns, so |Table| objects with multidimensional columns cannot be converted
  to :class:`~pandas.DataFrame`.

* Masked tables can be converted, but pandas uses a sentinel such as `numpy.nan` to
  represent missing values. Astropy uses a separate mask array to represent masked
  values, where the value "under the mask" is preserved. When converting any astropy
  masked column to a pandas DataFrame, the original values are lost.

* Tables with :ref:`mixin_columns` such as `~astropy.time.Time`,
  `~astropy.coordinates.SkyCoord`, and |Quantity| can be converted, but
  *with loss of information or fidelity*. For instance, `~astropy.time.Time` columns
  will be converted to a `pandas TimeSeries
  <https://pandas.pydata.org/docs/user_guide/timeseries.html>`_, but this object has
  only 64-bit precision and does not support leap seconds or time scales.

These issues are highlighted below in a more complex example with a table that includes
masked and mixin columns.

.. EXAMPLE START: Interfacing Tables with the Pandas Package (Complex Example)

First we create a table with a masked columns and a mixin column::

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

Now we convert this table to a :class:`~pandas.DataFrame`::

    >>> df = t.to_pandas()
    >>> df
          a    b    c         tm  sc.ra  sc.dec    q
    0     1  1.0  NaN 2021-01-01    1.0     4.0  1.0
    1  <NA>  2.0    b 2021-01-02    2.0     5.0  2.0
    2     3  NaN    c 2021-01-03    3.0     6.0  3.0

    >>> df.info()
    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 3 entries, 0 to 2
    Data columns (total 7 columns):
     #   Column  Non-Null Count  Dtype
    ---  ------  --------------  -----
     0   a       2 non-null      Int64
     1   b       2 non-null      float64
     2   c       2 non-null      object
     3   tm      3 non-null      datetime64[ns]
     4   sc.ra   3 non-null      float64
     5   sc.dec  3 non-null      float64
     6   q       3 non-null      float64
    dtypes: Int64(1), datetime64[ns](1), float64(4), object(1)
    memory usage: 303.0+ bytes

Notice a few things:

- The masked values in the original table are replaced with sentinel values
  in pandas. The integer column ``a`` is converted to a nullable integer column, and
  the string column ``c`` is converted to an ``object`` column.
- The `~astropy.time.Time` object is converted to a pandas TimeSeries using
  ``datetime64[ns]``.
- The `~astropy.coordinates.SkyCoord` object is converted to two float columns
  ``sc.ra`` and ``sc.dec``, and the unit is lost.
- The `~astropy.units.Quantity` object is converted to a float column and the unit is
  lost.

Now convert back to a table::

    >>> t_df = QTable.from_pandas(df)
    >>> t_df
    <QTable length=3>
      a      b     c              tm            sc.ra   sc.dec    q
    int64 float64 str1           Time          float64 float64 float64
    ----- ------- ---- ----------------------- ------- ------- -------
        1     1.0   -- 2021-01-01T00:00:00.000     1.0     4.0     1.0
       --     2.0    b 2021-01-02T00:00:00.000     2.0     5.0     2.0
        3      --    c 2021-01-03T00:00:00.000     3.0     6.0     3.0

The `~astropy.time.Time` column is restored (subject to the limitations discussed
previously), but the `~astropy.coordinates.SkyCoord` and `~astropy.units.Quantity`
columns are not restored as they were in the original table.

Finally see that the masked values in the original table are replaced with zero or "" in
the round-trip conversion::

    # Original data values
    >>> for nm in 'a', 'b', 'c':
    ...     print(t[nm].data.data)
    [1 2 3]
    [1. 2. 3.]
    ['a' 'b' 'c']

    # Data values after round-trip conversion
    >>> for nm in 'a', 'b', 'c':
    ...     print(t_df[nm].data.data)
    [1 0 3]
    [ 1.  2. nan]
    ['' 'b' 'c']

.. EXAMPLE END
