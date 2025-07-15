.. doctest-skip-all

.. _df_narwhals:

Interfacing with Pandas/Polars DataFrames via Narwhals
******************************************************

`Narwhals <https://narwhals-dev.github.io/narwhals/>`_ provides a unified translation layer for DataFrame-like objects in Python, including `pandas <https://pandas.pydata.org/>`_, `polars <https://pola.rs/>`__, and others. This enables seamless interoperability between :class:`~astropy.table.Table` and DataFrame objects using generic methods: :meth:`~astropy.table.Table.to_df` and :meth:`~astropy.table.Table.from_df`.

Narwhals supports multiple DataFrame libraries. The following is supported within Astropy:

* Importing from any Narwhals compatible DataFrame, such as `pandas`, ``polars``, ``pyarrow``, ``duckDB``, and others.
* Exporting to Narwhals compatible DataFrames with eager support, such as `pandas`, ``polars``, ``pyarrow``, ``modin``, or ``cudf``.

The following examples and caveats will focus on common DataFrame libraries, primarily `pandas`.

Basic Example
-------------

.. EXAMPLE START: Interfacing Tables with DataFrames via Narwhals

To demonstrate, we can create a minimal table::

    >>> from astropy.table import Table
    >>> t = Table()
    >>> t['a'] = [1, 2, 3, 4]
    >>> t['b'] = ['a', 'b', 'c', 'd']

Which we can then convert to a pandas DataFrame::

    >>> df = t.to_df("pandas")
    >>> df
       a  b
    0  1  a
    1  2  b
    2  3  c
    3  4  d
    >>> type(df)
    <class 'pandas.core.frame.DataFrame'>

To convert to a polars DataFrame, we can specify the polars backend. We can also pass the backend as a module directly::

    >>> import polars as pl
    >>> df_polars = t.to_df(pl)
    >>> type(df_polars)
    <class 'polars.dataframe.frame.DataFrame'>

It is also possible to create a table from a DataFrame::

    >>> t2 = Table.from_df(df)
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

Details and Limitations
-----------------------

Conversions to and from DataFrames are subject to several caveats:

* Multidimensional columns have varied support across DataFrame libraries. Pandas does not support multidimensional columns natively, while Polars can support Nd arrays. PyArrow can support 1D arrays, but higher dimensions are not supported. Exercise caution when converting tables with multidimensional columns.
* Mixed-type columns, i.e. object dtype columns, have varied support across DataFrame libraries. Where possible, we attempt to preserve the original data type, but this may not always be possible. For example, Pandas will preserve object-typed columns while Polars and PyArrow may fail to import the data at all.
* Masked tables are partially supported. DataFrames use sentinel values (e.g., `numpy.nan` or `None`) for missing data, while Astropy preserves the original value under the mask. Masked values may be lost or replaced during conversion.
* Tables with :ref:`mixin_columns` such as `~astropy.time.Time`, `~astropy.coordinates.SkyCoord`, and |Quantity| can be converted, but with *loss of information or fidelity*. These columns are broken down into primitive types (e.g., float, string), and units or metadata are not preserved.
* For time columns, pandas uses 64-bit precision and does not support leap seconds or astronomical time scales. Astropy uses 128-bit precision and supports these features. Polars stores times as strings and does not preserve time semantics.
* Sky coordinate columns are split into separate columns (e.g., ``ra``, ``dec``) and lose angular units and SkyCoord semantics.
* Quantity columns become plain float columns without unit tracking.

Complex Example
---------------

.. EXAMPLE START: Interfacing Tables with DataFrames via Narwhals (Complex Example)

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

Now we convert this table to a pandas DataFrame::

    >>> df = t.to_df(backend="pandas")
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

Or to a polars DataFrame::

    >>> df_polars = t.to_df(backend="polars")
    >>> df_polars
    shape: (3, 7)
    ┌──────┬──────┬──────┬─────────────────────┬───────┬────────┬─────┐
    │ a    ┆ b    ┆ c    ┆ tm                  ┆ sc.ra ┆ sc.dec ┆ q   │
    │ ---  ┆ ---  ┆ ---  ┆ ---                 ┆ ---   ┆ ---    ┆ --- │
    │ i64  ┆ f64  ┆ str  ┆ datetime[ns]        ┆ f64   ┆ f64    ┆ f64 │
    ╞══════╪══════╪══════╪═════════════════════╪═══════╪════════╪═════╡
    │ 1    ┆ 1.0  ┆ null ┆ 2021-01-01 00:00:00 ┆ 1.0   ┆ 4.0    ┆ 1.0 │
    │ null ┆ 2.0  ┆ b    ┆ 2021-01-02 00:00:00 ┆ 2.0   ┆ 5.0    ┆ 2.0 │
    │ 3    ┆ null ┆ c    ┆ 2021-01-03 00:00:00 ┆ 3.0   ┆ 6.0    ┆ 3.0 │
    └──────┴──────┴──────┴─────────────────────┴───────┴────────┴─────┘

Now convert back to a table::

    >>> t_df = QTable.from_df(df)
    >>> t_df
    <QTable length=3>
      a      b     c              tm            sc.ra   sc.dec    q
    int64 float64 str1           Time          float64 float64 float64
    ----- ------- ---- ----------------------- ------- ------- -------
        1     1.0   -- 2021-01-01T00:00:00.000     1.0     4.0     1.0
       --     2.0    b 2021-01-02T00:00:00.000     2.0     5.0     2.0
        3      --    c 2021-01-03T00:00:00.000     3.0     6.0     3.0

The `~astropy.time.Time` column is restored (subject to the limitations discussed previously), but the `~astropy.coordinates.SkyCoord` and `~astropy.units.Quantity` columns are not restored as they were in the original table.

Finally, see that masked values in the original table are replaced with zero, empty string, or ``nan`` in the round-trip conversion::

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

See the `Narwhals <https://narwhals-dev.github.io/narwhals/>`_ documentation for more details on Narwhals and supported backends.
