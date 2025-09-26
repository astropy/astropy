.. doctest-skip-all

.. _df_narwhals:

Interfacing with DataFrames
***************************

The :class:`~astropy.table.Table` class provides comprehensive support for interfacing with DataFrame libraries through two complementary approaches:

1. **Generic narwhals-based methods**: :meth:`~astropy.table.Table.to_df` and :meth:`~astropy.table.Table.from_df` for multi-backend DataFrame support
2. **Legacy pandas-specific methods**: :meth:`~astropy.table.Table.to_pandas` and :meth:`~astropy.table.Table.from_pandas` for direct pandas integration

The narwhals-based approach uses `Narwhals <https://narwhals-dev.github.io/narwhals/>`_ as a unified translation layer, enabling seamless interoperability with multiple DataFrame libraries including `pandas <https://pandas.pydata.org/>`_, `polars <https://pola.rs/>`__, `pyarrow <https://arrow.apache.org/docs/python/>`__, and others.

Generic Multi-Backend Methods
=============================

For users working with multiple DataFrame libraries or seeking broader compatibility, the generic methods provide a unified interface through narwhals.

Supported Backends
------------------

The generic methods support DataFrame libraries with eager execution, including:

* **pandas** - A popular and pioneering DataFrame library
* **polars** - High-performance DataFrame library with different handling of multidimensional data
* **pyarrow** - In-memory columnar format with good performance
* **modin** - Distributed pandas-compatible DataFrames
* **cudf** - GPU-accelerated DataFrames

.. note::
   **Testing and Support**: **pandas**, **polars**, and **pyarrow** are directly tested in the Astropy test suite. While other narwhals-compatible backends should work in principle, they may exhibit unexpected behavior or incompatibilities. If you encounter issues with any backend, please file a bug report on the `Astropy GitHub repository <https://github.com/astropy/astropy/issues>`_.

.. warning::
   **Backend Differences**: Different DataFrame libraries implement varying data models, type systems, and computational paradigms. These fundamental differences can lead to inconsistent behavior across backends, particularly with respect to data type handling, missing value representation, and memory layout. Users should verify that round-trip conversions preserve the expected data integrity for their specific use case and chosen backend.

Basic Multi-Backend Example
---------------------------

.. EXAMPLE START: Using Generic Multi-Backend Methods

Create a table and convert to different DataFrame backends::

    >>> from astropy.table import Table
    >>> t = Table()
    >>> t['a'] = [1, 2, 3, 4]
    >>> t['b'] = ['a', 'b', 'c', 'd']

    # Convert to pandas DataFrame
    >>> df_pandas = t.to_df("pandas")
    >>> type(df_pandas)
    <class 'pandas.core.frame.DataFrame'>

    # Convert to polars DataFrame
    >>> df_polars = t.to_df("polars")
    >>> type(df_polars)
    <class 'polars.dataframe.frame.DataFrame'>

Create a table from any supported DataFrame::

    >>> t2 = Table.from_df(df_pandas)  # From pandas
    >>> t3 = Table.from_df(df_polars)  # From polars

.. EXAMPLE END

Known Backend-Specific Differences
----------------------------------

Different DataFrame backends handle data differently:

**Multidimensional Columns:**
  - Pandas: Not supported, raises an error
  - Polars: Supported as Array type for arbitrary dimensions
  - PyArrow: Limited support for 1D arrays, currently unavailable.

**Index Support:**
  - Pandas: Full index support with :meth:`~astropy.table.Table.to_df`
  - Other backends: Index parameter raises an error (not supported)

**Missing Value Handling:**
  - All backends use sentinel values (NaN, null) rather than Astropy's mask arrays

Pandas-Specific Methods
=======================

For pandas users, Astropy provides dedicated methods that offer the most direct and feature-complete integration with pandas DataFrames.

Basic Pandas Example
--------------------

.. EXAMPLE START: Using Pandas-Specific Methods

To demonstrate, we can create a minimal table::

    >>> from astropy.table import Table
    >>> t = Table()
    >>> t['a'] = [1, 2, 3, 4]
    >>> t['b'] = ['a', 'b', 'c', 'd']

Convert to a pandas DataFrame using the pandas-specific method::

    >>> df = t.to_pandas()
    >>> df
       a  b
    0  1  a
    1  2  b
    2  3  c
    3  4  d
    >>> type(df)
    <class 'pandas.core.frame.DataFrame'>

Create a table from a pandas DataFrame::

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

Pandas Index Support
--------------------

The pandas-specific methods provide full support for DataFrame indexing, which is a unique pandas feature::

    >>> from astropy.time import Time
    >>> tm = Time([1998, 2002], format="jyear")
    >>> x = [1, 2]
    >>> t = Table([tm, x], names=["tm", "x"])

    # Use a column as the DataFrame index
    >>> df = t.to_pandas(index="tm")
    >>> df.index.name
    'tm'

    # Convert back including the index as a column
    >>> t_back = Table.from_pandas(df, index=True)
    >>> t_back.colnames
    ['tm', 'x']

When to Use Which Method
========================

**Use pandas-specific methods** (:meth:`~astropy.table.Table.to_pandas`, :meth:`~astropy.table.Table.from_pandas`) when:

* Working exclusively with pandas
* Need DataFrame index support
* Want the most battle-tested and feature-complete pandas integration
* Require the best performance for pandas-specific workflows

**Use generic methods** (:meth:`~astropy.table.Table.to_df`, :meth:`~astropy.table.Table.from_df`) when:

* Working with multiple DataFrame backends
* Need to support polars, pyarrow, or other backends
* Building library code that should work with various DataFrame types
* Want forward compatibility as new backends are added to narwhals

Conversion Details and Limitations
===================================

Both approaches share common limitations when converting between Tables and DataFrames:

Data Type Limitations
---------------------

* **Multidimensional columns**: Support varies by backend. At the time of writing, Pandas does not support them, while Polars can handle them as Array types.

* **Masked tables**: DataFrames typically use sentinel values (e.g., `numpy.nan` or `None`) for missing data, while Astropy preserves the original value under the mask. The original values under the mask are lost during conversion.

* **Mixed-type columns**: Object dtype columns have varied support. Pandas preserves them while other backends may fail to import such data.

Mixin Column Limitations
------------------------

Tables with :ref:`mixin_columns` such as `~astropy.time.Time`, `~astropy.coordinates.SkyCoord`, and |Quantity| can be converted, but **with loss of information**:

* **Time columns**: Converted to native datetime types with reduced precision and loss of astronomical time scale information
* **SkyCoord columns**: Split into separate coordinate component columns (e.g., ``ra``, ``dec``) with loss of units and coordinate frame information
* **Quantity columns**: Converted to plain numeric columns with complete loss of unit information

Complex Example with Both Methods
==================================

.. EXAMPLE START: Complex DataFrame Conversion Example

Create a table with masked and mixin columns::

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

Convert using the pandas-specific method::

    >>> df_pandas = t.to_pandas()
    >>> df_pandas
          a    b    c         tm  sc.ra  sc.dec    q
    0     1  1.0  NaN 2021-01-01    1.0     4.0  1.0
    1  <NA>  2.0    b 2021-01-02    2.0     5.0  2.0
    2     3  NaN    c 2021-01-03    3.0     6.0  3.0

Convert using the generic method to pandas::

    >>> df_generic = t.to_df("pandas")
    >>> # Results are identical to df_pandas

Convert to polars using the generic method::

    >>> df_polars = t.to_df("polars")
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

Convert back to tables::

    >>> t_from_pandas = QTable.from_pandas(df_pandas)  # Using pandas-specific method
    >>> t_from_generic = QTable.from_df(df_pandas)     # Using generic method
    >>> t_from_polars = QTable.from_df(df_polars)      # From polars DataFrame

Note the data transformations that occurred:

1. **Masked values**: Original values under the mask are lost and replaced with sentinel values
2. **Time columns**: Converted to basic datetime objects, preserving temporal information but losing astropy Time features
3. **SkyCoord columns**: Split into separate ``ra`` and ``dec`` columns, losing coordinate frame and unit information
4. **Quantity columns**: Converted to plain float columns, completely losing unit information

.. EXAMPLE END

Method Reference
================

Pandas-Specific Methods
-----------------------

- :meth:`~astropy.table.Table.to_pandas` - Convert Table to pandas DataFrame
- :meth:`~astropy.table.Table.from_pandas` - Create Table from pandas DataFrame

Generic Multi-Backend Methods
-----------------------------

- :meth:`~astropy.table.Table.to_df` - Convert Table to DataFrame using specified backend
- :meth:`~astropy.table.Table.from_df` - Create Table from any narwhals-compatible DataFrame

See the `Narwhals documentation <https://narwhals-dev.github.io/narwhals/>`_ for more information about supported backends and their capabilities.
