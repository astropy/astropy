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

Example
-------

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

The conversions to and from ``pandas`` are subject to the following caveats:

* The :class:`~pandas.DataFrame` structure does not support multidimensional
  columns, so |Table| objects with multidimensional columns cannot be converted
  to :class:`~pandas.DataFrame`.

* Masked tables can be converted, but in columns of ``float`` or string values
  the resulting :class:`~pandas.DataFrame` uses `numpy.nan` to indicate missing
  values. For ``float`` columns, the conversion therefore does not necessarily
  round-trip if converting back to an ``astropy`` table, because the
  distinction between `numpy.nan` and masked values is lost. This is not a
  problem for integer columns.

* Tables with :ref:`mixin_columns` can not be converted.
