.. doctest-skip-all

.. _pandas:

Interfacing with the pandas package
===================================

The `pandas <http://pandas.pydata.org/>`__ package is a package for high
performance data analysis of table-like structures that is complementary to the
:class:`~astropy.table.Table` class in Astropy.

In order to be able to easily exchange data between the :class:`~astropy.table.Table` class and the pandas `DataFrame`_ class (the main data structure in pandas), the :class:`~astropy.table.Table` class includes two methods, :meth:`~astropy.table.Table.to_pandas` and :meth:`~astropy.table.Table.from_pandas`.

To demonstrate these, we can create a simple table::

    >>> from astropy.table import Table
    >>> t = Table()
    >>> t['a'] = [1, 2, 3, 4]
    >>> t['b'] = ['a', 'b', 'c', 'd']
    
which we can then convert to a pandas `DataFrame`_::

    >>> df = t.to_pandas()
    >>> df
       a  b
    0  1  a
    1  2  b
    2  3  c
    3  4  d    
    >>> type(df)
    <class 'pandas.core.frame.DataFrame'>

It is also possible to create a table from a `DataFrame`_::

    >>> t2 = Table.from_pandas(df)
    >>> t2
    <Table masked=False length=4>
      a      b
    int64 string8
    ----- -------
        1       a
        2       b
        3       c
        4       d
        
The conversions to/from pandas are subject to the following caveats:

* The pandas `DataFrame`_ structure does not support multi-dimensional
  columns, so :class:`~astropy.table.Table` objects with multi-dimensional
  columns cannot be converted to `DataFrame`_.

* Masked tables can be converted, but `DataFrame`_ uses ``numpy.nan`` to
  indicate masked values, so all numerical columns (integer or float) are
  converted to ``numpy.float`` columns in `DataFrame`_, and string columns with
  missing values are converted to object columns with ``numpy.nan`` values to
  indicate missing values. For numerical columns, the conversion therefore does
  not necessarily round-trip if converting back to an Astropy table, because the
  distinction between ``numpy.nan`` and masked values is lost, and the different
  for example integer columns will be converted to floating-point.
  
* Tables with mixin columns can currently not be converted, but this may be
  implemented in the future.

.. _DataFrame: http://pandas.pydata.org/pandas-docs/dev/generated/pandas.DataFrame.html
