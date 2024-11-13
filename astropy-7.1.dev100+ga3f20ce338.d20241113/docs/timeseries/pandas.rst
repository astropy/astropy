.. _timeseries-pandas:

Interfacing with the Pandas Package
***********************************

The `astropy.timeseries` package is not the only package to provide
functionality related to time series. Another notable package is `pandas
<https://pandas.pydata.org/>`_, which provides a :class:`pandas.DataFrame`
class. The main benefits of `astropy.timeseries` in the context of astronomical
research are the following:

* The time column is a |Time| object that supports very high precision
  representation of times, and makes it easy to convert between different
  time scales and formats (e.g., ISO 8601 timestamps, Julian Dates, and so on).
* The data columns can include |Quantity| objects with units.
* The |BinnedTimeSeries| class includes variable-width time bins.
* There are built-in readers for common time series file formats, as well as
  the ability to define custom readers/writers.

Nevertheless, there are cases where using pandas :class:`~pandas.DataFrame`
objects might make sense, so we provide methods to convert to/from
:class:`~pandas.DataFrame` objects.

Example
-------

.. EXAMPLE START: Interfacing between Time Series and the Pandas DataFrame

Consider a concise example starting from a :class:`~pandas.DataFrame`:

.. doctest-requires:: pandas

    >>> import pandas
    >>> import numpy as np
    >>> from astropy.utils.introspection import minversion
    >>> df = pandas.DataFrame()
    >>> df['a'] = [1, 2, 3]
    >>> times = np.array(['2015-07-04', '2015-07-05', '2015-07-06'], dtype=np.datetime64)
    >>> df.set_index(pandas.DatetimeIndex(times), inplace=True)
    >>> df
        a
    2015-07-04  1
    2015-07-05  2
    2015-07-06  3

We can convert this to an ``astropy`` |TimeSeries| using
:meth:`~astropy.timeseries.TimeSeries.from_pandas`:

.. doctest-requires:: pandas

    >>> from astropy.timeseries import TimeSeries
    >>> ts = TimeSeries.from_pandas(df)
    >>> ts
    <TimeSeries length=3>
                 time               a
                 Time             int64
    ----------------------------- -----
    2015-07-04T00:00:00.000000000     1
    2015-07-05T00:00:00.000000000     2
    2015-07-06T00:00:00.000000000     3

Converting to :class:`~pandas.DataFrame` can also be done with
:meth:`~astropy.timeseries.TimeSeries.to_pandas`:

.. doctest-requires:: pandas

    >>> ts['b'] = [1.2, 3.4, 5.4]
    >>> df_new = ts.to_pandas()
    >>> df_new
                a    b
    time
    2015-07-04  1  1.2
    2015-07-05  2  3.4
    2015-07-06  3  5.4

Missing values in the time column are supported and correctly converted to a
pandas' NaT object:

.. doctest-requires:: pandas

    >>> ts.time[2] = np.nan
    >>> ts
    <TimeSeries length=3>
                 time               a      b
                 Time             int64 float64
    ----------------------------- ----- -------
    2015-07-04T00:00:00.000000000     1     1.2
    2015-07-05T00:00:00.000000000     2     3.4
                              ———     3     5.4
    >>> df_missing = ts.to_pandas()
    >>> df_missing
               a    b
    time
    2015-07-04  1  1.2
    2015-07-05  2  3.4
    NaT         3  5.4

.. EXAMPLE END
