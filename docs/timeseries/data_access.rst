.. _timeseries-data-access:

Accessing data in time series
*****************************

.. |Time| replace:: :class:`~astropy.time.Time`
.. |Table| replace:: :class:`~astropy.table.Table`
.. |QTable| replace:: :class:`~astropy.table.QTable`
.. |TimeSeries| replace:: :class:`~astropy.timeseries.TimeSeries`
.. |BinnedTimeSeries| replace:: :class:`~astropy.timeseries.BinnedTimeSeries`
.. |time_attr| replace:: :attr:`~astropy.timeseries.TimeSeries.time`
.. |time_bin_start| replace:: :attr:`~astropy.timeseries.BinnedTimeSeries.time_bin_start`
.. |time_bin_center| replace:: :attr:`~astropy.timeseries.BinnedTimeSeries.time_bin_center`
.. |time_bin_end| replace:: :attr:`~astropy.timeseries.BinnedTimeSeries.time_bin_end`
.. |time_bin_size| replace:: :attr:`~astropy.timeseries.BinnedTimeSeries.time_bin_size`

Accessing data
==============

For the examples in this page, we will consider a simple sampled time series
with two data columns - ``flux`` and ``temp``::

    >>> from astropy import units as u
    >>> from astropy.timeseries import TimeSeries
    >>> ts = TimeSeries(time_start='2016-03-22T12:30:31',
    ...                 time_delta=3 * u.s,
    ...                 data={'flux': [1., 4., 5., 3., 2.] * u.Jy,
    ...                       'temp': [40., 41., 39., 24., 20.] * u.K},
    ...                 names=('flux', 'temp'))

As for |Table|, columns can be accessed by name::

    >>> ts['flux']  # doctest: +FLOAT_CMP
    <Quantity [ 1., 4., 5., 3., 2.] Jy>
    >>> ts['time']
    <Time object: scale='utc' format='isot' value=['2016-03-22T12:30:31.000' '2016-03-22T12:30:34.000'
     '2016-03-22T12:30:37.000' '2016-03-22T12:30:40.000'
     '2016-03-22T12:30:43.000']>

and rows can be accessed by index::

    >>> ts[0]
    <Row index=0>
              time            flux    temp
                               Jy      K
             object         float64 float64
    ----------------------- ------- -------
    2016-03-22T12:30:31.000     1.0    40.0

Accessing individual values can then be done either by accessing a column then a
row, or vice-versa::

    >>> ts[0]['flux']  # doctest: +FLOAT_CMP
    <Quantity 1. Jy>

    >>> ts['temp'][2]  # doctest: +FLOAT_CMP
    <Quantity 39. K>

.. _timeseries-accessing-times:

Accessing times
===============

For |TimeSeries|, the ``time`` column can be accessed using the regular column
access notation, as shown in `Accessing data`_, but they can also be accessed
more conveniently using the |time_attr| attribute::

    >>> ts.time
    <Time object: scale='utc' format='isot' value=['2016-03-22T12:30:31.000' '2016-03-22T12:30:34.000'
     '2016-03-22T12:30:37.000' '2016-03-22T12:30:40.000'
     '2016-03-22T12:30:43.000']>

For |BinnedTimeSeries|, we provide three attributes: |time_bin_start|,
|time_bin_center|, and |time_bin_end|::

    >>> from astropy.timeseries import BinnedTimeSeries
    >>> bts = BinnedTimeSeries(time_bin_start='2016-03-22T12:30:31',
    ...                        time_bin_size=3 * u.s, n_bins=5)
    >>> bts.time_bin_start
    <Time object: scale='utc' format='isot' value=['2016-03-22T12:30:31.000' '2016-03-22T12:30:34.000'
     '2016-03-22T12:30:37.000' '2016-03-22T12:30:40.000'
     '2016-03-22T12:30:43.000']>
    >>> bts.time_bin_center
    <Time object: scale='utc' format='isot' value=['2016-03-22T12:30:32.500' '2016-03-22T12:30:35.500'
     '2016-03-22T12:30:38.500' '2016-03-22T12:30:41.500'
     '2016-03-22T12:30:44.500']>
    >>> bts.time_bin_end
    <Time object: scale='utc' format='isot' value=['2016-03-22T12:30:34.000' '2016-03-22T12:30:37.000'
     '2016-03-22T12:30:40.000' '2016-03-22T12:30:43.000'
     '2016-03-22T12:30:46.000']>

In addition, the |time_bin_size| attribute can be used to access the bin sizes::

    >>> bts.time_bin_size  # doctest: +SKIP
    <Quantity [3., 3., 3., 3., 3.] s>

Note that only |time_bin_start| and |time_bin_size| are available as actual
columns, and |time_bin_center| and |time_bin_end| are computed on-the-fly.

See :ref:`timeseries-times` for more information about changing between
different representations of time.

Extracting a subset of columns
==============================

We can create a new time series with just the ``flux`` column by doing::

   >>> ts['time', 'flux']
   <TimeSeries length=5>
             time            flux
                              Jy
            object         float64
   ----------------------- -------
   2016-03-22T12:30:31.000     1.0
   2016-03-22T12:30:34.000     4.0
   2016-03-22T12:30:37.000     5.0
   2016-03-22T12:30:40.000     3.0
   2016-03-22T12:30:43.000     2.0

Note that the new columns will be copies (not views) of the original columns.
We can also create a plain |QTable| by extracting just the ``flux`` and
``temp`` columns::

   >>> ts['flux', 'temp']
   <QTable length=5>
     flux    temp
       Jy      K
   float64 float64
   ------- -------
       1.0    40.0
       4.0    41.0
       5.0    39.0
       3.0    24.0
       2.0    20.0

Extracting a subset of rows
===========================

Time series objects can be sliced by rows, using the same syntax as for |Time|,
e.g.::

   >>> ts[0:2]
   <TimeSeries length=2>
             time            flux    temp
                              Jy      K
            object         float64 float64
   ----------------------- ------- -------
   2016-03-22T12:30:31.000     1.0    40.0
   2016-03-22T12:30:34.000     4.0    41.0

Time series objects are also automatically indexed using the functionality
described in :ref:`table-indexing`. This provides the ability to access rows and
subset of rows using the :attr:`~astropy.timeseries.TimeSeries.loc` and
:attr:`~astropy.timeseries.TimeSeries.iloc` attributes.

The :attr:`~astropy.timeseries.TimeSeries.loc` attribute can be used to slice
the time series by time. For example, the following can be used to extract all
entries for a given timestamp::

   >>> from astropy.time import Time
   >>> ts.loc[Time('2016-03-22T12:30:31.000')]  # doctest: +SKIP
   <Row index=0>
             time            flux    temp
                              Jy      K
            object         float64 float64
   ----------------------- ------- -------
   2016-03-22T12:30:31.000     1.0    40.0

or within a time range::

   >>> ts.loc['2016-03-22T12:30:30':'2016-03-22T12:30:41']
   <TimeSeries length=4>
             time            flux    temp
                              Jy      K
            object         float64 float64
   ----------------------- ------- -------
   2016-03-22T12:30:31.000     1.0    40.0
   2016-03-22T12:30:34.000     4.0    41.0
   2016-03-22T12:30:37.000     5.0    39.0
   2016-03-22T12:30:40.000     3.0    24.0

Note that in this case we didn't specify |Time| - this isn't needed if the
string is an ISO 8601 time string. Also, as for the |QTable| and |Table| class
``loc`` attribute, and to be consistent with `pandas
<https://pandas.pydata.org/>`_, the last item in the ``loc`` range is inclusive.

Note that the result will always be sorted by time. Similarly, the
:attr:`~astropy.timeseries.TimeSeries.iloc` attribute can be used to fetch
rows from the time series *sorted by time*, so for example the two first
entries (by time) can be accessed with::

   >>> ts.iloc[0:2]
   <TimeSeries length=2>
             time            flux    temp
                              Jy      K
            object         float64 float64
   ----------------------- ------- -------
   2016-03-22T12:30:31.000     1.0    40.0
   2016-03-22T12:30:34.000     4.0    41.0
