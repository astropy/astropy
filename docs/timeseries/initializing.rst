.. _timeseries-initializing:

Creating time series
***************************************

.. |Time| replace:: :class:`~astropy.time.Time`
.. |TimeDelta| replace:: :class:`~astropy.time.TimeDelta`
.. |Table| replace:: :class:`~astropy.table.Table`
.. |QTable| replace:: :class:`~astropy.table.Table`
.. |TimeSeries| replace:: :class:`~astropy.timeseries.TimeSeries`
.. |BinnedTimeSeries| replace:: :class:`~astropy.timeseries.BinnedTimeSeries`

Initializing a simple time series
=================================

The first type of time series that we will look at here is |TimeSeries|,
which can be used for a time series which samples a continuous variable at
discrete, instantaneous times. Initializing a |TimeSeries| can be done
in the same ways as initializing a |Table| object (see :ref:`Data Tables <astropy-table>`),
but additional arguments related to the times should be specified.

Evenly sampled time series
--------------------------

The easiest way to construct an evenly sampled time series is to specify the
start time, the time interval, and the number of samples, for evenly sampled
time series::

    >>> from astropy import units as u
    >>> from astropy.timeseries import TimeSeries
    >>> ts1 = TimeSeries(time_start='2016-03-22T12:30:31',
    ...                  time_delta=3 * u.s,
    ...                  n_samples=5)
    >>> ts1
    <TimeSeries length=5>
              time
             object
    -----------------------
    2016-03-22T12:30:31.000
    2016-03-22T12:30:34.000
    2016-03-22T12:30:37.000
    2016-03-22T12:30:40.000
    2016-03-22T12:30:43.000

The ``time`` keyword argument can be set to anything that can be passed to the
|Time| class (see also :ref:`Time and Dates <astropy-time>`) or |Time| objects
directly. Note that the ``n_samples`` argument is only needed if you are not
also passing in data during initialization (see `Passing data during
initialization`_).

Arbitrarily sampled time series
-------------------------------

To construct a sampled time series with samples at arbitrary times, you can
pass multiple times to the ``time`` argument::

    >>> ts2 = TimeSeries(time=['2016-03-22T12:30:31',
    ...                        '2016-03-22T12:30:38',
    ...                        '2016-03-22T12:34:40'])
    >>> ts2
    <TimeSeries length=3>
              time
             object
    -----------------------
    2016-03-22T12:30:31.000
    2016-03-22T12:30:38.000
    2016-03-22T12:34:40.000

You can also specify a vector |Time| object directly as the ``time=`` argument,
or a vector |TimeDelta| argument or a quantity array to the ``time_delta=``
argument.::

    >>> TimeSeries(time_start="2011-01-01T00:00:00",
    ...            time_delta=[0.1, 0.2, 0.1, 0.3, 0.2]*u.s)
    <TimeSeries length=5>
              time
            object
    -----------------------
    2011-01-01T00:00:00.000
    2011-01-01T00:00:00.100
    2011-01-01T00:00:00.300
    2011-01-01T00:00:00.400
    2011-01-01T00:00:00.700

Initializing a binned time series
=================================

The |BinnedTimeSeries| can be used to represent time series where each entry
corresponds to measurements taken over a range in time - for example a light
curve constructed by binning X-ray photon events. This class supports equal-size
or uneven bins, and contiguous and non-contiguous bins. As for
|TimeSeries|, initializing a |BinnedTimeSeries| can be done in the same
ways as initializing a |Table| object (see :ref:`Data Tables <astropy-table>`), but additional
arguments related to the times should be specified as described below.

Equal-sized contiguous bins
---------------------------

To create a binned time series with equal-size contiguous bins, it is sufficient
to specify a start time as well as a bin size::

    >>> from astropy.timeseries import BinnedTimeSeries
    >>> ts3 = BinnedTimeSeries(time_bin_start='2016-03-22T12:30:31',
    ...                        time_bin_size=3 * u.s, n_bins=10)
    >>> ts3
    <BinnedTimeSeries length=10>
        time_bin_start     time_bin_size
                                 s
            object            float64
    ----------------------- -------------
    2016-03-22T12:30:31.000           3.0
    2016-03-22T12:30:34.000           3.0
    2016-03-22T12:30:37.000           3.0
    2016-03-22T12:30:40.000           3.0
    2016-03-22T12:30:43.000           3.0
    2016-03-22T12:30:46.000           3.0
    2016-03-22T12:30:49.000           3.0
    2016-03-22T12:30:52.000           3.0
    2016-03-22T12:30:55.000           3.0
    2016-03-22T12:30:58.000           3.0

Note that the ``n_bins`` argument is only needed if you are not also passing in
data during initialization (see `Passing data during initialization`_).

Uneven contiguous bins
----------------------

Creating a binned time series with uneven contiguous bins, the bin size can be
changed to give multiple values (note that in this case ``n_bins`` is not
required)::

    >>> ts4 = BinnedTimeSeries(time_bin_start='2016-03-22T12:30:31',
    ...                        time_bin_size=[3, 3, 2, 3] * u.s)
    >>> ts4
    <BinnedTimeSeries length=4>
         time_bin_start     time_bin_size
                                  s
             object            float64
    ----------------------- -------------
    2016-03-22T12:30:31.000           3.0
    2016-03-22T12:30:34.000           3.0
    2016-03-22T12:30:37.000           2.0
    2016-03-22T12:30:39.000           3.0

Alternatively, you can create the same time series by giving an array of start
times as well as a single end time::

    >>> ts5 = BinnedTimeSeries(time_bin_start=['2016-03-22T12:30:31',
    ...                                        '2016-03-22T12:30:34',
    ...                                        '2016-03-22T12:30:37',
    ...                                        '2016-03-22T12:30:39'],
    ...                        time_bin_end='2016-03-22T12:30:42')
    >>> ts5  # doctest: +FLOAT_CMP
    <BinnedTimeSeries length=4>
        time_bin_start            time_bin_size
                                 s
          object              float64
    ----------------------- -----------------
    2016-03-22T12:30:31.000               3.0
    2016-03-22T12:30:34.000               3.0
    2016-03-22T12:30:37.000               2.0
    2016-03-22T12:30:39.000               3.0

Uneven non-contiguous bins
--------------------------

To create a binned time series with non-contiguous bins, you can either
specify an array of start times and bin widths::

    >>> ts6 = BinnedTimeSeries(time_bin_start=['2016-03-22T12:30:31',
    ...                                        '2016-03-22T12:30:38',
    ...                                        '2016-03-22T12:34:40'],
    ...                        time_bin_size=[5, 100, 2]*u.s)
    >>> ts6
    <BinnedTimeSeries length=3>
         time_bin_start     time_bin_size
                                  s
             object            float64
    ----------------------- -------------
    2016-03-22T12:30:31.000           5.0
    2016-03-22T12:30:38.000         100.0
    2016-03-22T12:34:40.000           2.0

Or in the most general case, you can also specify multiple times for
``time_bin_start`` and ``time_bin_end``::

    >>> ts7 = BinnedTimeSeries(time_bin_start=['2016-03-22T12:30:31',
    ...                                        '2016-03-22T12:30:33',
    ...                                        '2016-03-22T12:30:40'],
    ...                        time_bin_end=['2016-03-22T12:30:32',
    ...                                      '2016-03-22T12:30:35',
    ...                                      '2016-03-22T12:30:41'])
    >>> ts7  # doctest: +FLOAT_CMP
    <BinnedTimeSeries length=3>
        time_bin_start        time_bin_size
                                    s
             object              float64
    ----------------------- ------------------
    2016-03-22T12:30:31.000                1.0
    2016-03-22T12:30:33.000                2.0
    2016-03-22T12:30:40.000                1.0

Adding data to the time series
==============================

The above examples show how to initialize time series objects, but these don't
include any data aside from the times. There are different ways of adding data,
as for the |Table| class.

Adding data after initalization
-------------------------------

Once the time series is initialized, you can add columns/fields to it as you
would for a |Table| object::

    >>> from astropy import units as u
    >>> ts1['flux'] = [1., 4., 5., 6., 4.] * u.mJy
    >>> ts1
    <TimeSeries length=5>
              time            flux
                              mJy
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:34.000     4.0
    2016-03-22T12:30:37.000     5.0
    2016-03-22T12:30:40.000     6.0
    2016-03-22T12:30:43.000     4.0

Passing data during initialization
----------------------------------

It is also possible to pass the data during the initialization, as for
|Table|, e.g.::

    >>> ts8 = BinnedTimeSeries(time_bin_start=['2016-03-22T12:30:31',
    ...                                        '2016-03-22T12:30:34',
    ...                                        '2016-03-22T12:30:37',
    ...                                        '2016-03-22T12:30:39'],
    ...                        time_bin_end='2016-03-22T12:30:42',
    ...                        data={'flux': [1., 4., 5., 6.] * u.mJy})
    >>> ts8  # doctest: +FLOAT_CMP
    <BinnedTimeSeries length=4>
           time_bin_start            time_bin_size       flux
                                    s           mJy
             object              float64      float64
    ----------------------- ----------------- -------
    2016-03-22T12:30:31.000               3.0     1.0
    2016-03-22T12:30:34.000               3.0     4.0
    2016-03-22T12:30:37.000               2.0     5.0
    2016-03-22T12:30:39.000               3.0     6.0

Adding rows
-----------

Adding rows to |TimeSeries| or |BinnedTimeSeries| can be done using the
:meth:`~astropy.table.Table.add_row` method, as for |Table| and |QTable|. This
method takes a dictionary where the keys are column names::

    >>> ts8.add_row({'time_bin_start': '2016-03-22T12:30:44.000',
    ...              'time_bin_size': 2 * u.s,
    ...              'flux': 3 * u.mJy})
    >>> ts8  # doctest: +FLOAT_CMP
    <BinnedTimeSeries length=5>
        time_bin_start       time_bin_size      flux
                                    s           mJy
             object              float64      float64
    ----------------------- ----------------- -------
    2016-03-22T12:30:31.000               3.0     1.0
    2016-03-22T12:30:34.000               3.0     4.0
    2016-03-22T12:30:37.000               2.0     5.0
    2016-03-22T12:30:39.000               3.0     6.0
    2016-03-22T12:30:44.000               2.0     3.0

If you want to be able to skip some values when adding rows, you should make
sure that masking is enabled - see :ref:`timeseries-masking` for more details.
