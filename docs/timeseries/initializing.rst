.. _timeseries-initializing:

Creating sampled and binned time series
***************************************

.. |Time| replace:: :class:`~astropy.time.Time`
.. |Table| replace:: :class:`~astropy.table.Table`
.. |SampledTimeSeries| replace:: :class:`~astropy.timeseries.SampledTimeSeries`
.. |BinnedTimeSeries| replace:: :class:`~astropy.timeseries.BinnedTimeSeries`

Initializing a sampled time series
==================================

The first type of time series that we will look at here is |SampledTimeSeries|,
which can be used for a time series which samples a continuous variable at
discrete and instantaneous times. Initializing a |SampledTimeSeries| can be done
in the same ways as initializing a |Table| object (see :ref:`Data Tables <astropy-table>`),
but additional arguments related to the times should be specified.

Evenly sampled time series
--------------------------

The easiest way to construct an evenly sampled time series is to specify the
start time, the time interval, and the number of samples, for evenly sampled
time series::

    >>> from astropy import units as u
    >>> from astropy.timeseries import SampledTimeSeries
    >>> ts1 = SampledTimeSeries(time='2016-03-22T12:30:31',
    ...                         time_delta=3 * u.s,
    ...                         n_samples=10)

The ``time`` keyword argument can be set to anything that can be passed to the
|Time| class (see also :ref:`Time and Dates <astropy-time>`). Note that the
``n_samples`` argument is only needed if you are not also passing in data during
initialization (see `Passing data during initialization`_).

Arbitrarily sampled time series
-------------------------------

To construct a sampled time series with samples at arbitrary times, you can
pass multiple times to the ``time`` argument::

    >>> ts2 = SampledTimeSeries(time=['2016-03-22T12:30:31',
    ...                               '2016-03-22T12:30:38',
    ...                               '2016-03-22T12:34:40'])
    >>> ts2
    <SampledTimeSeries length=3>
              time
             object
    -----------------------
    2016-03-22T12:30:31.000
    2016-03-22T12:30:38.000
    2016-03-22T12:34:40.000

You can also specify a vector |Time| object directly.

Initializing a binned time series
=================================

The |BinnedTimeSeries| can be used to represent time series where each entry
corresponds to measurements taken over a range in time - for example a light
curve constructed by binning X-ray photon events. This class supports equal-size
or uneven bins, and contiguous and non-contiguous bins. As for
|SampledTimeSeries|, initializing a |BinnedTimeSeries| can be done in the same
ways as initializing a |Table| object (see :ref:`Data Tables <astropy-table>`), but additional
arguments related to the times should be specified as described below.

Equal-sized contiguous bins
---------------------------

To create a binned time series with equal-size contiguous bins, it is sufficient
to specify a start time as well as a bin size::

    >>> from astropy.timeseries import BinnedTimeSeries
    >>> ts3 = BinnedTimeSeries(start_time='2016-03-22T12:30:31',
    ...                        bin_size=3 * u.s, n_bins=10)
    >>> ts3
    <BinnedTimeSeries length=10>
        start_time       bin_size
                            s
          object         float64
    ----------------------- --------
    2016-03-22T12:30:31.000      3.0
    2016-03-22T12:30:34.000      3.0
    2016-03-22T12:30:37.000      3.0
    2016-03-22T12:30:40.000      3.0
    2016-03-22T12:30:43.000      3.0
    2016-03-22T12:30:46.000      3.0
    2016-03-22T12:30:49.000      3.0
    2016-03-22T12:30:52.000      3.0
    2016-03-22T12:30:55.000      3.0
    2016-03-22T12:30:58.000      3.0

Note that the ``n_bins`` argument is only needed if you are not also passing in
data during initialization (see `Passing data during initialization`_).

Uneven contiguous bins
----------------------

Creating a binned time series with uneven contiguous bins, the bin size can be
changed to give multiple values (note that in this case ``n_bins`` is not
required)::

    >>> ts4 = BinnedTimeSeries(start_time='2016-03-22T12:30:31',
    ...                        bin_size=[3, 3, 2, 3] * u.s)
    >>> ts4
    <BinnedTimeSeries length=4>
        start_time       bin_size
                            s
          object         float64
    ----------------------- --------
    2016-03-22T12:30:31.000      3.0
    2016-03-22T12:30:34.000      3.0
    2016-03-22T12:30:37.000      2.0
    2016-03-22T12:30:39.000      3.0

Alternatively, you can create the same time series by giving an array of start
times as well as a single end time::


    >>> ts5 = BinnedTimeSeries(start_time=['2016-03-22T12:30:31',
    ...                                    '2016-03-22T12:30:34',
    ...                                    '2016-03-22T12:30:37',
    ...                                    '2016-03-22T12:30:39'],
    ...                        end_time='2016-03-22T12:30:42')
    >>> ts5  # doctest: +FLOAT_CMP
    <BinnedTimeSeries length=4>
        start_time            bin_size
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

    >>> ts6 = BinnedTimeSeries(start_time=['2016-03-22T12:30:31',
    ...                                    '2016-03-22T12:30:38',
    ...                                    '2016-03-22T12:34:40'],
    ...                        bin_size=[5, 100, 2]*u.s)
    >>> ts6
    <BinnedTimeSeries length=3>
        start_time       bin_size
                            s
          object         float64
    ----------------------- --------
    2016-03-22T12:30:31.000      5.0
    2016-03-22T12:30:38.000    100.0
    2016-03-22T12:34:40.000      2.0


Or in the most general case, you can also specify multiple times for
``start_time`` and ``end_time``::

    >>> ts7 = BinnedTimeSeries(start_time=['2016-03-22T12:30:31',
    ...                                    '2016-03-22T12:30:33',
    ...                                    '2016-03-22T12:30:40'],
    ...                        end_time=['2016-03-22T12:30:32',
    ...                                  '2016-03-22T12:30:35',
    ...                                  '2016-03-22T12:30:41'])
    >>> ts7  # doctest: +FLOAT_CMP
    <BinnedTimeSeries length=3>
           start_time            bin_size
                                    s
             object              float64
    ----------------------- ------------------
    2016-03-22T12:30:31.000                1.0
    2016-03-22T12:30:33.000                2.0
    2016-03-22T12:30:40.000                1.0

You can also specify vector |Time| objects directly.

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
    >>> ts1['flux'] = [1., 4., 5., 6., 4., 5., 4., 3., 2., 3.] * u.mJy
    >>> ts1
    <SampledTimeSeries length=10>
              time            flux
                              mJy
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:34.000     4.0
    2016-03-22T12:30:37.000     5.0
    2016-03-22T12:30:40.000     6.0
    2016-03-22T12:30:43.000     4.0
    2016-03-22T12:30:46.000     5.0
    2016-03-22T12:30:49.000     4.0
    2016-03-22T12:30:52.000     3.0
    2016-03-22T12:30:55.000     2.0
    2016-03-22T12:30:58.000     3.0

Passing data during initialization
----------------------------------

It is also possible to pass the data during the initialization, as for
|Table|, e.g.::

    >>> ts8 = BinnedTimeSeries(start_time=['2016-03-22T12:30:31',
    ...                                    '2016-03-22T12:30:34',
    ...                                    '2016-03-22T12:30:37',
    ...                                    '2016-03-22T12:30:39'],
    ...                        end_time='2016-03-22T12:30:42',
    ...                        data={'flux': [1., 4., 5., 6.] * u.mJy})
    >>> ts8  # doctest: +FLOAT_CMP
    <BinnedTimeSeries length=4>
           start_time            bin_size       flux
                                    s           mJy
             object              float64      float64
    ----------------------- ----------------- -------
    2016-03-22T12:30:31.000               3.0     1.0
    2016-03-22T12:30:34.000               3.0     4.0
    2016-03-22T12:30:37.000               2.0     5.0
    2016-03-22T12:30:39.000               3.0     6.0

Adding rows
-----------

.. warning:: Doesn't work yet, see https://github.com/astropy/astropy/issues/7894
