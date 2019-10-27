.. _timeseries-times:

Converting between different time representations
*************************************************

.. |Time| replace:: :class:`~astropy.time.Time`
.. |TimeDelta| replace:: :class:`~astropy.time.TimeDelta`
.. |TimeSeries| replace:: :class:`~astropy.timeseries.TimeSeries`
.. |BinnedTimeSeries| replace:: :class:`~astropy.timeseries.BinnedTimeSeries`

In :ref:`timeseries-accessing-times`, we saw how to access the time
columns/attributes of the |TimeSeries| and |BinnedTimeSeries| classes. Here we
look in more detail at how to manipulate the resulting times.

Converting times
================

Since the time column in time series is always a |Time| object, it is possible to use the
usual attributes on |Time| to convert the time to different formats or scales.
For example, to get the times as modified Julian Dates from a simple time series::

    >>> from astropy import units as u
    >>> from astropy.timeseries import TimeSeries
    >>> ts = TimeSeries(time_start='2016-03-22T12:30:31', time_delta=3 * u.s,
    ...                 data={'flux': [1., 3., 4., 2., 4.]})
    >>> ts.time.mjd  # doctest: +FLOAT_CMP
    array([57469.52119213, 57469.52122685, 57469.52126157, 57469.5212963 ,
           57469.52133102])

or to convert the times to the Temps Atomique International (TAI) scale

    >>> ts.time.tai
    <Time object: scale='tai' format='isot' value=['2016-03-22T12:31:07.000' '2016-03-22T12:31:10.000'
     '2016-03-22T12:31:13.000' '2016-03-22T12:31:16.000'
     '2016-03-22T12:31:19.000']>

To find the current time scale of the data, you can do::

    >>> ts.time.scale
    'utc'

See :ref:`astropy-time` for more documentation on how to access and convert
times.

Formatting times
================

Since the various time columns are |Time| objects, the default format and scale
to use for the display of the time series can be changed using the ``format``
and ``scale`` attributes::

    >>> ts.time.format = 'isot'
    >>> ts
    <TimeSeries length=5>
              time            flux
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:34.000     3.0
    2016-03-22T12:30:37.000     4.0
    2016-03-22T12:30:40.000     2.0
    2016-03-22T12:30:43.000     4.0
    >>> ts.time.format = 'unix'
    >>> ts  # doctest: +FLOAT_CMP
    <TimeSeries length=5>
        time       flux
       object    float64
    ------------ -------
    1458649831.0     1.0
    1458649834.0     3.0
    1458649837.0     4.0
    1458649840.0     2.0
    1458649843.0     4.0

Times relative to other times
=============================

In some cases, it can be useful to use relative rather than absolute times.
This can be done by using the |TimeDelta| class instead of the |Time| class,
for example by subtracting a reference time from an existing time object::

    >>> ts_rel = TimeSeries(time=ts.time - ts.time[0])
    >>> ts_rel  # doctest: +FLOAT_CMP
    <TimeSeries length=5>
             time
            object
    ----------------------
                       0.0
     3.472222222222765e-05
      6.94444444444553e-05
    0.00010416666666657193
    0.00013888888888879958

The |TimeDelta| values can be converted to a different time unit (e.g., second) using::

    >>> ts_rel.time.to('second')
    <Quantity [ 0.,  3.,  6.,  9., 12.] s>
