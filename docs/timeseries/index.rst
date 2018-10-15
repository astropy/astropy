.. _astropy-timeseries:

****************************************************
Time series (`astropy.timeseries`)
****************************************************

.. |Time| replace:: :class:`~astropy.time.Time`
.. |Table| replace:: :class:`~astropy.table.Table`
.. |QTable| replace:: :class:`~astropy.table.QTable`
.. |Quantity| replace:: :class:`~astropy.units.Quantity`
.. |Longitude| replace:: :class:`~astropy.coordinates.Longitude`
.. |EarthLocation| replace:: :class:`~astropy.coordinates.EarthLocation`
.. |SkyCoord| replace:: :class:`~astropy.coordinates.SkyCoord`
.. |SampledTimeSeries| replace:: :class:`~astropy.timeseries.SampledTimeSeries`
.. |BinnedTimeSeries| replace:: :class:`~astropy.timeseries.BinnedTimeSeries`

Introduction
============

Many different areas of astrophysics have to deal with 1D time series data,
either sampling a continuous variable at fixed times or counting some events
binned into time windows. The `astropy.timeseries` package therefore provides
classes to represent and manipulate time series, including the following
functionality:

* Extending time series with extra rows
* Concatenating multiple time series objects
* Sorting
* Slicing / selecting time ranges (indexing)
* Re-binning and re-sampling time series
* Interpolating to different time stamps
* Masking
* Support for subtraction and addition (e.g. background)

While this functionality is also found in non-domain specific packages such as
`pandas <https://pandas.pydata.org/>`_, the time series classes here also
provide some functionality which is more specific to Astronomy:

* Converting between time systems
* Astropy unit support
* Support for variable width time bins.

The time series classes presented below are |QTable| objects that have special
columns to represent times using the |Time| class. Therefore, most of the
functionality described in :ref:`astropy-table` applies here.

Getting Started
===============

Initializing a sampled time series
----------------------------------

The first type of time series that we will look at here is |SampledTimeSeries|,
which can be used for a time series which samples a continuous variable at
discrete and instantaneous times. Initializing a |SampledTimeSeries| can be done
in the same ways as initializing a |Table| object (see :ref:`Data Tables <astropy-table>`),
but additional arguments related to the times should be specified.

Evenly sampled time series
^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
---------------------------------

The |BinnedTimeSeries| can be used to represent time series where each entry
corresponds to measurements taken over a range in time - for example a light
curve constructed by binning X-ray photon events. This class supports equal-size
or uneven bins, and contiguous and non-contiguous bins. As for
|SampledTimeSeries|, initializing a |BinnedTimeSeries| can be done in the same
ways as initializing a |Table| object (see :ref:`Data Tables <astropy-table>`), but additional
arguments related to the times should be specified as described below.

Equal-sized contiguous bins
^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^

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
------------------------------

The above examples show how to initialize time series objects, but these don't
include any data aside from the times. There are different ways of adding data,
as for the |Table| class.

Adding data after initalization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^

.. warning:: Doesn't work yet, see https://github.com/astropy/astropy/issues/7894

Accessing data
--------------

As for |Table|, columns can be accessed by name::

    >>> ts = SampledTimeSeries(time='2016-03-22T12:30:31',
    ...                        time_delta=3 * u.s,
    ...                        data={'flux': [1., 4., 5., 3., 2.]})
    >>> ts['flux']
    <Column name='flux' dtype='float64' length=5>
    1.0
    4.0
    5.0
    3.0
    2.0
    >>> ts['time']
    <Time object: scale='utc' format='isot' value=['2016-03-22T12:30:31.000' '2016-03-22T12:30:34.000'
     '2016-03-22T12:30:37.000' '2016-03-22T12:30:40.000'
     '2016-03-22T12:30:43.000']>

and rows can be accessed by index::

    >>> ts[0]
    <Row index=0>
              time            flux
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0

Accessing individual values can then be done either by accessing a column then a
row, or vice-versa::

    >>> ts[0]['time']
    <Time object: scale='utc' format='isot' value=2016-03-22T12:30:31.000>

    >>> ts['time'][0]
    <Time object: scale='utc' format='isot' value=2016-03-22T12:30:31.000>

Accessing times
---------------

The ``time`` column (for |SampledTimeSeries|) and the ``start_time``
and ``end_time`` columns (for |BinnedTimeSeries|) can be accessed using the
regular column access notation, as shown in `Accessing data`_, but they can also
be accessed with convenience attributes ``time``, ``start_time``, and
``end_time``::

    >>> ts.time
    <Time object: scale='utc' format='isot' value=['2016-03-22T12:30:31.000' '2016-03-22T12:30:34.000'
     '2016-03-22T12:30:37.000' '2016-03-22T12:30:40.000'
     '2016-03-22T12:30:43.000']>

Since these columns/attributes are |Time| arrays, it is possible to use the
usual attributes on |Time| to convert the time to different formats or scales::

    >>> ts.time.mjd  # doctest: +FLOAT_CMP
    array([57469.52119213, 57469.52122685, 57469.52126157, 57469.5212963 ,
           57469.52133102])

    >>> ts.time.tai
    <Time object: scale='tai' format='isot' value=['2016-03-22T12:31:07.000' '2016-03-22T12:31:10.000'
     '2016-03-22T12:31:13.000' '2016-03-22T12:31:16.000'
     '2016-03-22T12:31:19.000']>

Formatting times
----------------

Since the various time columns are |Time| objects, the default format and scale
to use for the display of the time series can be changed using the ``format``
and ``scale`` attributes::

    >>> ts.time.format = 'isot'
    >>> ts
    <SampledTimeSeries length=5>
              time            flux
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:34.000     4.0
    2016-03-22T12:30:37.000     5.0
    2016-03-22T12:30:40.000     3.0
    2016-03-22T12:30:43.000     2.0
    >>> ts.time.format = 'unix'
    >>> ts
    <SampledTimeSeries length=5>
           time          flux
          object       float64
    ------------------ -------
          1458649831.0     1.0
          1458649834.0     4.0
    1458649837.0000002     5.0
    1458649840.0000002     3.0
          1458649843.0     2.0

Combining time series
---------------------

The  :func:`~astropy.table.vstack`, and :func:`~astropy.table.hstack` functions
from the :mod:`astropy.table` module can be used to stack time series in
different ways.

Time series can be stacked 'vertically' or row-wise using the
:func:`~astropy.table.vstack` function (although note that sampled time
series cannot be combined with binned time series and vice-versa)::

    >>> from astropy.table import vstack
    >>> ts_a = SampledTimeSeries(time='2016-03-22T12:30:31',
    ...                          time_delta=3 * u.s,
    ...                          data={'flux': [1, 4, 5, 3, 2] * u.mJy})
    >>> ts_b = SampledTimeSeries(time='2016-03-22T12:50:31',
    ...                          time_delta=3 * u.s,
    ...                          data={'flux': [4, 3, 1, 2, 3] * u.mJy})
    >>> ts_ab = vstack([ts_a, ts_b])
    >>> ts_ab
    <SampledTimeSeries length=10>
              time            flux
                              mJy
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:34.000     4.0
    2016-03-22T12:30:37.000     5.0
    2016-03-22T12:30:40.000     3.0
    2016-03-22T12:30:43.000     2.0
    2016-03-22T12:50:31.000     4.0
    2016-03-22T12:50:34.000     3.0
    2016-03-22T12:50:37.000     1.0
    2016-03-22T12:50:40.000     2.0
    2016-03-22T12:50:43.000     3.0

Time series can also be combined 'horizontally' or column-wise with other tables
using the :func:`~astropy.table.hstack` function, though these should not be
time series (as having multiple time columns would be confusing)::

    >>> from astropy.table import Table, hstack
    >>> data = Table(data={'temperature': [40., 41., 40., 39., 30.] * u.K})
    >>> ts_a_data = hstack([ts_a, data])
    >>> ts_a_data
    <SampledTimeSeries length=5>
              time            flux  temperature
                              mJy          K
             object         float64    float64
    ----------------------- ------- -----------
    2016-03-22T12:30:31.000     1.0        40.0
    2016-03-22T12:30:34.000     4.0        41.0
    2016-03-22T12:30:37.000     5.0        40.0
    2016-03-22T12:30:40.000     3.0        39.0
    2016-03-22T12:30:43.000     2.0        30.0

Sorting time series
-------------------

Sorting time series in-place can be done using the
:meth:`~astropy.table.Table.sort` method, as for |Table|::

    >>> ts9 = SampledTimeSeries(time='2016-03-22T12:30:31',
    ...                         time_delta=3 * u.s,
    ...                         data={'flux': [1., 4., 5., 3., 2.]})
    >>> ts9
    <SampledTimeSeries length=5>
              time            flux
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:34.000     4.0
    2016-03-22T12:30:37.000     5.0
    2016-03-22T12:30:40.000     3.0
    2016-03-22T12:30:43.000     2.0
    >>> ts9.sort('flux')
    >>> ts9
    <SampledTimeSeries length=5>
              time            flux
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:43.000     2.0
    2016-03-22T12:30:40.000     3.0
    2016-03-22T12:30:34.000     4.0
    2016-03-22T12:30:37.000     5.0

Extracting a subset of columns
------------------------------

Let's consider a case where a time series has two data columns::

   >>> from collections import OrderedDict
   >>> ts = SampledTimeSeries(time='2016-03-22T12:30:31',
   ...                        time_delta=3 * u.s,
   ...                        data=OrderedDict(flux=[1., 4., 5., 3., 2.],
   ...                                         temp=[40., 41., 39., 24., 20.]))

We can create a new time series with just the flux column by doing::

    >>> ts['time', 'flux']
    <SampledTimeSeries length=5>
              time            flux
             object         float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:34.000     4.0
    2016-03-22T12:30:37.000     5.0
    2016-03-22T12:30:40.000     3.0
    2016-03-22T12:30:43.000     2.0

And we can also create a plain |QTable| by extracting just the ``flux`` and
``temp`` columns::

    >>> ts['flux', 'temp']
    <QTable length=5>
      flux    temp
    float64 float64
    ------- -------
        1.0    40.0
        4.0    41.0
        5.0    39.0
        3.0    24.0
        2.0    20.0

Extracting a subset of rows
---------------------------

Time series objects can be sliced by row index, using the same syntax as for
|Time|, e.g.::

    >>> ts[0:2]
    <SampledTimeSeries length=2>
              time            flux    temp
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
    >>> ts.loc[Time('2016-03-22T12:30:31.000')] # doctest: +SKIP
    <Row index=0>
              time            flux    temp
             object         float64 float64
    ----------------------- ------- -------
    2016-03-22T12:30:31.000     1.0    40.0

or within a time range::

    >>> ts.loc[Time('2016-03-22T12:30:31'):Time('2016-03-22T12:30:40')]
    <SampledTimeSeries length=4>
              time            flux    temp
             object         float64 float64
    ----------------------- ------- -------
    2016-03-22T12:30:31.000     1.0    40.0
    2016-03-22T12:30:34.000     4.0    41.0
    2016-03-22T12:30:37.000     5.0    39.0
    2016-03-22T12:30:40.000     3.0    24.0

.. TODO: make it so that Time() is not required above

Note that the result will always be sorted by time. Similarly, the
:attr:`~astropy.timeseries.TimeSeries.iloc` attribute can be used to fetch
rows from the time series *sorted by time*, so for example the two first
entries (by time) can be accessed with::

    >>> ts.iloc[0:2]
    <SampledTimeSeries length=2>
              time            flux    temp
             object         float64 float64
    ----------------------- ------- -------
    2016-03-22T12:30:31.000     1.0    40.0
    2016-03-22T12:30:34.000     4.0    41.0

Reading/writing time series
---------------------------

Since |SampledTimeSeries| and |BinnedTimeSeries| are sub-classes of |QTable|,
they have :meth:`~astropy.table.Table.read` and
:meth:`~astropy.table.Table.write` methods that can be used to read time series
from files. At the moment only a few formats are defined in astropy itself, but
it is easy for other packages and users to define their own readers/writers
with the :ref:`table_io`. Here is an example of using Kepler FITS time series:

.. plot::
   :include-source:
   :context:

   import matplotlib.pyplot as plt
   from astropy.timeseries import SampledTimeSeries
   from astropy.utils.data import get_pkg_data_filename

   example_data = get_pkg_data_filename('kepler/kplr010666592-2009131110544_slc.fits')
   kepler = SampledTimeSeries.read(example_data, format='kepler.fits')

   plt.plot(kepler.time.jd, kepler['sap_flux'])
   plt.xlabel('Barycentric Julian Date')
   plt.ylabel('SAP Flux (e-/s)')

Resampling
----------

The |SampledTimeSeries| class has a
:meth:`~astropy.timeseries.SampledTimeSeries.downsample` method that can be used
to bin values from the time series into bins of equal time, using a custom
function (mean, median, etc.). This operation returns a |BinnedTimeSeries|.
The following example shows how to use this to bin a light curve from the Kepler
mission into 20 minute bins using a median function:

.. plot::
   :context: reset
   :nofigs:

   import matplotlib.pyplot as plt
   from astropy.timeseries import SampledTimeSeries
   from astropy.utils.data import get_pkg_data_filename

   example_data = get_pkg_data_filename('kepler/kplr010666592-2009131110544_slc.fits')
   kepler = SampledTimeSeries.read(example_data, format='kepler.fits')

.. plot::
   :include-source:
   :context:

    import numpy as np
    from astropy import units as u

    kepler_binned = kepler.downsample(bin_size=20 * u.min, func=np.nanmedian)

    plt.plot(kepler.time.jd, kepler['sap_flux'], 'k.')
    plt.plot(kepler_binned.start_time.jd, kepler_binned['sap_flux'], 'r-', drawstyle='steps-pre')
    plt.xlabel('Barycentric Julian Date')
    plt.ylabel('SAP Flux (e-/s)')

Folding
-------

The |SampledTimeSeries| class has a
:meth:`~astropy.timeseries.SampledTimeSeries.fold` method that can be used to
return a new time series with a relative and folded time axis. This method
takes the period as a :class:`~astropy.units.Quantity`, and optionally takes
an epoch as a :class:`~astropy.time.Time`, which defines a zero time offset:

.. plot::
   :context: reset
   :nofigs:

   import numpy as np
   from astropy import units as u
   import matplotlib.pyplot as plt
   from astropy.timeseries import SampledTimeSeries
   from astropy.utils.data import get_pkg_data_filename

   example_data = get_pkg_data_filename('kepler/kplr010666592-2009131110544_slc.fits')
   kepler = SampledTimeSeries.read(example_data, format='kepler.fits')

.. plot::
   :include-source:
   :context:

    kepler_folded = kepler.fold(period=2.2 * u.day, midpoint_epoch='2009-05-02T20:53:40')

    plt.plot(kepler_folded.time.jd, kepler_folded['sap_flux'], 'k.', markersize=1)
    plt.xlabel('Time from midpoint epoch (days)')
    plt.ylabel('SAP Flux (e-/s)')

Arithmetic
----------

Since time series objects are sub-classes of |QTable|, they naturally support
arithmetic on any of the data columns. As an example, we can take the folded
Kepler time series we have seen in the examples above, and normalize it to the
sigma-clipped median value.

.. plot::
   :context: reset
   :nofigs:

   import numpy as np
   from astropy import units as u
   import matplotlib.pyplot as plt
   from astropy.timeseries import SampledTimeSeries
   from astropy.utils.data import get_pkg_data_filename

   example_data = get_pkg_data_filename('kepler/kplr010666592-2009131110544_slc.fits')
   kepler = SampledTimeSeries.read(example_data, format='kepler.fits')
   kepler_folded = kepler.fold(period=2.2 * u.day, midpoint_epoch='2009-05-02T20:53:40')

.. plot::
   :include-source:
   :context:

    from astropy.stats import sigma_clipped_stats

    mean, median, stddev = sigma_clipped_stats(kepler_folded['sap_flux'])

    kepler_folded['sap_flux_norm'] = kepler_folded['sap_flux'] / median

    plt.plot(kepler_folded.time.jd, kepler_folded['sap_flux_norm'], 'k.', markersize=1)
    plt.xlabel('Time from midpoint epoch (days)')
    plt.ylabel('Normalized flux')


Reference/API
=============

.. automodapi:: astropy.timeseries
   :inherited-members:
