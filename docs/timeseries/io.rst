.. _timeseries-io:

Reading/writing time series
***************************

.. |Table| replace:: :class:`~astropy.table.Table`
.. |TimeSeries| replace:: :class:`~astropy.timeseries.TimeSeries`
.. |BinnedTimeSeries| replace:: :class:`~astropy.timeseries.BinnedTimeSeries`

Build-in readers
================

Since |TimeSeries| and |BinnedTimeSeries| are sub-classes of |Table|,
they have :meth:`~astropy.table.Table.read` and
:meth:`~astropy.table.Table.write` methods that can be used to read time series
from files. We include a few readers for well-defined formats in `astropy.timeseries` -
for example we have readers for light curves in FITS format from the
`Kepler <https://www.nasa.gov/mission_pages/kepler/main/index.html>`_ and
`TESS <https://tess.gsfc.nasa.gov/>`_ missions.

Here is an example of using Kepler FITS time series - we start off by fetching an example
file:

.. plot::
   :include-source:
   :context: reset
   :nofigs:

   from astropy.utils.data import get_pkg_data_filename
   example_data = get_pkg_data_filename('timeseries/kplr010666592-2009131110544_slc.fits')

.. note::
    The light curve provided here is hand-picked for example purposes.  To get
    other Kepler light curves for science purposes using Python, see the
    `astroquery <https://astroquery.readthedocs.io>`_ affiliated package.

This will set ``example_data`` to the filename of the downloaded file (so you can
replace this by the filename for the file you want to read in). We can then read in
the time series using:

.. plot::
   :include-source:
   :context:
   :nofigs:

   from astropy.timeseries import TimeSeries
   kepler = TimeSeries.read(example_data, format='kepler.fits')

Let's check that the time series has been read in correctly:

.. plot::
   :include-source:
   :context:

   import matplotlib.pyplot as plt

   plt.plot(kepler.time.jd, kepler['sap_flux'], 'k.', markersize=1)
   plt.xlabel('Julian Date')
   plt.ylabel('SAP Flux (e-/s)')

Reading other formats
=====================

At the moment only a few formats are defined in ``astropy`` itself, in part because
there are not many well documented formats for storing time series. So in many
cases, you will likely have to first read in your files using the more
generic |Table| class (see :ref:`read_write_tables`). In fact, the
:meth:`TimeSeries.read <astropy.timeseries.TimeSeries.read>` and
:meth:`BinnedTimeSeries.read <astropy.timeseries.BinnedTimeSeries.read>` methods
can do this behind the scenes - if the table cannot be read by any of the time
series readers, these methods will try and use some of the default :class:`~astropy.table.Table`
readers and then require users to specify the name of the important columns.

For example, if you are reading in a file called :download:`sampled.csv
<sampled.csv>` where the time column is called ``Date`` and is an ISO string,
you can do::

    >>> from astropy.timeseries import TimeSeries
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> sampled_filename = get_pkg_data_filename('data/sampled.csv',
    ...                                          package='astropy.timeseries.tests')
    >>> ts = TimeSeries.read(sampled_filename, format='ascii.csv',
    ...                      time_column='Date')
    >>> ts[:3]
    <TimeSeries length=3>
              time             A       B       C       D       E       F       G
             object         float64 float64 float64 float64 float64 float64 float64
    ----------------------- ------- ------- ------- ------- ------- ------- -------
    2008-03-18 00:00:00.000   24.68  164.93  114.73   26.27   19.21   28.87   63.44
    2008-03-19 00:00:00.000   24.18  164.89  114.75   26.22   19.07   27.76   59.98
    2008-03-20 00:00:00.000   23.99  164.63  115.04   25.78   19.01   27.04   59.61

If you are reading in a binned time series from a file called
:download:`binned.csv <binned.csv>` and with a column ``time_start`` giving the start time
and ``bin_size`` giving the size of each bin, you can do::

    >>> from astropy import units as u
    >>> from astropy.timeseries import BinnedTimeSeries
    >>> binned_filename = get_pkg_data_filename('data/binned.csv',
    ...                                          package='astropy.timeseries.tests')
    >>> ts = BinnedTimeSeries.read(binned_filename, format='ascii.csv',
    ...                            time_bin_start_column='time_start',
    ...                            time_bin_size_column='bin_size',
    ...                            time_bin_size_unit=u.s)
    >>> ts[:3]
    <BinnedTimeSeries length=3>
         time_bin_start     time_bin_size ...    E       F
                                  s       ...
             object            float64    ... float64 float64
    ----------------------- ------------- ... ------- -------
    2016-03-22T12:30:31.000           3.0 ...   28.87   63.44
    2016-03-22T12:30:34.000           3.0 ...   27.76   59.98
    2016-03-22T12:30:37.000           3.0 ...   27.04   59.61

See the documentation for :meth:`TimeSeries.read
<astropy.timeseries.TimeSeries.read>` and :meth:`BinnedTimeSeries.read
<astropy.timeseries.BinnedTimeSeries.read>` for more details.

Alternatively, you can read in the table using your own code then construct the
time series object as described in :ref:`timeseries-initializing`, although then
you cannot write out another time series in the same format.

If you have written a reader/writer for a commonly used format, please feel free
to contribute it to ``astropy``!
