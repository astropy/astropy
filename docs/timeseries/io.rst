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
`Kepler <https://www.nasa.gov/mission_pages/timeseries/main/index.html>`_ and
`TESS <https://tess.gsfc.nasa.gov/>`_ missions.

Here is an example of using Kepler FITS time series - we start off by fetching an example
file:

.. plot::
   :include-source:
   :context: reset
   :nofigs:

   from astropy.utils.data import get_pkg_data_filename
   example_data = get_pkg_data_filename('timeseries/kplr010666592-2009131110544_slc.fits')

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
   plt.xlabel('Barycentric Julian Date')
   plt.ylabel('SAP Flux (e-/s)')

Reading other formats
=====================

At the moment only a few formats are defined in astropy itself, in part because
there are not many well documented formats for storing time series. So in many cases,
you will likely have to first read in your files using the more generic |Table|
class (see :ref:`read_write_tables`), and then construct the time series object as
described in :ref:`timeseries-initializing`.

If you have written a reader/writer for a commonly used format, and it is well
documented, please feel free to contribute it to astropy!
