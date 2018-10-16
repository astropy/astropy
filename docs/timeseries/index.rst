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
classes to represent and manipulate time series

The time series classes presented below are |QTable| objects that have special
columns to represent times using the |Time| class. Therefore, much of the
functionality described in :ref:`astropy-table` applies here.

Getting Started
===============

In this section, we take a quick look at how to read in a time series, access
the data, and carry out some basic analysis. For more details about creating and
using time series, see the full documentation in :ref:`using-timeseries`.

To start off, we retrieve a FITS file containing a Kepler light curve for a
source::

    >>> from astropy.utils.data import get_pkg_data_filename
    >>> filename = get_pkg_data_filename('kepler/kplr010666592-2009131110544_slc.fits')

We can then use the |SampledTimeSeries| class to read in this file::

    >>> from astropy.timeseries import SampledTimeSeries
    >>> ts = SampledTimeSeries.read(filename, format='kepler.fits')

Time series are specialized kinds of astropy |Table| objects, and in fact if we
look at the ``ts``, it will look like a normal table::

    >>> ts
    <SampledTimeSeries length=14280>
        time             timecorr   ...   pos_corr1      pos_corr2
                            d       ...     pixels         pixels
       object            float32    ...    float32        float32
    ----------------------- ------------- ... -------------- --------------
    2009-05-02T00:41:40.338  6.630610e-04 ...  1.5822421e-03 -1.4463664e-03
    2009-05-02T00:42:39.187  6.630857e-04 ...  1.5743829e-03 -1.4540013e-03
    2009-05-02T00:43:38.045  6.631103e-04 ...  1.5665225e-03 -1.4616371e-03
    2009-05-02T00:44:36.894  6.631350e-04 ...  1.5586632e-03 -1.4692718e-03
    2009-05-02T00:45:35.752  6.631597e-04 ...  1.5508028e-03 -1.4769078e-03
    2009-05-02T00:46:34.601  6.631844e-04 ...  1.5429436e-03 -1.4845425e-03
    2009-05-02T00:47:33.451  6.632091e-04 ...  1.5350844e-03 -1.4921773e-03
    2009-05-02T00:48:32.291  6.632337e-04 ...  1.5272264e-03 -1.4998110e-03
    2009-05-02T00:49:31.149  6.632584e-04 ...  1.5193661e-03 -1.5074468e-03
                  ...           ... ...            ...            ...
    2009-05-11T17:58:22.526  1.014493e-03 ...  3.6121816e-03  3.1950327e-03
    2009-05-11T17:59:21.376  1.014518e-03 ...  3.6102540e-03  3.1872767e-03
    2009-05-11T18:00:20.225  1.014542e-03 ...  3.6083264e-03  3.1795206e-03
    2009-05-11T18:01:19.065  1.014567e-03 ...  3.6063993e-03  3.1717657e-03
    2009-05-11T18:02:17.923  1.014591e-03 ...  3.6044715e-03  3.1640085e-03
    2009-05-11T18:03:16.772  1.014615e-03 ...  3.6025438e-03  3.1562524e-03
    2009-05-11T18:04:15.630  1.014640e-03 ...  3.6006160e-03  3.1484952e-03
    2009-05-11T18:05:14.479  1.014664e-03 ...  3.5986886e-03  3.1407392e-03
    2009-05-11T18:06:13.328  1.014689e-03 ...  3.5967610e-03  3.1329831e-03
    2009-05-11T18:07:12.186  1.014713e-03 ...  3.5948332e-03  3.1252259e-03

As for |Table|, the various columns and rows can be accessed and sliced using
index notation::

    >>> ts['sap_flux']
    <Quantity [1027045.06, 1027184.44, 1027076.25, ..., 1025451.56, 1025468.5 ,
               1025930.9 ] electron / s>

    >>> ts['time', 'sap_flux']
    <SampledTimeSeries length=14280>
              time             sap_flux
                             electron / s
             object            float32
    ----------------------- --------------
    2009-05-02T00:41:40.338  1.0270451e+06
    2009-05-02T00:42:39.187  1.0271844e+06
    2009-05-02T00:43:38.045  1.0270762e+06
    2009-05-02T00:44:36.894  1.0271414e+06
    2009-05-02T00:45:35.752  1.0271569e+06
    2009-05-02T00:46:34.601  1.0272296e+06
    2009-05-02T00:47:33.451  1.0273199e+06
    2009-05-02T00:48:32.291  1.0271497e+06
    2009-05-02T00:49:31.149  1.0271755e+06
                        ...            ...
    2009-05-11T17:58:22.526  1.0234769e+06
    2009-05-11T17:59:21.376  1.0234574e+06
    2009-05-11T18:00:20.225  1.0238128e+06
    2009-05-11T18:01:19.065  1.0243234e+06
    2009-05-11T18:02:17.923  1.0244257e+06
    2009-05-11T18:03:16.772  1.0248654e+06
    2009-05-11T18:04:15.630  1.0250156e+06
    2009-05-11T18:05:14.479  1.0254516e+06
    2009-05-11T18:06:13.328  1.0254685e+06
    2009-05-11T18:07:12.186  1.0259309e+06

    >>> ts[0:10]
    <SampledTimeSeries length=10>
              time             timecorr   ...   pos_corr1      pos_corr2
                                  d       ...     pixels         pixels
             object            float32    ...    float32        float32
    ----------------------- ------------- ... -------------- --------------
    2009-05-02T00:41:40.338  6.630610e-04 ...  1.5822421e-03 -1.4463664e-03
    2009-05-02T00:42:39.187  6.630857e-04 ...  1.5743829e-03 -1.4540013e-03
    2009-05-02T00:43:38.045  6.631103e-04 ...  1.5665225e-03 -1.4616371e-03
    2009-05-02T00:44:36.894  6.631350e-04 ...  1.5586632e-03 -1.4692718e-03
    2009-05-02T00:45:35.752  6.631597e-04 ...  1.5508028e-03 -1.4769078e-03
    2009-05-02T00:46:34.601  6.631844e-04 ...  1.5429436e-03 -1.4845425e-03
    2009-05-02T00:47:33.451  6.632091e-04 ...  1.5350844e-03 -1.4921773e-03
    2009-05-02T00:48:32.291  6.632337e-04 ...  1.5272264e-03 -1.4998110e-03
    2009-05-02T00:49:31.149  6.632584e-04 ...  1.5193661e-03 -1.5074468e-03
    2009-05-02T00:50:29.998  6.632830e-04 ...  1.5115069e-03 -1.5150816e-03

All |SampledTimeSeries| objects have a ``time`` column, which is always the
first column::

    >>> ts['time']
    <Time object: scale='tcb' format='isot' value=['2009-05-02T00:41:40.338' '2009-05-02T00:42:39.187'
     '2009-05-02T00:43:38.045' ... '2009-05-11T18:05:14.479'
     '2009-05-11T18:06:13.328' '2009-05-11T18:07:12.186']>

This column can also be accessed using the ``.time`` attribute::

    >>> ts.time
    <Time object: scale='tcb' format='isot' value=['2009-05-02T00:41:40.338' '2009-05-02T00:42:39.187'
     '2009-05-02T00:43:38.045' ... '2009-05-11T18:05:14.479'
     '2009-05-11T18:06:13.328' '2009-05-11T18:07:12.186']>

and is always a |Time| object (see :ref:`astropy-time`), which therefore
supports the ability to convert to different time scales and formats::

    >>> ts.time.mjd
    array([54953.0289391 , 54953.02962023, 54953.03030145, ...,
           54962.7536398 , 54962.75432093, 54962.75500215])
    >>> ts.time.unix
    array([1.24122482e+09, 1.24122488e+09, 1.24122494e+09, ...,
           1.24206503e+09, 1.24206509e+09, 1.24206515e+09])

Let's use what we've seen so far to make a plot

.. plot::
   :context: reset
   :nofigs:

    from astropy.utils.data import get_pkg_data_filename
    filename = get_pkg_data_filename('kepler/kplr010666592-2009131110544_slc.fits')
    from astropy.timeseries import SampledTimeSeries
    ts = SampledTimeSeries.read(filename, format='kepler.fits')

.. plot::
   :include-source:
   :context:

   >>> import matplotlib.pyplot as plt
   >>> plt.plot(ts.time.jd, ts['sap_flux'])  # doctest: +SKIP
   >>> plt.xlabel('Barycentric Julian Date')  # doctest: +SKIP
   >>> plt.ylabel('SAP Flux (e-/s)')  # doctest: +SKIP

It looks like there are a few transits! Let's use the :ref:`stats-bls`
functionality to estimate the period, using a box with a duration of 0.2 days:

.. plot::
   :context:
   :include-source:
   :nofigs:

   >>> import numpy as np
   >>> from astropy import units as u
   >>> from astropy.stats import BoxLeastSquares
   >>> keep = ~np.isnan(ts['sap_flux'])
   >>> periodogram = BoxLeastSquares(ts.time.jd[keep] * u.day, ts['sap_flux'][keep]).autopower(0.2 * u.day)
   >>> period = periodogram.period[np.argmax(periodogram.power)]
   >>> period
   <Quantity 2.21584977 d>

.. plot::
   :context:
   :nofigs:

   >>> plt.clf()

We can also take a look at the periodogram:

.. plot::
   :context:
   :include-source:

   >>> plt.plot(periodogram.period, periodogram.power)  # doctest: +SKIP
   >>> plt.xlabel('Period (days)')  # doctest: +SKIP
   >>> plt.ylabel('Power')  # doctest: +SKIP

We can now fold the time series using the period we've found above using the
:meth:`~astropy.timeseries.SampledTimeSeries.fold` method:

.. plot::
   :context:
   :nofigs:

   >>> plt.clf()

.. plot::
   :context:
   :include-source:

   >>> ts_folded = ts.fold(period=period)
   >>> plt.plot(ts_folded.time.jd, ts_folded['sap_flux'], '.', markersize=1)  # doctest: +SKIP
   >>> plt.xlabel('Time (days)')  # doctest: +SKIP
   >>> plt.ylabel('SAP Flux (e-/s)')  # doctest: +SKIP

Using the :ref:`astropy-stats` module, we can normalize the flux by sigma-clipping
the data to determine the baseline flux:

.. plot::
   :context:
   :nofigs:

   >>> plt.clf()

.. plot::
   :context:
   :include-source:

   >>> from astropy.stats import sigma_clipped_stats
   >>> mean, median, stddev = sigma_clipped_stats(ts_folded['sap_flux'])
   >>> ts_folded['sap_flux_norm'] = ts_folded['sap_flux'] / median
   >>> plt.plot(ts_folded.time.jd, ts_folded['sap_flux_norm'], '.', markersize=1)  # doctest: +SKIP
   >>> plt.xlabel('Time (days)')  # doctest: +SKIP
   >>> plt.ylabel('Normalized flux')  # doctest: +SKIP

And we can downsample the time series by binning the points into bins of equal
time - this returns a |BinnedTimeSeries|:

.. plot::
   :context:
   :nofigs:

   >>> plt.clf()

.. plot::
   :context:
   :include-source:

   >>> ts_binned = ts_folded.downsample(0.03 * u.day)
   >>> plt.plot(ts_folded.time.jd, ts_folded['sap_flux_norm'], '.', markersize=1)  # doctest: +SKIP
   >>> plt.plot(ts_binned.start_time.jd, ts_binned['sap_flux_norm'], drawstyle='steps-post')  # doctest: +SKIP
   >>> plt.xlabel('Time (days)')  # doctest: +SKIP
   >>> plt.ylabel('Normalized flux')  # doctest: +SKIP

.. _using-timeseries:

Using ``timeseries``
====================

The details of using `astropy.timeseries` are provided in the following sections:

Initializing and reading in time series
---------------------------------------

.. toctree::
   :maxdepth: 2

   initializing.rst
   io.rst

Accessing data and manipulating time series
-------------------------------------------

.. toctree::
   :maxdepth: 2

   data_access.rst
   times.rst
   analyzing.rst

Comparison to other packages
----------------------------

The `astropy.timeseries` package is not the only package to provide
functionality related to time series. For example, another notable package is
`pandas <https://pandas.pydata.org/>`_, which provides a :class:`pandas.Series`
class. The main benefits of `astropy.timeseries` in the context of astronomical
research are the following:

* The time column is a |Time| object that supports very high precision
  representation of times, and makes it easy to convert between different
  time scales and formats (e.g. ISO 8601 timestamps, Julian Dates, and so on).
* The data columns can include |Quantity| objects with units
* The |BinnedTimeSeries| class includes variable width time bins
* There are built-in readers for common time series file formats, as well as
  the ability to define custom readers/writers.

Reference/API
=============

.. automodapi:: astropy.timeseries
   :inherited-members:
