.. _utils-iers:

************************************************
IERS data access (`astropy.utils.iers`)
************************************************

Introduction
============

The `~astropy.utils.iers` package provides access to the tables provided by
the International Earth Rotation and Reference Systems (IERS) service, in
particular files allowing interpolation of published UT1-UTC and polar motion
values for given times.  The UT1-UTC values are used in `astropy.time` to
provide UT1 values, and the polar motions are used in `astropy.coordinates` to
determine Earth orientation for celestial-to-terrestrial coordinate
transformations.

.. note:: The package also provides machinery to track leap seconds.  Since it
          generally should not be necessary to deal with those by hand, this
          is not discussed below.  For details, see the documentation of
          `~astropy.utils.iers.LeapSeconds`.

Getting started
===============

Starting with astropy 1.2, the latest IERS values (which include approximately
one year of predictive values) are automatically downloaded from the IERS
service when required.  This happens when a time or coordinate transformation
needs a value which is not already available via the download cache.  In most
cases there is no need for invoking the `~astropy.utils.iers` classes oneself,
but it is useful to understand the situations when a download will occur
and how this can be controlled.

Basic usage
-----------

By default, the IERS data are managed via instances of the
:class:`~astropy.utils.iers.IERS_Auto` class.  These instances are created
internally within the relevant time and coordinate objects during
transformations.  If the astropy data cache does not have the required IERS
data file then astropy will request the file from the IERS service.  This will
occur the first time such a transform is done for a new setup or on a new
machine.  Here is an example that shows the typical download progress bar::

  >>> from astropy.time import Time
  >>> t = Time('2016:001')
  >>> t.ut1  # doctest: +SKIP
  Downloading https://maia.usno.navy.mil/ser7/finals2000A.all
  |==================================================================| 3.0M/3.0M (100.00%)         6s
  <Time object: scale='ut1' format='yday' value=2016:001:00:00:00.082>

Note that you can forcibly clear the download cache as follows::

  >>> from astropy.utils.data import clear_download_cache
  >>> clear_download_cache()  # doctest: +SKIP

The default IERS data used automatically is updated by the service every 7 days
and includes transforms dating back to 1973-01-01.

.. note:: The :class:`~astropy.utils.iers.IERS_Auto` class contains machinery
    to ensure that the IERS table is kept up to date by auto-downloading the
    latest version as needed.  This means that the IERS table is assured of
    having the state-of-the-art definitive and predictive values for Earth
    rotation.  As a user it is **your responsibility** to understand the
    accuracy of IERS predictions if your science depends on that.  If you
    request ``UT1-UTC`` or polar motions for times beyond the range of IERS
    table data then the nearest available values will be provided.


Configuration parameters
------------------------

There are three configuration parameters that control the behavior
of the automatic IERS downloading:

  auto_download:
    Enable auto-downloading of the latest IERS data.  If set to ``False`` then
    the local IERS-B file will be used by default (even if the full IERS file
    with predictions was already downloaded and cached).  This replicates the
    behavior prior to astropy 1.2.  (default=True)

  auto_max_age:
    Maximum age of predictive data before auto-downloading (days).  See
    next section for details. (default=30)

  iers_auto_url:
    URL for auto-downloading IERS file data

  iers_auto_url_mirror:
    Mirror URL for auto-downloading IERS file data.

  remote_timeout:
    Remote timeout downloading IERS file data (seconds)

Auto refresh behavior
---------------------

The first time that one attempts a time or coordinate transformation that
requires IERS data, the latest version of the IERS table (from 1973 through
one year into the future) will be downloaded and stored in the astropy cache.

Transformations will then use the cached data file if possible.  However, the
``IERS_Auto`` table is automatically updated in place from the network if the
following two conditions a met when the table is queried for ``UT1-UTC`` or
polar motion values:

- Any of the requested IERS values are *predictive*, meaning that they have
  been extrapolated into the future with a model that is fit to measured data.
  The IERS table contains approximately one year of predictive data from the
  time it is created.
- The first predictive values in the table are at least ``conf.auto_max_age
  days`` old relative to the current actual time (i.e. ``Time.now()``).  This
  means that the IERS table is out of date and a newer version can be found on
  the IERS service.

The IERS Service provides the default online table
(set by ``astropy.utils.iers.IERS_A_URL``) and updates the content
once each 7 days.  The default value of ``auto_max_age`` is 30 days to avoid
unnecessary network access, but one can reduce this to as low as 10 days.

.. _iers-working-offline:

Working offline
---------------

If you are working without an internet connection and doing transformations
that require IERS data, there are a couple of options.

**Disable auto downloading**

Here you can do::

  >>> from astropy.utils import iers
  >>> iers.conf.auto_download = False  # doctest: +SKIP

In this case any transforms will use the bundled IERS-B data which covers
the time range from 1962 to just before the astropy release date.  Any
transforms outside of this range will not be allowed.

**Set the auto-download max age parameter**

*Only do this if you understand what you are doing, THIS CAN GIVE INACCURATE
ANSWERS!* Assuming you have previously been connected to the internet and have
downloaded and cached the IERS auto values previously, then do the following::

  >>> iers.conf.auto_max_age = None  # doctest: +SKIP

This disables the check of whether the IERS values are sufficiently recent, and
all the transformations (even those outside the time range of available IERS
data) will succeed with at most warnings.

Direct table access
-------------------

In most cases the automatic interface will suffice, but you may need to
directly load and manipulate IERS tables.  IERS-B values are provided as part
of astropy and can be used to calculate time offsets and polar motion
directly, or set up for internal use in further time and coordinate
transformations.  For example::

  >>> from astropy.utils import iers
  >>> t = Time('2010:001')
  >>> iers_b = iers.IERS_B.open()
  >>> iers_b.ut1_utc(t)  # doctest: +FLOAT_CMP
  <Quantity 0.114033 s>
  >>> iers.earth_orientation_table.set(iers_b)
  <ScienceState earth_orientation_table: <IERS_B length=...>...>
  >>> t.ut1.iso
  '2010-01-01 00:00:00.114'

Instead of local copies of IERS files, one can also download them, using
``iers.IERS_A_URL`` (or ``iers.IERS_A_URL_MIRROR``) and ``iers.IERS_B_URL``,
and then use those for future time and coordinate transformations (in this
example, just for a single calculation, by using
`~astropy.utils.iers.earth_orientation_table` as a context manager)::

  >>> iers_a = iers.IERS_A.open(iers.IERS_A_URL)  # doctest: +SKIP
  >>> with iers.earth_orientation_table.set(iers_a):  # doctest: +SKIP
  ...     print(t.ut1.iso)
  2010-01-01 00:00:00.114

To reset to the default, pass in `None` (which is equivalent to passing in
``iers.IERS_Auto.open()``)::

  >>> iers.earth_orientation_table.set(None)  # doctest: +REMOTE_DATA
  <ScienceState earth_orientation_table: <IERS_Auto length=...>...>

To see the internal IERS data that gets used in astropy you can do the
following::

  >>> dat = iers.earth_orientation_table.get()  # doctest: +REMOTE_DATA
  >>> type(dat)  # doctest: +REMOTE_DATA
  <class 'astropy.utils.iers.iers.IERS_Auto'>
  >>> dat  # doctest: +SKIP
  <IERS_Auto length=16196>
   year month  day    MJD   PolPMFlag_A ... UT1Flag    PM_x     PM_y   PolPMFlag
                       d                ...           arcsec   arcsec
  int64 int64 int64 float64     str1    ... unicode1 float64  float64   unicode1
  ----- ----- ----- ------- ----------- ... -------- -------- -------- ---------
     73     1     2 41684.0           I ...        B    0.143    0.137         B
     73     1     3 41685.0           I ...        B    0.141    0.134         B
     73     1     4 41686.0           I ...        B    0.139    0.131         B
     73     1     5 41687.0           I ...        B    0.137    0.128         B
    ...   ...   ...     ...         ... ...      ...      ...      ...       ...
     17     5     2 57875.0           P ...        P 0.007211  0.44884         P
     17     5     3 57876.0           P ...        P 0.008757 0.450321         P
     17     5     4 57877.0           P ...        P 0.010328 0.451777         P
     17     5     5 57878.0           P ...        P 0.011924 0.453209         P
     17     5     6 57879.0           P ...        P 0.013544 0.454617         P

The explanation for most of the columns can be found in the file named
``iers.IERS_A_README``.  The important columns of this table are MJD, UT1_UTC,
UT1Flag, PM_x, PM_y, PolPMFlag::

  >>> dat['MJD', 'UT1_UTC', 'UT1Flag', 'PM_x', 'PM_y', 'PolPMFlag']  # doctest: +SKIP
  <IERS_Auto length=16196>
    MJD    UT1_UTC   UT1Flag    PM_x     PM_y   PolPMFlag
     d        s                arcsec   arcsec
  float64  float64   unicode1 float64  float64   unicode1
  ------- ---------- -------- -------- -------- ---------
  41684.0     0.8075        B    0.143    0.137         B
  41685.0     0.8044        B    0.141    0.134         B
  41686.0     0.8012        B    0.139    0.131         B
  41687.0     0.7981        B    0.137    0.128         B
      ...        ...      ...      ...      ...       ...
  57875.0 -0.6545408        P 0.007211  0.44884         P
  57876.0 -0.6559528        P 0.008757 0.450321         P
  57877.0 -0.6573705        P 0.010328 0.451777         P
  57878.0 -0.6587712        P 0.011924 0.453209         P
  57879.0  -0.660187        P 0.013544 0.454617         P
