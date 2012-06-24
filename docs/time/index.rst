.. _astropy_time:

.. include:: references.txt

****************************************************
Time and Dates (`astropy.time`)
****************************************************

.. note:: This documentation is largely aimed at the developer community right now.

Introduction
============

The `astropy.time` package provides functionality for manipulating times and
dates.  Specific emphasis is placed on supporting time scales (e.g. UTC, TAI, UT1) and
time representations (e.g. JD, MJD, ISO 8601) that are used in astronomy.
It uses Cython to wrap the C language `SOFA`_ time and calendar
routines.  All time scale conversions are done by Cython vectorized versions
of the `SOFA`_ routines and are fast and memory efficient.

Other parts of the current implementation are pure Python but with a goal of
using Cython routines where possible.

Following the SOFA implementation, the internal representation of time is a
pair of doubles that sum up to the time JD in the current scale.  The SOFA
routines take care throughout to maintain overall precision of the double pair.
The user is free to choose the way in which total JD is distributed between the
two values.  MOST IMPORTANTLY the user free to entirely ignore the whole issue
and just supply time values as strings or doubles and not even know what is
happening underneath.  After working with SOFA for a week I am convinced the
pair-of-doubles strategy is quite sound and useful.


Getting Started
===============

The basic way to use `astropy.time` is to create a |Time|
object by supplying one or more input time values as well as the format and
scale of those values.  The input time(s) can either be a single scalar like
`"2010-01-01 00:00:00"` or a list (or `numpy` array) of values as shown below.
In general any output values have the same shape (scalar or array) as the input.

  >>> import astropy.time as astrotime

  >>> times = ['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00']
  >>> t = astrotime.Time(times, format='iso', scale='utc')
  >>> t
  <Time object: scale='utc' format='iso' vals=['1999-01-01 00:00:00.123' '2010-01-01 00:00:00.000']>
  >>> t.jd1
  array([ 2451179.5,  2455197.5])
  >>> t.jd2
  array([  1.42889802e-06,   0.00000000e+00])

Terminology
------------

In this context of `astropy.time` a "time" is a single instant of time which is
independent of the way the time is represented (the "format") and the time
"scale" which specifies the offset and scaling relation of the unit of time.

Format
^^^^^^^^^^^

The time format specifies how an instant of time is represented.  The currently
available formats are can be found in the ``astrotime.TIME_FORMATS`` dict and are:

=========  ===================================================
Format            Class
=========  ===================================================
byear      :class:`~astropy.time.astrotime.TimeBesselianEpoch`
cxcsec     :class:`~astropy.time.astrotime.TimeCxcSec`
iso        :class:`~astropy.time.astrotime.TimeISO`
isot       :class:`~astropy.time.astrotime.TimeISOT`
jd         :class:`~astropy.time.astrotime.TimeJD`
jyear      :class:`~astropy.time.astrotime.TimeJulianEpoch`
mjd        :class:`~astropy.time.astrotime.TimeMJD`
unix       :class:`~astropy.time.astrotime.TimeUnix`
=========  ===================================================

Time Scale
^^^^^^^^^^

The time scale (or `time standard
<http://en.wikipedia.org/wiki/Time_standard>`_) is "a specification for
measuring time: either the rate at which time passes; or points in time; or
both" [#]_. See also [#]_ and [#]_. 

  >>> astrotime.TIME_SCALES
  ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')

====== =================================
Scale        Description
====== =================================
tai    International Atomic Time   (TAI)
tcb    Barycentric Coordinate Time (TCB)
tcg    Geocentric Coordinate Time  (TCG)
tdb    Barycentric Dynamical Time  (TDB)
tt     Terrestrial Time            (TT)
ut1    Universal Time              (UT1)
utc    Coordinated Universal Time  (UTC)
====== =================================

.. [#] Wikipedia `time standard <http://en.wikipedia.org/wiki/Time_standard>`_ article
.. [#] SOFA Time Scale and Calendar Tools 
       `(PDF) <http://www.iausofa.org/2012_0301_C/sofa/sofa_ts_c.pdf>`_
.. [#] `<http://www.ucolick.org/~sla/leapsecs/timescales.html>`_

Examples
-----------

Set time scale to TAI::

  >>> t2 = t.tai
  <Time object: scale='tai' format='iso' vals=['1999-01-01 00:00:32.123' '2010-01-01 00:00:34.000']>
  >>> t2.jd1
  array([ 2451179.5,  2455197.5])
  >>> t2.jd2
  array([ 0.0003718 ,  0.00039352])

Get a new ``Time`` object which is referenced to the TT scale (internal JD1 and JD1 are
now with respect to the TT time scale)::

  >>> t2.tt  # scale property returns a new Time object
  <Time object: scale='tt' format='iso' vals=['1999-01-01 00:01:04.307' '2010-01-01 00:01:06.184']>

Get the representation of the ``Time`` object in a particular format (in this
case seconds since 1998.0).  This returns either a scalar or array, depending
on whether the input was a scalar or array::

  >>> t2.cxcsec  # format property returns an array or scalar of that representation
  array([  3.15360643e+07,   3.78691266e+08])


Use properties to convert scales and formats.  Note that the UT1 to UTC
transformation requires a supplementary value (``delta_ut1_utc``) that can be
obtained by interpolating from a table supplied by IERS.  This will be included
in the package later.
::

  >>> t = astrotime.Time('2010-01-01 00:00:00', format='iso', scale='utc')
  >>> t.set_delta_ut1_utc(0.3341)  # Explicitly set one part of the transformation
  >>> t.jd        # JD representation of time in current scale (UTC)
  2455197.5
  >>> t.iso       # ISO representation of time in current scale (UTC)
  '2010-01-01 00:00:00.000'
  >>> t.tt.iso    # ISO representation of time in TT scale
  '2010-01-01 00:01:06.184'
  >>> t.tai.iso   # ISO representation of time in TAI scale
  '2010-01-01 00:00:34.000'
  >>> t.utc.jd    # JD representation of time in UTC scale
  2455197.5
  >>> t.ut1.jd    # JD representation of time in UT1 scale
  2455197.500003867
  >>> t.tcg.isot  # ISO time with a "T" in the middle
  '2010-01-01T00:00:00.000'
  >>> t.unix      # seconds since 1970.0 (utc) excluding leapseconds
  1262304000.0
  >>> t.cxcsec    # SI seconds since 1998.0 (tt)
  378691266.184

Set the output precision which is used for some formats::

  >>> t.set_opt(precision=9)
  >>> t.iso
  '2010-01-01 00:00:00.000000000'

Transform from UTC to all supported time scales (TAI, TCB, TCG, TDB, TT, UT1,
UTC).  This requires auxilliary information (latitude and longitude).
::

  >>> lat = 19.48125
  >>> lon = -155.933222
  >>> t = astrotime.Time('2006-01-15 21:24:37.5', format='iso', scale='utc',
  ...                    lat=lat, lon=lon)
  >>> t.set_opt(precision=6)
  >>> t.set_delta_ut1_utc(0.3341)  # Explicitly set one part of the transformation
  >>> t.utc.iso
  '2006-01-15 21:24:37.500000'
  >>> t.ut1.iso
  '2006-01-15 21:24:37.834100'
  >>> t.tai.iso
  '2006-01-15 21:25:10.500000'
  >>> t.tt.iso
  '2006-01-15 21:25:42.684000'
  >>> t.tcg.iso
  '2006-01-15 21:25:43.322690'
  >>> t.tdb.iso
  '2006-01-15 21:25:42.683799'
  >>> t.tcb.iso
  '2006-01-15 21:25:56.893378'


See Also
===================

Other time references TBD.


Reference/API
=============

.. automodapi:: astropy.time


Acknowledgments and Licenses
=======================================

This package makes use of the `SOFA Software
<http://www.iausofa.org/index.html>`_ ANSI C library.  The copyright of the SOFA
Software belongs to the Standards Of Fundamental Astronomy Board of the
International Astronomical Union.  This library is made available under the
terms of the `SOFA license <http://www.iausofa.org/tandc.html>`_.
