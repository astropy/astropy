"""
========================
Time class for astronomy
========================

Summary
---------

This package is a candidate implementation of the ``astropy.time``
sub-package.  It uses Cython to wrap the C language SOFA time and calendar
routines.  All time system conversions are done by Cython vectorized versions
of the SOFA routines.  These transformations will be fast and memory efficient
(no temporary arrays).  

Other parts of the current implementation are pure Python but with a goal of
using Cython routines where possible.

Following the SOFA implementation, the internal representation of time is a
pair of doubles that sum up to the time JD in the current system.  The SOFA
routines take care throughout to maintain overall precision of the double pair.
The user is free to choose the way in which total JD is distributed between the
two values.  MOST IMPORTANTLY the user free to entirely ignore the whole issue
and just supply time values as strings or doubles and not even know what is
happening underneath.  After working with SOFA for a week I am convinced the
pair-of-doubles strategy is quite sound and useful.

No formal test suite is available yet but I have confirmed that the
transformations in the examples below match expected values based on the
examples in the `SOFA Time Scale and Calendar Tools
<http://www.iausofa.org/2012_0301_C/sofa/sofa_ts_c.pdf>`_ tutorial.

Note: SOFA and most official references talk about time "scales" not time
"systems".  I find time "scale" confusing but we could easily change to that
terminology.

Build for testing
------------------

The current commit is not integrated with the astropy build system.  Instead
build inplace with (requires Cython)::  

  % python setup.py build_ext --inplace

Examples
--------

::

  >>> import astropy.time as astrotime

  >>> times = ['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00']
  >>> t = astrotime.Time(times, format='iso', system='utc')
  >>> t
  <Time object: system='utc' format='iso' vals=['1999-01-01 00:00:00.123' '2010-01-01 00:00:00.000']>
  >>> t.jd1
  array([ 2451179.5,  2455197.5])
  >>> t.jd2
  array([  1.42889802e-06,   0.00000000e+00])

Set system to TAI::

  >>> t.set_system('tai')
  >>> t
  <Time object: system='tai' format='iso' vals=['1999-01-01 00:00:32.123' '2010-01-01 00:00:34.000']>
  >>> t.jd1
  array([ 2451179.5,  2455197.5])
  >>> t.jd2
  array([ 0.0003718 ,  0.00039352])

Get a new ``Time`` object which is referenced to the TT system (internal JD1 and JD1 are
now with respect to TT system)::

  >>> t.tt  # system property returns a new Time object
  <Time object: system='tt' format='iso' vals=['1999-01-01 00:01:04.307' '2010-01-01 00:01:06.184']>

Get the representation of the ``Time`` object in a particular format (in this
case seconds since 1998.0).  This returns either a scalar or array, depending
on whether the input was a scalar or array::

  >>> t.cxcsec  # format property returns an array or scalar of that representation
  array([  3.15360643e+07,   3.78691266e+08])


Use properties to convert systems and formats.  Note that the UT1 to UTC
transformation requires a supplementary value (``delta_ut1_utc``) that can be
obtained by interpolating from a table supplied by IERS.  This will be included
in the package later.
::

  >>> t = astrotime.Time('2010-01-01 00:00:00', format='iso', system='utc')
  >>> t.set_delta_ut1_utc(0.3341)  # Explicitly set one part of the transformation
  >>> t.jd        # JD representation of time in current system (UTC)
  2455197.5
  >>> t.iso       # ISO representation of time in current system (UTC)
  '2010-01-01 00:00:00.000'
  >>> t.tt.iso    # ISO representation of time in TT system
  '2010-01-01 00:01:06.184'
  >>> t.tai.iso   # ISO representation of time in TAI system
  '2010-01-01 00:00:34.000'
  >>> t.utc.jd    # JD representation of time in UTC system
  2455197.5
  >>> t.ut1.jd    # JD representation of time in UT1 system
  2455197.500003867
  >>> t.tcg.isot  # ISO time with a "T" in the middle
  '2010-01-01T00:00:00.000'
  >>> t.unix      # seconds since 1970.0 (utc) excluding leapseconds
  1262304000.0
  >>> t.cxcsec    # SI seconds since 1998.0 (tt)
  378691266.184

Set the output precision which is used for some formats::

  >>> t.precision = 9
  >>> t.iso
  '2010-01-01 00:00:00.000000000'

Transform from UTC to all supported time systems (TAI, TCB, TCG, TDB, TT, UT1,
UTC).  This requires auxilliary information (latitude and longitude).
::

  >>> lat = 19.48125
  >>> lon = -155.933222
  >>> t = astrotime.Time('2006-01-15 21:24:37.5', format='iso', system='utc',
  ...                    precision=6, lat=lat, lon=lon)
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
"""
import sys
import time
from itertools import izip
import numpy as np
try:
    import sofa_time
except ImportError:
    pass


SECS_PER_DAY = 86400
TIME_SYSTEMS = ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')
MULTI_HOPS = {('tai', 'tcb'): ('tt', 'tdb'),
              ('tai', 'tcg'): ('tt',),
              ('tai', 'ut1'): ('utc',),
              ('tai', 'tdb'): ('tt',),
              ('tcb', 'tcg'): ('tdb', 'tt'),
              ('tcb', 'tt'): ('tdb',),
              ('tcb', 'ut1'): ('tdb', 'tt', 'tai', 'utc'),
              ('tcb', 'utc'): ('tdb', 'tt', 'tai'),
              ('tcg', 'tdb'): ('tt', 'tdb'),
              ('tcg', 'ut1'): ('tt', 'tai', 'utc'),
              ('tcg', 'utc'): ('tt', 'tai'),
              ('tdb', 'ut1'): ('tt', 'tai', 'utc'),
              ('tdb', 'utc'): ('tt', 'tai'),
              ('tt', 'ut1'): ('tai', 'utc'),
              ('tt', 'utc'): ('tai',),
              }


class Time(object):
    def __init__(self, val, val2=None, format=None, system=None,
                 precision=3, lat=0.0, lon=0.0):
        self.lat = lat
        self.lon = lon
        self.precision = precision

        # Coerce val into a 1-d array
        val, val_ndim = _make_1d_array(val)

        # If val2 is None then replace with zeros of the same length
        if val2 is None:
            val2 = np.zeros(len(val), dtype=np.double)
            val2_ndim = val_ndim
        else:
            val2, val2_ndim = _make_1d_array(val2)

        # Consistency checks
        if len(val) != len(val2):
            raise ValueError('Input val and val2 must match in length')

        self.is_scalar = (val_ndim == 0)
        if val_ndim != val2_ndim:
            raise ValueError('Input val and val2 must have same dimensions')

        # To Do: auto determine format and system if not supplied.  For now
        # raise exception.
        if format not in TIME_FORMATS:
            raise ValueError('Must supply valid format in {0}'
                             .format(sorted(TIME_FORMATS)))

        if system is not None and system not in TIME_SYSTEMS:
            raise ValueError('System {0} is not in the allowed systems {1}'
                             .format(system, sorted(TIME_SYSTEMS)))

        self._time = TIME_FORMATS[format](val, val2, system)
        self._time.precision = self.precision
        # XXX setting precision here and other places is awkward.  Improve?
        self._format = format
        self._system = system

    def get_format(self):
        return self._format

    def set_format(self, format):
        NewFormat = TIME_FORMATS[format]
        # If the new format class has a "system" class attr then that system is
        # required and the input jd1,2 has to be converted first.
        if hasattr(NewFormat, 'system'):
            system = getattr(NewFormat, 'system')
            new = getattr(self, system)  # self JDs converted to system
            self._time = NewFormat(new.jd1, new.jd2, system, from_jd=True)
        else:
            self._time = NewFormat(self.jd1, self.jd2,
                                   self.system, from_jd=True)

        self._time.precision = self.precision
        self._format = format

    def __repr__(self):
        return ("<Time object: system='%s' format='%s' vals=%s>" % (
                self.system, self.format, getattr(self, self.format)))

    def __str__(self):
        return str(getattr(self, self.format))

    format = property(get_format, set_format)

    def get_system(self):
        return self._system

    def set_system(self, system):
        if system == self._system:
            return
        if system not in TIME_SYSTEMS:
            raise ValueError('System {0} is not in the allowed systems {1}'
                             .format(system, sorted(TIME_SYSTEMS)))

        # Determine the chain of system transformations to get from the current
        # system to the new system.  MULTI_HOPS contains a dict of all
        # transformations (xforms) that require intermediate xforms.
        # The MULTI_HOPS dict is keyed by (sys1, sys2) in alphabetical order.
        xform = (self._system, system)
        xform_sort = tuple(sorted(xform))
        multi = MULTI_HOPS.get(xform_sort, ())
        xforms = xform_sort[:1] + multi + xform_sort[-1:]
        # If we made the reverse xform then reverse it now.
        if xform_sort != xform:
            xforms = tuple(reversed(xforms))

        # Transform the jd1,2 pairs through the chain of system xforms.
        jd1, jd2 = self.jd1, self.jd2
        for sys1, sys2 in izip(xforms[:-1], xforms[1:]):
            # Some xforms require an additional delta_ argument that is
            # provided through Time methods.  These values may be supplied by
            # the user or computed based on available approximations.  The
            # get_delta_ methods are available for only one combination of
            # sys1, sys2 though the property applies for both xform directions.
            args = [jd1, jd2]
            for sys12 in ((sys1, sys2), (sys2, sys1)):
                dt_method = 'get_delta_{0}_{1}'.format(*sys12)
                try:
                    get_dt = getattr(self, dt_method)
                except AttributeError:
                    pass
                else:
                    args.append(get_dt(jd1, jd2))
                    break

            conv_func = getattr(sofa_time, sys1 + '_' + sys2)
            jd1, jd2 = conv_func(*args)
        self._time = TIME_FORMATS[self.format](jd1, jd2, system,
                                              from_jd=True)
        self._time.precision = self.precision
        self._system = system

    system = property(get_system, set_system)

    @property
    def jd1(self):
        return self._time.jd1

    @property
    def jd2(self):
        return self._time.jd2

    @property
    def vals(self):
        return self._time.vals

    def _get_time_object(self, format):
        tm = Time(self.jd1, self.jd2,
                  format='jd', system=self.system)
        tm.set_format(format)
        attrs = ('precision', 'is_scalar',
                 '_delta_ut1_utc', '_delta_tdb_tt',
                 'lat', 'lon')
        for attr in attrs:
            try:
                setattr(tm, attr, getattr(self, attr))
            except AttributeError:
                pass
        tm._time.precision = tm.precision
        return tm

    def __getattr__(self, attr):
        if attr in TIME_SYSTEMS:
            tm = self._get_time_object(format=self.format)
            tm.set_system(attr)
            return tm

        elif attr in TIME_FORMATS:
            tm = self._get_time_object(format=attr)
            return (tm.vals[0].tolist() if self.is_scalar else tm.vals)

        else:
            # Should raise AttributeError
            return self.__getattribute__(attr)

    def _match_len(self, val):
        """Ensure that `val` is matched to length of self.
        If val has length 1 then broadcast, otherwise cast to double
        and make sure length matches.
        """
        val, ndim = _make_1d_array(val)
        if len(val) == 1:
            oval = val
            val = np.empty(len(self), dtype=np.double)
            val[:] = oval
        elif len(val) != len(self):
            raise ValueError('Attribute length must match Time object length')
        return val

    # SOFA DUT arg = UT1 - UTC
    def get_delta_ut1_utc(self, jd1, jd2):
        """
        Sec. 4.3.1: the arg DUT is the quantity delta_UT1 = UT1 - UTC in
        seconds. It can be obtained from tables published by the IERS.
        XXX - get that table when needed and interpolate or whatever.
        """
        if not hasattr(self, '_delta_ut1_utc'):
            self._delta_ut1_utc = np.zeros(len(self), dtype=np.double)

        return self._delta_ut1_utc

    def set_delta_ut1_utc(self, val):
        self._delta_ut1_utc = self._match_len(val)

    # SOFA DTR arg = TDB - TT
    def get_delta_tdb_tt(self, jd1, jd2):
        if not hasattr(self, '_delta_tdb_tt'):
            # First go from the current input time (which is either
            # TDB or TT) to an approximate UTC.  Since TT and TDB are
            # pretty close (few msec?), assume TT.
            njd1, njd2 = sofa_time.tt_tai(jd1, jd2)
            njd1, njd2 = sofa_time.tai_utc(njd1, njd2)
            # XXX actually need to go to UT1 which needs DUT.
            ut = njd1 + njd2

            # Compute geodetic params needed for d_tdb_tt()
            phi = np.radians(self.lat)
            elon = np.radians(self.lon)
            xyz = sofa_time.iau_gd2gc(1, elon, phi, 0.0)
            u = np.sqrt(xyz[0] ** 2 + xyz[1] ** 2)
            v = xyz[2]

            self._delta_tdb_tt = sofa_time.d_tdb_tt(jd1, jd2, ut, elon, u, v)

        return self._delta_tdb_tt

    def set_delta_tdb_tt(self, val):
        self._delta_tdb_tt = self._match_len(val)

    def __len__(self):
        return len(self.jd1)


class TimeFormat(object):
    """
    Base class for time representations.
    """
    def __init__(self, val1, val2, system, from_jd=False):
        if hasattr(self.__class__, 'system'):
            # This format class has a required time system
            cls_system = getattr(self.__class__, 'system')
            if (system is not None and system != cls_system):
                raise ValueError('Class {0} requires system={1} or None'
                                 .format(self.__class__.__name__, cls_system))
        else:
            self.system = system
        self.n_times = len(val1)
        if len(val1) != len(val2):
            raise ValueError('Input val1 and val2 must match in length')

        if from_jd:
            self.jd1 = val1
            self.jd2 = val2
        else:
            self.set_jds(val1, val2)


class TimeJD(TimeFormat):
    name = 'jd'

    def set_jds(self, val1, val2):
        self.jd1 = val1
        self.jd2 = val2

    @property
    def vals(self):
        return self.jd1 + self.jd2


class TimeFromEpoch(TimeFormat):
    """Base class for times that represent the interval from a particular
    epoch as a floating point multiple of a unit time interval (e.g. seconds
    or days).
    """
    def __init__(self, val1, val2, system, from_jd=False):
        epoch = Time(self.epoch_val, self.epoch_val2, system=self.epoch_system,
                     format=self.epoch_format)
        self.epoch = getattr(epoch, self.system)
        super(TimeFromEpoch, self).__init__(val1, val2, system, from_jd)

    def set_jds(self, val1, val2):
        self.jd1 = self.epoch.jd1 + val2 * self.unit
        self.jd2 = self.epoch.jd2 + val1 * self.unit

    @property
    def vals(self):
        return ((self.jd1 - self.epoch.jd1) +
                (self.jd2 - self.epoch.jd2)) / self.unit


class TimeUnix(TimeFromEpoch):
    """Unix time: seconds from 1970-01-01 00:00:00 UTC.

    NOTE: this quantity is not exactly unix time and differs from the strict
    POSIX definition by up to 1 second on days with a leap second.  POSIX
    unix time actually jumps backward by 1 second at midnight on leap second
    days while this class value is monotonically increasing at 86400 seconds
    per UTC day.
    """
    name = 'unix'
    unit = 1.0 / SECS_PER_DAY  # in days (1 day == 86400 seconds)
    epoch_val = '1970-01-01 00:00:00'
    epoch_val2 = None
    epoch_system = 'utc'
    epoch_format = 'iso'
    system = 'utc'


class TimeCxcSec(TimeFromEpoch):
    """Chandra X-ray Center seconds from 1998-01-01 00:00:00 TT.
    """
    name = 'cxcsec'
    unit = 1.0 / SECS_PER_DAY  # in days (1 day == 86400 seconds)
    epoch_val = '1998-01-01 00:00:00'
    epoch_val2 = None
    epoch_system = 'tt'
    epoch_format = 'iso'
    system = 'tai'


class TimeString(TimeFormat):
    """Base class for string-like time represetations.
    """
    precision = 3  # default fractional second precision (millisecs)

    def set_jds(self, val1, val2):
        """
        Parse the time strings contained in val1 and set jd1, jd2.
        """
        iy = np.empty(self.n_times, dtype=np.intc)
        im = np.empty(self.n_times, dtype=np.intc)
        id = np.empty(self.n_times, dtype=np.intc)
        ihr = np.empty(self.n_times, dtype=np.intc)
        imin = np.empty(self.n_times, dtype=np.intc)
        dsec = np.empty(self.n_times, dtype=np.double)

        if self.strptime_fmt.endswith('%S'):
            ends_with_secs = True

        for i, timestr in enumerate(val1):
            if ends_with_secs:
                try:
                    idot = timestr.rindex('.')
                except:
                    fracsec = 0.0
                else:
                    timestr, fracsec = timestr[:idot], timestr[idot:]
                    fracsec = float(fracsec)

            tm = time.strptime(timestr, self.strptime_fmt)
            iy[i] = tm.tm_year
            im[i] = tm.tm_mon
            id[i] = tm.tm_mday
            ihr[i] = tm.tm_hour
            imin[i] = tm.tm_min
            dsec[i] = tm.tm_sec + fracsec

        self.jd1, self.jd2 = sofa_time.dtf_jd(self.system,
                                              iy, im, id, ihr, imin, dsec)

    def str_kwargs(self):
        iys, ims, ids, ihmsfs = sofa_time.jd_dtf(self.system.upper(),
                                                 self.precision,
                                                 self.jd1, self.jd2)
        for iy, im, id, ihmsf in izip(iys, ims, ids, ihmsfs):
            ihr, imin, isec, ifracsec = ihmsf
            yield {'year': int(iy), 'mon': int(im), 'day': int(id),
                   'hour': int(ihr), 'min': int(imin), 'sec': int(isec),
                   'fracsec': int(ifracsec)}

    @property
    def vals(self):
        str_fmt = self.str_fmt
        if self.precision > 0:
            str_fmt += '.{fracsec:0' + str(self.precision) + 'd}'

        # Try to optimize this later.  Can't pre-allocate because length of
        # output could change, e.g. year rolls from 999 to 1000.
        outs = []
        for kwargs in self.str_kwargs():
            outs.append(str_fmt.format(**kwargs))

        return np.array(outs)


class TimeISO(TimeString):
    name = 'iso'
    strptime_fmt = '%Y-%m-%d %H:%M:%S'
    str_fmt = '{year:d}-{mon:02d}-{day:02d} {hour:02d}:{min:02d}:{sec:02d}'


class TimeISOT(TimeString):
    name = 'isot'
    system = 'utc'
    strptime_fmt = '%Y-%m-%dT%H:%M:%S'
    str_fmt = '{year:d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}:{sec:02d}'


# Set module constant with names of all available time formats
TIME_FORMATS = {}
module = sys.modules[__name__]
for name in dir(module):
    val = getattr(module, name)
    try:
        ok = issubclass(val, TimeFormat)
    except:
        pass
    else:
        if ok and hasattr(val, 'name'):
            TIME_FORMATS[val.name] = val


def _make_1d_array(val):
    val = np.asarray(val)
    val_ndim = val.ndim  # remember original ndim
    if val.ndim == 0:
        val = np.asarray([val])
    elif val_ndim > 1:
        # Maybe lift this restriction later to allow multi-dim in/out?
        raise ValueError('Input val must be zero or one dimensional')

    # Allow only string or float arrays as input (XXX datetime later...)
    if val.dtype.kind == 'i':
        val = np.asarray(val, dtype=np.float64)

    return val, val_ndim
