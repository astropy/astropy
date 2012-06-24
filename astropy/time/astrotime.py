"""
The astropy.time package provides functionality for manipulating times and
dates. Specific emphasis is placed on supporting time scales (e.g. UTC, TAI,
UT1) and time representations (e.g. JD, MJD, ISO 8601) that are used in
astronomy.
"""
import sys
import time
import itertools
import numpy as np

try:
    from . import sofa_time
except ImportError:
    pass

MJD_ZERO = 2400000.5
SECS_PER_DAY = 86400
TIME_SCALES = ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')
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

_global_opt = {'precision': 3}


def set_opt(**kwargs):
    """Set default options that affect TimeFormat class behavior.

    Example::

      set_opt(precision=6)

    Parameters
    ----------
    precision : int between 0 and 9 inclusive
        Decimal precision when outputting seconds as floating point
    """
    _global_opt.update(**kwargs)


class OptDict(dict):
    """Only a limited set of keys are allowed and the values are
    validated.  If a special value ``__base__`` is set to a dictionary
    then that supplies defaults.
    """
    def __setitem__(self, item, val):
        if item == '__base__':
            pass
        elif item == 'precision':
            if not (isinstance(val, int) and val >= 0 and val < 10):
                raise ValueError('Precision option must be int '
                                 'between 0 and 9 inclusive')
        else:
            raise KeyError('{0} is not a valid option key'.format(item))
        super(OptDict, self).__setitem__(item, val)

    def __getitem__(self, item):
        if '__base__' in self:
            base = super(OptDict, self).__getitem__('__base__').copy()
        else:
            base = {}
        base.update(self)
        return base[item]


class Time(object):
    """Represent and manipulate times and dates for astronomy.

    Specific emphasis is placed on supporting time scales (e.g. UTC, TAI, UT1)
    and time representations (e.g. JD, MJD, ISO 8601) that are used in
    astronomy.

    A Time object is initialized with one or more times in the ``val``
    argument.  The input times in ``val`` must conform to the specified
    ``format`` and must correspond to the specified time ``scale``.  The
    optional ``val2`` time input should be supplied only for numeric input
    formats (e.g. JD) where very high precision (better than 64-bit precision)
    is required.

    Parameters
    ----------
    val : numpy ndarray, list, str, or number
        Data to initialize table.
    val2 : numpy ndarray, list, str, or number; optional
        Data to initialize table.
    format : str
        Format of input value(s)
    scale : str
        Time scale of input value(s)
    opt : dict, optional
        options
    lat : float, optional
        Earth latitude of observer
    lon : float, optional
        Earth longitude of observer
    """
    def __init__(self, val, val2=None, format=None, scale=None,
                 opt={}, lat=0.0, lon=0.0):
        if 'astropy.time.sofa_time' not in sys.modules:
            raise ImportError('Failed to import astropy.time.sofa_time '
                              'extension module (check installation)')
        self.lat = lat
        self.lon = lon
        opt['__base__'] = _global_opt
        self.opt = OptDict(opt)

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

        # To Do: auto determine format and scale if not supplied.  For now
        # raise exception.
        if format not in TIME_FORMATS:
            raise ValueError('Must supply valid format in {0}'
                             .format(sorted(TIME_FORMATS)))

        if scale is not None and scale not in TIME_SCALES:
            raise ValueError('Scale {0} is not in the allowed scales {1}'
                             .format(scale, sorted(TIME_SCALES)))

        self._time = TIME_FORMATS[format](val, val2, scale, self.opt)
        self._format = format
        self._scale = scale

    def _get_format(self):
        return self._format

    def _set_format(self, format):
        NewFormat = TIME_FORMATS[format]
        # If the new format class has a "scale" class attr then that scale is
        # required and the input jd1,2 has to be converted first.
        if hasattr(NewFormat, 'scale'):
            scale = getattr(NewFormat, 'scale')
            new = getattr(self, scale)  # self JDs converted to scale
            self._time = NewFormat(new.jd1, new.jd2, scale, self.opt,
                                   from_jd=True)
        else:
            self._time = NewFormat(self._time.jd1, self._time.jd2,
                                   self.scale, self.opt, from_jd=True)

        self._format = format

    def __repr__(self):
        return ("<Time object: scale='%s' format='%s' vals=%s>" % (
                self.scale, self.format, getattr(self, self.format)))

    def __str__(self):
        return str(getattr(self, self.format))

    format = property(_get_format)

    def _get_scale(self):
        return self._scale

    def _set_scale(self, scale):
        if scale == self._scale:
            return
        if scale not in TIME_SCALES:
            raise ValueError('Scale {0} is not in the allowed scales {1}'
                             .format(scale, sorted(TIME_SCALES)))

        # Determine the chain of scale transformations to get from the current
        # scale to the new scale.  MULTI_HOPS contains a dict of all
        # transformations (xforms) that require intermediate xforms.
        # The MULTI_HOPS dict is keyed by (sys1, sys2) in alphabetical order.
        xform = (self._scale, scale)
        xform_sort = tuple(sorted(xform))
        multi = MULTI_HOPS.get(xform_sort, ())
        xforms = xform_sort[:1] + multi + xform_sort[-1:]
        # If we made the reverse xform then reverse it now.
        if xform_sort != xform:
            xforms = tuple(reversed(xforms))

        # Transform the jd1,2 pairs through the chain of scale xforms.
        jd1, jd2 = self._time.jd1, self._time.jd2
        for sys1, sys2 in itertools.izip(xforms[:-1], xforms[1:]):
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
        self._time = TIME_FORMATS[self.format](jd1, jd2, scale, self.opt,
                                              from_jd=True)
        self._scale = scale

    scale = property(_get_scale)

    def set_opt(self, **kwargs):
        """Set options that affect TimeFormat class behavior for this Time
        instance.

        Example::

          t = Time(100.0, format='mjd', scale='utc')
          t.set_opt(precision=1)

        Parameters
        ----------
        precision : int between 0 and 9 inclusive
            Decimal precision when outputting seconds as floating point
        """
        self.opt.update(**kwargs)

    @property
    def jd1(self):
        vals = self._time.jd1
        return (vals[0].tolist() if self.is_scalar else vals)

    @property
    def jd2(self):
        vals = self._time.jd2
        return (vals[0].tolist() if self.is_scalar else vals)

    @property
    def vals(self):
        return self._time.vals

    def _get_time_object(self, format):
        """Turn this into copy??"""
        tm = Time(self._time.jd1, self._time.jd2,
                  format='jd', scale=self.scale, opt=self.opt)
        tm._set_format(format)
        attrs = ('is_scalar',
                 '_delta_ut1_utc', '_delta_tdb_tt',
                 'lat', 'lon')
        for attr in attrs:
            try:
                setattr(tm, attr, getattr(self, attr))
            except AttributeError:
                pass
        return tm

    def __getattr__(self, attr):
        if attr in TIME_SCALES:
            tm = self._get_time_object(format=self.format)
            tm._set_scale(attr)
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
        return len(self._time.jd1)


class TimeFormat(object):
    """
    Base class for time representations.
    """
    def __init__(self, val1, val2, scale, opt, from_jd=False):
        if hasattr(self.__class__, 'scale'):
            # This format class has a required time scale
            cls_scale = getattr(self.__class__, 'scale')
            if (scale is not None and scale != cls_scale):
                raise ValueError('Class {0} requires scale={1} or None'
                                 .format(self.__class__.__name__, cls_scale))
        else:
            self.scale = scale
        self.opt = OptDict(opt)
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


class TimeMJD(TimeFormat):
    name = 'mjd'

    def set_jds(self, val1, val2):
        # XXX - this routine and vals should be Cythonized to follow the SOFA
        # convention of preserving precision by adding to the larger of the two
        # values in a vectorized operation.  But in most practical cases the
        # first one is probably biggest.
        self.jd1 = val1 + MJD_ZERO
        self.jd2 = val2

    @property
    def vals(self):
        return (self.jd1 - MJD_ZERO) + self.jd2


class TimeFromEpoch(TimeFormat):
    """Base class for times that represent the interval from a particular
    epoch as a floating point multiple of a unit time interval (e.g. seconds
    or days).
    """
    def __init__(self, val1, val2, scale, opt, from_jd=False):
        epoch = Time(self.epoch_val, self.epoch_val2, scale=self.epoch_scale,
                     format=self.epoch_format, opt=opt)
        self.epoch = getattr(epoch, self.scale)
        super(TimeFromEpoch, self).__init__(val1, val2, scale, opt, from_jd)

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
    epoch_scale = 'utc'
    epoch_format = 'iso'
    scale = 'utc'


class TimeCxcSec(TimeFromEpoch):
    """Chandra X-ray Center seconds from 1998-01-01 00:00:00 TT.
    """
    name = 'cxcsec'
    unit = 1.0 / SECS_PER_DAY  # in days (1 day == 86400 seconds)
    epoch_val = '1998-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'tt'
    epoch_format = 'iso'
    scale = 'tai'


class TimeString(TimeFormat):
    """Base class for string-like time represetations.
    """

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

        self.jd1, self.jd2 = sofa_time.dtf_jd(self.scale,
                                              iy, im, id, ihr, imin, dsec)

    def str_kwargs(self):
        iys, ims, ids, ihmsfs = sofa_time.jd_dtf(self.scale.upper(),
                                                 self.opt['precision'],
                                                 self.jd1, self.jd2)
        for iy, im, id, ihmsf in itertools.izip(iys, ims, ids, ihmsfs):
            ihr, imin, isec, ifracsec = ihmsf
            yield {'year': int(iy), 'mon': int(im), 'day': int(id),
                   'hour': int(ihr), 'min': int(imin), 'sec': int(isec),
                   'fracsec': int(ifracsec)}

    @property
    def vals(self):
        str_fmt = self.str_fmt
        if self.opt['precision'] > 0:
            str_fmt += '.{fracsec:0' + str(self.opt['precision']) + 'd}'

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
    scale = 'utc'
    strptime_fmt = '%Y-%m-%dT%H:%M:%S'
    str_fmt = '{year:d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}:{sec:02d}'


class TimeEpochDate(TimeFormat):
    """Base class for support Besselian and Julian epoch dates (e.g.
    B1950.0 or J2000.0 etc).
    """
    def set_jds(self, val1, val2):
        epoch_to_jd = getattr(sofa_time, self.epoch_to_jd)
        self.jd1, self.jd2 = epoch_to_jd(val1 + val2)

    @property
    def vals(self):
        jd_to_epoch = getattr(sofa_time, self.jd_to_epoch)
        return jd_to_epoch(self.jd1, self.jd2)


class TimeBesselianEpoch(TimeEpochDate):
    """Besselian Epoch year"""
    name = 'byear'
    epoch_to_jd = 'besselian_epoch_jd'
    jd_to_epoch = 'jd_besselian_epoch'


class TimeJulianEpoch(TimeEpochDate):
    """Julian Epoch year"""
    name = 'jyear'
    epoch_to_jd = 'julian_epoch_jd'
    jd_to_epoch = 'jd_julian_epoch'


# Set module constant with names of all available time formats
TIME_FORMATS = {}
_module = sys.modules[__name__]
for name in dir(_module):
    val = getattr(_module, name)
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
