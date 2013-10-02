# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The astropy.time package provides functionality for manipulating times and
dates. Specific emphasis is placed on supporting time scales (e.g. UTC, TAI,
UT1) and time representations (e.g. JD, MJD, ISO 8601) that are used in
astronomy.
"""
from datetime import datetime
import time
import itertools
import numpy as np

from ..utils import deprecated, deprecated_attribute
from ..utils.compat.misc import override__dir__

__all__ = ['Time', 'TimeDelta', 'TimeFormat', 'TimeJD', 'TimeMJD',
           'TimeFromEpoch', 'TimeUnix', 'TimeCxcSec', 'TimeGPS', 'TimePlotDate',
           'TimeDatetime',
           'TimeString', 'TimeISO', 'TimeISOT', 'TimeYearDayTime', 'TimeEpochDate',
           'TimeBesselianEpoch', 'TimeJulianEpoch', 'TimeDeltaFormat',
           'TimeDeltaSec', 'TimeDeltaJD', 'ScaleValueError',
           'OperandTypeError', 'TimeEpochDateString',
           'TimeBesselianEpochString', 'TimeJulianEpochString',
           'TIME_FORMATS', 'TIME_DELTA_FORMATS', 'TIME_SCALES',
           'TIME_DELTA_SCALES']

__doctest_skip__ = ['TimePlotDate']

try:
    # Not guaranteed available at setup time
    from . import erfa_time
except ImportError:
    if not _ASTROPY_SETUP_:
        raise

MJD_ZERO = 2400000.5
SECS_PER_DAY = 86400

# These both get filled in at end after TimeFormat subclasses defined
TIME_FORMATS = {}
TIME_DELTA_FORMATS = {}

TIME_SCALES = ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')
TIME_DELTA_SCALES = ('tai',)

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
    """
    Represent and manipulate times and dates for astronomy.

    A `Time` object is initialized with one or more times in the ``val``
    argument.  The input times in ``val`` must conform to the specified
    ``format`` and must correspond to the specified time ``scale``.  The
    optional ``val2`` time input should be supplied only for numeric input
    formats (e.g. JD) where very high precision (better than 64-bit precision)
    is required.

    Parameters
    ----------
    val : sequence, str, number, or `~astropy.time.Time` object
        Value(s) to initialize the time or times.
    val2 : sequence, str, or number; optional
        Value(s) to initialize the time or times.
    format : str, optional
        Format of input value(s)
    scale : str, optional
        Time scale of input value(s)
    lat : float, optional
        Earth latitude of observer (decimal degrees)
    lon : float, optional
        Earth longitude of observer (decimal degrees)
    copy : bool, optional
        Make a copy of the input values
    """

    is_scalar = deprecated_attribute(name='is_scalar', since='0.3',
                                     alternative='isscalar')

    _precision = 3  # Precision when for seconds as floating point
    _in_subfmt = '*'  # Select subformat for inputting string times
    _out_subfmt = '*'  # Select subformat for outputting string times

    SCALES = TIME_SCALES
    """List of time scales"""

    FORMATS = TIME_FORMATS
    """Dict of time formats"""

    def __new__(cls, val, val2=None, format=None, scale=None,
                precision=None, in_subfmt=None, out_subfmt=None,
                lat=0.0, lon=0.0, copy=False):

        if isinstance(val, cls):
            self = val.replicate(format=format, copy=copy)
        else:
            self = super(Time, cls).__new__(cls)
        return self

    def __init__(self, val, val2=None, format=None, scale=None,
                 precision=None, in_subfmt=None, out_subfmt=None,
                 lat=0.0, lon=0.0, copy=False):

        self.lat = lat
        self.lon = lon
        if precision is not None:
            self.precision = precision
        if in_subfmt is not None:
            self.in_subfmt = in_subfmt
        if out_subfmt is not None:
            self.out_subfmt = out_subfmt

        if isinstance(val, self.__class__):
            if scale is not None:
                self._set_scale(scale)
        else:
            self._init_from_vals(val, val2, format, scale, copy)

    def _init_from_vals(self, val, val2, format, scale, copy):
        """
        Set the internal _format, scale, and _time attrs from user
        inputs.  This handles coercion into the correct shapes and
        some basic input validation.
        """

        # check whether input is some form of list of Time objects,
        # since these should be treated separately
        try:
            val0 = val[0]
        except:
            isiterable_of_times = False
        else:
            isiterable_of_times = isinstance(val0, self.__class__)

        if isiterable_of_times:
            if val2 is not None:
                raise ValueError(
                    'non-None second value for list of {0!r} objects'
                    .format(self.__class__.__name__))
            self.isscalar = False
        else:
            # Coerce val into a 1-d array
            val, val_ndim = _make_1d_array(val, copy)

            # If val2 is None then replace with zeros of the same length
            if val2 is None:
                val2 = np.zeros(len(val), dtype=np.double)
                val2_ndim = val_ndim
            else:
                val2, val2_ndim = _make_1d_array(val2, copy)

            # Consistency checks
            if len(val) != len(val2):
                raise ValueError('Input val and val2 must match in length')

            self.isscalar = (val_ndim == 0)
            if val_ndim != val2_ndim:
                raise ValueError('Input val and val2 must have same dimensions')

        if scale is not None:
            if not (isinstance(scale, basestring) and
                    scale.lower() in self.SCALES):
                raise ScaleValueError("Scale {0} is not in the allowed scales "
                                      "{1}".format(repr(scale),
                                                   sorted(self.SCALES)))

        # Parse / convert input values into internal jd1, jd2 based on format
        self._time = self._get_time_fmt(val, val2, format, scale)
        self._format = self._time.name

    def _get_time_fmt(self, val, val2, format, scale):
        """
        Given the supplied val, val2, format and scale try to instantiate
        the corresponding TimeFormat class to convert the input values into
        the internal jd1 and jd2.

        If format is None and the input is a string-type or object array then guess
        available formats and stop when one matches.
        """

        if format is None and (isinstance(val[0], self.__class__) or
                               val.dtype.kind in ('S', 'U', 'O')):
            formats = [(name, cls) for name, cls in self.FORMATS.items()
                       if issubclass(cls, TimeUnique)]
            err_msg = 'any of the formats where the format keyword is optional {0}'.format(
                [name for name, cls in formats])
        elif not (isinstance(format, basestring) and
                  format.lower() in self.FORMATS):
            if format is None:
                raise ValueError("No time format was given, and the input is "
                                 "not unique")
            else:
                raise ValueError("Format {0} is not one of the allowed "
                                 "formats {1}".format(repr(format),
                                                      sorted(self.FORMATS)))
        else:
            formats = [(format, self.FORMATS[format])]
            err_msg = 'the format class {0}'.format(format)

        for format, FormatClass in formats:
            try:
                return FormatClass(val, val2, scale, self.precision,
                                   self.in_subfmt, self.out_subfmt)
            except (ValueError, TypeError):
                pass
        else:
            raise ValueError('Input values did not match {0}'.format(err_msg))

    @classmethod
    def now(cls):
        """
        Creates a new object corresponding to the instant in time this
        method is called.

        .. note::
            "Now" is determined using the `datetime.utcnow` function, so
            its accuracy and precision is determined by that function.
            Generally that means it is set by the accuracy of your
            system clock.

        Returns
        -------
        nowtime
            A new `Time` object (or a subclass of `Time` if this is called from
            such a subclass) at the current time.
        """
        # call `utcnow` immediately to be sure it's ASAP
        dtnow = datetime.utcnow()
        return cls(val=dtnow, format='datetime', scale='utc')

    @property
    def format(self):
        """Time format"""
        return self._format

    def __repr__(self):
        return ("<{0} object: scale='{1}' format='{2}' value={3}>"
                .format(self.__class__.__name__, self.scale, self.format,
                        getattr(self, self.format)))

    def __str__(self):
        return str(getattr(self, self.format))

    @property
    def scale(self):
        """Time scale"""
        return self._time.scale

    def _set_scale(self, scale):
        """
        This is the key routine that actually does time scale conversions.
        This is not public and not connected to the read-only scale property.
        """

        if scale == self.scale:
            return
        if scale not in self.SCALES:
            raise ValueError("Scale {0} is not in the allowed scales {1}"
                             .format(repr(scale), sorted(self.SCALES)))

        # Determine the chain of scale transformations to get from the current
        # scale to the new scale.  MULTI_HOPS contains a dict of all
        # transformations (xforms) that require intermediate xforms.
        # The MULTI_HOPS dict is keyed by (sys1, sys2) in alphabetical order.
        xform = (self.scale, scale)
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
                dt_method = '_get_delta_{0}_{1}'.format(*sys12)
                try:
                    get_dt = getattr(self, dt_method)
                except AttributeError:
                    pass
                else:
                    args.append(get_dt(jd1, jd2))
                    break

            conv_func = getattr(erfa_time, sys1 + '_' + sys2)
            jd1, jd2 = conv_func(*args)
        self._time = self.FORMATS[self.format](jd1, jd2, scale, self.precision,
                                               self.in_subfmt, self.out_subfmt,
                                               from_jd=True)

    @property
    def precision(self):
        """
        Decimal precision when outputting seconds as floating point (int
        value between 0 and 9 inclusive).
        """
        return self._precision

    @precision.setter
    def precision(self, val):
        if not isinstance(val, int) or val < 0 or val > 9:
            raise ValueError('precision attribute must be an int between '
                             '0 and 9')
        self._precision = val

    @property
    def in_subfmt(self):
        """
        Unix wildcard pattern to select subformats for parsing string input
        times
        """
        return self._in_subfmt

    @in_subfmt.setter
    def in_subfmt(self, val):
        if not isinstance(val, basestring):
            raise ValueError('in_subfmt attribute must be a string')
        self._in_subfmt = val

    @property
    def out_subfmt(self):
        """
        Unix wildcard pattern to select subformats for outputting times
        """
        return self._out_subfmt

    @out_subfmt.setter
    def out_subfmt(self, val):
        if not isinstance(val, basestring):
            raise ValueError('out_subfmt attribute must be a string')
        self._out_subfmt = val

    def _shaped_like_input(self, values):
        if self.isscalar:
            value0 = values[0]
            try:
                return value0.tolist()
            except AttributeError:
                return value0
        else:
            return values

    @property
    def jd1(self):
        """
        First of the two doubles that internally store time value(s) in JD
        """
        return self._shaped_like_input(self._time.jd1)

    @property
    def jd2(self):
        """
        Second of the two doubles that internally store time value(s) in JD
        """
        return self._shaped_like_input(self._time.jd2)

    @property
    def value(self):
        """Time value(s) in current format"""
        return self._shaped_like_input(self._time.value)

    @property
    @deprecated("0.3", name="val", alternative="value")
    def val(self):
        return self.value

    @property
    @deprecated("0.3", name="vals", alternative="value")
    def vals(self):
        """Time values in current format as a numpy array"""
        return self.value

    def copy(self, format=None):
        """
        Return a fully independent copy the Time object, optionally changing
        the format.

        If ``format`` is supplied then the time format of the returned Time
        object will be set accordingly, otherwise it will be unchanged from the
        original.

        In this method a full copy of the internal time arrays will be made.
        The internal time arrays are normally not changeable by the user so in
        most cases the ``replicate()`` method should be used.

        Parameters
        ----------
        format : str, optional
            Time format of the copy.

        Returns
        -------
        tm: Time object
            Copy of this object
        """
        return self.replicate(format, copy=True)

    def replicate(self, format=None, copy=False):
        """
        Return a replica of the Time object, optionally changing the format.

        If ``format`` is supplied then the time format of the returned Time
        object will be set accordingly, otherwise it will be unchanged from the
        original.

        If ``copy`` is set to True then a full copy of the internal time arrays
        will be made.  By default the replica will use a reference to the
        original arrays when possible to save memory.  The internal time arrays
        are normally not changeable by the user so in most cases it should not
        be necessary to set ``copy`` to True.

        The convenience method copy() is available in which ``copy`` is True
        by default.

        Parameters
        ----------
        format : str, optional
            Time format of the replica.
        copy : bool, optional
            Return a true copy instead of using references where possible.

        Returns
        -------
        tm: Time object
            Replica of this object
        """
        # To avoid recalculating integer day + fraction, no longer just
        # instantiate a new class instance, but rather do the steps by hand.
        # This also avoids quite a bit of unnecessary work in __init__
        ###  tm = self.__class__(self._time.jd1, self._time.jd2,
        ###                      format='jd', scale=self.scale, copy=copy)
        tm = super(Time, self.__class__).__new__(self.__class__)
        tm._time = TimeJD(self._time.jd1.copy() if copy else self._time.jd1,
                          self._time.jd2.copy() if copy else self._time.jd2,
                          self.scale, self.precision,
                          self.in_subfmt, self.out_subfmt, from_jd=True)
        # Optional or non-arg attributes
        attrs = ('isscalar', '_delta_ut1_utc', '_delta_tdb_tt',
                 'lat', 'lon', 'precision', 'in_subfmt', 'out_subfmt')
        for attr in attrs:
            try:
                setattr(tm, attr, getattr(self, attr))
            except AttributeError:
                pass

        if format is None:
            format = self.format

        # Make the new internal _time object corresponding to the format
        # in the copy.  If the format is unchanged this process is lightweight
        # and does not create any new arrays.

        NewFormat = tm.FORMATS[format]
        # If the new format class has a "scale" class attr then that scale is
        # required and the input jd1,2 has to be converted first.
        if hasattr(NewFormat, 'required_scale'):
            scale = NewFormat.required_scale
            new = getattr(tm, scale)  # self JDs converted to scale
            tm._time = NewFormat(new._time.jd1, new._time.jd2, scale,
                                 tm.precision,
                                 tm.in_subfmt, tm.out_subfmt,
                                 from_jd=True)
        else:
            tm._time = NewFormat(tm._time.jd1, tm._time.jd2,
                                 tm.scale, tm.precision,
                                 tm.in_subfmt, tm.out_subfmt,
                                 from_jd=True)
        tm._format = format

        return tm

    def __copy__(self):
        """
        Overrides the default behavior of the `copy.copy` function in
        the python stdlib to behave like `Time.copy`. Does *not* make a
        copy of the JD arrays - only copies by reference.
        """
        return self.replicate()

    def __deepcopy__(self, memo):
        """
        Overrides the default behavior of the `copy.deepcopy` function
        in the python stdlib to behave like `Time.copy`. Does make a
        copy of the JD arrays.
        """
        return self.copy()

    def __getitem__(self, item):
        if self.isscalar:
            raise TypeError('scalar {0!r} object is not subscriptable.'.format(
                self.__class__.__name__))
        tm = self.replicate()
        jd1 = self._time.jd1[item]
        tm.isscalar = jd1.ndim == 0

        def keepasarray(x, isscalar=tm.isscalar):
            return np.array([x]) if isscalar else x
        tm._time.jd1 = keepasarray(jd1)
        tm._time.jd2 = keepasarray(self._time.jd2[item])
        attrs = ('_delta_ut1_utc', '_delta_tdb_tt')
        for attr in attrs:
            if hasattr(self, attr):
                val = getattr(self, attr)
                setattr(tm, attr, keepasarray(val[item]))
        return tm

    def __getattr__(self, attr):
        """
        Get dynamic attributes to output format or do timescale conversion.
        """
        if attr in self.SCALES:
            tm = self.replicate()
            tm._set_scale(attr)
            return tm

        elif attr in self.FORMATS:
            tm = self.replicate(format=attr)
            if self.isscalar:
                out = tm._time.value[0]
                # convert to native python for non-object dtypes
                if tm._time.value.dtype.kind != 'O':
                    out = out.tolist()
            else:
                out = tm.value
            return out

        else:
            # Should raise AttributeError
            return self.__getattribute__(attr)

    @override__dir__
    def __dir__(self):
        return set(list(self.SCALES) + self.FORMATS.keys())

    def _match_len(self, val):
        """
        Ensure that `val` is matched to length of self.  If val has length 1
        then broadcast, otherwise cast to double and make sure length matches.
        """
        val, ndim = _make_1d_array(val, copy=True)  # be conservative and copy
        if len(val) == 1:
            oval = val
            val = np.empty(len(self), dtype=np.double)
            val[:] = oval
        elif len(val) != len(self):
            raise ValueError('Attribute length must match Time object length')
        return val

    def get_delta_ut1_utc(self, iers_table=None, return_status=False):
        """Find UT1 - UTC differences by interpolating in IERS Table.

        Parameters
        ----------
        iers_table: `~astropy.time.iers.IERS` table, optional
            Table containing UT1-UTC differences from IERS Bulletins A and/or B
            If None, use default version (see `~astropy.time.iers`)
        return_status : bool
            Whether to return status values.  If `False` (default),
            iers raises `IndexError` if any time is out of the range covered
            by the IERS table.

        Returns
        -------
        ut1_utc: float or float array
            UT1-UTC, interpolated in IERS Table
        status: int or int array
            Status values (if `return_status`=`True`)::
            `~astropy.time.iers.FROM_IERS_B`
            `~astropy.time.iers.FROM_IERS_A`
            `~astropy.time.iers.FROM_IERS_A_PREDICTION`
            `~astropy.time.iers.TIME_BEFORE_IERS_RANGE`
            `~astropy.time.iers.TIME_BEYOND_IERS_RANGE`

        Notes
        -----
        In normal usage, UT1-UTC differences are calculated automatically
        on the first instance ut1 is needed.

        Examples
        --------
        To check in code whether any times are before the IERS table range::

            >>> from astropy.time.iers import TIME_BEFORE_IERS_RANGE
            >>> t = Time(['1961-01-01', '2000-01-01'], scale='utc')
            >>> delta, status = t.get_delta_ut1_utc(return_status=True)
            >>> status == TIME_BEFORE_IERS_RANGE
            array([ True, False], dtype=bool)

        To use an updated IERS A bulletin to calculate UT1-UTC
        (see also `~astropy.time.iers`)::

            >>> from astropy.time.iers import IERS_A, IERS_A_URL
            >>> from astropy.utils.data import download_file
            >>> iers_a_file = download_file(IERS_A_URL,
            ...                             cache=True)        # doctest: +SKIP
            >>> iers_a = IERS_A.open(iers_a_file)              # doctest: +SKIP
            >>> t.delta_ut1_utc = t.get_delta_ut1_utc(iers_a)  # doctest: +SKIP

        The delta_ut1_utc property will be used to calculate t.ut1;
        raises IndexError if any of the times is out of range.

        """
        if iers_table is None:
            from .iers import IERS
            iers_table = IERS.open()

        return iers_table.ut1_utc(self.utc, return_status=return_status)

    # Property for ERFA DUT arg = UT1 - UTC
    def _get_delta_ut1_utc(self, jd1=None, jd2=None):
        """
        Get ERFA DUT arg = UT1 - UTC.  This getter takes optional jd1 and
        jd2 args because it gets called that way when converting time scales.
        If delta_ut1_utc is not yet set, this will interpolate them from the
        the IERS table.
        """

        # Sec. 4.3.1: the arg DUT is the quantity delta_UT1 = UT1 - UTC in
        # seconds. It is obtained from tables published by the IERS.
        if not hasattr(self, '_delta_ut1_utc'):
            from .iers import IERS
            iers_table = IERS.open()
            self._set_delta_ut1_utc(iers_table.ut1_utc(jd1, jd2))

        return self._delta_ut1_utc

    def _set_delta_ut1_utc(self, val):
        self._delta_ut1_utc = self._match_len(val)

    # Note can't use @property because _get_delta_tdb_tt is explicitly
    # called with the optional jd1 and jd2 args.
    delta_ut1_utc = property(_get_delta_ut1_utc, _set_delta_ut1_utc)
    """UT1 - UTC time scale offset"""

    # Property for ERFA DTR arg = TDB - TT
    def _get_delta_tdb_tt(self, jd1=None, jd2=None):
        if not hasattr(self, '_delta_tdb_tt'):
            # If jd1 and jd2 are not provided (which is the case for property
            # attribute access) then require that the time scale is TT or TDB.
            # Otherwise the computations here are not correct.
            if jd1 is None or jd2 is None:
                if self.scale not in ('tt', 'tdb'):
                    raise ValueError('Accessing the delta_tdb_tt attribute '
                                     'is only possible for TT or TDB time '
                                     'scales')
                else:
                    jd1 = self._time.jd1
                    jd2 = self._time.jd2

            # First go from the current input time (which is either
            # TDB or TT) to an approximate UTC.  Since TT and TDB are
            # pretty close (few msec?), assume TT.
            njd1, njd2 = erfa_time.tt_tai(jd1, jd2)
            njd1, njd2 = erfa_time.tai_utc(njd1, njd2)
            # TODO: actually need to go to UT1 which needs DUT.
            ut = njd1 + njd2

            # Compute geodetic params needed for d_tdb_tt()
            phi = np.radians(self.lat)
            elon = np.radians(self.lon)
            xyz = erfa_time.era_gd2gc(1, elon, phi, 0.0)
            u = np.sqrt(xyz[0] ** 2 + xyz[1] ** 2)
            v = xyz[2]

            self._delta_tdb_tt = erfa_time.d_tdb_tt(jd1, jd2, ut, elon, u, v)

        return self._delta_tdb_tt

    def _set_delta_tdb_tt(self, val):
        self._delta_tdb_tt = self._match_len(val)

    # Note can't use @property because _get_delta_tdb_tt is explicitly
    # called with the optional jd1 and jd2 args.
    delta_tdb_tt = property(_get_delta_tdb_tt, _set_delta_tdb_tt)
    """TDB - TT time scale offset"""

    def __len__(self):
        return len(self._time)

    def __sub__(self, other):
        self_tai = self.tai
        if not isinstance(other, Time):
            raise OperandTypeError(self, other)

        other_tai = other.tai
        # Note: jd1 is exact, and jd2 carry-over is done in
        # Time/TimeDelta initialisation
        jd1 = self_tai.jd1 - other_tai.jd1
        jd2 = self_tai.jd2 - other_tai.jd2

        # T      - Tdelta = T
        # Tdelta - Tdelta = Tdelta
        # T      - T      = Tdelta
        # Tdelta - T      = error
        self_delta = isinstance(self, TimeDelta)
        other_delta = isinstance(other, TimeDelta)
        self_time = not self_delta  # only 2 possibilities
        other_time = not other_delta
        if (self_delta and other_delta) or (self_time and other_time):
            out = TimeDelta(jd1, jd2, format='jd')
            if self_delta:
                out = out.replicate(format=self.format)
            return out
        elif (self_time and other_delta):
            tai = Time(jd1, jd2, format='jd', scale='tai', copy=False)
            return getattr(tai.replicate(format=self.format), self.scale)
        else:
            raise OperandTypeError(self, other)

    def __add__(self, other):
        self_tai = self.tai
        if not isinstance(other, Time):
            raise OperandTypeError(self, other)

        other_tai = other.tai
        # Note: jd1 is exact, and jd2 carry-over is done in
        # Time/TimeDelta initialisation
        jd1 = self_tai.jd1 + other_tai.jd1
        jd2 = self_tai.jd2 + other_tai.jd2

        # T      + Tdelta = T
        # Tdelta + Tdelta = Tdelta
        # T      + T      = error
        # Tdelta + T      = T
        self_delta = isinstance(self, TimeDelta)
        other_delta = isinstance(other, TimeDelta)
        self_time = not self_delta  # only 2 possibilities
        other_time = not other_delta
        if (self_delta and other_delta):
            out = TimeDelta(jd1, jd2, format='jd')
            return out.replicate(format=self.format)
        elif (self_time and other_delta) or (self_delta and other_time):
            format = self.format if self_time else other.format
            scale = self.scale if self_time else other.scale
            tai = Time(jd1, jd2, format='jd', scale='tai', copy=False)
            return getattr(tai.replicate(format=format), scale)
        else:
            raise OperandTypeError(self, other)

    def _tai_difference(self, other):
        """If other is of same class as self, return difference in TAI.
        Otherwise, raise OperandTypeError.
        """
        if other.__class__ is not self.__class__:
            raise OperandTypeError(self, other)
        self_tai = self.tai
        other_tai = other.tai
        return (self_tai.jd1 - other_tai.jd1) + (self_tai.jd2 - other_tai.jd2)

    def __lt__(self, other):
        return self._tai_difference(other) < 0.

    def __le__(self, other):
        return self._tai_difference(other) <= 0.

    def __eq__(self, other):
        return self._tai_difference(other) == 0.

    def __ne__(self, other):
        return self._tai_difference(other) != 0.

    def __gt__(self, other):
        return self._tai_difference(other) > 0.

    def __ge__(self, other):
        return self._tai_difference(other) >= 0.


class TimeDelta(Time):
    """
    Represent the time difference between two times.

    A TimeDelta object is initialized with one or more times in the ``val``
    argument.  The input times in ``val`` must conform to the specified
    ``format``.  The optional ``val2`` time input should be supplied only for
    numeric input formats (e.g. JD) where very high precision (better than
    64-bit precision) is required.

    Parameters
    ----------
    val : numpy ndarray, list, str, number, or `~astropy.time.TimeDelta` object
        Data to initialize table.
    val2 : numpy ndarray, list, str, or number; optional
        Data to initialize table.
    format : str, optional
        Format of input value(s)
    copy : bool, optional
        Make a copy of the input values
    """
    SCALES = TIME_DELTA_SCALES
    """List of time delta scales"""

    FORMATS = TIME_DELTA_FORMATS
    """Dict of time delta formats"""

    def __init__(self, val, val2=None, format=None, scale=None, copy=False):
        # Note: scale is not used but is needed because of the inheritance
        # from Time.
        if not isinstance(val, self.__class__):
            self._init_from_vals(val, val2, format, 'tai', copy)

    def __neg__(self):
        """Negation of a `TimeDelta` object."""
        new = self.copy()
        new._time.jd1 = -self._time.jd1
        new._time.jd2 = -self._time.jd2
        return new

    def __abs__(self):
        """Absolute value of a `TimeDelta` object."""
        jd1, jd2 = self._time.jd1, self._time.jd2
        negative = jd1 + jd2 < 0
        new = self.copy()
        new._time.jd1 = np.where(negative, -jd1, jd1)
        new._time.jd2 = np.where(negative, -jd2, jd2)
        return new

    def __mul__(self, other):
        """Multiplication of `TimeDelta` objects by numbers/arrays."""
        # check needed since otherwise the self.jd1 * other multiplication
        # would enter here again (via __rmul__)
        if isinstance(other, Time):
            raise OperandTypeError(self, other)

        jd1, jd2 = day_frac(self.jd1, self.jd2, factor=other)
        out = TimeDelta(jd1, jd2, format='jd')
        if self.format != 'jd':
            out = out.replicate(format=self.format)
        return out

    def __rmul__(self, other):
        """Multiplication of numbers/arrays with `TimeDelta` objects."""
        return self.__mul__(other)

    def __div__(self, other):
        """Division of `TimeDelta` objects by numbers/arrays."""
        return self.__truediv__(other)

    def __truediv__(self, other):
        """Division of `TimeDelta` objects by numbers/arrays."""
        # cannot do __mul__(1./other) as that looses precision
        jd1, jd2 = day_frac(self.jd1, self.jd2, divisor=other)
        out = TimeDelta(jd1, jd2, format='jd')
        if self.format != 'jd':
            out = out.replicate(format=self.format)
        return out


class TimeFormat(object):
    """
    Base class for time representations.

    Parameters
    ----------
    val1 : numpy ndarray, list, str, or number
        Data to initialize table.
    val2 : numpy ndarray, list, str, or number; optional
        Data to initialize table.
    scale : str
        Time scale of input value(s)
    precision : int
        Precision for seconds as floating point
    in_subfmt : str
        Select subformat for inputting string times
    out_subfmt : str
        Select subformat for outputting string times
    from_jd : bool
        If true then val1, val2 are jd1, jd2
    """
    def __init__(self, val1, val2, scale, precision,
                 in_subfmt, out_subfmt, from_jd=False):
        self.scale = scale  # validation of scale done later with _check_scale
        self.precision = precision
        self.in_subfmt = in_subfmt
        self.out_subfmt = out_subfmt
        if len(val1) != len(val2):
            raise ValueError('Input val1 and val2 must match in length')

        if from_jd:
            self.jd1 = val1
            self.jd2 = val2
        else:
            self._check_val_type(val1, val2)
            self.set_jds(val1, val2)

    def __len__(self):
        return len(self.jd1)

    @property
    def scale(self):
        """Time scale"""
        self._scale = self._check_scale(self._scale)
        return self._scale

    @scale.setter
    def scale(self, val):
        self._scale = val

    def _check_val_type(self, val1, val2):
        """Input value validation, typically overridden by derived classes"""
        if val1.dtype.type != np.double or val2.dtype.type != np.double:
            raise TypeError('Input values for {0} class must be doubles'
                            .format(self.name))

    def _check_scale(self, scale):
        """
        Return a validated scale value.

        If there is a class attribute 'scale' then that defines the default /
        required time scale for this format.  In this case if a scale value was
        provided that needs to match the class default, otherwise return
        the class default.

        Otherwise just make sure that scale is in the allowed list of
        scales.  Provide a different error message if None (no value) was
        supplied.
        """
        if hasattr(self.__class__, 'epoch_scale') and scale is None:
            scale = self.__class__.epoch_scale

        if scale not in TIME_SCALES:
            if scale is None:
                raise ScaleValueError("No scale value supplied but it is "
                                      "required for class {0}"
                                      .format(self.__class__.__name__))
            raise ScaleValueError("Scale value '{0}' not in "
                                  "allowed values {1}"
                                  .format(scale, TIME_SCALES))

        return scale

    def set_jds(self, val1, val2):
        """
        Set internal jd1 and jd2 from val1 and val2.  Must be provided
        by derived classes.
        """
        raise NotImplementedError

    @property
    def value(self):
        """
        Return time representation from internal jd1 and jd2.  Must be
        provided by derived classes.
        """
        raise NotImplementedError


class TimeJD(TimeFormat):
    """Julian Date time format"""
    name = 'jd'

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2)

    @property
    def value(self):
        return self.jd1 + self.jd2


class TimeMJD(TimeFormat):
    """Modified Julian Date time format"""
    name = 'mjd'

    def set_jds(self, val1, val2):
        # TODO - this routine and vals should be Cythonized to follow the ERFA
        # convention of preserving precision by adding to the larger of the two
        # values in a vectorized operation.  But in most practical cases the
        # first one is probably biggest.
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2)
        self.jd1 += MJD_ZERO

    @property
    def value(self):
        return (self.jd1 - MJD_ZERO) + self.jd2


class TimeFromEpoch(TimeFormat):
    """
    Base class for times that represent the interval from a particular
    epoch as a floating point multiple of a unit time interval (e.g. seconds
    or days).
    """
    def __init__(self, val1, val2, scale, precision,
                 in_subfmt, out_subfmt, from_jd=False):
        self.scale = scale
        # Initialize the reference epoch which is a single time defined in subclasses
        epoch = Time(self.epoch_val, self.epoch_val2, scale=self.epoch_scale,
                     format=self.epoch_format)
        self.epoch = epoch

        # Now create the TimeFormat object as normal
        super(TimeFromEpoch, self).__init__(val1, val2, scale, precision,
                                            in_subfmt, out_subfmt, from_jd)

    def set_jds(self, val1, val2):
        """
        Initialize the internal jd1 and jd2 attributes given val1 and val2.  For an
        TimeFromEpoch subclass like TimeUnix these will be floats giving the effective
        seconds since an epoch time (e.g. 1970-01-01 00:00:00).
        """
        # Form new JDs based on epoch time + time from epoch (converted to JD).
        # One subtlety that might not be obvious is that 1.000 Julian days in UTC
        # can be 86400 or 86401 seconds.  For the TimeUnix format the assumption
        # is that every day is exactly 86400 seconds, so in principle this
        # is doing the math incorrectly, *except* that it matches the definition
        # of Unix time which does not include leap seconds.

        # note: use divisor=1./self.unit, since this is either 1 or 1/86400,
        # and 1/86400 is not exactly representable as a float64, so multiplying
        # by that will cause rounding errors. (But inverting it as a float64
        # recovers the exact number)
        day, frac = day_frac(val1, val2, divisor=1. / self.unit)

        jd1 = self.epoch.jd1 + day
        jd2 = self.epoch.jd2 + frac

        # Create a temporary Time object corresponding to the current new (jd1, jd2) in
        # the epoch scale (e.g. UTC for TimeUnix) then convert that to the desired time
        # scale for this object.
        #
        # A known limitation is that the transform from self.epoch_scale to self.scale
        # cannot involve any metadata like lat or lon.
        tm = getattr(Time(jd1, jd2, scale=self.epoch_scale, format='jd'), self.scale)

        self.jd1 = tm.jd1
        self.jd2 = tm.jd2

    @property
    def value(self):
        # Create a temporary Time object corresponding to the parent Time, then transform
        # that to the epoch scale.  From there do the simple math to compute delta time
        # from the epoch in Julian days, and then in the desired output units
        # (e.g. seconds).
        #
        # A known limitation is that the transform from self.scale to self.epoch_scale
        # cannot involve any metadata like lat or lon.
        tm = getattr(Time(self.jd1, self.jd2, scale=self.scale, format='jd'), self.epoch_scale)
        time_from_epoch = ((tm.jd1 - self.epoch.jd1) +
                           (tm.jd2 - self.epoch.jd2)) / self.unit
        return time_from_epoch


class TimeUnix(TimeFromEpoch):
    """
    Unix time: seconds from 1970-01-01 00:00:00 UTC.

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


class TimeCxcSec(TimeFromEpoch):
    """Chandra X-ray Center seconds from 1998-01-01 00:00:00 TT"""
    name = 'cxcsec'
    unit = 1.0 / SECS_PER_DAY  # in days (1 day == 86400 seconds)
    epoch_val = '1998-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'tt'
    epoch_format = 'iso'


class TimeGPS(TimeFromEpoch):
    """GPS time: seconds from 1980-01-06 00:00:00 UTC

    GPS time includes leap seconds.  For details, see
    http://tycho.usno.navy.mil/gpstt.html
    """
    name = 'gps'
    unit = 1.0 / SECS_PER_DAY  # in days (1 day == 86400 seconds)
    epoch_val = '1980-01-06 00:00:19'
    # above epoch is the same as Time('1980-01-06 00:00:00', scale='utc').tai
    epoch_val2 = None
    epoch_scale = 'tai'
    epoch_format = 'iso'


class TimePlotDate(TimeFromEpoch):
    """
    Matplotlib `~matplotlib.pyplot.plot_date` input: 1 + number of days from 0001-01-01 00:00:00 UTC

    This can be used directly in the matplotlib `~matplotlib.pyplot.plot_date` function::

      >>> import matplotlib.pyplot as plt
      >>> jyear = np.linspace(2000, 2001, 20)
      >>> t = Time(jyear, format='jyear', scale='utc')
      >>> plt.plot_date(t.plot_date, jyear)
      >>> plt.gcf().autofmt_xdate()  # orient date labels at a slant
      >>> plt.draw()
    """
    # This corresponds to the zero reference time for matplotlib plot_date().
    # Note that TAI and UTC are equivalent at the reference time, but
    # specifying epoch_scale = 'utc' here generates WARNINGS when the
    # class is first used.  Just use 'tai' instead.
    name = 'plot_date'
    unit = 1.0
    epoch_val = 1721424.5  # Time('0001-01-01 00:00:00', scale='tai').jd - 1
    epoch_val2 = None
    epoch_scale = 'tai'
    epoch_format = 'jd'


class TimeUnique(TimeFormat):
    """
    Base class for time formats that can uniquely create a time object
    without requiring an explicit format specifier.  This class does
    nothing but provide inheritance to identify a class as unique.
    """
    pass


class TimeAstropyTime(TimeUnique):
    """
    Instantiate date from an Astropy Time object (or list thereof).

    This is purely for instantiating from a Time object.  The output
    format is the same as the first time instance.
    """
    name = 'astropy_time'

    def __new__(cls, val1, val2, scale, precision,
                in_subfmt, out_subfmt, from_jd=False):
        """
        Use __new__ instead of __init__ to output a class instance that
        is the same as the class of the first Time object in the list.
        """

        if not all(isinstance(val, Time) for val in val1):
            raise TypeError('Input values for {0} class must be datetime objects'
                            .format(cls.name))

        if scale is None:
            scale = val1[0].scale
        jd1 = np.concatenate([getattr(val, scale)._time.jd1 for val in val1])
        jd2 = np.concatenate([getattr(val, scale)._time.jd2 for val in val1])
        OutTimeFormat = val1[0]._time.__class__
        self = OutTimeFormat(jd1, jd2, scale, precision, in_subfmt, out_subfmt, from_jd=True)

        return self


class TimeDatetime(TimeUnique):
    """
    Represent date as Python standard library `~datetime.datetime` object

    Example::

      >>> from datetime import datetime
      >>> t = Time(datetime(2000, 1, 2, 12, 0, 0), scale='utc')
      >>> t.iso
      '2000-01-02 12:00:00.000'
      >>> t.tt.datetime
      datetime.datetime(2000, 1, 2, 12, 1, 4, 184000)
    """
    name = 'datetime'

    def _check_val_type(self, val1, val2):
        if not all(isinstance(val, datetime) for val in val1):
            raise TypeError('Input values for {0} class must be datetime objects'
                            .format(self.name))
            # Note: don't care about val2 for this classes

    def set_jds(self, val1, val2):
        """Convert datetime object contained in val1 to jd1, jd2"""
        n_times = len(val1)
        iy = np.empty(n_times, dtype=np.intc)
        im = np.empty(n_times, dtype=np.intc)
        id = np.empty(n_times, dtype=np.intc)
        ihr = np.empty(n_times, dtype=np.intc)
        imin = np.empty(n_times, dtype=np.intc)
        dsec = np.empty(n_times, dtype=np.double)

        # Iterate through the datetime objects
        for i, val in enumerate(val1):
            iy[i] = val.year
            im[i] = val.month
            id[i] = val.day
            ihr[i] = val.hour
            imin[i] = val.minute
            dsec[i] = val.second + val.microsecond / 1e6

        self.jd1, self.jd2 = erfa_time.dtf_jd(self.scale.upper().encode('utf8'),
                                              iy, im, id, ihr, imin, dsec)

    @property
    def value(self):
        iys, ims, ids, ihmsfs = erfa_time.jd_dtf(self.scale.upper()
                                                 .encode('utf8'),
                                                 6,  # precision = 6 for microseconds
                                                 self.jd1, self.jd2)

        out = np.empty(len(self), dtype=np.object)
        idxs = itertools.count()
        for idx, iy, im, id, ihmsf in itertools.izip(idxs, iys, ims, ids, ihmsfs):
            ihr, imin, isec, ifracsec = ihmsf
            out[idx] = datetime(int(iy), int(im), int(id),
                                int(ihr), int(imin), int(isec), int(ifracsec))

        return out


class TimeString(TimeUnique):
    """
    Base class for string-like time represetations.

    This class assumes that anything following the last decimal point to the
    right is a fraction of a second.

    This is a reference implementation can be made much faster with effort.
    """
    def _check_val_type(self, val1, val2):
        if val1.dtype.kind not in ('S', 'U'):
            raise TypeError('Input values for {0} class must be strings'
                            .format(self.name))
            # Note: don't care about val2 for these classes

    def set_jds(self, val1, val2):
        """Parse the time strings contained in val1 and set jd1, jd2"""
        n_times = len(val1)  # val1,2 already checked to have same len
        iy = np.empty(n_times, dtype=np.intc)
        im = np.empty(n_times, dtype=np.intc)
        id = np.empty(n_times, dtype=np.intc)
        ihr = np.empty(n_times, dtype=np.intc)
        imin = np.empty(n_times, dtype=np.intc)
        dsec = np.empty(n_times, dtype=np.double)

        # Select subformats based on current self.in_subfmt
        subfmts = self._select_subfmts(self.in_subfmt)

        for i, timestr in enumerate(val1):
            # Assume that anything following "." on the right side is a
            # floating fraction of a second.
            try:
                idot = timestr.rindex('.')
            except:
                fracsec = 0.0
            else:
                timestr, fracsec = timestr[:idot], timestr[idot:]
                fracsec = float(fracsec)

            for _, strptime_fmt, _ in subfmts:
                try:
                    tm = time.strptime(timestr, strptime_fmt)
                except ValueError:
                    pass
                else:
                    iy[i] = tm.tm_year
                    im[i] = tm.tm_mon
                    id[i] = tm.tm_mday
                    ihr[i] = tm.tm_hour
                    imin[i] = tm.tm_min
                    dsec[i] = tm.tm_sec + fracsec
                    break
            else:
                raise ValueError('Time {0} does not match {1} format'
                                 .format(timestr, self.name))

        self.jd1, self.jd2 = erfa_time.dtf_jd(self.scale.upper().encode('utf8'),
                                              iy, im, id, ihr, imin, dsec)

    def str_kwargs(self):
        """
        Generator that yields a dict of values corresponding to the
        calendar date and time for the internal JD values.
        """
        iys, ims, ids, ihmsfs = erfa_time.jd_dtf(self.scale.upper()
                                                 .encode('utf8'),
                                                 self.precision,
                                                 self.jd1, self.jd2)

        # Get the str_fmt element of the first allowed output subformat
        _, _, str_fmt = self._select_subfmts(self.out_subfmt)[0]

        if '{yday:' in str_fmt:
            has_yday = True
        else:
            has_yday = False
            yday = None

        for iy, im, id, ihmsf in itertools.izip(iys, ims, ids, ihmsfs):
            ihr, imin, isec, ifracsec = ihmsf
            if has_yday:
                yday = datetime(iy, im, id).timetuple().tm_yday

            yield {'year': int(iy), 'mon': int(im), 'day': int(id),
                   'hour': int(ihr), 'min': int(imin), 'sec': int(isec),
                   'fracsec': int(ifracsec), 'yday': yday}

    @property
    def value(self):
        # Select the first available subformat based on current
        # self.out_subfmt
        subfmts = self._select_subfmts(self.out_subfmt)
        _, _, str_fmt = subfmts[0]

        # TODO: fix this ugly hack
        if self.precision > 0 and str_fmt.endswith('{sec:02d}'):
            str_fmt += '.{fracsec:0' + str(self.precision) + 'd}'

        # Try to optimize this later.  Can't pre-allocate because length of
        # output could change, e.g. year rolls from 999 to 1000.
        outs = []
        for kwargs in self.str_kwargs():
            outs.append(str_fmt.format(**kwargs))

        return np.array(outs)

    def _select_subfmts(self, pattern):
        """
        Return a list of subformats where name matches ``pattern`` using
        fnmatch.
        """
        from fnmatch import fnmatchcase
        subfmts = [x for x in self.subfmts if fnmatchcase(x[0], pattern)]
        if len(subfmts) == 0:
            raise ValueError('No subformats match {0}'.format(pattern))
        return subfmts


class TimeISO(TimeString):
    """
    ISO 8601 compliant date-time format "YYYY-MM-DD HH:MM:SS.sss...".

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """

    name = 'iso'
    subfmts = (('date_hms',
                '%Y-%m-%d %H:%M:%S',
                # XXX To Do - use strftime for output ??
                '{year:d}-{mon:02d}-{day:02d} {hour:02d}:{min:02d}:{sec:02d}'),
               ('date_hm',
                '%Y-%m-%d %H:%M',
                '{year:d}-{mon:02d}-{day:02d} {hour:02d}:{min:02d}'),
               ('date',
                '%Y-%m-%d',
                '{year:d}-{mon:02d}-{day:02d}'))


class TimeISOT(TimeString):
    """
    ISO 8601 compliant date-time format "YYYY-MM-DDTHH:MM:SS.sss...".
    This is the same as TimeISO except for a "T" instead of space between
    the date and time.

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """

    name = 'isot'
    subfmts = (('date_hms',
                '%Y-%m-%dT%H:%M:%S',
                '{year:d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}:{sec:02d}'),
               ('date_hm',
                '%Y-%m-%dT%H:%M',
                '{year:d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}'),
               ('date',
                '%Y-%m-%d',
                '{year:d}-{mon:02d}-{day:02d}'))


class TimeYearDayTime(TimeString):
    """
    Year, day-of-year and time as "YYYY:DOY:HH:MM:SS.sss...".  The
    day-of-year (DOY) goes from 001 to 365 (366 in leap years).

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """

    name = 'yday'
    subfmts = (('date_hms',
                '%Y:%j:%H:%M:%S',
                '{year:d}:{yday:03d}:{hour:02d}:{min:02d}:{sec:02d}'),
               ('date_hm',
                '%Y:%j:%H:%M',
                '{year:d}:{yday:03d}:{hour:02d}:{min:02d}'),
               ('date',
                '%Y:%j',
                '{year:d}:{yday:03d}'))


class TimeEpochDate(TimeFormat):
    """
    Base class for support floating point Besselian and Julian epoch dates
    """
    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # validate scale.
        epoch_to_jd = getattr(erfa_time, self.epoch_to_jd)
        jd1, jd2 = epoch_to_jd(val1 + val2)
        self.jd1, self.jd2 = day_frac(jd1, jd2)

    @property
    def value(self):
        jd_to_epoch = getattr(erfa_time, self.jd_to_epoch)
        return jd_to_epoch(self.jd1, self.jd2)


class TimeBesselianEpoch(TimeEpochDate):
    """Besselian Epoch year as floating point value(s) like 1950.0"""
    name = 'byear'
    epoch_to_jd = 'besselian_epoch_jd'
    jd_to_epoch = 'jd_besselian_epoch'


class TimeJulianEpoch(TimeEpochDate):
    """Julian Epoch year as floating point value(s) like 2000.0"""
    name = 'jyear'
    epoch_to_jd = 'julian_epoch_jd'
    jd_to_epoch = 'jd_julian_epoch'


class TimeEpochDateString(TimeString):
    """
    Base class to support string Besselian and Julian epoch dates
    such as 'B1950.0' or 'J2000.0' respectively.
    """
    def set_jds(self, val1, val2):
        years = np.empty(len(val1), dtype=np.double)
        epoch_prefix = self.epoch_prefix

        for i, time_str in enumerate(val1):
            try:
                epoch_type, year_str = time_str[0], time_str[1:]
                year = float(year_str)
                if epoch_type.upper() != epoch_prefix:
                    raise ValueError
            except (IndexError, ValueError) as err:
                raise ValueError('Time {0} does not match {1} format'
                                 .format(time_str, self.name))
            else:
                years[i] = year

        self._check_scale(self._scale)  # validate scale.
        epoch_to_jd = getattr(erfa_time, self.epoch_to_jd)
        self.jd1, self.jd2 = epoch_to_jd(years)

    @property
    def value(self):
        jd_to_epoch = getattr(erfa_time, self.jd_to_epoch)
        years = jd_to_epoch(self.jd1, self.jd2)
        # Use old-style format since it is a factor of 2 faster
        str_fmt = self.epoch_prefix + '%.' + str(self.precision) + 'f'
        outs = [str_fmt % year for year in years]
        return np.array(outs)


class TimeBesselianEpochString(TimeEpochDateString):
    """Besselian Epoch year as string value(s) like 'B1950.0'"""
    name = 'byear_str'
    epoch_to_jd = 'besselian_epoch_jd'
    jd_to_epoch = 'jd_besselian_epoch'
    epoch_prefix = 'B'


class TimeJulianEpochString(TimeEpochDateString):
    """Julian Epoch year as string value(s) like 'J2000.0'"""
    name = 'jyear_str'
    epoch_to_jd = 'julian_epoch_jd'
    jd_to_epoch = 'jd_julian_epoch'
    epoch_prefix = 'J'


class TimeDeltaFormat(TimeFormat):
    """Base class for time delta representations"""
    pass


class TimeDeltaSec(TimeDeltaFormat):
    """Time delta in SI seconds"""
    name = 'sec'

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2, divisor=SECS_PER_DAY)

    @property
    def value(self):
        return (self.jd1 + self.jd2) * SECS_PER_DAY


class TimeDeltaJD(TimeDeltaFormat):
    """Time delta in Julian days (86400 SI seconds)"""
    name = 'jd'

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2)

    @property
    def value(self):
        return self.jd1 + self.jd2


class ScaleValueError(Exception):
    pass


# Set module constant with names of all available time formats
for name, val in locals().items():
    try:
        is_timeformat = issubclass(val, TimeFormat)
        is_timedeltaformat = issubclass(val, TimeDeltaFormat)
    except:
        pass
    else:
        if hasattr(val, 'name'):
            if is_timedeltaformat:
                TIME_DELTA_FORMATS[val.name] = val
            elif is_timeformat:
                TIME_FORMATS[val.name] = val


def _make_1d_array(val, copy=False):
    """
    Take ``val`` and convert/reshape to a 1-d array.  If ``copy`` is True
    then copy input values.

    Returns
    -------
    val, val_ndim: ndarray, int
        Array version of ``val`` and the number of dims in original.
    """
    val = np.array(val, copy=copy)
    val_ndim = val.ndim  # remember original ndim
    if val.ndim == 0:
        val = (val.reshape(1) if val.dtype.kind == 'O' else np.array([val]))
    elif val_ndim > 1:
        # Maybe lift this restriction later to allow multi-dim in/out?
        raise TypeError('Input val must be zero or one dimensional')

    # Allow only string or float arrays as input (XXX datetime later...)
    if val.dtype.kind == 'i':
        val = np.asarray(val, dtype=np.float64)

    return val, val_ndim


def day_frac(val1, val2, factor=1., divisor=1.):
    """
    Return the sum of ``val1`` and ``val2`` as two float64s, an integer part and the
    fractional remainder.  If ``factor`` is not 1.0 then multiply the sum by ``factor``.
    If ``divisor`` is not 1.0 then divide the sum by ``divisor``.

    The arithmetic is all done with exact floating point operations so no precision is
    lost to rounding error.  This routine assumes the sum is less than about 1e16,
    otherwise the ``frac`` part will be greater than 1.0.

    Returns
    -------
    day, frac: float64
        Integer and fractional part of val1 + val2.
    """

    # Add val1 and val2 exactly, returning the result as two float64s.  The first is the
    # approximate sum (with some floating point error) and the second is the error of the
    # float64 sum.
    sum12, err12 = two_sum(val1, val2)

    if np.any(factor != 1.):
        sum12, carry = two_product(sum12, factor)
        carry += err12 * factor
        sum12, err12 = two_sum(sum12, carry)

    if np.any(divisor != 1.):
        q1 = sum12 / divisor
        p1, p2 = two_product(q1, divisor)
        d1, d2 = two_sum(sum12, -p1)
        d2 += err12
        d2 -= p2
        q2 = (d1 + d2) / divisor  # 3-part float fine here; nothing can be lost
        sum12, err12 = two_sum(q1, q2)

    # get integer fraction
    day = np.round(sum12)
    extra, frac = two_sum(sum12, -day)
    frac += extra + err12
    return day, frac


def two_sum(a, b):
    """
    Add ``a`` and ``b`` exactly, returning the result as two float64s.  The first is the
    approximate sum (with some floating point error) and the second is the error of the
    float64 sum.

    Using the procedure of Shewchuk, 1997,
    Discrete & Computational Geometry 18(3):305-363
    http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

    Returns
    -------
    sum, err: float64
        Approximate sum of a + b and the exact floating point error
    """
    x = a + b
    eb = x - a
    eb = b - eb
    ea = x - b
    ea = a - ea
    return x, ea + eb


def two_product(a, b):
    """
    Multiple ``a`` and ``b`` exactly, returning the result as two float64s.  The first is
    the approximate prodcut (with some floating point error) and the second is the error
    of the float64 product.

    Uses the procedure of Shewchuk, 1997,
    Discrete & Computational Geometry 18(3):305-363
    http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

    Returns
    -------
    prod, err: float64
        Approximate product a * b and the exact floating point error
    """
    x = a * b
    ah, al = split(a)
    bh, bl = split(b)
    y1 = ah * bh
    y = x - y1
    y2 = al * bh
    y -= y2
    y3 = ah * bl
    y -= y3
    y4 = al * bl
    y = y4 - y
    return x, y


def split(a):
    """
    Split float64 in two aligned parts.

    Uses the procedure of Shewchuk, 1997,
    Discrete & Computational Geometry 18(3):305-363
    http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

    """
    c = 134217729. * a  # 2**27+1.
    abig = c - a
    ah = c - abig
    al = a - ah
    return ah, al


class OperandTypeError(TypeError):
    def __init__(self, left, right):
        self.value = ("unsupported operand type(s) for -: "
                      "'{0}' and '{1}'".format(left.__class__.__name__,
                                               right.__class__.__name__))
