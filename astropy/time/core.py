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

__all__ = ['Time', 'TimeDelta', 'TimeFormat', 'TimeJD', 'TimeMJD',
           'TimeFromEpoch', 'TimeUnix', 'TimeCxcSec', 'TimePlotDate',
           'TimeDatetime',
           'TimeString', 'TimeISO', 'TimeISOT', 'TimeYearDayTime', 'TimeEpochDate',
           'TimeBesselianEpoch', 'TimeJulianEpoch', 'TimeDeltaFormat',
           'TimeDeltaSec', 'TimeDeltaJD', 'ScaleValueError',
           'OperandTypeError', 'TimeEpochDateString',
           'TimeBesselianEpochString', 'TimeJulianEpochString',
           'TIME_FORMATS', 'TIME_DELTA_FORMATS', 'TIME_SCALES',
           'TIME_DELTA_SCALES']

try:
    # Not guaranteed available at setup time
    from . import sofa_time
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
    opt : dict, optional
        options
    lat : float, optional
        Earth latitude of observer (decimal degrees)
    lon : float, optional
        Earth longitude of observer (decimal degrees)
    copy : bool, optional
        Make a copy of the input values
    """

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

        self.is_scalar = (val_ndim == 0)
        if val_ndim != val2_ndim:
            raise ValueError('Input val and val2 must have same dimensions')

        if scale is not None and scale not in self.SCALES:
            raise ScaleValueError("Scale {0} is not in the allowed scales {1}"
                                  .format(repr(scale), sorted(self.SCALES)))

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

        if format is None and val.dtype.kind in ('S', 'U', 'O'):
            formats = [(name, cls) for name, cls in self.FORMATS.items()
                       if issubclass(cls, TimeUnique)]
            err_msg = 'any of the formats where the format keyword is optional {0}'.format(
                [name for name, cls in formats])
        elif format not in self.FORMATS:
            if format is None:
                raise ValueError("No time format was given, and the input is not unique")
            else:
                raise ValueError("Format {0} is not one of the allowed "
                                 "formats {1}".format(repr(format), sorted(self.FORMATS)))
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
        #call `utcnow` immediately to be sure it's ASAP
        dtnow = datetime.utcnow()
        return cls(val=dtnow, format='datetime', scale='utc')

    @property
    def format(self):
        """Time format"""
        return self._format

    def __repr__(self):
        return ("<{0} object: scale='{1}' format='{2}' vals={3}>"
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

            conv_func = getattr(sofa_time, sys1 + '_' + sys2)
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

    def _shaped_like_input(self, vals):
        return (vals[0].tolist() if self.is_scalar else vals)

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
    def val(self):
        """Time value(s) in current format"""
        return self._shaped_like_input(self._time.vals)

    @property
    def vals(self):
        """Time values in current format as a numpy array"""
        return self._time.vals

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
        tm = self.__class__(self._time.jd1, self._time.jd2,
                            format='jd', scale=self.scale, copy=copy)

        # Optional or non-arg attributes
        attrs = ('is_scalar', '_delta_ut1_utc', '_delta_tdb_tt',
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

    def _getAttributeNames(self):
        """
        Add dynamic attribute names for IPython completer.
        """
        return list(self.SCALES) + self.FORMATS.keys()

    def __getattr__(self, attr):
        """
        Get dynamic attributes to output format or do timescale conversion.
        """
        # The following is needed for the IPython completer
        if attr == 'trait_names':
            return []

        if attr in self.SCALES:
            tm = self.replicate()
            tm._set_scale(attr)
            return tm

        elif attr in self.FORMATS:
            tm = self.replicate(format=attr)
            if self.is_scalar:
                out = tm.vals[0]
                # convert to native python for non-object dtypes
                if tm.vals.dtype.kind != 'O':
                    out = out.tolist()
            else:
                out = tm.vals
            return out

        else:
            # Should raise AttributeError
            return self.__getattribute__(attr)

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

    # Property for SOFA DUT arg = UT1 - UTC
    def _get_delta_ut1_utc(self, jd1=None, jd2=None):
        """
        Get SOFA DUT arg = UT1 - UTC.  This getter takes optional jd1 and
        jd2 args because it gets called that way when converting time scales.
        The current code ignores these, but when the IERS table is interpolated
        by this module they will be used.
        """

        # Sec. 4.3.1: the arg DUT is the quantity delta_UT1 = UT1 - UTC in
        # seconds. It can be obtained from tables published by the IERS.
        # TODO - get that table when needed and interpolate or whatever.
        if not hasattr(self, '_delta_ut1_utc'):
            # Exception until the IERS table is available to the package
            raise ValueError('Must set the delta_ut1_utc attribute in '
                             'order to do the scale transform.')

        return self._delta_ut1_utc

    def _set_delta_ut1_utc(self, val):
        self._delta_ut1_utc = self._match_len(val)

    # Note can't use @property because _get_delta_tdb_tt is explicitly
    # called with the optional jd1 and jd2 args.
    delta_ut1_utc = property(_get_delta_ut1_utc, _set_delta_ut1_utc)
    """UT1 - UTC time scale offset"""

    # Property for SOFA DTR arg = TDB - TT
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
            njd1, njd2 = sofa_time.tt_tai(jd1, jd2)
            njd1, njd2 = sofa_time.tai_utc(njd1, njd2)
            # TODO: actually need to go to UT1 which needs DUT.
            ut = njd1 + njd2

            # Compute geodetic params needed for d_tdb_tt()
            phi = np.radians(self.lat)
            elon = np.radians(self.lon)
            xyz = sofa_time.iau_gd2gc(1, elon, phi, 0.0)
            u = np.sqrt(xyz[0] ** 2 + xyz[1] ** 2)
            v = xyz[2]

            self._delta_tdb_tt = sofa_time.d_tdb_tt(jd1, jd2, ut, elon, u, v)

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
            return TimeDelta(jd1, jd2, format='jd')
        elif (self_time and other_delta):
            self_tai._time.jd1 = jd1
            self_tai._time.jd2 = jd2
            return getattr(self_tai, self.scale)
        else:
            raise OperandTypeError(self, other)

    def __add__(self, other):
        self_tai = self.tai
        if not isinstance(other, Time):
            raise OperandTypeError(self, other)

        other_tai = other.tai
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
            return TimeDelta(jd1, jd2, format='jd')
        elif (self_time and other_delta) or (self_delta and other_time):
            tai = self_tai if self_time else other_tai
            scale = self.scale if self_time else other.scale
            tai._time.jd1 = jd1
            tai._time.jd2 = jd2
            return getattr(tai, scale)
        else:
            raise OperandTypeError(self, other)


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
    def vals(self):
        """
        Return time representation from internal jd1 and jd2.  Must be
        provided by by derived classes.
        """
        raise NotImplementedError


class TimeJD(TimeFormat):
    """Julian Date time format"""
    name = 'jd'

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1 = val1
        self.jd2 = val2

    @property
    def vals(self):
        return self.jd1 + self.jd2


class TimeMJD(TimeFormat):
    """Modified Julian Date time format"""
    name = 'mjd'

    def set_jds(self, val1, val2):
        # TODO - this routine and vals should be Cythonized to follow the SOFA
        # convention of preserving precision by adding to the larger of the two
        # values in a vectorized operation.  But in most practical cases the
        # first one is probably biggest.
        self._check_scale(self._scale)  # Validate scale.
        self.jd1 = val1 + MJD_ZERO
        self.jd2 = val2

    @property
    def vals(self):
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
        # Convert the epoch time to TAI so that calculations involving adding/subtracting
        # JDs are always linear and well-defined (a JD in TAI is always 86400 uniform
        # seconds, unlike other scales, in particular UTC where one JD isn't always the
        # same duration).
        self.epoch = epoch
        super(TimeFromEpoch, self).__init__(val1, val2, scale, precision,
                                            in_subfmt, out_subfmt, from_jd)

    def set_jds(self, val1, val2):
        # Form new JDs based on epoch time (TAI) + time from epoch (converted to JD)
        jd1 = self.epoch.jd1 + val2 * self.unit
        jd2 = self.epoch.jd2 + val1 * self.unit

        # Create a Time object corresponding to the new (jd1, jd2) in the epoch scale,
        # then convert that to the time scale for this object.
        tm = getattr(Time(jd1, jd2, scale=self.epoch_scale, format='jd'), self.scale)

        self.jd1 = tm.jd1
        self.jd2 = tm.jd2

    @property
    def vals(self):
        # Convert current JDs to TAI
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


class TimePlotDate(TimeFromEpoch):
    """
    Matplotlib `~matplotlib.pyplot.plot_date` input: 1 + number of days from 0001-01-01 00:00:00 UTC

    This can be used directly in the matplotlib `~matplotlib.pyplot.plot_date` function::

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

        self.jd1, self.jd2 = sofa_time.dtf_jd(self.scale.upper().encode('utf8'),
                                              iy, im, id, ihr, imin, dsec)

    @property
    def vals(self):
        iys, ims, ids, ihmsfs = sofa_time.jd_dtf(self.scale.upper()
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

        self.jd1, self.jd2 = sofa_time.dtf_jd(self.scale.upper().encode('utf8'),
                                              iy, im, id, ihr, imin, dsec)

    def str_kwargs(self):
        """
        Generator that yields a dict of values corresponding to the
        calendar date and time for the internal JD values.
        """
        iys, ims, ids, ihmsfs = sofa_time.jd_dtf(self.scale.upper()
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
    def vals(self):
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
        epoch_to_jd = getattr(sofa_time, self.epoch_to_jd)
        self.jd1, self.jd2 = epoch_to_jd(val1 + val2)

    @property
    def vals(self):
        jd_to_epoch = getattr(sofa_time, self.jd_to_epoch)
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
        epoch_to_jd = getattr(sofa_time, self.epoch_to_jd)
        self.jd1, self.jd2 = epoch_to_jd(years)

    @property
    def vals(self):
        jd_to_epoch = getattr(sofa_time, self.jd_to_epoch)
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
        self.jd1 = val1 / SECS_PER_DAY
        self.jd2 = val2 / SECS_PER_DAY

    @property
    def vals(self):
        return (self.jd1 + self.jd2) * SECS_PER_DAY


class TimeDeltaJD(TimeDeltaFormat):
    """Time delta in Julian days (86400 SI seconds)"""
    name = 'jd'

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1 = val1
        self.jd2 = val2

    @property
    def vals(self):
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


class OperandTypeError(TypeError):
    def __init__(self, left, right):
        self.value = ("unsupported operand type(s) for -: "
                      "'{0}' and '{1}'".format(left.__class__.__name__,
                                               right.__class__.__name__))
