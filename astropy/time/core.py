# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The astropy.time package provides functionality for manipulating times and
dates. Specific emphasis is placed on supporting time scales (e.g. UTC, TAI,
UT1) and time representations (e.g. JD, MJD, ISO 8601) that are used in
astronomy.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import fnmatch
import time
import re

from datetime import datetime

import numpy as np

from .. import units as u
from .. import _erfa as erfa
from ..units import UnitConversionError
from ..utils.compat.odict import OrderedDict
from ..utils.compat.misc import override__dir__
from ..utils.data_info import InfoDescriptor, DataInfo, data_info_factory
from ..utils.compat.numpy import broadcast_to
from ..extern import six


__all__ = ['Time', 'TimeDelta', 'TimeFormat', 'TimeJD', 'TimeMJD',
           'TimeFromEpoch', 'TimeUnix', 'TimeCxcSec', 'TimeGPS',
           'TimeDecimalYear',
           'TimePlotDate', 'TimeDatetime', 'TimeString',
           'TimeISO', 'TimeISOT', 'TimeFITS', 'TimeYearDayTime',
           'TimeEpochDate', 'TimeBesselianEpoch', 'TimeJulianEpoch',
           'TimeDeltaFormat', 'TimeDeltaSec', 'TimeDeltaJD',
           'ScaleValueError', 'OperandTypeError', 'TimeEpochDateString',
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

# These both get filled in at end after TimeFormat subclasses defined.
# Use an OrderedDict to fix the order in which formats are tried.
# This ensures, e.g., that 'isot' gets tried before 'fits'.
TIME_FORMATS = OrderedDict()
TIME_DELTA_FORMATS = OrderedDict()

TIME_SCALES = ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')
MULTI_HOPS = {('tai', 'tcb'): ('tt', 'tdb'),
              ('tai', 'tcg'): ('tt',),
              ('tai', 'ut1'): ('utc',),
              ('tai', 'tdb'): ('tt',),
              ('tcb', 'tcg'): ('tdb', 'tt'),
              ('tcb', 'tt'): ('tdb',),
              ('tcb', 'ut1'): ('tdb', 'tt', 'tai', 'utc'),
              ('tcb', 'utc'): ('tdb', 'tt', 'tai'),
              ('tcg', 'tdb'): ('tt',),
              ('tcg', 'ut1'): ('tt', 'tai', 'utc'),
              ('tcg', 'utc'): ('tt', 'tai'),
              ('tdb', 'ut1'): ('tt', 'tai', 'utc'),
              ('tdb', 'utc'): ('tt', 'tai'),
              ('tt', 'ut1'): ('tai', 'utc'),
              ('tt', 'utc'): ('tai',),
              }
GEOCENTRIC_SCALES = ('tai', 'tt', 'tcg')
BARYCENTRIC_SCALES = ('tcb', 'tdb')
ROTATIONAL_SCALES = ('ut1',)
TIME_DELTA_TYPES = dict((scale, scales)
                        for scales in (GEOCENTRIC_SCALES, BARYCENTRIC_SCALES,
                                       ROTATIONAL_SCALES) for scale in scales)
TIME_DELTA_SCALES = TIME_DELTA_TYPES.keys()
# For time scale changes, we need L_G and L_B, which are stored in erfam.h as
#   /* L_G = 1 - d(TT)/d(TCG) */
#   define ERFA_ELG (6.969290134e-10)
#   /* L_B = 1 - d(TDB)/d(TCB), and TDB (s) at TAI 1977/1/1.0 */
#   define ERFA_ELB (1.550519768e-8)
# These are exposed in erfa as erfa.ELG and erfa.ELB.
# Implied: d(TT)/d(TCG) = 1-L_G
# and      d(TCG)/d(TT) = 1/(1-L_G) = 1 + (1-(1-L_G))/(1-L_G) = 1 + L_G/(1-L_G)
# scale offsets as second = first + first * scale_offset[(first,second)]
SCALE_OFFSETS = {('tt', 'tai'): None,
                 ('tai', 'tt'): None,
                 ('tcg', 'tt'): -erfa.ELG,
                 ('tt', 'tcg'): erfa.ELG / (1. - erfa.ELG),
                 ('tcg', 'tai'): -erfa.ELG,
                 ('tai', 'tcg'): erfa.ELG / (1. - erfa.ELG),
                 ('tcb', 'tdb'): -erfa.ELB,
                 ('tdb', 'tcb'): erfa.ELB / (1. - erfa.ELB)}

# Translations between deprecated FITS timescales defined by
# Rots et al. 2015, A&A 574:A36, and timescales used here.
FITS_DEPRECATED_SCALES = {'TDT': 'tt', 'TDT': 'tt', 'ET': 'tt',
                          'GMT': 'utc', 'UT': 'utc', 'IAT': 'tai'}

# triple-level dictionary, yay!
SIDEREAL_TIME_MODELS = {
    'mean': {
        'IAU2006': {'function': erfa_time.gmst06, 'scales': ('ut1', 'tt')},
        'IAU2000': {'function': erfa_time.gmst00, 'scales': ('ut1', 'tt')},
        'IAU1982': {'function': erfa_time.gmst82, 'scales': ('ut1',)}},
    'apparent': {
        'IAU2006A': {'function': erfa_time.gst06a, 'scales': ('ut1', 'tt')},
        'IAU2000A': {'function': erfa_time.gst00a, 'scales': ('ut1', 'tt')},
        'IAU2000B': {'function': erfa_time.gst00b, 'scales': ('ut1',)},
        'IAU1994': {'function': erfa_time.gst94, 'scales': ('ut1',)}}}


class TimeInfo(DataInfo):
    """
    Container for meta information like name, description, format.  This is
    required when the object is used as a mixin column within a table, but can
    be used as a general way to store meta information.
    """
    attrs_from_parent = set(['unit'])  # unit is read-only and None

    @property
    def unit(self):
        return None

    info_summary_stats = staticmethod(
        data_info_factory(names=DataInfo._stats,
                          funcs=[getattr(np, stat) for stat in DataInfo._stats]))
    # When Time has mean, std, min, max methods:
    # funcs = [lambda x: getattr(x, stat)() for stat_name in DataInfo._stats])

class Time(object):
    """
    Represent and manipulate times and dates for astronomy.

    A `Time` object is initialized with one or more times in the ``val``
    argument.  The input times in ``val`` must conform to the specified
    ``format`` and must correspond to the specified time ``scale``.  The
    optional ``val2`` time input should be supplied only for numeric input
    formats (e.g. JD) where very high precision (better than 64-bit precision)
    is required.

    The allowed values for ``format`` can be listed with::

      >>> list(Time.FORMATS)
      ['jd', 'mjd', 'decimalyear', 'unix', 'cxcsec', 'gps', 'plot_date',
       'datetime', 'iso', 'isot', 'yday', 'fits', 'byear', 'jyear', 'byear_str',
       'jyear_str']

    Parameters
    ----------
    val : sequence, str, number, or `~astropy.time.Time` object
        Value(s) to initialize the time or times.
    val2 : sequence, str, or number; optional
        Value(s) to initialize the time or times.
    format : str, optional
        Format of input value(s)
    scale : str, optional
        Time scale of input value(s), must be one of the following:
        ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')
    precision : int, optional
        Digits of precision in string representation of time
    in_subfmt : str, optional
        Subformat for inputting string times
    out_subfmt : str, optional
        Subformat for outputting string times
    location : `~astropy.coordinates.EarthLocation` or tuple, optional
        If given as an tuple, it should be able to initialize an
        an EarthLocation instance, i.e., either contain 3 items with units of
        length for geocentric coordinates, or contain a longitude, latitude,
        and an optional height for geodetic coordinates.
        Can be a single location, or one for each input time.
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

    # Make sure that reverse arithmetic (e.g., TimeDelta.__rmul__)
    # gets called over the __mul__ of Numpy arrays.
    __array_priority__ = 20000

    # Declare that Time can be used as a Table column by defining the
    # attribute where column attributes will be stored.
    _astropy_column_attrs = None

    def __new__(cls, val, val2=None, format=None, scale=None,
                precision=None, in_subfmt=None, out_subfmt=None,
                location=None, copy=False):

        if isinstance(val, cls):
            self = val.replicate(format=format, copy=copy)
        else:
            self = super(Time, cls).__new__(cls)

        return self

    def __getnewargs__(self):
        return (self._time,)

    def __init__(self, val, val2=None, format=None, scale=None,
                 precision=None, in_subfmt=None, out_subfmt=None,
                 location=None, copy=False):

        if location is not None:
            from ..coordinates import EarthLocation
            if isinstance(location, EarthLocation):
                self.location = location
            else:
                self.location = EarthLocation(*location)
        else:
            self.location = None

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

        if self.location and (self.location.size > 1
                              and self.location.shape != self.shape):
            try:
                # check the location can be broadcast to self's shape.
                self.location = broadcast_to(self.location, self.shape,
                                             subok=True)
            except:
                raise ValueError('The location with shape {0} cannot be '
                                 'broadcast against time with shape {1}. '
                                 'Typically, either give a single location or '
                                 'one for each time.'
                                 .format(self.location.shape, self.shape))

    def _init_from_vals(self, val, val2, format, scale, copy):
        """
        Set the internal _format, scale, and _time attrs from user
        inputs.  This handles coercion into the correct shapes and
        some basic input validation.
        """

        # Coerce val into an array
        val = _make_array(val, copy)

        # If val2 is not None, ensure consistency
        if val2 is not None:
            val2 = _make_array(val2, copy)
            try:
                np.broadcast(val, val2)
            except ValueError:
                raise ValueError('Input val and val2 have inconsistent shape; '
                                 'they cannot be broadcast together.')

        if scale is not None:
            if not (isinstance(scale, six.string_types) and
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

        If format is `None` and the input is a string-type or object array then
        guess available formats and stop when one matches.
        """

        if format is None and val.dtype.kind in ('S', 'U', 'O'):
            formats = [(name, cls) for name, cls in self.FORMATS.items()
                       if issubclass(cls, TimeUnique)]
            err_msg = ('any of the formats where the format keyword is '
                       'optional {0}'.format([name for name, cls in formats]))
            # AstropyTime is a pseudo-format that isn't in the TIME_FORMATS registry,
            # but try to guess it at the end.
            formats.append(('astropy_time', TimeAstropyTime))

        elif not (isinstance(format, six.string_types) and
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
            except UnitConversionError:
                raise
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
            "Now" is determined using the `~datetime.datetime.utcnow`
            function, so its accuracy and precision is determined by that
            function.  Generally that means it is set by the accuracy of
            your system clock.

        Returns
        -------
        nowtime
            A new `Time` object (or a subclass of `Time` if this is called from
            such a subclass) at the current time.
        """
        # call `utcnow` immediately to be sure it's ASAP
        dtnow = datetime.utcnow()
        return cls(val=dtnow, format='datetime', scale='utc')

    info = InfoDescriptor(TimeInfo)

    @property
    def format(self):
        """
        Get or set time format.

        The format defines the way times are represented when accessed via the
        ``.value`` attribute.  By default it is the same as the format used for
        initializing the `Time` instance, but it can be set to any other value
        that could be used for initialization.  These can be listed with::

          >>> list(Time.FORMATS)
          ['jd', 'mjd', 'decimalyear', 'unix', 'cxcsec', 'gps', 'plot_date',
           'datetime', 'iso', 'isot', 'yday', 'fits', 'byear', 'jyear', 'byear_str',
           'jyear_str']
        """
        return self._format

    @format.setter
    def format(self, format):
        """Set time format"""
        if format not in self.FORMATS:
            raise ValueError('format must be one of {0}'
                             .format(list(self.FORMATS)))
        format_cls = self.FORMATS[format]

        # If current output subformat is not in the new format then replace
        # with default '*'
        if hasattr(format_cls, 'subfmts'):
            subfmt_names = [subfmt[0] for subfmt in format_cls.subfmts]
            if self.out_subfmt not in subfmt_names:
                self.out_subfmt = '*'

        self._time = format_cls(self._time.jd1, self._time.jd2,
                                self._time._scale, self.precision,
                                in_subfmt=self.in_subfmt,
                                out_subfmt=self.out_subfmt,
                                from_jd=True)
        self._format = format

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
        for sys1, sys2 in six.moves.zip(xforms[:-1], xforms[1:]):
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
        times.
        """
        return self._in_subfmt

    @in_subfmt.setter
    def in_subfmt(self, val):
        if not isinstance(val, six.string_types):
            raise ValueError('in_subfmt attribute must be a string')
        self._in_subfmt = val

    @property
    def out_subfmt(self):
        """
        Unix wildcard pattern to select subformats for outputting times.
        """
        return self._out_subfmt

    @out_subfmt.setter
    def out_subfmt(self, val):
        if not isinstance(val, six.string_types):
            raise ValueError('out_subfmt attribute must be a string')
        self._out_subfmt = val

    @property
    def ndim(self):
        return self._time.jd1.ndim

    @property
    def shape(self):
        """The shape of the time instances.

        Like `~numpy.ndarray.shape`, can be set to a new shape by assigning a
        tuple.

        Raises
        ------
        AttributeError: if the shape of the ``jd1``, ``jd2``, ``location``,
        ``delta_ut1_utc``, or ``delta_tdb_tt`` attributes cannot be changed
        without the arrays being copied.  For these cases, use the
        `Time.reshape` method.
        """
        return self._time.jd1.shape

    @shape.setter
    def shape(self, shape):
        self._time.jd1.shape = shape
        self._time.jd2.shape = shape
        for attr in ('_delta_ut1_utc', '_delta_tdb_tt', 'location'):
            val = getattr(self, attr, None)
            if val is not None and val.size > 1:
                val.shape = shape

    @property
    def size(self):
        return self._time.jd1.size

    def __bool__(self):
        """Any time should evaluate to True, except when it is empty."""
        return self.size > 0

    # In python2, __bool__ is not defined.
    __nonzero__ = __bool__

    @property
    def isscalar(self):
        return self.shape == ()

    def _shaped_like_input(self, value):
        return value if self._time.jd1.shape else value.item()

    @property
    def jd1(self):
        """
        First of the two doubles that internally store time value(s) in JD.
        """
        return self._shaped_like_input(self._time.jd1)

    @property
    def jd2(self):
        """
        Second of the two doubles that internally store time value(s) in JD.
        """
        return self._shaped_like_input(self._time.jd2)

    @property
    def value(self):
        """Time value(s) in current format"""
        return self._shaped_like_input(self._time.value)

    def sidereal_time(self, kind, longitude=None, model=None):
        """Calculate sidereal time.

        Parameters
        ---------------
        kind : str
            ``'mean'`` or ``'apparent'``, i.e., accounting for precession
            only, or also for nutation.
        longitude : `~astropy.units.Quantity`, `str`, or `None`; optional
            The longitude on the Earth at which to compute the sidereal time.
            Can be given as a `~astropy.units.Quantity` with angular units
            (or an `~astropy.coordinates.Angle` or
            `~astropy.coordinates.Longitude`), or as a name of an
            observatory (currently, only ``'greenwich'`` is supported,
            equivalent to 0 deg).  If `None` (default), the ``lon`` attribute of
            the Time object is used.
        model : str or `None`; optional
            Precession (and nutation) model to use.  The available ones are:
            - {0}: {1}
            - {2}: {3}
            If `None` (default), the last (most recent) one from the appropriate
            list above is used.

        Returns
        -------
        sidereal time : `~astropy.coordinates.Longitude`
            Sidereal time as a quantity with units of hourangle
        """  # docstring is formatted below

        from ..coordinates import Longitude

        if kind.lower() not in SIDEREAL_TIME_MODELS.keys():
            raise ValueError('The kind of sidereal time has to be {0}'.format(
                ' or '.join(sorted(SIDEREAL_TIME_MODELS.keys()))))

        available_models = SIDEREAL_TIME_MODELS[kind.lower()]

        if model is None:
            model = sorted(available_models.keys())[-1]
        else:
            if model.upper() not in available_models:
                raise ValueError(
                    'Model {0} not implemented for {1} sidereal time; '
                    'available models are {2}'
                    .format(model, kind, sorted(available_models.keys())))

        if longitude is None:
            if self.location is None:
                raise ValueError('No longitude is given but the location for '
                                 'the Time object is not set.')
            longitude = self.location.longitude
        elif longitude == 'greenwich':
            longitude = Longitude(0., u.degree,
                                  wrap_angle=180.*u.degree)
        else:
            # sanity check on input
            longitude = Longitude(longitude, u.degree,
                                  wrap_angle=180.*u.degree)

        gst = self._erfa_sidereal_time(available_models[model.upper()])
        return Longitude(gst + longitude, u.hourangle)
    sidereal_time.__doc__ = sidereal_time.__doc__.format(
        'apparent', sorted(SIDEREAL_TIME_MODELS['apparent'].keys()),
        'mean', sorted(SIDEREAL_TIME_MODELS['mean'].keys()))

    def _erfa_sidereal_time(self, model):
        """Calculate a sidereal time using a IAU precession/nutation model."""

        from ..coordinates import Longitude

        erfa_function = model['function']
        erfa_parameters = [getattr(getattr(self, scale)._time, jd_part)
                           for scale in model['scales']
                           for jd_part in ('jd1', 'jd2')]

        sidereal_time = erfa_function(*erfa_parameters)

        return Longitude(sidereal_time, u.radian).to(u.hourangle)

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
        tm : Time object
            Copy of this object
        """
        return self._replicate(format=format, method='copy')

    def replicate(self, format=None, copy=False):
        """
        Return a replica of the Time object, optionally changing the format.

        If ``format`` is supplied then the time format of the returned Time
        object will be set accordingly, otherwise it will be unchanged from the
        original.

        If ``copy`` is set to `True` then a full copy of the internal time arrays
        will be made.  By default the replica will use a reference to the
        original arrays when possible to save memory.  The internal time arrays
        are normally not changeable by the user so in most cases it should not
        be necessary to set ``copy`` to `True`.

        The convenience method copy() is available in which ``copy`` is `True`
        by default.

        Parameters
        ----------
        format : str, optional
            Time format of the replica.
        copy : bool, optional
            Return a true copy instead of using references where possible.

        Returns
        -------
        tm : Time object
            Replica of this object
        """
        # To avoid recalculating integer day + fraction, no longer just
        # instantiate a new class instance, but rather do the steps by hand.
        # This also avoids quite a bit of unnecessary work in __init__
        ###  tm = self.__class__(self._time.jd1, self._time.jd2,
        ###                      format='jd', scale=self.scale, copy=copy)
        return self._replicate(format=format, method='copy' if copy else None)


    def _replicate(self, method=None, *args, **kwargs):
        """Replicate a time object, possibly applying a method to the arrays.

        Parameters
        ----------
        method : str, optional
            If given, the method is applied to the internal ``jd1`` and ``jd2``
            arrays, as well as to possible ``location``, ``_delta_ut1_utc``,
            and ``_delta_tdb_tt`` arrays, broadcasting the latter as required.
            Example methods: ``copy``, ``__getitem__``, ``reshape``.
        args : tuple
            Any positional arguments for ``method``.
        kwargs : dict
            Any keyword arguments for ``method``.  If the ``format`` keyword
            argument is present, this will be used as the Time format of the
            replica.
        """
        new_format = kwargs.pop('format', None)
        if new_format is None:
            new_format = self.format

        jd1, jd2 = self._time.jd1, self._time.jd2
        if method is not None:
            jd1 = getattr(jd1, method)(*args, **kwargs)
            jd2 = getattr(jd2, method)(*args, **kwargs)

        tm = super(Time, self.__class__).__new__(self.__class__)
        tm._time = TimeJD(jd1, jd2, self.scale, self.precision,
                          self.in_subfmt, self.out_subfmt, from_jd=True)
        # Optional ndarray attributes.
        for attr in ('_delta_ut1_utc', '_delta_tdb_tt', 'location'):
            try:
                val = getattr(self, attr)
            except AttributeError:
                continue

            # Appy the method to any value arrays (though skip if there is only
            # a single element and the method would return a view, since in
            # that case nothing would change).
            if method is not None and val is not None:
                if method == 'copy' or method == 'flatten' and val.size == 1:
                    val = val.copy()

                elif val.size > 1:
                    val = getattr(val, method)(*args, **kwargs)

            setattr(tm, attr, val)

        for attr in ('precision', 'in_subfmt', 'out_subfmt'):
            val = getattr(self, attr)
            if method in ('copy', 'flatten') and hasattr(val, 'copy'):
                val = val.copy()

            setattr(tm, attr, val)

        # Make the new internal _time object corresponding to the format
        # in the copy.  If the format is unchanged this process is lightweight
        # and does not create any new arrays.
        if new_format not in tm.FORMATS:
            raise ValueError('format must be one of {0}'
                             .format(list(tm.FORMATS)))

        NewFormat = tm.FORMATS[new_format]
        tm._time = NewFormat(tm._time.jd1, tm._time.jd2,
                             tm._time._scale, tm.precision,
                             tm.in_subfmt, tm.out_subfmt,
                             from_jd=True)
        tm._format = new_format

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
        tm._time.jd1 = self._time.jd1[item]
        tm._time.jd2 = self._time.jd2[item]
        attrs = ('_delta_ut1_utc', '_delta_tdb_tt', 'location')
        for attr in attrs:
            val = getattr(self, attr, None)
            if val is not None:
                try:
                    setattr(tm, attr, val[item])
                except IndexError:  # location may be scalar (same for all)
                    continue

        return tm

    def reshape(self, *args, **kwargs):
        """Returns a time instance containing the same data with a new shape.

        Parameters are as for :meth:`~numpy.ndarray.reshape`.  Note that it is
        not always possible to change the shape of an array without copying the
        data. If you want an error to be raise if the data is copied, you
        should assign the new shape to the shape attribute.
        """
        return self._replicate('reshape', *args, **kwargs)

    def ravel(self, *args, **kwargs):
        """Return an instance with the time array collapsed into one dimension.

        Parameters are as for :meth:`~numpy.ndarray.ravel`. Note that it is
        not always possible to unravel an array without copying the data.
        If you want an error to be raise if the data is copied, you should
        should assign shape ``(-1,)`` to the shape attribute.
        """
        return self._replicate('ravel', *args, **kwargs)

    def flatten(self, *args, **kwargs):
        """Return a copy with the time array collapsed into one dimension.

        Parameters are as for :meth:`~numpy.ndarray.flatten`.
        """
        return self._replicate('flatten', *args, **kwargs)

    def transpose(self, *args, **kwargs):
        """Return a time instance with the data transposed.

        Parameters are as for :meth:`~numpy.ndarray.transpose`.  All internal
        data are views of the data of the original.
        """
        return self._replicate('transpose', *args, **kwargs)

    @property
    def T(self):
        """Return a time instance with the data transposed.

        Parameters are as for :attr:`~numpy.ndarray.T`.  All internal
        data are views of the data of the original.
        """
        if self.ndim < 2:
            return self
        else:
            return self.transpose()

    def swapaxes(self, *args, **kwargs):
        """Return a time instance with the given axes interchanged.

        Parameters are as for :meth:`~numpy.ndarray.swapaxes`.  All internal
        data are views of the data of the original.
        """
        return self._replicate('swapaxes', *args, **kwargs)

    def diagonal(self, *args, **kwargs):
        """Return a time instance with the specified diagonals.

        Parameters are as for :meth:`~numpy.ndarray.diagonal`.  All internal
        data are views of the data of the original.
        """
        return self._replicate('diagonal', *args, **kwargs)

    def squeeze(self, *args, **kwargs):
        """Return a time instance with single-dimensional shape entries removed

        Parameters are as for :meth:`~numpy.ndarray.squeeze`.  All internal
        data are views of the data of the original.
        """
        return self._replicate('squeeze', *args, **kwargs)

    def take(self, indices, axis=None, mode='raise'):
        """Return a Time object formed from the elements the given indices.

        Parameters are as for :meth:`~numpy.ndarray.take`, except that,
        obviously, no output array can be given.
        """
        return self._replicate('take', indices, axis=axis, mode=mode)

    def __getattr__(self, attr):
        """
        Get dynamic attributes to output format or do timescale conversion.
        """
        if attr in self.SCALES and self.scale is not None:
            tm = self.replicate()
            tm._set_scale(attr)
            return tm

        elif attr in self.FORMATS:
            tm = self.replicate(format=attr)
            if hasattr(self.FORMATS[attr], 'epoch_scale'):
                tm._set_scale(self.FORMATS[attr].epoch_scale)
            return tm.value

        elif attr in TIME_SCALES:  # allowed ones done above (self.SCALES)
            if self.scale is None:
                raise ScaleValueError("Cannot convert TimeDelta with "
                                      "undefined scale to any defined scale.")
            else:
                raise ScaleValueError("Cannot convert {0} with scale "
                                      "'{1}' to scale '{2}'"
                                      .format(self.__class__.__name__,
                                              self.scale, attr))

        else:
            # Should raise AttributeError
            return self.__getattribute__(attr)

    @override__dir__
    def __dir__(self):
        result = set(self.SCALES)
        result.update(self.FORMATS)
        return result

    def _match_shape(self, val):
        """
        Ensure that `val` is matched to length of self.  If val has length 1
        then broadcast, otherwise cast to double and make sure shape matches.
        """
        val = _make_array(val, copy=True)  # be conservative and copy
        if val.size > 1 and val.shape != self.shape:
            try:
                # check the value can be broadcast to the shape of self.
                val = broadcast_to(val, self.shape, subok=True)
            except:
                raise ValueError('Attribute shape must match or be '
                                 'broadcastable to that of Time object. '
                                 'Typically, give either a single value or '
                                 'one for each time.')

        return val

    def get_delta_ut1_utc(self, iers_table=None, return_status=False):
        """Find UT1 - UTC differences by interpolating in IERS Table.

        Parameters
        ----------
        iers_table : ``astropy.utils.iers.IERS`` table, optional
            Table containing UT1-UTC differences from IERS Bulletins A
            and/or B.  If `None`, use default version (see
            ``astropy.utils.iers``)
        return_status : bool
            Whether to return status values.  If `False` (default), iers
            raises `~.exceptions.IndexError` if any time is out of the range
            covered by the IERS table.

        Returns
        -------
        ut1_utc : float or float array
            UT1-UTC, interpolated in IERS Table
        status : int or int array
            Status values (if ``return_status=`True```)::
            ``astropy.utils.iers.FROM_IERS_B``
            ``astropy.utils.iers.FROM_IERS_A``
            ``astropy.utils.iers.FROM_IERS_A_PREDICTION``
            ``astropy.utils.iers.TIME_BEFORE_IERS_RANGE``
            ``astropy.utils.iers.TIME_BEYOND_IERS_RANGE``

        Notes
        -----
        In normal usage, UT1-UTC differences are calculated automatically
        on the first instance ut1 is needed.

        Examples
        --------
        To check in code whether any times are before the IERS table range::

            >>> from astropy.utils.iers import TIME_BEFORE_IERS_RANGE
            >>> t = Time(['1961-01-01', '2000-01-01'], scale='utc')
            >>> delta, status = t.get_delta_ut1_utc(return_status=True)
            >>> status == TIME_BEFORE_IERS_RANGE
            array([ True, False], dtype=bool)

        To use an updated IERS A bulletin to calculate UT1-UTC
        (see also ``astropy.utils.iers``)::

            >>> from astropy.utils.iers import IERS_A, IERS_A_URL
            >>> from astropy.utils.data import download_file
            >>> t = Time(['1974-01-01', '2000-01-01'], scale='utc')
            >>> iers_a_file = download_file(IERS_A_URL,
            ...                             cache=True)        # doctest: +REMOTE_DATA
            Downloading ... [Done]
            >>> iers_a = IERS_A.open(iers_a_file)              # doctest: +REMOTE_DATA
            >>> t.delta_ut1_utc = t.get_delta_ut1_utc(iers_a)  # doctest: +REMOTE_DATA

        The delta_ut1_utc property will be used to calculate t.ut1;
        raises IndexError if any of the times is out of range.

        """
        if iers_table is None:
            from ..utils.iers import IERS
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
            from ..utils.iers import IERS
            iers_table = IERS.open()
            # jd1, jd2 are normally set (see above), except if delta_ut1_utc
            # is access directly; ensure we behave as expected for that case
            if jd1 is None:
                self_utc = self.utc
                jd1, jd2 = self_utc.jd1, self_utc.jd2
                scale = 'utc'
            else:
                scale = self.scale
            # interpolate UT1-UTC in IERS table
            delta = iers_table.ut1_utc(jd1, jd2)
            # if we interpolated using UT1 jds, we may be off by one
            # second near leap seconds (and very slightly off elsewhere)
            if scale == 'ut1':
                # calculate UTC using the offset we got; the ERFA routine
                # is tolerant of leap seconds, so will do this right
                jd1_utc, jd2_utc = erfa_time.ut1_utc(jd1, jd2, delta)
                # calculate a better estimate using the nearly correct UTC
                delta = iers_table.ut1_utc(jd1_utc, jd2_utc)

            self._set_delta_ut1_utc(delta)

        return self._delta_ut1_utc

    def _set_delta_ut1_utc(self, val):
        val = self._match_shape(val)
        if hasattr(val, 'to'):  # Matches Quantity but also TimeDelta.
            val = val.to(u.second).value
        self._delta_ut1_utc = val

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
            # TDB or TT) to an approximate UT1.  Since TT and TDB are
            # pretty close (few msec?), assume TT.  Similarly, since the
            # UT1 terms are very small, use UTC instead of UT1.
            njd1, njd2 = erfa_time.tt_tai(jd1, jd2)
            njd1, njd2 = erfa_time.tai_utc(njd1, njd2)
            # subtract 0.5, so UT is fraction of the day from midnight
            ut = day_frac(njd1 - 0.5, njd2)[1]

            if self.location is None:
                from ..coordinates import EarthLocation
                location = EarthLocation.from_geodetic(0., 0., 0.)
            else:
                location = self.location
            # Geodetic params needed for d_tdb_tt()
            lon = location.longitude
            rxy = np.hypot(location.x, location.y)
            z = location.z
            self._delta_tdb_tt = erfa_time.d_tdb_tt(
                jd1, jd2, ut, lon.to(u.radian).value,
                rxy.to(u.km).value, z.to(u.km).value)

        return self._delta_tdb_tt

    def _set_delta_tdb_tt(self, val):
        val = self._match_shape(val)
        if hasattr(val, 'to'):  # Matches Quantity but also TimeDelta.
            val = val.to(u.second).value
        self._delta_tdb_tt = val

    # Note can't use @property because _get_delta_tdb_tt is explicitly
    # called with the optional jd1 and jd2 args.
    delta_tdb_tt = property(_get_delta_tdb_tt, _set_delta_tdb_tt)
    """TDB - TT time scale offset"""

    def __len__(self):
        if self.isscalar:
            raise TypeError("Scalar {0} object has no len()"
                            .format(self.__class__.__name__))
        return len(self.jd1)

    def __sub__(self, other):
        if not isinstance(other, Time):
            try:
                other = TimeDelta(other)
            except:
                raise OperandTypeError(self, other, '-')

        # Tdelta - something is dealt with in TimeDelta, so we have
        # T      - Tdelta = T
        # T      - T      = Tdelta
        other_is_delta = isinstance(other, TimeDelta)

        # we need a constant scale to calculate, which is guaranteed for
        # TimeDelta, but not for Time (which can be UTC)
        if other_is_delta:  # T - Tdelta
            if self.scale in other.SCALES:
                out = self.replicate()
                if other.scale not in (out.scale, None):
                    other = getattr(other, out.scale)
            else:
                out = getattr(self, (other.scale if other.scale is not None
                                     else 'tai'))
            # remove attributes that are invalidated by changing time
            for attr in ('_delta_ut1_utc', '_delta_tdb_tt'):
                if hasattr(out, attr):
                    delattr(out, attr)

        else:  # T - T
            self_time = (self._time if self.scale in TIME_DELTA_SCALES
                         else self.tai._time)
            # set up TimeDelta, subtraction to be done shortly
            out = TimeDelta(self_time.jd1, self_time.jd2, format='jd',
                            scale=self_time.scale)

            if other.scale != out.scale:
                other = getattr(other, out.scale)

        jd1 = out._time.jd1 - other._time.jd1
        jd2 = out._time.jd2 - other._time.jd2

        out._time.jd1, out._time.jd2 = day_frac(jd1, jd2)

        if other_is_delta and out.scale != self.scale:
            return getattr(out, self.scale)

        return out

    def __add__(self, other):
        if not isinstance(other, Time):
            try:
                other = TimeDelta(other)
            except:
                raise OperandTypeError(self, other, '+')

        # Tdelta + something is dealt with in TimeDelta, so we have
        # T      + Tdelta = T
        # T      + T      = error

        if not isinstance(other, TimeDelta):
            raise OperandTypeError(self, other, '+')

        # ideally, we calculate in the scale of the Time item, since that is
        # what we want the output in, but this may not be possible, since
        # TimeDelta cannot be converted arbitrarily
        if self.scale in other.SCALES:
            out = self.replicate()
            if other.scale not in (out.scale, None):
                other = getattr(other, out.scale)
        else:
            out = getattr(self, (other.scale if other.scale is not None
                                 else 'tai'))

        # remove attributes that are invalidated by changing time
        for attr in ('_delta_ut1_utc', '_delta_tdb_tt'):
            if hasattr(out, attr):
                delattr(out, attr)

        jd1 = out._time.jd1 + other._time.jd1
        jd2 = out._time.jd2 + other._time.jd2

        out._time.jd1, out._time.jd2 = day_frac(jd1, jd2)

        if out.scale != self.scale:
            return getattr(out, self.scale)

        return out

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        out = self.__sub__(other)
        return -out

    def _time_difference(self, other, op=None):
        """If other is of same class as self, return difference in self.scale.
        Otherwise, raise OperandTypeError.
        """
        if other.__class__ is not self.__class__:
            try:
                other = self.__class__(other, scale=self.scale)
            except:
                raise OperandTypeError(self, other, op)

        if(self.scale is not None and self.scale not in other.SCALES or
           other.scale is not None and other.scale not in self.SCALES):
            raise TypeError("Cannot compare TimeDelta instances with scales "
                            "'{0}' and '{1}'".format(self.scale, other.scale))

        if self.scale is not None and other.scale is not None:
            other = getattr(other, self.scale)

        return (self.jd1 - other.jd1) + (self.jd2 - other.jd2)

    def __lt__(self, other):
        return self._time_difference(other, '<') < 0.

    def __le__(self, other):
        return self._time_difference(other, '<=') <= 0.

    def __eq__(self, other):
        """
        If other is an incompatible object for comparison, return `False`.
        Otherwise, return `True` if the time difference between self and
        other is zero.
        """
        try:
            diff = self._time_difference(other)
        except OperandTypeError:
            return False
        return diff == 0.

    def __ne__(self, other):
        """
        If other is an incompatible object for comparison, return `True`.
        Otherwise, return `False` if the time difference between self and
        other is zero.
        """
        try:
            diff = self._time_difference(other)
        except OperandTypeError:
            return True
        return diff != 0.

    def __gt__(self, other):
        return self._time_difference(other, '>') > 0.

    def __ge__(self, other):
        return self._time_difference(other, '>=') >= 0.


class TimeDelta(Time):
    """
    Represent the time difference between two times.

    A TimeDelta object is initialized with one or more times in the ``val``
    argument.  The input times in ``val`` must conform to the specified
    ``format``.  The optional ``val2`` time input should be supplied only for
    numeric input formats (e.g. JD) where very high precision (better than
    64-bit precision) is required.

    The allowed values for ``format`` can be listed with::

      >>> list(TimeDelta.FORMATS)
      ['sec', 'jd']

    Note that for time differences, the scale can be among three groups:
    geocentric ('tai', 'tt', 'tcg'), barycentric ('tcb', 'tdb'), and rotational
    ('ut1'). Within each of these, the scales for time differences are the
    same. Conversion between geocentric and barycentric is possible, as there
    is only a scale factor change, but one cannot convert to or from 'ut1', as
    this requires knowledge of the actual times, not just their difference. For
    a similar reason, 'utc' is not a valid scale for a time difference: a UTC
    day is not always 86400 seconds.

    Parameters
    ----------
    val : numpy ndarray, list, str, number, or `~astropy.time.TimeDelta` object
        Data to initialize table.
    val2 : numpy ndarray, list, str, or number; optional
        Data to initialize table.
    format : str, optional
        Format of input value(s)
    scale : str, optional
        Time scale of input value(s), must be one of the following values:
        ('tdb', 'tt', 'ut1', 'tcg', 'tcb', 'tai'). If not given (or
        ``None``), the scale is arbitrary; when added or subtracted from a
        ``Time`` instance, it will be used without conversion.
    copy : bool, optional
        Make a copy of the input values
    """
    SCALES = TIME_DELTA_SCALES
    """List of time delta scales."""

    FORMATS = TIME_DELTA_FORMATS
    """Dict of time delta formats."""

    def __init__(self, val, val2=None, format=None, scale=None, copy=False):
        if isinstance(val, self.__class__):
            if scale is not None:
                self._set_scale(scale)
        else:
            if format is None:
                try:
                    val = val.to(u.day)
                    if val2 is not None:
                        val2 = val2.to(u.day)
                except:
                    raise ValueError('Only Quantities with Time units can '
                                     'be used to initiate {0} instances .'
                                     .format(self.__class__.__name__))
                format = 'jd'

            self._init_from_vals(val, val2, format, scale, copy)

            if scale is not None:
                self.SCALES = TIME_DELTA_TYPES[scale]

    def replicate(self, *args, **kwargs):
        out = super(TimeDelta, self).replicate(*args, **kwargs)
        out.SCALES = self.SCALES
        return out

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

        # For TimeDelta, there can only be a change in scale factor,
        # which is written as time2 - time1 = scale_offset * time1
        scale_offset = SCALE_OFFSETS[(self.scale, scale)]
        if scale_offset is None:
            self._time.scale = scale
        else:
            jd1, jd2 = self._time.jd1, self._time.jd2
            offset1, offset2 = day_frac(jd1, jd2, factor=scale_offset)
            self._time = self.FORMATS[self.format](
                jd1 + offset1, jd2 + offset2, scale,
                self.precision, self.in_subfmt,
                self.out_subfmt, from_jd=True)

    def __add__(self, other):
        # only deal with TimeDelta + TimeDelta
        if isinstance(other, Time):
            if not isinstance(other, TimeDelta):
                return other.__add__(self)
        else:
            try:
                other = TimeDelta(other)
            except:
                raise OperandTypeError(self, other, '+')

        # the scales should be compatible (e.g., cannot convert TDB to TAI)
        if(self.scale is not None and self.scale not in other.SCALES or
           other.scale is not None and other.scale not in self.SCALES):
            raise TypeError("Cannot add TimeDelta instances with scales "
                            "'{0}' and '{1}'".format(self.scale, other.scale))

        # adjust the scale of other if the scale of self is set (or no scales)
        if self.scale is not None or other.scale is None:
            out = self.replicate()
            if other.scale is not None:
                other = getattr(other, self.scale)
        else:
            out = other.replicate()

        jd1 = self._time.jd1 + other._time.jd1
        jd2 = self._time.jd2 + other._time.jd2

        out._time.jd1, out._time.jd2 = day_frac(jd1, jd2)

        return out

    def __sub__(self, other):
        # only deal with TimeDelta - TimeDelta
        if isinstance(other, Time):
            if not isinstance(other, TimeDelta):
                raise OperandTypeError(self, other, '-')
        else:
            try:
                other = TimeDelta(other)
            except:
                raise OperandTypeError(self, other, '-')

        # the scales should be compatible (e.g., cannot convert TDB to TAI)
        if(self.scale is not None and self.scale not in other.SCALES or
           other.scale is not None and other.scale not in self.SCALES):
            raise TypeError("Cannot subtract TimeDelta instances with scales "
                            "'{0}' and '{1}'".format(self.scale, other.scale))

        # adjust the scale of other if the scale of self is set (or no scales)
        if self.scale is not None or other.scale is None:
            out = self.replicate()
            if other.scale is not None:
                other = getattr(other, self.scale)
        else:
            out = other.replicate()

        jd1 = self._time.jd1 - other._time.jd1
        jd2 = self._time.jd2 - other._time.jd2

        out._time.jd1, out._time.jd2 = day_frac(jd1, jd2)

        return out

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
            raise OperandTypeError(self, other, '*')

        try:   # convert to straight float if dimensionless quantity
            other = other.to(1)
        except:
            pass

        try:
            jd1, jd2 = day_frac(self.jd1, self.jd2, factor=other)
            out = TimeDelta(jd1, jd2, format='jd', scale=self.scale)
        except Exception as err:  # try downgrading self to a quantity
            try:
                return self.to(u.day) * other
            except:
                raise err

        if self.format != 'jd':
            out = out.replicate(format=self.format)
        return out

    def __rmul__(self, other):
        """Multiplication of numbers/arrays with `TimeDelta` objects."""
        return self.__mul__(other)

    def __div__(self, other):
        """Division of `TimeDelta` objects by numbers/arrays."""
        return self.__truediv__(other)

    def __rdiv__(self, other):
        """Division by `TimeDelta` objects of numbers/arrays."""
        return self.__rtruediv__(other)

    def __truediv__(self, other):
        """Division of `TimeDelta` objects by numbers/arrays."""
        # cannot do __mul__(1./other) as that looses precision
        try:
            other = other.to(1)
        except:
            pass

        try:   # convert to straight float if dimensionless quantity
            jd1, jd2 = day_frac(self.jd1, self.jd2, divisor=other)
            out = TimeDelta(jd1, jd2, format='jd', scale=self.scale)
        except Exception as err:  # try downgrading self to a quantity
            try:
                return self.to(u.day) / other
            except:
                raise err

        if self.format != 'jd':
            out = out.replicate(format=self.format)
        return out

    def __rtruediv__(self, other):
        """Division by `TimeDelta` objects of numbers/arrays."""
        return other / self.to(u.day)

    def to(self, *args, **kwargs):
        return u.Quantity(self._time.jd1 + self._time.jd2,
                          u.day).to(*args, **kwargs)


def _regexify_subfmts(subfmts):
    """
    Iterate through each of the sub-formats and try substituting simple
    regular expressions for the strptime codes for year, month, day-of-month,
    hour, minute, second.  If no % characters remain then turn the final string
    into a compiled regex.  This assumes time formats do not have a % in them.

    This is done both to speed up parsing of strings and to allow mixed formats
    where strptime does not quite work well enough.
    """
    new_subfmts = []
    for subfmt_tuple in subfmts:
        subfmt_in = subfmt_tuple[1]
        for strptime_code, regex in (('%Y', r'(?P<year>\d\d\d\d)'),
                                     ('%m', r'(?P<mon>\d{1,2})'),
                                     ('%d', r'(?P<mday>\d{1,2})'),
                                     ('%H', r'(?P<hour>\d{1,2})'),
                                     ('%M', r'(?P<min>\d{1,2})'),
                                     ('%S', r'(?P<sec>\d{1,2})')):
            subfmt_in = re.sub(strptime_code, regex, subfmt_in)

        if '%' not in subfmt_in:
            subfmt_tuple = (subfmt_tuple[0],
                            re.compile(subfmt_in + '$'),
                            subfmt_tuple[2])
        new_subfmts.append(subfmt_tuple)

    return tuple(new_subfmts)


class TimeFormatMeta(type):
    """
    Metaclass that adds `TimeFormat` and `TimeDeltaFormat` to the
    `TIME_FORMATS` and `TIME_DELTA_FORMATS` registries, respectively.
    """

    _registry = TIME_FORMATS

    def __new__(mcls, name, bases, members):
        cls = super(TimeFormatMeta, mcls).__new__(mcls, name, bases, members)

        # Register time formats that have a name, but leave out astropy_time since
        # it is not a user-accessible format and is only used for initialization into
        # a different format.
        if 'name' in members and cls.name != 'astropy_time':
            mcls._registry[cls.name] = cls

        if 'subfmts' in members:
            cls.subfmts = _regexify_subfmts(members['subfmts'])

        return cls


@six.add_metaclass(TimeFormatMeta)
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

        if from_jd:
            self.jd1 = val1
            self.jd2 = val2
        else:
            val1, val2 = self._check_val_type(val1, val2)
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
        if not (val1.dtype == np.double and np.all(np.isfinite(val1)) and
                (val2 is None or
                 val2.dtype == np.double and np.all(np.isfinite(val2)))):
            raise TypeError('Input values for {0} class must be finite doubles'
                            .format(self.name))

        if getattr(val1, 'unit', None) is not None:
            # set possibly scaled unit any quantities should be converted to
            _unit = u.CompositeUnit(getattr(self, 'unit', 1.), [u.day], [1])
            val1 = val1.to(_unit).value
            if val2 is not None:
                val2 = val2.to(_unit).value
        elif getattr(val2, 'unit', None) is not None:
            raise TypeError('Cannot mix float and Quantity inputs')

        if val2 is None:
            val2 = np.zeros_like(val1)

        return val1, val2

    def _check_scale(self, scale):
        """
        Return a validated scale value.

        If there is a class attribute 'scale' then that defines the default /
        required time scale for this format.  In this case if a scale value was
        provided that needs to match the class default, otherwise return
        the class default.

        Otherwise just make sure that scale is in the allowed list of
        scales.  Provide a different error message if `None` (no value) was
        supplied.
        """
        if hasattr(self.__class__, 'epoch_scale') and scale is None:
            scale = self.__class__.epoch_scale

        if scale is None:
            scale = 'utc'  # Default scale as of astropy 0.4

        if scale not in TIME_SCALES:
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
    """
    Julian Date time format.
    This represents the number of days since the beginning of
    the Julian Period.
    For example, 2451544.5 in JD is midnight on January 1, 2000.
    """
    name = 'jd'

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2)

    @property
    def value(self):
        return self.jd1 + self.jd2


class TimeMJD(TimeFormat):
    """
    Modified Julian Date time format.
    This represents the number of days since midnight on November 17, 1858.
    For example, 51544.0 in MJD is midnight on January 1, 2000.
    """
    name = 'mjd'

    def set_jds(self, val1, val2):
        # TODO - this routine and vals should be Cythonized to follow the ERFA
        # convention of preserving precision by adding to the larger of the two
        # values in a vectorized operation.  But in most practical cases the
        # first one is probably biggest.
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2)
        self.jd1 += erfa.DJM0  # erfa.DJM0=2400000.5 (from erfam.h)

    @property
    def value(self):
        return (self.jd1 - erfa.DJM0) + self.jd2


class TimeDecimalYear(TimeFormat):
    """
    Time as a decimal year, with integer values corresponding to midnight
    of the first day of each year.  For example 2000.5 corresponds to the
    ISO time '2000-07-02 00:00:00'.
    """
    name = 'decimalyear'

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.

        sum12, err12 = two_sum(val1, val2)
        iy_start = np.trunc(sum12).astype(np.int)
        extra, y_frac = two_sum(sum12, -iy_start)
        y_frac += extra + err12

        val = (val1 + val2).astype(np.double)
        iy_start = np.trunc(val).astype(np.int)

        imon = np.ones_like(iy_start)
        iday = np.ones_like(iy_start)
        ihr = np.zeros_like(iy_start)
        imin = np.zeros_like(iy_start)
        isec = np.zeros_like(y_frac)

        # Possible enhancement: use np.unique to only compute start, stop
        # for unique values of iy_start.
        jd1_start, jd2_start = erfa_time.dtf_jd(
            self.scale.upper().encode('utf8'), iy_start, imon, iday, ihr, imin, isec)
        jd1_end, jd2_end = erfa_time.dtf_jd(
            self.scale.upper().encode('utf8'), iy_start + 1, imon, iday, ihr, imin, isec)

        t_start = Time(jd1_start, jd2_start, scale=self.scale, format='jd')
        t_end = Time(jd1_end, jd2_end, scale=self.scale, format='jd')
        t_frac = t_start + (t_end - t_start) * y_frac

        self.jd1, self.jd2 = day_frac(t_frac.jd1, t_frac.jd2)

    @property
    def value(self):
        iy_start, ims, ids, ihmsfs = erfa_time.jd_dtf(self.scale.upper().encode('utf8'),
                                                 0,  # precision=0
                                                 self.jd1, self.jd2)
        imon = np.ones_like(iy_start)
        iday = np.ones_like(iy_start)
        ihr = np.zeros_like(iy_start)
        imin = np.zeros_like(iy_start)
        isec = np.zeros_like(self.jd1)

        # Possible enhancement: use np.unique to only compute start, stop
        # for unique values of iy_start.
        jd1_start, jd2_start = erfa_time.dtf_jd(
            self.scale.upper().encode('utf8'), iy_start, imon, iday, ihr, imin, isec)
        jd1_end, jd2_end = erfa_time.dtf_jd(
            self.scale.upper().encode('utf8'), iy_start + 1, imon, iday, ihr, imin, isec)

        dt = (self.jd1 - jd1_start) + (self.jd2 - jd2_start)
        dt_end = (jd1_end - jd1_start) + (jd2_end - jd2_start)
        decimalyear = iy_start + dt / dt_end

        return decimalyear


class TimeFromEpoch(TimeFormat):
    """
    Base class for times that represent the interval from a particular
    epoch as a floating point multiple of a unit time interval (e.g. seconds
    or days).
    """
    def __init__(self, val1, val2, scale, precision,
                 in_subfmt, out_subfmt, from_jd=False):
        self.scale = scale
        # Initialize the reference epoch (a single time defined in subclasses)
        epoch = Time(self.epoch_val, self.epoch_val2, scale=self.epoch_scale,
                     format=self.epoch_format)
        self.epoch = epoch

        # Now create the TimeFormat object as normal
        super(TimeFromEpoch, self).__init__(val1, val2, scale, precision,
                                            in_subfmt, out_subfmt, from_jd)

    def set_jds(self, val1, val2):
        """
        Initialize the internal jd1 and jd2 attributes given val1 and val2.
        For an TimeFromEpoch subclass like TimeUnix these will be floats giving
        the effective seconds since an epoch time (e.g. 1970-01-01 00:00:00).
        """
        # Form new JDs based on epoch time + time from epoch (converted to JD).
        # One subtlety that might not be obvious is that 1.000 Julian days in
        # UTC can be 86400 or 86401 seconds.  For the TimeUnix format the
        # assumption is that every day is exactly 86400 seconds, so this is, in
        # principle, doing the math incorrectly, *except* that it matches the
        # definition of Unix time which does not include leap seconds.

        # note: use divisor=1./self.unit, since this is either 1 or 1/86400,
        # and 1/86400 is not exactly representable as a float64, so multiplying
        # by that will cause rounding errors. (But inverting it as a float64
        # recovers the exact number)
        day, frac = day_frac(val1, val2, divisor=1. / self.unit)

        jd1 = self.epoch.jd1 + day
        jd2 = self.epoch.jd2 + frac

        # Create a temporary Time object corresponding to the new (jd1, jd2) in
        # the epoch scale (e.g. UTC for TimeUnix) then convert that to the
        # desired time scale for this object.
        #
        # A known limitation is that the transform from self.epoch_scale to
        # self.scale cannot involve any metadata like lat or lon.
        try:
            tm = getattr(Time(jd1, jd2, scale=self.epoch_scale,
                              format='jd'), self.scale)
        except Exception as err:
            raise ScaleValueError("Cannot convert from '{0}' epoch scale '{1}'"
                                  "to specified scale '{2}', got error:\n{3}"
                                  .format(self.name, self.epoch_scale,
                                          self.scale, err))

        self.jd1 = tm._time.jd1
        self.jd2 = tm._time.jd2

    @property
    def value(self):
        # when we get here, getattr will already have ensured our scale
        # equals epoch scale, so we can just subtract the epoch and convert
        time_from_epoch = ((self.jd1 - self.epoch.jd1) +
                           (self.jd2 - self.epoch.jd2)) / self.unit
        return time_from_epoch


class TimeUnix(TimeFromEpoch):
    """
    Unix time: seconds from 1970-01-01 00:00:00 UTC.
    For example, 946684800.0 in Unix time is midnight on January 1, 2000.

    NOTE: this quantity is not exactly unix time and differs from the strict
    POSIX definition by up to 1 second on days with a leap second.  POSIX
    unix time actually jumps backward by 1 second at midnight on leap second
    days while this class value is monotonically increasing at 86400 seconds
    per UTC day.
    """
    name = 'unix'
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = '1970-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'utc'
    epoch_format = 'iso'


class TimeCxcSec(TimeFromEpoch):
    """
    Chandra X-ray Center seconds from 1998-01-01 00:00:00 TT.
    For example, 63072064.184 is midnight on January 1, 2000.
    """
    name = 'cxcsec'
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = '1998-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'tt'
    epoch_format = 'iso'


class TimeGPS(TimeFromEpoch):
    """GPS time: seconds from 1980-01-06 00:00:00 UTC
    For example, 630720013.0 is midnight on January 1, 2000.

    Notes
    =====
    This implementation is strictly a representation of the number of seconds
    (including leap seconds) since midnight UTC on 1980-01-06.  GPS can also be
    considered as a time scale which is ahead of TAI by a fixed offset
    (to within about 100 nanoseconds).

    For details, see http://tycho.usno.navy.mil/gpstt.html
    """
    name = 'gps'
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = '1980-01-06 00:00:19'
    # above epoch is the same as Time('1980-01-06 00:00:00', scale='utc').tai
    epoch_val2 = None
    epoch_scale = 'tai'
    epoch_format = 'iso'


class TimePlotDate(TimeFromEpoch):
    """
    Matplotlib `~matplotlib.pyplot.plot_date` input:
    1 + number of days from 0001-01-01 00:00:00 UTC

    This can be used directly in the matplotlib `~matplotlib.pyplot.plot_date`
    function::

      >>> import matplotlib.pyplot as plt
      >>> jyear = np.linspace(2000, 2001, 20)
      >>> t = Time(jyear, format='jyear', scale='utc')
      >>> plt.plot_date(t.plot_date, jyear)
      >>> plt.gcf().autofmt_xdate()  # orient date labels at a slant
      >>> plt.draw()

    For example, 730120.0003703703 is midnight on January 1, 2000.
    """
    # This corresponds to the zero reference time for matplotlib plot_date().
    # Note that TAI and UTC are equivalent at the reference time.
    name = 'plot_date'
    unit = 1.0
    epoch_val = 1721424.5  # Time('0001-01-01 00:00:00', scale='tai').jd - 1
    epoch_val2 = None
    epoch_scale = 'utc'
    epoch_format = 'jd'


class TimeUnique(TimeFormat):
    """
    Base class for time formats that can uniquely create a time object
    without requiring an explicit format specifier.  This class does
    nothing but provide inheritance to identify a class as unique.
    """


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
        val1_0 = val1.flat[0]
        if not (isinstance(val1_0, Time) and all(type(val) is type(val1_0)
                                                 for val in val1.flat)):
            raise TypeError('Input values for {0} class must all be same '
                            'astropy Time type.'.format(cls.name))

        if scale is None:
            scale = val1_0.scale
        if val1.shape:
            vals = [getattr(val, scale)._time for val in val1]
            jd1 = np.concatenate([np.atleast_1d(val.jd1) for val in vals])
            jd2 = np.concatenate([np.atleast_1d(val.jd2) for val in vals])
        else:
            val = getattr(val1_0, scale)._time
            jd1, jd2 = val.jd1, val.jd2

        OutTimeFormat = val1_0._time.__class__
        self = OutTimeFormat(jd1, jd2, scale, precision, in_subfmt, out_subfmt,
                             from_jd=True)

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
        # Note: don't care about val2 for this class
        if not all(isinstance(val, datetime) for val in val1.flat):
            raise TypeError('Input values for {0} class must be '
                            'datetime objects'.format(self.name))
        return val1, None

    def set_jds(self, val1, val2):
        """Convert datetime object contained in val1 to jd1, jd2"""
        # Iterate through the datetime objects, getting year, month, etc.
        iterator = np.nditer([val1, None, None, None, None, None, None],
                             flags=['refs_ok'],
                             op_dtypes=[np.object] + 5*[np.intc] + [np.double])
        for val, iy, im, id, ihr, imin, dsec in iterator:
            dt = val.item()

            if dt.tzinfo is not None:
                dt = (dt - dt.utcoffset()).replace(tzinfo=None)

            iy[...] = dt.year
            im[...] = dt.month
            id[...] = dt.day
            ihr[...] = dt.hour
            imin[...] = dt.minute
            dsec[...] = dt.second + dt.microsecond / 1e6

        self.jd1, self.jd2 = erfa_time.dtf_jd(
            self.scale.upper().encode('utf8'), *iterator.operands[1:])

    @property
    def value(self):
        iys, ims, ids, ihmsfs = erfa_time.jd_dtf(self.scale.upper()
                                                 .encode('utf8'),
                                                 6,  # precision=6 for microsec
                                                 self.jd1, self.jd2)
        ihrs = ihmsfs[..., 0]
        imins = ihmsfs[..., 1]
        isecs = ihmsfs[..., 2]
        ifracs = ihmsfs[..., 3]
        iterator = np.nditer([iys, ims, ids, ihrs, imins, isecs, ifracs, None],
                             flags=['refs_ok'],
                             op_dtypes=7*[iys.dtype] + [np.object])
        for iy, im, id, ihr, imin, isec, ifracsec, out in iterator:
            out[...] = datetime(iy, im, id, ihr, imin, isec, ifracsec)

        return iterator.operands[-1]


class TimeString(TimeUnique):
    """
    Base class for string-like time representations.

    This class assumes that anything following the last decimal point to the
    right is a fraction of a second.

    This is a reference implementation can be made much faster with effort.
    """
    def _check_val_type(self, val1, val2):
        # Note: don't care about val2 for these classes
        if val1.dtype.kind not in ('S', 'U'):
            raise TypeError('Input values for {0} class must be strings'
                            .format(self.name))
        return val1, None

    def parse_string(self, timestr, subfmts):
        """Read time from a single string, using a set of possible formats."""
        # Datetime components required for conversion to JD by ERFA, along
        # with the default values.
        components = ('year', 'mon', 'mday', 'hour', 'min', 'sec')
        defaults = (None, 1, 1, 0, 0, 0)
        # Assume that anything following "." on the right side is a
        # floating fraction of a second.
        try:
            idot = timestr.rindex('.')
        except:
            fracsec = 0.0
        else:
            timestr, fracsec = timestr[:idot], timestr[idot:]
            fracsec = float(fracsec)

        for _, strptime_fmt_or_regex, _ in subfmts:
            if isinstance(strptime_fmt_or_regex, six.string_types):
                try:
                    tm = time.strptime(timestr, strptime_fmt_or_regex)
                except ValueError:
                    continue
                else:
                    vals = [getattr(tm, 'tm_' + component)
                            for component in components]

            else:
                tm = re.match(strptime_fmt_or_regex, timestr)
                if tm is None:
                    continue
                tm = tm.groupdict()
                vals = [int(tm.get(component, default)) for component, default
                        in six.moves.zip(components, defaults)]

            # Add fractional seconds
            vals[-1] = vals[-1] + fracsec
            return vals
        else:
            raise ValueError('Time {0} does not match {1} format'
                             .format(timestr, self.name))

    def set_jds(self, val1, val2):
        """Parse the time strings contained in val1 and set jd1, jd2"""
        # Select subformats based on current self.in_subfmt
        subfmts = self._select_subfmts(self.in_subfmt)

        iterator = np.nditer([val1, None, None, None, None, None, None],
                             op_dtypes=[val1.dtype] + 5*[np.intc] + [np.double])


        for val, iy, im, id, ihr, imin, dsec in iterator:
            iy[...], im[...], id[...], ihr[...], imin[...], dsec[...] = (
                self.parse_string(val.item(), subfmts))

        self.jd1, self.jd2 = erfa_time.dtf_jd(
            self.scale.upper().encode('utf8'), *iterator.operands[1:])

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

        ihrs = ihmsfs[..., 0]
        imins = ihmsfs[..., 1]
        isecs = ihmsfs[..., 2]
        ifracs = ihmsfs[..., 3]
        for iy, im, id, ihr, imin, isec, ifracsec in np.nditer(
                [iys, ims, ids, ihrs, imins, isecs, ifracs]):
            if has_yday:
                yday = datetime(iy, im, id).timetuple().tm_yday

            yield {'year': int(iy), 'mon': int(im), 'day': int(id),
                   'hour': int(ihr), 'min': int(imin), 'sec': int(isec),
                   'fracsec': int(ifracsec), 'yday': yday}

    def format_string(self, str_fmt, **kwargs):
        """Write time to a string using a given format.

        By default, just interprets str_fmt as a format string,
        but subclasses can add to this.
        """
        return str_fmt.format(**kwargs)

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
            outs.append(str(self.format_string(str_fmt, **kwargs)))

        return np.array(outs).reshape(self.jd1.shape)

    def _select_subfmts(self, pattern):
        """
        Return a list of subformats where name matches ``pattern`` using
        fnmatch.
        """

        fnmatchcase = fnmatch.fnmatchcase
        subfmts = [x for x in self.subfmts if fnmatchcase(x[0], pattern)]
        if len(subfmts) == 0:
            raise ValueError('No subformats match {0}'.format(pattern))
        return subfmts


class TimeISO(TimeString):
    """
    ISO 8601 compliant date-time format "YYYY-MM-DD HH:MM:SS.sss...".
    For example, 2000-01-01 00:00:00.000 is midnight on January 1, 2000.

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

    def parse_string(self, timestr, subfmts):
        # Handle trailing 'Z' for UTC time
        if timestr.endswith('Z'):
            if self.scale != 'utc':
                raise ValueError("Time input terminating in 'Z' must have "
                                 "scale='UTC'")
            timestr = timestr[:-1]
        return super(TimeISO, self).parse_string(timestr, subfmts)


class TimeISOT(TimeISO):
    """
    ISO 8601 compliant date-time format "YYYY-MM-DDTHH:MM:SS.sss...".
    This is the same as TimeISO except for a "T" instead of space between
    the date and time.
    For example, 2000-01-01T00:00:00.000 is midnight on January 1, 2000.

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


class TimeYearDayTime(TimeISO):
    """
    Year, day-of-year and time as "YYYY:DOY:HH:MM:SS.sss...".
    The day-of-year (DOY) goes from 001 to 365 (366 in leap years).
    For example, 2000:001:00:00:00.000 is midnight on January 1, 2000.

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


class TimeFITS(TimeString):
    """
    FITS format: "[±Y]YYYY-MM-DD[THH:MM:SS[.sss]][(SCALE[(REALIZATION)])]".

    ISOT with two extensions:
    - Can give signed five-digit year (mostly for negative years);
    - A possible time scale (and realization) appended in parentheses.

    Note: FITS supports some deprecated names for timescales; these are
    translated to the formal names upon initialization.  Furthermore, any
    specific realization information is stored only as long as the time scale
    is not changed.

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date': date
    - 'longdate_hms': as 'date_hms', but with signed 5-digit year
    - 'longdate': as 'date', but with signed 5-digit year

    See Rots et al., 2015, A&A 574:A36 (arXiv:1409.7583).
    """
    name = 'fits'
    subfmts = (
        ('date_hms',
         (r'(?P<year>\d{4})-(?P<mon>\d\d)-(?P<mday>\d\d)T'
          r'(?P<hour>\d\d):(?P<min>\d\d):(?P<sec>\d\d(\.\d*)?)'),
         '{year:04d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}:{sec:02d}'),
        ('date',
         r'(?P<year>\d{4})-(?P<mon>\d\d)-(?P<mday>\d\d)',
         '{year:04d}-{mon:02d}-{day:02d}'),
        ('longdate_hms',
         (r'(?P<year>[+-]\d{5})-(?P<mon>\d\d)-(?P<mday>\d\d)T'
          r'(?P<hour>\d\d):(?P<min>\d\d):(?P<sec>\d\d(\.\d*)?)'),
         '{year:+06d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}:{sec:02d}'),
        ('longdate',
         r'(?P<year>[+-]\d{5})-(?P<mon>\d\d)-(?P<mday>\d\d)',
         '{year:+06d}-{mon:02d}-{day:02d}'))
    # Add the regex that parses the scale and possible realization.
    subfmts = tuple(
        (subfmt[0],
         subfmt[1] + r'(\((?P<scale>\w+)(\((?P<realization>\w+)\))?\))?',
         subfmt[2]) for subfmt in subfmts)
    _fits_scale = None
    _fits_realization = None

    def parse_string(self, timestr, subfmts):
        """Read time and set scale according to trailing scale codes."""
        # Try parsing with any of the allowed sub-formats.
        for _, regex, _ in subfmts:
            tm = re.match(regex, timestr)
            if tm:
                break
        else:
            raise ValueError('Time {0} does not match {1} format'
                             .format(timestr, self.name))
        tm = tm.groupdict()
        if tm['scale'] is not None:
            # If a scale was given, translate from a possible deprecated
            # timescale identifier to the scale used by Time.
            fits_scale = tm['scale'].upper()
            scale = FITS_DEPRECATED_SCALES.get(fits_scale, fits_scale.lower())
            if scale not in TIME_SCALES:
                raise ValueError("Scale {0} is not in the allowed scales {1}"
                                 .format(repr(scale), sorted(TIME_SCALES)))
            # If no scale was given in the initialiser, set the scale to
            # that given in the string.  Also store a possible realization,
            # so we can round-trip (as long as no scale changes are made).
            fits_realization = (tm['realization'].upper()
                                if tm['realization'] else None)
            if self._scale is None:
                self._scale = scale
                self._fits_scale = fits_scale
                self._fits_realization = fits_realization
            elif (scale != self.scale or fits_scale != self._fits_scale or
                  fits_realization != self._fits_realization):
                raise ValueError("Input strings for {0} class must all "
                                 "have consistent time scales."
                                 .format(self.name))
        return [int(tm['year']), int(tm['mon']), int(tm['mday']),
                int(tm.get('hour', 0)), int(tm.get('min', 0)),
                float(tm.get('sec', 0.))]

    def format_string(self, str_fmt, **kwargs):
        """Format time-string: append the scale to the normal ISOT format."""
        time_str = super(TimeFITS, self).format_string(str_fmt, **kwargs)
        if self._fits_scale and self._fits_realization:
            return '{0}({1}({2}))'.format(time_str, self._fits_scale,
                                          self._fits_realization)
        else:
            return '{0}({1})'.format(time_str, self._scale.upper())

    @property
    def value(self):
        """Convert times to strings, using signed 5 digit if necessary."""
        if 'long' not in self.out_subfmt:
            # If we have times before year 0 or after year 9999, we can
            # output only in a "long" format, using signed 5-digit years.
            jd = self.jd1 + self.jd2
            if jd.min() < 1721425.5 or jd.max() >= 5373484.5:
                self.out_subfmt = 'long' + self.out_subfmt
        return super(TimeFITS, self).value


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

    def _check_val_type(self, val1, val2):
        """Input value validation, typically overridden by derived classes"""
        if hasattr(val1, 'to') and hasattr(val1, 'unit'):
            raise ValueError("Cannot use Quantities for 'byear' format, "
                             "as the interpretation would be ambiguous. "
                             "Use float with Besselian year instead. ")

        return super(TimeBesselianEpoch, self)._check_val_type(val1, val2)


class TimeJulianEpoch(TimeEpochDate):
    """Julian Epoch year as floating point value(s) like 2000.0"""
    name = 'jyear'
    unit = erfa.DJY  # 365.25, the Julian year, for conversion to quantities
    epoch_to_jd = 'julian_epoch_jd'
    jd_to_epoch = 'jd_julian_epoch'


class TimeEpochDateString(TimeString):
    """
    Base class to support string Besselian and Julian epoch dates
    such as 'B1950.0' or 'J2000.0' respectively.
    """
    def set_jds(self, val1, val2):
        epoch_prefix = self.epoch_prefix
        iterator = np.nditer([val1, None], op_dtypes=[val1.dtype, np.double])
        for val, years in iterator:
            time_str = val.item()
            try:
                epoch_type, year_str = time_str[0], time_str[1:]
                year = float(year_str)
                if epoch_type.upper() != epoch_prefix:
                    raise ValueError
            except (IndexError, ValueError):
                raise ValueError('Time {0} does not match {1} format'
                                 .format(time_str, self.name))
            else:
                years[...] = year

        self._check_scale(self._scale)  # validate scale.
        epoch_to_jd = getattr(erfa_time, self.epoch_to_jd)
        self.jd1, self.jd2 = epoch_to_jd(iterator.operands[-1])

    @property
    def value(self):
        jd_to_epoch = getattr(erfa_time, self.jd_to_epoch)
        years = jd_to_epoch(self.jd1, self.jd2)
        # Use old-style format since it is a factor of 2 faster
        str_fmt = self.epoch_prefix + '%.' + str(self.precision) + 'f'
        outs = [str_fmt % year for year in years.flat]
        return np.array(outs).reshape(self.jd1.shape)


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


class TimeDeltaFormatMeta(TimeFormatMeta):
    _registry = TIME_DELTA_FORMATS


@six.add_metaclass(TimeDeltaFormatMeta)
class TimeDeltaFormat(TimeFormat):
    """Base class for time delta representations"""

    def _check_scale(self, scale):
        """
        Check that the scale is in the allowed list of scales, or is `None`
        """
        if scale is not None and scale not in TIME_DELTA_SCALES:
            raise ScaleValueError("Scale value '{0}' not in "
                                  "allowed values {1}"
                                  .format(scale, TIME_DELTA_SCALES))

        return scale

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2, divisor=1./self.unit)

    @property
    def value(self):
        return (self.jd1 + self.jd2) / self.unit


class TimeDeltaSec(TimeDeltaFormat):
    """Time delta in SI seconds"""
    name = 'sec'
    unit = 1. / erfa.DAYSEC  # for quantity input


class TimeDeltaJD(TimeDeltaFormat):
    """Time delta in Julian days (86400 SI seconds)"""
    name = 'jd'
    unit = 1.


class ScaleValueError(Exception):
    pass


def _make_array(val, copy=False):
    """
    Take ``val`` and convert/reshape to an array.  If ``copy`` is `True`
    then copy input values.

    Returns
    -------
    val : ndarray
        Array version of ``val``.
    """
    val = np.array(val, copy=copy, subok=True)
    # Allow only float64, string or object arrays as input
    # (object is for datetime, maybe add more specific test later?)
    # This also ensures the right byteorder for float64 (closes #2942).
    if not (val.dtype == np.float64 or val.dtype.kind in 'OSUa'):
        val = np.asanyarray(val, dtype=np.float64)

    return val


def day_frac(val1, val2, factor=1., divisor=1.):
    """
    Return the sum of ``val1`` and ``val2`` as two float64s, an integer part
    and the fractional remainder.  If ``factor`` is not 1.0 then multiply the
    sum by ``factor``.  If ``divisor`` is not 1.0 then divide the sum by
    ``divisor``.

    The arithmetic is all done with exact floating point operations so no
    precision is lost to rounding error.  This routine assumes the sum is less
    than about 1e16, otherwise the ``frac`` part will be greater than 1.0.

    Returns
    -------
    day, frac : float64
        Integer and fractional part of val1 + val2.
    """
    # Add val1 and val2 exactly, returning the result as two float64s.
    # The first is the approximate sum (with some floating point error)
    # and the second is the error of the float64 sum.
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
    Add ``a`` and ``b`` exactly, returning the result as two float64s.
    The first is the approximate sum (with some floating point error)
    and the second is the error of the float64 sum.

    Using the procedure of Shewchuk, 1997,
    Discrete & Computational Geometry 18(3):305-363
    http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

    Returns
    -------
    sum, err : float64
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
    Multiple ``a`` and ``b`` exactly, returning the result as two float64s.
    The first is the approximate product (with some floating point error)
    and the second is the error of the float64 product.

    Uses the procedure of Shewchuk, 1997,
    Discrete & Computational Geometry 18(3):305-363
    http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

    Returns
    -------
    prod, err : float64
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
    def __init__(self, left, right, op=None):
        op_string = '' if op is None else ' for {0}'.format(op)
        super(OperandTypeError, self).__init__(
            "Unsupported operand type(s){0}: "
            "'{1}' and '{2}'".format(op_string,
                                     left.__class__.__name__,
                                     right.__class__.__name__))
