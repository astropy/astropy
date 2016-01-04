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

import itertools
from datetime import datetime
from collections import defaultdict

import numpy as np

from .. import units as u
from .. import _erfa as erfa
from ..units import UnitConversionError
from ..utils.decorators import lazyproperty
from ..utils.compat.numpycompat import NUMPY_LT_1_7
from ..utils.compat.misc import override__dir__
from ..utils.data_info import MixinInfo, data_info_factory
from ..utils.compat.numpy import broadcast_to
from ..extern import six
from .utils import day_frac
from .formats import (TIME_FORMATS, TIME_DELTA_FORMATS,
                      TimeJD, TimeUnique, TimeAstropyTime, TimeDatetime)
# Import TimeFromEpoch to avoid breaking code that followed the old example of
# making a custom timescale in the documentation.
from .formats import TimeFromEpoch


__all__ = ['Time', 'TimeDelta', 'TIME_SCALES', 'TIME_DELTA_SCALES',
           'ScaleValueError', 'OperandTypeError']


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

# triple-level dictionary, yay!
SIDEREAL_TIME_MODELS = {
    'mean': {
        'IAU2006': {'function': erfa.gmst06, 'scales': ('ut1', 'tt')},
        'IAU2000': {'function': erfa.gmst00, 'scales': ('ut1', 'tt')},
        'IAU1982': {'function': erfa.gmst82, 'scales': ('ut1',)}},
    'apparent': {
        'IAU2006A': {'function': erfa.gst06a, 'scales': ('ut1', 'tt')},
        'IAU2000A': {'function': erfa.gst00a, 'scales': ('ut1', 'tt')},
        'IAU2000B': {'function': erfa.gst00b, 'scales': ('ut1',)},
        'IAU1994': {'function': erfa.gst94, 'scales': ('ut1',)}}}


class TimeInfo(MixinInfo):
    """
    Container for meta information like name, description, format.  This is
    required when the object is used as a mixin column within a table, but can
    be used as a general way to store meta information.
    """
    attrs_from_parent = set(['unit'])  # unit is read-only and None
    _supports_indexing = True

    @property
    def unit(self):
        return None

    info_summary_stats = staticmethod(
        data_info_factory(names=MixinInfo._stats,
                          funcs=[getattr(np, stat) for stat in MixinInfo._stats]))
    # When Time has mean, std, min, max methods:
    # funcs = [lambda x: getattr(x, stat)() for stat_name in MixinInfo._stats])

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

        if isinstance(val, self.__class__):
            # Update _time formatting parameters if explicitly specified
            if precision is not None:
                self._time.precision = precision
            if in_subfmt is not None:
                self._time.in_subfmt = in_subfmt
            if out_subfmt is not None:
                self._time.out_subfmt = out_subfmt

            if scale is not None:
                self._set_scale(scale)
        else:
            self._init_from_vals(val, val2, format, scale, copy,
                                 precision, in_subfmt, out_subfmt)

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

    def _init_from_vals(self, val, val2, format, scale, copy,
                        precision=None, in_subfmt=None, out_subfmt=None):
        """
        Set the internal _format, scale, and _time attrs from user
        inputs.  This handles coercion into the correct shapes and
        some basic input validation.
        """
        if precision is None:
            precision = 3
        if in_subfmt is None:
            in_subfmt = '*'
        if out_subfmt is None:
            out_subfmt = '*'

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
        self._time = self._get_time_fmt(val, val2, format, scale,
                                        precision, in_subfmt, out_subfmt)
        self._format = self._time.name

    def _get_time_fmt(self, val, val2, format, scale,
                      precision, in_subfmt, out_subfmt):
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
                return FormatClass(val, val2, scale, precision, in_subfmt, out_subfmt)
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

    info = TimeInfo()

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

            conv_func = getattr(erfa, sys1 + sys2)
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
        return self._time.precision

    @precision.setter
    def precision(self, val):
        if not isinstance(val, int) or val < 0 or val > 9:
            raise ValueError('precision attribute must be an int between '
                             '0 and 9')
        self._time.precision = val
        del self.cache

    @property
    def in_subfmt(self):
        """
        Unix wildcard pattern to select subformats for parsing string input
        times.
        """
        return self._time.in_subfmt

    @in_subfmt.setter
    def in_subfmt(self, val):
        if not isinstance(val, six.string_types):
            raise ValueError('in_subfmt attribute must be a string')
        self._time.in_subfmt = val
        del self.cache

    @property
    def out_subfmt(self):
        """
        Unix wildcard pattern to select subformats for outputting times.
        """
        return self._time.out_subfmt

    @out_subfmt.setter
    def out_subfmt(self, val):
        if not isinstance(val, six.string_types):
            raise ValueError('out_subfmt attribute must be a string')
        self._time.out_subfmt = val
        del self.cache

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
        # The underlying way to get the time values for the current format is:
        #     self._shaped_like_input(self._time.to_value(parent=self))
        # This is done in __getattr__.  By calling getattr(self, self.format)
        # the ``value`` attribute is cached.
        return getattr(self, self.format)

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

    if isinstance(sidereal_time.__doc__, six.string_types):
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

            # Apply the method to any value arrays (though skip if there is only
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

        # Copy other 'info' attr only if it has actually been defined.
        # See PR #3898 for further explanation and justification, along
        # with Quantity.__array_finalize__
        if 'info' in self.__dict__:
            tm.info = self.info

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

    def __iter__(self):
        if self.isscalar:
            raise TypeError('scalar {0!r} object is not iterable.'.format(
                self.__class__.__name__))

        def time_iter():
            try:
                for idx in itertools.count():
                    yield self[idx]
            except IndexError:
                # Results in StopIteration
                pass

        return time_iter()

    def __getitem__(self, item):
        if self.isscalar:
            raise TypeError('scalar {0!r} object is not subscriptable.'.format(
                self.__class__.__name__))
        return self._replicate('__getitem__', item)

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

    def _advanced_index(self, indices, axis=None, keepdims=False):
        """Turn argmin, argmax output into an advanced index.

        Argmin, argmax output contains indices along a given axis in an array
        shaped like the other dimensions.  To use this to get values at the
        correct location, a list is constructed in which the other axes are
        indexed sequentially.  For ``keepdims`` is ``True``, the net result is
        the same as constructing an index grid with ``np.ogrid`` and then
        replacing the ``axis`` item with ``indices`` with its shaped expanded
        at ``axis``. For ``keepdims`` is ``False``, the result is the same but
        with the ``axis`` dimension removed from all list entries.

        For ``axis`` is ``None``, this calls :func:`~numpy.unravel_index`.

        Parameters
        ----------
        indices : array
            Output of argmin or argmax.
        axis : int or None
            axis along which argmin or argmax was used.
        keepdims : bool
            Whether to construct indices that keep or remove the axis along
            which argmin or argmax was used.  Default: ``False``.

        Returns
        -------
        advanced_index : list of arrays
            Suitable for use as an advanced index.
        """
        if axis is None:
            return np.unravel_index(indices, self.shape)

        ndim = self.ndim
        if axis < 0:
            axis = axis + ndim

        if keepdims and indices.ndim < self.ndim:
            indices = np.expand_dims(indices, axis)
        return [(indices if i == axis else np.arange(s).reshape(
            (1,)*(i if keepdims or i < axis else i-1) + (s,) +
            (1,)*(ndim-i-(1 if keepdims or i > axis else 2))))
                for i, s in enumerate(self.shape)]

    def argmin(self, axis=None, out=None):
        """Return indices of the minimum values along the given axis.

        This is similar to :meth:`~numpy.ndarray.argmin`, but adapted to ensure
        that the full precision given by the two doubles ``jd1`` and ``jd2``
        is used.  See :func:`~numpy.argmin` for detailed documentation.
        """
        # first get the minimum at normal precision.
        jd = self.jd1 + self.jd2
        if NUMPY_LT_1_7:
            approx = jd.min(axis)
            if axis is not None:
                approx = np.expand_dims(approx, axis)
        else:
            approx = jd.min(axis, keepdims=True)

        # Approx is very close to the true minimum, and by subtracting it at
        # full precision, all numbers near 0 can be represented correctly,
        # so we can be sure we get the true minimum.
        # The below is effectively what would be done for
        # dt = (self - self.__class__(approx, format='jd')).jd
        # which translates to:
        # approx_jd1, approx_jd2 = day_frac(approx, 0.)
        # dt = (self.jd1 - approx_jd1) + (self.jd2 - approx_jd2)
        dt = (self.jd1 - approx) + self.jd2
        return dt.argmin(axis, out)

    def argmax(self, axis=None, out=None):
        """Return indices of the maximum values along the given axis.

        This is similar to :meth:`~numpy.ndarray.argmax`, but adapted to ensure
        that the full precision given by the two doubles ``jd1`` and ``jd2``
        is used.  See :func:`~numpy.argmax` for detailed documentation.
        """
        # For procedure, see comment on argmin.
        jd = self.jd1 + self.jd2
        if NUMPY_LT_1_7:
            approx = jd.max(axis)
            if axis is not None:
                approx = np.expand_dims(approx, axis)
        else:
            approx = jd.max(axis, keepdims=True)

        dt = (self.jd1 - approx) + self.jd2
        return dt.argmax(axis, out)

    def argsort(self, axis=-1):
        """Returns the indices that would sort the time array.

        This is similar to :meth:`~numpy.ndarray.argsort`, but adapted to ensure
        that the full precision given by the two doubles ``jd1`` and ``jd2``
        is used, and that corresponding attributes are copied.  Internally,
        it uses :func:`~numpy.lexsort`, and hence no sort method can be chosen.
        """
        jd_approx = self.jd
        jd_remainder = (self - self.__class__(jd_approx, format='jd')).jd
        if axis is None:
            return np.lexsort((jd_remainder.ravel(), jd_approx.ravel()))
        else:
            return np.lexsort(keys=(jd_remainder, jd_approx), axis=axis)

    def min(self, axis=None, out=None, keepdims=False):
        """Minimum along a given axis.

        This is similar to :meth:`~numpy.ndarray.min`, but adapted to ensure
        that the full precision given by the two doubles ``jd1`` and ``jd2``
        is used, and that corresponding attributes are copied.

        Note that the ``out`` argument is present only for compatibility with
        ``np.min``; since `Time` instances are immutable, it is not possible
        to have an actual ``out`` to store the result in.
        """
        if out is not None:
            raise ValueError("Since `Time` instances are immutable, ``out`` "
                             "cannot be set to anything but ``None``.")
        return self[self._advanced_index(self.argmin(axis), axis, keepdims)]

    def max(self, axis=None, out=None, keepdims=False):
        """Maximum along a given axis.

        This is similar to :meth:`~numpy.ndarray.max`, but adapted to ensure
        that the full precision given by the two doubles ``jd1`` and ``jd2``
        is used, and that corresponding attributes are copied.

        Note that the ``out`` argument is present only for compatibility with
        ``np.max``; since `Time` instances are immutable, it is not possible
        to have an actual ``out`` to store the result in.
        """
        if out is not None:
            raise ValueError("Since `Time` instances are immutable, ``out`` "
                             "cannot be set to anything but ``None``.")
        return self[self._advanced_index(self.argmax(axis), axis, keepdims)]

    def ptp(self, axis=None, out=None, keepdims=False):
        """Peak to peak (maximum - minimum) along a given axis.

        This is similar to :meth:`~numpy.ndarray.ptp`, but adapted to ensure
        that the full precision given by the two doubles ``jd1`` and ``jd2``
        is used.

        Note that the ``out`` argument is present only for compatibility with
        `~numpy.ptp`; since `Time` instances are immutable, it is not possible
        to have an actual ``out`` to store the result in.
        """
        if out is not None:
            raise ValueError("Since `Time` instances are immutable, ``out`` "
                             "cannot be set to anything but ``None``.")
        return (self.max(axis, keepdims=keepdims) -
                self.min(axis, keepdims=keepdims))

    def sort(self, axis=-1):
        """Return a copy sorted along the specified axis.

        This is similar to :meth:`~numpy.ndarray.sort`, but internally uses
        indexing with :func:`~numpy.lexsort` to ensure that the full precision
        given by the two doubles ``jd1`` and ``jd2`` is kept, and that
        corresponding attributes are properly sorted and copied as well.

        Parameters
        ----------
        axis : int or None
            Axis to be sorted.  If ``None``, the flattened array is sorted.
            By default, sort over the last axis.
        """
        return self[self._advanced_index(self.argsort(axis), axis,
                                         keepdims=True)]

    @lazyproperty
    def cache(self):
        """
        Return the cache associated with this instance.
        """
        return defaultdict(dict)

    def __getattr__(self, attr):
        """
        Get dynamic attributes to output format or do timescale conversion.
        """
        if attr in self.SCALES and self.scale is not None:
            cache = self.cache['scale']
            if attr not in cache:
                if attr == self.scale:
                    tm = self
                else:
                    tm = self.replicate()
                    tm._set_scale(attr)
                cache[attr] = tm
            return cache[attr]

        elif attr in self.FORMATS:
            cache = self.cache['format']
            if attr not in cache:
                if attr == self.format:
                    tm = self
                else:
                    tm = self.replicate(format=attr)
                value = tm._shaped_like_input(tm._time.to_value(parent=tm))
                cache[attr] = value
            return cache[attr]

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
                jd1_utc, jd2_utc = erfa.ut1utc(jd1, jd2, delta)
                # calculate a better estimate using the nearly correct UTC
                delta = iers_table.ut1_utc(jd1_utc, jd2_utc)

            self._set_delta_ut1_utc(delta)

        return self._delta_ut1_utc

    def _set_delta_ut1_utc(self, val):
        val = self._match_shape(val)
        if hasattr(val, 'to'):  # Matches Quantity but also TimeDelta.
            val = val.to(u.second).value
        self._delta_ut1_utc = val
        del self.cache

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
            njd1, njd2 = erfa.tttai(jd1, jd2)
            njd1, njd2 = erfa.taiutc(njd1, njd2)
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
            self._delta_tdb_tt = erfa.dtdb(
                jd1, jd2, ut, lon.to(u.radian).value,
                rxy.to(u.km).value, z.to(u.km).value)

        return self._delta_tdb_tt

    def _set_delta_tdb_tt(self, val):
        val = self._match_shape(val)
        if hasattr(val, 'to'):  # Matches Quantity but also TimeDelta.
            val = val.to(u.second).value
        self._delta_tdb_tt = val
        del self.cache

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
            out = self.replicate()
            if self.scale in other.SCALES:
                if other.scale not in (out.scale, None):
                    other = getattr(other, out.scale)
            else:
                out._set_scale(other.scale if other.scale is not None
                               else 'tai')
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

        if other_is_delta:
            # Go back to left-side scale if needed
            out._set_scale(self.scale)

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
        out = self.replicate()
        if self.scale in other.SCALES:
            if other.scale not in (out.scale, None):
                other = getattr(other, out.scale)
        else:
            out._set_scale(other.scale if other.scale is not None else 'tai')

        # remove attributes that are invalidated by changing time
        for attr in ('_delta_ut1_utc', '_delta_tdb_tt'):
            if hasattr(out, attr):
                delattr(out, attr)

        jd1 = out._time.jd1 + other._time.jd1
        jd2 = out._time.jd2 + other._time.jd2

        out._time.jd1, out._time.jd2 = day_frac(jd1, jd2)

        # Go back to left-side scale if needed
        out._set_scale(self.scale)

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

    def to_datetime(self, timezone=None):
        tm = self.replicate(format='datetime')
        return tm._shaped_like_input(tm._time.to_value(timezone))

    to_datetime.__doc__ = TimeDatetime.to_value.__doc__



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


class OperandTypeError(TypeError):
    def __init__(self, left, right, op=None):
        op_string = '' if op is None else ' for {0}'.format(op)
        super(OperandTypeError, self).__init__(
            "Unsupported operand type(s){0}: "
            "'{1}' and '{2}'".format(op_string,
                                     left.__class__.__name__,
                                     right.__class__.__name__))
