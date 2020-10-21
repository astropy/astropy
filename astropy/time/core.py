# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The astropy.time package provides functionality for manipulating times and
dates. Specific emphasis is placed on supporting time scales (e.g. UTC, TAI,
UT1) and time representations (e.g. JD, MJD, ISO 8601) that are used in
astronomy.
"""

import os
import copy
import operator
from datetime import datetime, date, timedelta
from time import strftime
from warnings import warn

import numpy as np
import erfa

from astropy import units as u, constants as const
from astropy.units import UnitConversionError
from astropy.utils import ShapedLikeNDArray
from astropy.utils.compat.misc import override__dir__
from astropy.utils.data_info import MixinInfo, data_info_factory
from astropy.utils.exceptions import AstropyWarning
from .utils import day_frac
from .formats import (TIME_FORMATS, TIME_DELTA_FORMATS,
                      TimeJD, TimeUnique, TimeAstropyTime, TimeDatetime)
# Import TimeFromEpoch to avoid breaking code that followed the old example of
# making a custom timescale in the documentation.
from .formats import TimeFromEpoch  # noqa

from astropy.extern import _strptime

__all__ = ['TimeBase', 'Time', 'TimeDelta', 'TimeInfo', 'update_leap_seconds',
           'TIME_SCALES', 'STANDARD_TIME_SCALES', 'TIME_DELTA_SCALES',
           'ScaleValueError', 'OperandTypeError']


STANDARD_TIME_SCALES = ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')
LOCAL_SCALES = ('local',)
TIME_TYPES = dict((scale, scales) for scales in (STANDARD_TIME_SCALES, LOCAL_SCALES)
                  for scale in scales)
TIME_SCALES = STANDARD_TIME_SCALES + LOCAL_SCALES
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
                                       ROTATIONAL_SCALES, LOCAL_SCALES) for scale in scales)
TIME_DELTA_SCALES = GEOCENTRIC_SCALES + BARYCENTRIC_SCALES + ROTATIONAL_SCALES + LOCAL_SCALES
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


_LEAP_SECONDS_CHECKED = False


class TimeInfo(MixinInfo):
    """
    Container for meta information like name, description, format.  This is
    required when the object is used as a mixin column within a table, but can
    be used as a general way to store meta information.
    """
    attr_names = MixinInfo.attr_names | {'serialize_method'}
    _supports_indexing = True

    # The usual tuple of attributes needed for serialization is replaced
    # by a property, since Time can be serialized different ways.
    _represent_as_dict_extra_attrs = ('format', 'scale', 'precision',
                                      'in_subfmt', 'out_subfmt', 'location',
                                      '_delta_ut1_utc', '_delta_tdb_tt')

    # When serializing, write out the `value` attribute using the column name.
    _represent_as_dict_primary_data = 'value'

    mask_val = np.ma.masked

    @property
    def _represent_as_dict_attrs(self):
        method = self.serialize_method[self._serialize_context]
        if method == 'formatted_value':
            out = ('value',)
        elif method == 'jd1_jd2':
            out = ('jd1', 'jd2')
        else:
            raise ValueError("serialize method must be 'formatted_value' or 'jd1_jd2'")

        return out + self._represent_as_dict_extra_attrs

    def __init__(self, bound=False):
        super().__init__(bound)

        # If bound to a data object instance then create the dict of attributes
        # which stores the info attribute values.
        if bound:
            # Specify how to serialize this object depending on context.
            # If ``True`` for a context, then use formatted ``value`` attribute
            # (e.g. the ISO time string).  If ``False`` then use float jd1 and jd2.
            self.serialize_method = {'fits': 'jd1_jd2',
                                     'ecsv': 'formatted_value',
                                     'hdf5': 'jd1_jd2',
                                     'yaml': 'jd1_jd2',
                                     None: 'jd1_jd2'}

    def get_sortable_arrays(self):
        """
        Return a list of arrays which can be lexically sorted to represent
        the order of the parent column.

        Returns
        -------
        arrays : list of ndarray
        """
        parent = self._parent
        jd_approx = parent.jd
        jd_remainder = (parent - parent.__class__(jd_approx, format='jd')).jd
        return [jd_approx, jd_remainder]

    @property
    def unit(self):
        return None

    info_summary_stats = staticmethod(
        data_info_factory(names=MixinInfo._stats,
                          funcs=[getattr(np, stat) for stat in MixinInfo._stats]))
    # When Time has mean, std, min, max methods:
    # funcs = [lambda x: getattr(x, stat)() for stat_name in MixinInfo._stats])

    def _construct_from_dict_base(self, map):
        if 'jd1' in map and 'jd2' in map:
            # Initialize as JD but revert to desired format and out_subfmt (if needed)
            format = map.pop('format')
            out_subfmt = map.pop('out_subfmt', None)
            map['format'] = 'jd'
            map['val'] = map.pop('jd1')
            map['val2'] = map.pop('jd2')
            out = self._parent_cls(**map)
            out.format = format
            if out_subfmt is not None:
                out.out_subfmt = out_subfmt

        else:
            map['val'] = map.pop('value')
            out = self._parent_cls(**map)

        return out

    def _construct_from_dict(self, map):
        delta_ut1_utc = map.pop('_delta_ut1_utc', None)
        delta_tdb_tt = map.pop('_delta_tdb_tt', None)

        out = self._construct_from_dict_base(map)

        if delta_ut1_utc is not None:
            out._delta_ut1_utc = delta_ut1_utc
        if delta_tdb_tt is not None:
            out._delta_tdb_tt = delta_tdb_tt

        return out

    def new_like(self, cols, length, metadata_conflicts='warn', name=None):
        """
        Return a new Time instance which is consistent with the input Time objects
        ``cols`` and has ``length`` rows.

        This is intended for creating an empty Time instance whose elements can
        be set in-place for table operations like join or vstack.  It checks
        that the input locations and attributes are consistent.  This is used
        when a Time object is used as a mixin column in an astropy Table.

        Parameters
        ----------
        cols : list
            List of input columns (Time objects)
        length : int
            Length of the output column object
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output column name

        Returns
        -------
        col : Time (or subclass)
            Empty instance of this class consistent with ``cols``

        """
        # Get merged info attributes like shape, dtype, format, description, etc.
        attrs = self.merge_cols_attributes(cols, metadata_conflicts, name,
                                           ('meta', 'description'))
        attrs.pop('dtype')  # Not relevant for Time
        col0 = cols[0]

        # Check that location is consistent for all Time objects
        for col in cols[1:]:
            # This is the method used by __setitem__ to ensure that the right side
            # has a consistent location (and coerce data if necessary, but that does
            # not happen in this case since `col` is already a Time object).  If this
            # passes then any subsequent table operations via setitem will work.
            try:
                col0._make_value_equivalent(slice(None), col)
            except ValueError:
                raise ValueError('input columns have inconsistent locations')

        # Make a new Time object with the desired shape and attributes
        shape = (length,) + attrs.pop('shape')
        jd2000 = 2451544.5  # Arbitrary JD value J2000.0 that will work with ERFA
        jd1 = np.full(shape, jd2000, dtype='f8')
        jd2 = np.zeros(shape, dtype='f8')
        tm_attrs = {attr: getattr(col0, attr)
                    for attr in ('scale', 'location',
                                 'precision', 'in_subfmt', 'out_subfmt')}
        out = self._parent_cls(jd1, jd2, format='jd', **tm_attrs)
        out.format = col0.format

        # Set remaining info attributes
        for attr, value in attrs.items():
            setattr(out.info, attr, value)

        return out


class TimeDeltaInfo(TimeInfo):
    _represent_as_dict_extra_attrs = ('format', 'scale')

    def _construct_from_dict(self, map):
        return self._construct_from_dict_base(map)

    def new_like(self, cols, length, metadata_conflicts='warn', name=None):
        """
        Return a new TimeDelta instance which is consistent with the input Time objects
        ``cols`` and has ``length`` rows.

        This is intended for creating an empty Time instance whose elements can
        be set in-place for table operations like join or vstack.  It checks
        that the input locations and attributes are consistent.  This is used
        when a Time object is used as a mixin column in an astropy Table.

        Parameters
        ----------
        cols : list
            List of input columns (Time objects)
        length : int
            Length of the output column object
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output column name

        Returns
        -------
        col : Time (or subclass)
            Empty instance of this class consistent with ``cols``

        """
        # Get merged info attributes like shape, dtype, format, description, etc.
        attrs = self.merge_cols_attributes(cols, metadata_conflicts, name,
                                           ('meta', 'description'))
        attrs.pop('dtype')  # Not relevant for Time
        col0 = cols[0]

        # Make a new Time object with the desired shape and attributes
        shape = (length,) + attrs.pop('shape')
        jd1 = np.zeros(shape, dtype='f8')
        jd2 = np.zeros(shape, dtype='f8')
        out = self._parent_cls(jd1, jd2, format='jd', scale=col0.scale)
        out.format = col0.format

        # Set remaining info attributes
        for attr, value in attrs.items():
            setattr(out.info, attr, value)

        return out


class TimeBase(ShapedLikeNDArray):
    """Base time class from which Time and TimeDelta inherit."""

    # Make sure that reverse arithmetic (e.g., TimeDelta.__rmul__)
    # gets called over the __mul__ of Numpy arrays.
    __array_priority__ = 20000

    # Declare that Time can be used as a Table column by defining the
    # attribute where column attributes will be stored.
    _astropy_column_attrs = None

    def __getnewargs__(self):
        return (self._time,)

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
            if not (isinstance(scale, str)
                    and scale.lower() in self.SCALES):
                raise ScaleValueError("Scale {!r} is not in the allowed scales "
                                      "{}".format(scale,
                                                  sorted(self.SCALES)))

        # If either of the input val, val2 are masked arrays then
        # find the masked elements and fill them.
        mask, val, val2 = _check_for_masked_and_fill(val, val2)

        # Parse / convert input values into internal jd1, jd2 based on format
        self._time = self._get_time_fmt(val, val2, format, scale,
                                        precision, in_subfmt, out_subfmt)
        self._format = self._time.name

        # Hack from #9969 to allow passing the location value that has been
        # collected by the TimeAstropyTime format class up to the Time level.
        # TODO: find a nicer way.
        if hasattr(self._time, '_location'):
            self.location = self._time._location
            del self._time._location

        # If any inputs were masked then masked jd2 accordingly.  From above
        # routine ``mask`` must be either Python bool False or an bool ndarray
        # with shape broadcastable to jd2.
        if mask is not False:
            mask = np.broadcast_to(mask, self._time.jd2.shape)
            self._time.jd1[mask] = 2451544.5  # Set to JD for 2000-01-01
            self._time.jd2[mask] = np.nan

    def _get_time_fmt(self, val, val2, format, scale,
                      precision, in_subfmt, out_subfmt):
        """
        Given the supplied val, val2, format and scale try to instantiate
        the corresponding TimeFormat class to convert the input values into
        the internal jd1 and jd2.

        If format is `None` and the input is a string-type or object array then
        guess available formats and stop when one matches.
        """

        if (format is None
                and (val.dtype.kind in ('S', 'U', 'O', 'M') or val.dtype.names)):
            # Input is a string, object, datetime, or a table-like ndarray
            # (structured array, recarray). These input types can be
            # uniquely identified by the format classes.
            formats = [(name, cls) for name, cls in self.FORMATS.items()
                       if issubclass(cls, TimeUnique)]

            # AstropyTime is a pseudo-format that isn't in the TIME_FORMATS registry,
            # but try to guess it at the end.
            formats.append(('astropy_time', TimeAstropyTime))

        elif not (isinstance(format, str)
                  and format.lower() in self.FORMATS):
            if format is None:
                raise ValueError("No time format was given, and the input is "
                                 "not unique")
            else:
                raise ValueError("Format {!r} is not one of the allowed "
                                 "formats {}".format(format,
                                                     sorted(self.FORMATS)))
        else:
            formats = [(format, self.FORMATS[format])]

        assert formats
        problems = {}
        for name, cls in formats:
            try:
                return cls(val, val2, scale, precision, in_subfmt, out_subfmt)
            except UnitConversionError:
                raise
            except (ValueError, TypeError) as err:
                # If ``format`` specified then there is only one possibility, so raise
                # immediately and include the upstream exception message to make it
                # easier for user to see what is wrong.
                if len(formats) == 1:
                    raise ValueError(
                        f'Input values did not match the format class {format}:'
                        + os.linesep
                        + f'{err.__class__.__name__}: {err}'
                    ) from err
                else:
                    problems[name] = err
        else:
            raise ValueError(f'Input values did not match any of the formats '
                             f'where the format keyword is optional: '
                             f'{problems}') from problems[formats[0][0]]

    @property
    def writeable(self):
        return self._time.jd1.flags.writeable & self._time.jd2.flags.writeable

    @writeable.setter
    def writeable(self, value):
        self._time.jd1.flags.writeable = value
        self._time.jd2.flags.writeable = value

    @property
    def format(self):
        """
        Get or set time format.

        The format defines the way times are represented when accessed via the
        ``.value`` attribute.  By default it is the same as the format used for
        initializing the `Time` instance, but it can be set to any other value
        that could be used for initialization.  These can be listed with::

          >>> list(Time.FORMATS)
          ['jd', 'mjd', 'decimalyear', 'unix', 'unix_tai', 'cxcsec', 'gps', 'plot_date',
           'stardate', 'datetime', 'ymdhms', 'iso', 'isot', 'yday', 'datetime64',
           'fits', 'byear', 'jyear', 'byear_str', 'jyear_str']
        """
        return self._format

    @format.setter
    def format(self, format):
        """Set time format"""
        if format not in self.FORMATS:
            raise ValueError('format must be one of {}'
                             .format(list(self.FORMATS)))
        format_cls = self.FORMATS[format]

        # Get the new TimeFormat object to contain time in new format.  Possibly
        # coerce in/out_subfmt to '*' (default) if existing subfmt values are
        # not valid in the new format.
        self._time = format_cls(
            self._time.jd1, self._time.jd2,
            self._time._scale, self.precision,
            in_subfmt=format_cls._get_allowed_subfmt(self.in_subfmt),
            out_subfmt=format_cls._get_allowed_subfmt(self.out_subfmt),
            from_jd=True)

        self._format = format

    def __repr__(self):
        return ("<{} object: scale='{}' format='{}' value={}>"
                .format(self.__class__.__name__, self.scale, self.format,
                        getattr(self, self.format)))

    def __str__(self):
        return str(getattr(self, self.format))

    def __hash__(self):

        try:
            loc = getattr(self, 'location', None)
            if loc is not None:
                loc = loc.x.to_value(u.m), loc.y.to_value(u.m), loc.z.to_value(u.m)

            return hash((self.jd1, self.jd2, self.scale, loc))

        except TypeError:
            if self.ndim != 0:
                reason = '(must be scalar)'
            elif self.masked:
                reason = '(value is masked)'
            else:
                raise

            raise TypeError("unhashable type: '{}' {}"
                            .format(self.__class__.__name__, reason))

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
            raise ValueError("Scale {!r} is not in the allowed scales {}"
                             .format(scale, sorted(self.SCALES)))

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
        jd1, jd2 = self._time.jd1, self._time.jd2_filled
        for sys1, sys2 in zip(xforms[:-1], xforms[1:]):
            # Some xforms require an additional delta_ argument that is
            # provided through Time methods.  These values may be supplied by
            # the user or computed based on available approximations.  The
            # get_delta_ methods are available for only one combination of
            # sys1, sys2 though the property applies for both xform directions.
            args = [jd1, jd2]
            for sys12 in ((sys1, sys2), (sys2, sys1)):
                dt_method = '_get_delta_{}_{}'.format(*sys12)
                try:
                    get_dt = getattr(self, dt_method)
                except AttributeError:
                    pass
                else:
                    args.append(get_dt(jd1, jd2))
                    break

            conv_func = getattr(erfa, sys1 + sys2)
            jd1, jd2 = conv_func(*args)

        jd1, jd2 = day_frac(jd1, jd2)
        if self.masked:
            jd2[self.mask] = np.nan

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
        del self.cache
        if not isinstance(val, int) or val < 0 or val > 9:
            raise ValueError('precision attribute must be an int between '
                             '0 and 9')
        self._time.precision = val

    @property
    def in_subfmt(self):
        """
        Unix wildcard pattern to select subformats for parsing string input
        times.
        """
        return self._time.in_subfmt

    @in_subfmt.setter
    def in_subfmt(self, val):
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
        # Setting the out_subfmt property here does validation of ``val``
        self._time.out_subfmt = val
        del self.cache

    @property
    def shape(self):
        """The shape of the time instances.

        Like `~numpy.ndarray.shape`, can be set to a new shape by assigning a
        tuple.  Note that if different instances share some but not all
        underlying data, setting the shape of one instance can make the other
        instance unusable.  Hence, it is strongly recommended to get new,
        reshaped instances with the ``reshape`` method.

        Raises
        ------
        ValueError
            If the new shape has the wrong total number of elements.
        AttributeError
            If the shape of the ``jd1``, ``jd2``, ``location``,
            ``delta_ut1_utc``, or ``delta_tdb_tt`` attributes cannot be changed
            without the arrays being copied.  For these cases, use the
            `Time.reshape` method (which copies any arrays that cannot be
            reshaped in-place).
        """
        return self._time.jd1.shape

    @shape.setter
    def shape(self, shape):
        del self.cache

        # We have to keep track of arrays that were already reshaped,
        # since we may have to return those to their original shape if a later
        # shape-setting fails.
        reshaped = []
        oldshape = self.shape

        # In-place reshape of data/attributes.  Need to access _time.jd1/2 not
        # self.jd1/2 because the latter are not guaranteed to be the actual
        # data, and in fact should not be directly changeable from the public
        # API.
        for obj, attr in ((self._time, 'jd1'),
                          (self._time, 'jd2'),
                          (self, '_delta_ut1_utc'),
                          (self, '_delta_tdb_tt'),
                          (self, 'location')):
            val = getattr(obj, attr, None)
            if val is not None and val.size > 1:
                try:
                    val.shape = shape
                except Exception:
                    for val2 in reshaped:
                        val2.shape = oldshape
                    raise
                else:
                    reshaped.append(val)

    def _shaped_like_input(self, value):
        if self._time.jd1.shape:
            if isinstance(value, np.ndarray):
                return value
            else:
                raise TypeError(
                    f"JD is an array ({self._time.jd1!r}) but value "
                    f"is not ({value!r})")
        else:
            # zero-dimensional array, is it safe to unbox?
            if (isinstance(value, np.ndarray)
                    and not value.shape
                    and not np.ma.is_masked(value)):
                if value.dtype.kind == 'M':
                    # existing test doesn't want datetime64 converted
                    return value[()]
                elif value.dtype.fields:
                    # Unpack but keep field names; .item() doesn't
                    # Still don't get python types in the fields
                    return value[()]
                else:
                    return value.item()
            else:
                return value

    @property
    def jd1(self):
        """
        First of the two doubles that internally store time value(s) in JD.
        """
        jd1 = self._time.mask_if_needed(self._time.jd1)
        return self._shaped_like_input(jd1)

    @property
    def jd2(self):
        """
        Second of the two doubles that internally store time value(s) in JD.
        """
        jd2 = self._time.mask_if_needed(self._time.jd2)
        return self._shaped_like_input(jd2)

    def to_value(self, format, subfmt='*'):
        """Get time values expressed in specified output format.

        This method allows representing the ``Time`` object in the desired
        output ``format`` and optional sub-format ``subfmt``.  Available
        built-in formats include ``jd``, ``mjd``, ``iso``, and so forth. Each
        format can have its own sub-formats

        For built-in numerical formats like ``jd`` or ``unix``, ``subfmt`` can
        be one of 'float', 'long', 'decimal', 'str', or 'bytes'.  Here, 'long'
        uses ``numpy.longdouble`` for somewhat enhanced precision (with
        the enhancement depending on platform), and 'decimal'
        :class:`decimal.Decimal` for full precision.  For 'str' and 'bytes', the
        number of digits is also chosen such that time values are represented
        accurately.

        For built-in date-like string formats, one of 'date_hms', 'date_hm', or
        'date' (or 'longdate_hms', etc., for 5-digit years in
        `~astropy.time.TimeFITS`).  For sub-formats including seconds, the
        number of digits used for the fractional seconds is as set by
        `~astropy.time.Time.precision`.

        Parameters
        ----------
        format : str
            The format in which one wants the time values. Default: the current
            format.
        subfmt : str or `None`, optional
            Value or wildcard pattern to select the sub-format in which the
            values should be given.  The default of '*' picks the first
            available for a given format, i.e., 'float' or 'date_hms'.
            If `None`, use the instance's ``out_subfmt``.

        """
        # TODO: add a precision argument (but ensure it is keyword argument
        # only, to make life easier for TimeDelta.to_value()).
        if format not in self.FORMATS:
            raise ValueError('format must be one of {}'
                             .format(list(self.FORMATS)))

        cache = self.cache['format']
        # Try to keep cache behaviour like it was in astropy < 4.0.
        key = format if subfmt is None else (format, subfmt)
        if key not in cache:
            if format == self.format:
                tm = self
            else:
                tm = self.replicate(format=format)

            # Some TimeFormat subclasses may not be able to handle being passes
            # on a out_subfmt. This includes some core classes like
            # TimeBesselianEpochString that do not have any allowed subfmts. But
            # those do deal with `self.out_subfmt` internally, so if subfmt is
            # the same, we do not pass it on.
            kwargs = {}
            if subfmt is not None and subfmt != tm.out_subfmt:
                kwargs['out_subfmt'] = subfmt
            try:
                value = tm._time.to_value(parent=tm, **kwargs)
            except TypeError as exc:
                # Try validating subfmt, e.g. for formats like 'jyear_str' that
                # do not implement out_subfmt in to_value() (because there are
                # no allowed subformats).  If subfmt is not valid this gives the
                # same exception as would have occurred if the call to
                # `to_value()` had succeeded.
                tm._time._select_subfmts(subfmt)

                # Subfmt was valid, so fall back to the original exception to see
                # if it was lack of support for out_subfmt as a call arg.
                if "unexpected keyword argument 'out_subfmt'" in str(exc):
                    raise ValueError(
                        f"to_value() method for format {format!r} does not "
                        f"support passing a 'subfmt' argument") from None
                else:
                    # Some unforeseen exception so raise.
                    raise

            value = tm._shaped_like_input(value)
            cache[key] = value
        return cache[key]

    @property
    def value(self):
        """Time value(s) in current format"""
        return self.to_value(self.format, None)

    @property
    def masked(self):
        return self._time.masked

    @property
    def mask(self):
        return self._time.mask

    def insert(self, obj, values, axis=0):
        """
        Insert values before the given indices in the column and return
        a new `~astropy.time.Time` or  `~astropy.time.TimeDelta` object.

        The values to be inserted must conform to the rules for in-place setting
        of ``Time`` objects (see ``Get and set values`` in the ``Time``
        documentation).

        The API signature matches the ``np.insert`` API, but is more limited.
        The specification of insert index ``obj`` must be a single integer,
        and the ``axis`` must be ``0`` for simple row insertion before the
        index.

        Parameters
        ----------
        obj : int
            Integer index before which ``values`` is inserted.
        values : array_like
            Value(s) to insert.  If the type of ``values`` is different
            from that of quantity, ``values`` is converted to the matching type.
        axis : int, optional
            Axis along which to insert ``values``.  Default is 0, which is the
            only allowed value and will insert a row.

        Returns
        -------
        out : `~astropy.time.Time` subclass
            New time object with inserted value(s)

        """
        # Validate inputs: obj arg is integer, axis=0, self is not a scalar, and
        # input index is in bounds.
        try:
            idx0 = operator.index(obj)
        except TypeError:
            raise TypeError('obj arg must be an integer')

        if axis != 0:
            raise ValueError('axis must be 0')

        if not self.shape:
            raise TypeError('cannot insert into scalar {} object'
                            .format(self.__class__.__name__))

        if abs(idx0) > len(self):
            raise IndexError('index {} is out of bounds for axis 0 with size {}'
                             .format(idx0, len(self)))

        # Turn negative index into positive
        if idx0 < 0:
            idx0 = len(self) + idx0

        # For non-Time object, use numpy to help figure out the length.  (Note annoying
        # case of a string input that has a length which is not the length we want).
        if not isinstance(values, self.__class__):
            values = np.asarray(values)
        n_values = len(values) if values.shape else 1

        # Finally make the new object with the correct length and set values for the
        # three sections, before insert, the insert, and after the insert.
        out = self.__class__.info.new_like([self], len(self) + n_values, name=self.info.name)

        out._time.jd1[:idx0] = self._time.jd1[:idx0]
        out._time.jd2[:idx0] = self._time.jd2[:idx0]

        # This uses the Time setting machinery to coerce and validate as necessary.
        out[idx0:idx0 + n_values] = values

        out._time.jd1[idx0 + n_values:] = self._time.jd1[idx0:]
        out._time.jd2[idx0 + n_values:] = self._time.jd2[idx0:]

        return out

    def __setitem__(self, item, value):
        if not self.writeable:
            if self.shape:
                raise ValueError('{} object is read-only. Make a '
                                 'copy() or set "writeable" attribute to True.'
                                 .format(self.__class__.__name__))
            else:
                raise ValueError('scalar {} object is read-only.'
                                 .format(self.__class__.__name__))

        # Any use of setitem results in immediate cache invalidation
        del self.cache

        # Setting invalidates transform deltas
        for attr in ('_delta_tdb_tt', '_delta_ut1_utc'):
            if hasattr(self, attr):
                delattr(self, attr)

        if value is np.ma.masked or value is np.nan:
            self._time.jd2[item] = np.nan
            return

        value = self._make_value_equivalent(item, value)

        # Finally directly set the jd1/2 values.  Locations are known to match.
        if self.scale is not None:
            value = getattr(value, self.scale)
        self._time.jd1[item] = value._time.jd1
        self._time.jd2[item] = value._time.jd2

    def isclose(self, other, atol=None):
        """Returns a boolean or boolean array where two Time objects are
        element-wise equal within a time tolerance.

        This evaluates the expression below::

          abs(self - other) <= atol

        Parameters
        ----------
        other : `~astropy.time.Time`
            Time object for comparison.
        atol : `~astropy.units.Quantity` or `~astropy.time.TimeDelta`
            Absoute tolerance for equality with units of time (e.g. ``u.s`` or
            ``u.day``). Default is two bits in the 128-bit JD time representation,
            equivalent to about 40 picosecs.
        """
        if atol is None:
            # Note: use 2 bits instead of 1 bit based on experience in precision
            # tests, since taking the difference with a UTC time means one has
            # to do a scale change.
            atol = 2 * np.finfo(float).eps * u.day

        if not isinstance(atol, (u.Quantity, TimeDelta)):
            raise TypeError("'atol' argument must be a Quantity or TimeDelta instance, got "
                            f'{atol.__class__.__name__} instead')

        try:
            # Separate these out so user sees where the problem is
            dt = self - other
            dt = abs(dt)
            out = dt <= atol
        except Exception as err:
            raise TypeError("'other' argument must support subtraction with Time "
                            f"and return a value that supports comparison with "
                            f"{atol.__class__.__name__}: {err}")

        return out

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
        return self._apply('copy', format=format)

    def replicate(self, format=None, copy=False, cls=None):
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
        return self._apply('copy' if copy else 'replicate', format=format, cls=cls)

    def _apply(self, method, *args, format=None, cls=None, **kwargs):
        """Create a new time object, possibly applying a method to the arrays.

        Parameters
        ----------
        method : str or callable
            If string, can be 'replicate'  or the name of a relevant
            `~numpy.ndarray` method. In the former case, a new time instance
            with unchanged internal data is created, while in the latter the
            method is applied to the internal ``jd1`` and ``jd2`` arrays, as
            well as to possible ``location``, ``_delta_ut1_utc``, and
            ``_delta_tdb_tt`` arrays.
            If a callable, it is directly applied to the above arrays.
            Examples: 'copy', '__getitem__', 'reshape', `~numpy.broadcast_to`.
        args : tuple
            Any positional arguments for ``method``.
        kwargs : dict
            Any keyword arguments for ``method``.  If the ``format`` keyword
            argument is present, this will be used as the Time format of the
            replica.

        Examples
        --------
        Some ways this is used internally::

            copy : ``_apply('copy')``
            replicate : ``_apply('replicate')``
            reshape : ``_apply('reshape', new_shape)``
            index or slice : ``_apply('__getitem__', item)``
            broadcast : ``_apply(np.broadcast, shape=new_shape)``
        """
        new_format = self.format if format is None else format

        if callable(method):
            apply_method = lambda array: method(array, *args, **kwargs)

        else:
            if method == 'replicate':
                apply_method = None
            else:
                apply_method = operator.methodcaller(method, *args, **kwargs)

        jd1, jd2 = self._time.jd1, self._time.jd2
        if apply_method:
            jd1 = apply_method(jd1)
            jd2 = apply_method(jd2)

        # Get a new instance of our class and set its attributes directly.
        tm = super().__new__(cls or self.__class__)
        tm._time = TimeJD(jd1, jd2, self.scale, precision=0,
                          in_subfmt='*', out_subfmt='*', from_jd=True)

        # Optional ndarray attributes.
        for attr in ('_delta_ut1_utc', '_delta_tdb_tt', 'location'):
            try:
                val = getattr(self, attr)
            except AttributeError:
                continue

            if apply_method:
                # Apply the method to any value arrays (though skip if there is
                # only an array scalar and the method would return a view,
                # since in that case nothing would change).
                if getattr(val, 'shape', ()):
                    val = apply_method(val)
                elif method == 'copy' or method == 'flatten':
                    # flatten should copy also for a single element array, but
                    # we cannot use it directly for array scalars, since it
                    # always returns a one-dimensional array. So, just copy.
                    val = copy.copy(val)

            setattr(tm, attr, val)

        # Copy other 'info' attr only if it has actually been defined and the
        # time object is not a scalar (issue #10688).
        # See PR #3898 for further explanation and justification, along
        # with Quantity.__array_finalize__
        if 'info' in self.__dict__:
            tm.info = self.info

        # Make the new internal _time object corresponding to the format
        # in the copy.  If the format is unchanged this process is lightweight
        # and does not create any new arrays.
        if new_format not in tm.FORMATS:
            raise ValueError('format must be one of {}'
                             .format(list(tm.FORMATS)))

        NewFormat = tm.FORMATS[new_format]

        tm._time = NewFormat(
            tm._time.jd1, tm._time.jd2,
            tm._time._scale,
            precision=self.precision,
            in_subfmt=NewFormat._get_allowed_subfmt(self.in_subfmt),
            out_subfmt=NewFormat._get_allowed_subfmt(self.out_subfmt),
            from_jd=True)
        tm._format = new_format
        tm.SCALES = self.SCALES

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

        index = [indices
                 if i == axis
                 else np.arange(s).reshape(
                     (1,) * (i if keepdims or i < axis else i - 1)
                     + (s,)
                     + (1,) * (ndim - i - (1 if keepdims or i > axis else 2))
                 )
                 for i, s in enumerate(self.shape)]

        return tuple(index)

    def argmin(self, axis=None, out=None):
        """Return indices of the minimum values along the given axis.

        This is similar to :meth:`~numpy.ndarray.argmin`, but adapted to ensure
        that the full precision given by the two doubles ``jd1`` and ``jd2``
        is used.  See :func:`~numpy.argmin` for detailed documentation.
        """
        # first get the minimum at normal precision.
        jd = self.jd1 + self.jd2

        approx = np.min(jd, axis, keepdims=True)

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

        approx = np.max(jd, axis, keepdims=True)

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
        jd_remainder = (self - self.__class__(jd_approx, format='jd', scale=self.scale)).jd
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
        return (self.max(axis, keepdims=keepdims)
                - self.min(axis, keepdims=keepdims))

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

    @property
    def cache(self):
        """
        Return the cache associated with this instance.
        """
        return self._time.cache

    @cache.deleter
    def cache(self):
        del self._time.cache

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
                    if tm.shape:
                        # Prevent future modification of cached array-like object
                        tm.writeable = False
                cache[attr] = tm
            return cache[attr]

        elif attr in self.FORMATS:
            return self.to_value(attr, subfmt=None)

        elif attr in TIME_SCALES:  # allowed ones done above (self.SCALES)
            if self.scale is None:
                raise ScaleValueError("Cannot convert TimeDelta with "
                                      "undefined scale to any defined scale.")
            else:
                raise ScaleValueError("Cannot convert {} with scale "
                                      "'{}' to scale '{}'"
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
                val = np.broadcast_to(val, self.shape, subok=True)
            except Exception:
                raise ValueError('Attribute shape must match or be '
                                 'broadcastable to that of Time object. '
                                 'Typically, give either a single value or '
                                 'one for each time.')

        return val

    def _time_comparison(self, other, op):
        """If other is of same class as self, compare difference in self.scale.
        Otherwise, return NotImplemented
        """
        if other.__class__ is not self.__class__:
            try:
                other = self.__class__(other, scale=self.scale)
            except Exception:
                # Let other have a go.
                return NotImplemented

        if(self.scale is not None and self.scale not in other.SCALES
           or other.scale is not None and other.scale not in self.SCALES):
            # Other will also not be able to do it, so raise a TypeError
            # immediately, allowing us to explain why it doesn't work.
            raise TypeError("Cannot compare {} instances with scales "
                            "'{}' and '{}'".format(self.__class__.__name__,
                                                   self.scale, other.scale))

        if self.scale is not None and other.scale is not None:
            other = getattr(other, self.scale)

        return op((self.jd1 - other.jd1) + (self.jd2 - other.jd2), 0.)

    def __lt__(self, other):
        return self._time_comparison(other, operator.lt)

    def __le__(self, other):
        return self._time_comparison(other, operator.le)

    def __eq__(self, other):
        """
        If other is an incompatible object for comparison, return `False`.
        Otherwise, return `True` if the time difference between self and
        other is zero.
        """
        return self._time_comparison(other, operator.eq)

    def __ne__(self, other):
        """
        If other is an incompatible object for comparison, return `True`.
        Otherwise, return `False` if the time difference between self and
        other is zero.
        """
        return self._time_comparison(other, operator.ne)

    def __gt__(self, other):
        return self._time_comparison(other, operator.gt)

    def __ge__(self, other):
        return self._time_comparison(other, operator.ge)


class Time(TimeBase):
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
      ['jd', 'mjd', 'decimalyear', 'unix', 'unix_tai', 'cxcsec', 'gps', 'plot_date',
       'stardate', 'datetime', 'ymdhms', 'iso', 'isot', 'yday', 'datetime64',
       'fits', 'byear', 'jyear', 'byear_str', 'jyear_str']

    See also: http://docs.astropy.org/en/stable/time/

    Parameters
    ----------
    val : sequence, ndarray, number, str, bytes, or `~astropy.time.Time` object
        Value(s) to initialize the time or times.  Bytes are decoded as ascii.
    val2 : sequence, ndarray, or number; optional
        Value(s) to initialize the time or times.  Only used for numerical
        input, to help preserve precision.
    format : str, optional
        Format of input value(s)
    scale : str, optional
        Time scale of input value(s), must be one of the following:
        ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')
    precision : int, optional
        Digits of precision in string representation of time
    in_subfmt : str, optional
        Unix glob to select subformats for parsing input times
    out_subfmt : str, optional
        Unix glob to select subformat for outputting times
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

    def __new__(cls, val, val2=None, format=None, scale=None,
                precision=None, in_subfmt=None, out_subfmt=None,
                location=None, copy=False):

        # Because of import problems, this can only be done on
        # first call of Time.
        global _LEAP_SECONDS_CHECKED
        if not _LEAP_SECONDS_CHECKED:
            # *Must* set to True first as update_leap_seconds uses Time.
            # In principle, this may cause wrong leap seconds in
            # update_leap_seconds itself, but since expiration is in
            # units of days, that is fine.
            _LEAP_SECONDS_CHECKED = True
            update_leap_seconds()

        if isinstance(val, Time):
            self = val.replicate(format=format, copy=copy, cls=cls)
        else:
            self = super().__new__(cls)

        return self

    def __init__(self, val, val2=None, format=None, scale=None,
                 precision=None, in_subfmt=None, out_subfmt=None,
                 location=None, copy=False):

        if location is not None:
            from astropy.coordinates import EarthLocation
            if isinstance(location, EarthLocation):
                self.location = location
            else:
                self.location = EarthLocation(*location)
            if self.location.size == 1:
                self.location = self.location.squeeze()
        else:
            if not hasattr(self, 'location'):
                self.location = None

        if isinstance(val, Time):
            # Update _time formatting parameters if explicitly specified
            if precision is not None:
                self._time.precision = precision
            if in_subfmt is not None:
                self._time.in_subfmt = in_subfmt
            if out_subfmt is not None:
                self._time.out_subfmt = out_subfmt
            self.SCALES = TIME_TYPES[self.scale]
            if scale is not None:
                self._set_scale(scale)
        else:
            self._init_from_vals(val, val2, format, scale, copy,
                                 precision, in_subfmt, out_subfmt)
            self.SCALES = TIME_TYPES[self.scale]

        if self.location is not None and (self.location.size > 1
                                          and self.location.shape != self.shape):
            try:
                # check the location can be broadcast to self's shape.
                self.location = np.broadcast_to(self.location, self.shape,
                                                subok=True)
            except Exception as err:
                raise ValueError('The location with shape {} cannot be '
                                 'broadcast against time with shape {}. '
                                 'Typically, either give a single location or '
                                 'one for each time.'
                                 .format(self.location.shape, self.shape)) from err

    def _make_value_equivalent(self, item, value):
        """Coerce setitem value into an equivalent Time object"""

        # If there is a vector location then broadcast to the Time shape
        # and then select with ``item``
        if self.location is not None and self.location.shape:
            self_location = np.broadcast_to(self.location, self.shape, subok=True)[item]
        else:
            self_location = self.location

        if isinstance(value, Time):
            # Make sure locations are compatible.  Location can be either None or
            # a Location object.
            if self_location is None and value.location is None:
                match = True
            elif ((self_location is None and value.location is not None)
                  or (self_location is not None and value.location is None)):
                match = False
            else:
                match = np.all(self_location == value.location)
            if not match:
                raise ValueError('cannot set to Time with different location: '
                                 'expected location={} and '
                                 'got location={}'
                                 .format(self_location, value.location))
        else:
            try:
                value = self.__class__(value, scale=self.scale, location=self_location)
            except Exception:
                try:
                    value = self.__class__(value, scale=self.scale, format=self.format,
                                           location=self_location)
                except Exception as err:
                    raise ValueError('cannot convert value to a compatible Time object: {}'
                                     .format(err))
        return value

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

    @classmethod
    def strptime(cls, time_string, format_string, **kwargs):
        """
        Parse a string to a Time according to a format specification.
        See `time.strptime` documentation for format specification.

        >>> Time.strptime('2012-Jun-30 23:59:60', '%Y-%b-%d %H:%M:%S')
        <Time object: scale='utc' format='isot' value=2012-06-30T23:59:60.000>

        Parameters
        ----------
        time_string : str, sequence, or ndarray
            Objects containing time data of type string
        format_string : str
            String specifying format of time_string.
        kwargs : dict
            Any keyword arguments for ``Time``.  If the ``format`` keyword
            argument is present, this will be used as the Time format.

        Returns
        -------
        time_obj : `~astropy.time.Time`
            A new `~astropy.time.Time` object corresponding to the input
            ``time_string``.

        """
        time_array = np.asarray(time_string)

        if time_array.dtype.kind not in ('U', 'S'):
            err = "Expected type is string, a bytes-like object or a sequence"\
                  " of these. Got dtype '{}'".format(time_array.dtype.kind)
            raise TypeError(err)

        to_string = (str if time_array.dtype.kind == 'U' else
                     lambda x: str(x.item(), encoding='ascii'))
        iterator = np.nditer([time_array, None],
                             op_dtypes=[time_array.dtype, 'U30'])

        for time, formatted in iterator:
            tt, fraction = _strptime._strptime(to_string(time), format_string)
            time_tuple = tt[:6] + (fraction,)
            formatted[...] = '{:04}-{:02}-{:02}T{:02}:{:02}:{:02}.{:06}'\
                .format(*time_tuple)

        format = kwargs.pop('format', None)
        out = cls(*iterator.operands[1:], format='isot', **kwargs)
        if format is not None:
            out.format = format

        return out

    def strftime(self, format_spec):
        """
        Convert Time to a string or a numpy.array of strings according to a
        format specification.
        See `time.strftime` documentation for format specification.

        Parameters
        ----------
        format_spec : str
            Format definition of return string.

        Returns
        -------
        formatted : str or numpy.array
            String or numpy.array of strings formatted according to the given
            format string.

        """
        formatted_strings = []
        for sk in self.replicate('iso')._time.str_kwargs():
            date_tuple = date(sk['year'], sk['mon'], sk['day']).timetuple()
            datetime_tuple = (sk['year'], sk['mon'], sk['day'],
                              sk['hour'], sk['min'], sk['sec'],
                              date_tuple[6], date_tuple[7], -1)
            fmtd_str = format_spec
            if '%f' in fmtd_str:
                fmtd_str = fmtd_str.replace('%f', '{frac:0{precision}}'.format(
                    frac=sk['fracsec'], precision=self.precision))
            fmtd_str = strftime(fmtd_str, datetime_tuple)
            formatted_strings.append(fmtd_str)

        if self.isscalar:
            return formatted_strings[0]
        else:
            return np.array(formatted_strings).reshape(self.shape)

    def light_travel_time(self, skycoord, kind='barycentric', location=None, ephemeris=None):
        """Light travel time correction to the barycentre or heliocentre.

        The frame transformations used to calculate the location of the solar
        system barycentre and the heliocentre rely on the erfa routine epv00,
        which is consistent with the JPL DE405 ephemeris to an accuracy of
        11.2 km, corresponding to a light travel time of 4 microseconds.

        The routine assumes the source(s) are at large distance, i.e., neglects
        finite-distance effects.

        Parameters
        ----------
        skycoord : `~astropy.coordinates.SkyCoord`
            The sky location to calculate the correction for.
        kind : str, optional
            ``'barycentric'`` (default) or ``'heliocentric'``
        location : `~astropy.coordinates.EarthLocation`, optional
            The location of the observatory to calculate the correction for.
            If no location is given, the ``location`` attribute of the Time
            object is used
        ephemeris : str, optional
            Solar system ephemeris to use (e.g., 'builtin', 'jpl'). By default,
            use the one set with ``astropy.coordinates.solar_system_ephemeris.set``.
            For more information, see `~astropy.coordinates.solar_system_ephemeris`.

        Returns
        -------
        time_offset : `~astropy.time.TimeDelta`
            The time offset between the barycentre or Heliocentre and Earth,
            in TDB seconds.  Should be added to the original time to get the
            time in the Solar system barycentre or the Heliocentre.
            Also, the time conversion to BJD will then include the relativistic correction as well.
        """

        if kind.lower() not in ('barycentric', 'heliocentric'):
            raise ValueError("'kind' parameter must be one of 'heliocentric' "
                             "or 'barycentric'")

        if location is None:
            if self.location is None:
                raise ValueError('An EarthLocation needs to be set or passed '
                                 'in to calculate bary- or heliocentric '
                                 'corrections')
            location = self.location

        from astropy.coordinates import (UnitSphericalRepresentation, CartesianRepresentation,
                                         HCRS, ICRS, GCRS, solar_system_ephemeris)

        # ensure sky location is ICRS compatible
        if not skycoord.is_transformable_to(ICRS()):
            raise ValueError("Given skycoord is not transformable to the ICRS")

        # get location of observatory in ITRS coordinates at this Time
        try:
            itrs = location.get_itrs(obstime=self)
        except Exception:
            raise ValueError("Supplied location does not have a valid `get_itrs` method")

        with solar_system_ephemeris.set(ephemeris):
            if kind.lower() == 'heliocentric':
                # convert to heliocentric coordinates, aligned with ICRS
                cpos = itrs.transform_to(HCRS(obstime=self)).cartesian.xyz
            else:
                # first we need to convert to GCRS coordinates with the correct
                # obstime, since ICRS coordinates have no frame time
                gcrs_coo = itrs.transform_to(GCRS(obstime=self))
                # convert to barycentric (BCRS) coordinates, aligned with ICRS
                cpos = gcrs_coo.transform_to(ICRS()).cartesian.xyz

        # get unit ICRS vector to star
        spos = (skycoord.icrs.represent_as(UnitSphericalRepresentation).
                represent_as(CartesianRepresentation).xyz)

        # Move X,Y,Z to last dimension, to enable possible broadcasting below.
        cpos = np.rollaxis(cpos, 0, cpos.ndim)
        spos = np.rollaxis(spos, 0, spos.ndim)

        # calculate light travel time correction
        tcor_val = (spos * cpos).sum(axis=-1) / const.c
        return TimeDelta(tcor_val, scale='tdb')

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

        from astropy.coordinates import Longitude

        if kind.lower() not in SIDEREAL_TIME_MODELS.keys():
            raise ValueError('The kind of sidereal time has to be {}'.format(
                ' or '.join(sorted(SIDEREAL_TIME_MODELS.keys()))))

        available_models = SIDEREAL_TIME_MODELS[kind.lower()]

        if model is None:
            model = sorted(available_models.keys())[-1]
        else:
            if model.upper() not in available_models:
                raise ValueError(
                    'Model {} not implemented for {} sidereal time; '
                    'available models are {}'
                    .format(model, kind, sorted(available_models.keys())))

        if longitude is None:
            if self.location is None:
                raise ValueError('No longitude is given but the location for '
                                 'the Time object is not set.')
            longitude = self.location.lon
        elif longitude == 'greenwich':
            longitude = Longitude(0., u.degree,
                                  wrap_angle=180. * u.degree)
        else:
            # sanity check on input
            longitude = Longitude(longitude, u.degree,
                                  wrap_angle=180. * u.degree)

        gst = self._erfa_sidereal_time(available_models[model.upper()])
        return Longitude(gst + longitude, u.hourangle)

    if isinstance(sidereal_time.__doc__, str):
        sidereal_time.__doc__ = sidereal_time.__doc__.format(
            'apparent', sorted(SIDEREAL_TIME_MODELS['apparent'].keys()),
            'mean', sorted(SIDEREAL_TIME_MODELS['mean'].keys()))

    def _erfa_sidereal_time(self, model):
        """Calculate a sidereal time using a IAU precession/nutation model."""

        from astropy.coordinates import Longitude

        erfa_function = model['function']
        erfa_parameters = [getattr(getattr(self, scale)._time, jd_part)
                           for scale in model['scales']
                           for jd_part in ('jd1', 'jd2_filled')]

        sidereal_time = erfa_function(*erfa_parameters)

        if self.masked:
            sidereal_time[self.mask] = np.nan

        return Longitude(sidereal_time, u.radian).to(u.hourangle)

    def get_delta_ut1_utc(self, iers_table=None, return_status=False):
        """Find UT1 - UTC differences by interpolating in IERS Table.

        Parameters
        ----------
        iers_table : `~astropy.utils.iers.IERS` table, optional
            Table containing UT1-UTC differences from IERS Bulletins A
            and/or B.  Default: `~astropy.utils.iers.earth_orientation_table`
            (which in turn defaults to the combined version provided by
            `~astropy.utils.iers.IERS_Auto`).
        return_status : bool
            Whether to return status values.  If `False` (default), iers
            raises `IndexError` if any time is out of the range
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
            >>> delta, status = t.get_delta_ut1_utc(return_status=True)  # doctest: +REMOTE_DATA
            >>> status == TIME_BEFORE_IERS_RANGE  # doctest: +REMOTE_DATA
            array([ True, False]...)
        """
        if iers_table is None:
            from astropy.utils.iers import earth_orientation_table
            iers_table = earth_orientation_table.get()

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
            from astropy.utils.iers import earth_orientation_table
            iers_table = earth_orientation_table.get()
            # jd1, jd2 are normally set (see above), except if delta_ut1_utc
            # is access directly; ensure we behave as expected for that case
            if jd1 is None:
                self_utc = self.utc
                jd1, jd2 = self_utc._time.jd1, self_utc._time.jd2_filled
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
                jd1_utc, jd2_utc = erfa.ut1utc(jd1, jd2, delta.to_value(u.s))
                # calculate a better estimate using the nearly correct UTC
                delta = iers_table.ut1_utc(jd1_utc, jd2_utc)

            self._set_delta_ut1_utc(delta)

        return self._delta_ut1_utc

    def _set_delta_ut1_utc(self, val):
        del self.cache
        if hasattr(val, 'to'):  # Matches Quantity but also TimeDelta.
            val = val.to(u.second).value
        val = self._match_shape(val)
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
                    jd2 = self._time.jd2_filled

            # First go from the current input time (which is either
            # TDB or TT) to an approximate UT1.  Since TT and TDB are
            # pretty close (few msec?), assume TT.  Similarly, since the
            # UT1 terms are very small, use UTC instead of UT1.
            njd1, njd2 = erfa.tttai(jd1, jd2)
            njd1, njd2 = erfa.taiutc(njd1, njd2)
            # subtract 0.5, so UT is fraction of the day from midnight
            ut = day_frac(njd1 - 0.5, njd2)[1]

            if self.location is None:
                from astropy.coordinates import EarthLocation
                location = EarthLocation.from_geodetic(0., 0., 0.)
            else:
                location = self.location
            # Geodetic params needed for d_tdb_tt()
            lon = location.lon
            rxy = np.hypot(location.x, location.y)
            z = location.z
            self._delta_tdb_tt = erfa.dtdb(
                jd1, jd2, ut, lon.to_value(u.radian),
                rxy.to_value(u.km), z.to_value(u.km))

        return self._delta_tdb_tt

    def _set_delta_tdb_tt(self, val):
        del self.cache
        if hasattr(val, 'to'):  # Matches Quantity but also TimeDelta.
            val = val.to(u.second).value
        val = self._match_shape(val)
        self._delta_tdb_tt = val

    # Note can't use @property because _get_delta_tdb_tt is explicitly
    # called with the optional jd1 and jd2 args.
    delta_tdb_tt = property(_get_delta_tdb_tt, _set_delta_tdb_tt)
    """TDB - TT time scale offset"""

    def __sub__(self, other):
        # T      - Tdelta = T
        # T      - T      = Tdelta
        other_is_delta = not isinstance(other, Time)
        if other_is_delta:  # T - Tdelta
            # Check other is really a TimeDelta or something that can initialize.
            if not isinstance(other, TimeDelta):
                try:
                    other = TimeDelta(other)
                except Exception:
                    return NotImplemented

            # we need a constant scale to calculate, which is guaranteed for
            # TimeDelta, but not for Time (which can be UTC)
            out = self.replicate()
            if self.scale in other.SCALES:
                if other.scale not in (out.scale, None):
                    other = getattr(other, out.scale)
            else:
                if other.scale is None:
                    out._set_scale('tai')
                else:
                    if self.scale not in TIME_TYPES[other.scale]:
                        raise TypeError("Cannot subtract Time and TimeDelta instances "
                                        "with scales '{}' and '{}'"
                                        .format(self.scale, other.scale))
                    out._set_scale(other.scale)
            # remove attributes that are invalidated by changing time
            for attr in ('_delta_ut1_utc', '_delta_tdb_tt'):
                if hasattr(out, attr):
                    delattr(out, attr)

        else:  # T - T
            # the scales should be compatible (e.g., cannot convert TDB to LOCAL)
            if other.scale not in self.SCALES:
                raise TypeError("Cannot subtract Time instances "
                                "with scales '{}' and '{}'"
                                .format(self.scale, other.scale))
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
        # T      + Tdelta = T
        # T      + T      = error
        if isinstance(other, Time):
            raise OperandTypeError(self, other, '+')

        # Check other is really a TimeDelta or something that can initialize.
        if not isinstance(other, TimeDelta):
            try:
                other = TimeDelta(other)
            except Exception:
                return NotImplemented

        # ideally, we calculate in the scale of the Time item, since that is
        # what we want the output in, but this may not be possible, since
        # TimeDelta cannot be converted arbitrarily
        out = self.replicate()
        if self.scale in other.SCALES:
            if other.scale not in (out.scale, None):
                other = getattr(other, out.scale)
        else:
            if other.scale is None:
                out._set_scale('tai')
            else:
                if self.scale not in TIME_TYPES[other.scale]:
                    raise TypeError("Cannot add Time and TimeDelta instances "
                                    "with scales '{}' and '{}'"
                                    .format(self.scale, other.scale))
                out._set_scale(other.scale)
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

    # Reverse addition is possible: <something-Tdelta-ish> + T
    # but there is no case of <something> - T, so no __rsub__.
    def __radd__(self, other):
        return self.__add__(other)

    def to_datetime(self, timezone=None):
        # TODO: this could likely go through to_value, as long as that
        # had an **kwargs part that was just passed on to _time.
        tm = self.replicate(format='datetime')
        return tm._shaped_like_input(tm._time.to_value(timezone))

    to_datetime.__doc__ = TimeDatetime.to_value.__doc__


class TimeDelta(TimeBase):
    """
    Represent the time difference between two times.

    A TimeDelta object is initialized with one or more times in the ``val``
    argument.  The input times in ``val`` must conform to the specified
    ``format``.  The optional ``val2`` time input should be supplied only for
    numeric input formats (e.g. JD) where very high precision (better than
    64-bit precision) is required.

    The allowed values for ``format`` can be listed with::

      >>> list(TimeDelta.FORMATS)
      ['sec', 'jd', 'datetime']

    Note that for time differences, the scale can be among three groups:
    geocentric ('tai', 'tt', 'tcg'), barycentric ('tcb', 'tdb'), and rotational
    ('ut1'). Within each of these, the scales for time differences are the
    same. Conversion between geocentric and barycentric is possible, as there
    is only a scale factor change, but one cannot convert to or from 'ut1', as
    this requires knowledge of the actual times, not just their difference. For
    a similar reason, 'utc' is not a valid scale for a time difference: a UTC
    day is not always 86400 seconds.

    See also:

    - https://docs.astropy.org/en/stable/time/
    - https://docs.astropy.org/en/stable/time/index.html#time-deltas

    Parameters
    ----------
    val : sequence, ndarray, number, `~astropy.units.Quantity` or `~astropy.time.TimeDelta` object
        Value(s) to initialize the time difference(s). Any quantities will
        be converted appropriately (with care taken to avoid rounding
        errors for regular time units).
    val2 : sequence, ndarray, number, or `~astropy.units.Quantity`; optional
        Additional values, as needed to preserve precision.
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

    info = TimeDeltaInfo()

    def __new__(cls, val, val2=None, format=None, scale=None,
                precision=None, in_subfmt=None, out_subfmt=None,
                location=None, copy=False):

        if isinstance(val, TimeDelta):
            self = val.replicate(format=format, copy=copy, cls=cls)
        else:
            self = super().__new__(cls)

        return self

    def __init__(self, val, val2=None, format=None, scale=None, copy=False):
        if isinstance(val, TimeDelta):
            if scale is not None:
                self._set_scale(scale)
        else:
            if format is None:
                format = 'datetime' if isinstance(val, timedelta) else 'jd'

            self._init_from_vals(val, val2, format, scale, copy)

            if scale is not None:
                self.SCALES = TIME_DELTA_TYPES[scale]

    def replicate(self, *args, **kwargs):
        out = super().replicate(*args, **kwargs)
        out.SCALES = self.SCALES
        return out

    def to_datetime(self):
        """
        Convert to ``datetime.timedelta`` object.
        """
        tm = self.replicate(format='datetime')
        return tm._shaped_like_input(tm._time.value)

    def _set_scale(self, scale):
        """
        This is the key routine that actually does time scale conversions.
        This is not public and not connected to the read-only scale property.
        """

        if scale == self.scale:
            return
        if scale not in self.SCALES:
            raise ValueError("Scale {!r} is not in the allowed scales {}"
                             .format(scale, sorted(self.SCALES)))

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

    def _add_sub(self, other, op):
        """Perform common elements of addition / subtraction for two delta times"""
        # If not a TimeDelta then see if it can be turned into a TimeDelta.
        if not isinstance(other, TimeDelta):
            try:
                other = TimeDelta(other)
            except Exception:
                return NotImplemented

        # the scales should be compatible (e.g., cannot convert TDB to TAI)
        if(self.scale is not None and self.scale not in other.SCALES
           or other.scale is not None and other.scale not in self.SCALES):
            raise TypeError("Cannot add TimeDelta instances with scales "
                            "'{}' and '{}'".format(self.scale, other.scale))

        # adjust the scale of other if the scale of self is set (or no scales)
        if self.scale is not None or other.scale is None:
            out = self.replicate()
            if other.scale is not None:
                other = getattr(other, self.scale)
        else:
            out = other.replicate()

        jd1 = op(self._time.jd1, other._time.jd1)
        jd2 = op(self._time.jd2, other._time.jd2)

        out._time.jd1, out._time.jd2 = day_frac(jd1, jd2)

        return out

    def __add__(self, other):
        # If other is a Time then use Time.__add__ to do the calculation.
        if isinstance(other, Time):
            return other.__add__(self)

        return self._add_sub(other, operator.add)

    def __sub__(self, other):
        # TimeDelta - Time is an error
        if isinstance(other, Time):
            raise OperandTypeError(self, other, '-')

        return self._add_sub(other, operator.sub)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        out = self.__sub__(other)
        return -out

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
        # Check needed since otherwise the self.jd1 * other multiplication
        # would enter here again (via __rmul__)
        if isinstance(other, Time):
            raise OperandTypeError(self, other, '*')
        elif ((isinstance(other, u.UnitBase)
               and other == u.dimensionless_unscaled)
                or (isinstance(other, str) and other == '')):
            return self.copy()

        # If other is something consistent with a dimensionless quantity
        # (could just be a float or an array), then we can just multiple in.
        try:
            other = u.Quantity(other, u.dimensionless_unscaled, copy=False)
        except Exception:
            # If not consistent with a dimensionless quantity, try downgrading
            # self to a quantity and see if things work.
            try:
                return self.to(u.day) * other
            except Exception:
                # The various ways we could multiply all failed;
                # returning NotImplemented to give other a final chance.
                return NotImplemented

        jd1, jd2 = day_frac(self.jd1, self.jd2, factor=other.value)
        out = TimeDelta(jd1, jd2, format='jd', scale=self.scale)

        if self.format != 'jd':
            out = out.replicate(format=self.format)
        return out

    def __rmul__(self, other):
        """Multiplication of numbers/arrays with `TimeDelta` objects."""
        return self.__mul__(other)

    def __truediv__(self, other):
        """Division of `TimeDelta` objects by numbers/arrays."""
        # Cannot do __mul__(1./other) as that looses precision
        if ((isinstance(other, u.UnitBase)
             and other == u.dimensionless_unscaled)
                or (isinstance(other, str) and other == '')):
            return self.copy()

        # If other is something consistent with a dimensionless quantity
        # (could just be a float or an array), then we can just divide in.
        try:
            other = u.Quantity(other, u.dimensionless_unscaled, copy=False)
        except Exception:
            # If not consistent with a dimensionless quantity, try downgrading
            # self to a quantity and see if things work.
            try:
                return self.to(u.day) / other
            except Exception:
                # The various ways we could divide all failed;
                # returning NotImplemented to give other a final chance.
                return NotImplemented

        jd1, jd2 = day_frac(self.jd1, self.jd2, divisor=other.value)
        out = TimeDelta(jd1, jd2, format='jd', scale=self.scale)

        if self.format != 'jd':
            out = out.replicate(format=self.format)
        return out

    def __rtruediv__(self, other):
        """Division by `TimeDelta` objects of numbers/arrays."""
        # Here, we do not have to worry about returning NotImplemented,
        # since other has already had a chance to look at us.
        return other / self.to(u.day)

    def to(self, unit, equivalencies=[]):
        """
        Convert to a quantity in the specified unit.

        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            The unit to convert to.
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not directly
            convertible (see :ref:`unit_equivalencies`). If `None`, no
            equivalencies will be applied at all, not even any set globallyq
            or within a context.

        Returns
        -------
        quantity : `~astropy.units.Quantity`
            The quantity in the units specified.

        See also
        --------
        to_value : get the numerical value in a given unit.
        """
        return u.Quantity(self._time.jd1 + self._time.jd2,
                          u.day).to(unit, equivalencies=equivalencies)

    def to_value(self, *args, **kwargs):
        """Get time delta values expressed in specified output format or unit.

        This method is flexible and handles both conversion to a specified
        ``TimeDelta`` format / sub-format AND conversion to a specified unit.
        If positional argument(s) are provided then the first one is checked
        to see if it is a valid ``TimeDelta`` format, and next it is checked
        to see if it is a valid unit or unit string.

        To convert to a ``TimeDelta`` format and optional sub-format the options
        are::

          tm = TimeDelta(1.0 * u.s)
          tm.to_value('jd')  # equivalent of tm.jd
          tm.to_value('jd', 'decimal')  # convert to 'jd' as a Decimal object
          tm.to_value('jd', subfmt='decimal')
          tm.to_value(format='jd', subfmt='decimal')

        To convert to a unit with optional equivalencies, the options are::

          tm.to_value('hr')  # convert to u.hr (hours)
          tm.to_value('hr', [])  # specify equivalencies as a positional arg
          tm.to_value('hr', equivalencies=[])
          tm.to_value(unit='hr', equivalencies=[])

        The built-in `~astropy.time.TimeDelta` options for ``format`` are:
        {'jd', 'sec', 'datetime'}.

        For the two numerical formats 'jd' and 'sec', the available ``subfmt``
        options are: {'float', 'long', 'decimal', 'str', 'bytes'}. Here, 'long'
        uses ``numpy.longdouble`` for somewhat enhanced precision (with the
        enhancement depending on platform), and 'decimal' instances of
        :class:`decimal.Decimal` for full precision.  For the 'str' and 'bytes'
        sub-formats, the number of digits is also chosen such that time values
        are represented accurately.  Default: as set by ``out_subfmt`` (which by
        default picks the first available for a given format, i.e., 'float').

        Parameters
        ----------
        format : str, optional
            The format in which one wants the `~astropy.time.TimeDelta` values.
            Default: the current format.
        subfmt : str, optional
            Possible sub-format in which the values should be given. Default: as
            set by ``out_subfmt`` (which by default picks the first available
            for a given format, i.e., 'float' or 'date_hms').
        unit : `~astropy.units.UnitBase` instance or str, optional
            The unit in which the value should be given.
        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not directly
            convertible (see :ref:`unit_equivalencies`). If `None`, no
            equivalencies will be applied at all, not even any set globally or
            within a context.

        Returns
        -------
        value : `~numpy.ndarray` or scalar
            The value in the format or units specified.

        See also
        --------
        to : Convert to a `~astropy.units.Quantity` instance in a given unit.
        value : The time value in the current format.

        """
        if not (args or kwargs):
            raise TypeError('to_value() missing required format or unit argument')

        # TODO: maybe allow 'subfmt' also for units, keeping full precision
        # (effectively, by doing the reverse of quantity_day_frac)?
        # This way, only equivalencies could lead to possible precision loss.
        if ('format' in kwargs
                or (args != () and (args[0] is None or args[0] in self.FORMATS))):
            # Super-class will error with duplicate arguments, etc.
            return super().to_value(*args, **kwargs)

        # With positional arguments, we try parsing the first one as a unit,
        # so that on failure we can give a more informative exception.
        if args:
            try:
                unit = u.Unit(args[0])
            except ValueError as exc:
                raise ValueError("first argument is not one of the known "
                                 "formats ({}) and failed to parse as a unit."
                                 .format(list(self.FORMATS))) from exc
            args = (unit,) + args[1:]

        return u.Quantity(self._time.jd1 + self._time.jd2,
                          u.day).to_value(*args, **kwargs)

    def _make_value_equivalent(self, item, value):
        """Coerce setitem value into an equivalent TimeDelta object"""
        if not isinstance(value, TimeDelta):
            try:
                value = self.__class__(value, scale=self.scale, format=self.format)
            except Exception as err:
                raise ValueError('cannot convert value to a compatible TimeDelta '
                                 'object: {}'.format(err))
        return value

    def isclose(self, other, atol=None, rtol=0.0):
        """Returns a boolean or boolean array where two TimeDelta objects are
        element-wise equal within a time tolerance.

        This effectively evaluates the expression below::

          abs(self - other) <= atol + rtol * abs(other)

        Parameters
        ----------
        other : `~astropy.units.Quantity` or `~astropy.time.TimeDelta`
            Quantity or TimeDelta object for comparison.
        atol : `~astropy.units.Quantity` or `~astropy.time.TimeDelta`
            Absolute tolerance for equality with units of time (e.g. ``u.s`` or
            ``u.day``). Default is one bit in the 128-bit JD time representation,
            equivalent to about 20 picosecs.
        rtol : float
            Relative tolerance for equality
        """
        try:
            other_day = other.to_value(u.day)
        except Exception as err:
            raise TypeError(f"'other' argument must support conversion to days: {err}")

        if atol is None:
            atol = np.finfo(float).eps * u.day

        if not isinstance(atol, (u.Quantity, TimeDelta)):
            raise TypeError("'atol' argument must be a Quantity or TimeDelta instance, got "
                            f'{atol.__class__.__name__} instead')

        return np.isclose(self.to_value(u.day), other_day,
                          rtol=rtol, atol=atol.to_value(u.day))


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
    if isinstance(val, (tuple, list)) and len(val) > 0 and isinstance(val[0], Time):
        dtype = object
    else:
        dtype = None

    val = np.array(val, copy=copy, subok=True, dtype=dtype)

    # Allow only float64, string or object arrays as input
    # (object is for datetime, maybe add more specific test later?)
    # This also ensures the right byteorder for float64 (closes #2942).
    if val.dtype.kind == "f" and val.dtype.itemsize >= np.dtype(np.float64).itemsize:
        pass
    elif val.dtype.kind in 'OSUMaV':
        pass
    else:
        val = np.asanyarray(val, dtype=np.float64)

    return val


def _check_for_masked_and_fill(val, val2):
    """
    If ``val`` or ``val2`` are masked arrays then fill them and cast
    to ndarray.

    Returns a mask corresponding to the logical-or of masked elements
    in ``val`` and ``val2``.  If neither is masked then the return ``mask``
    is ``None``.

    If either ``val`` or ``val2`` are masked then they are replaced
    with filled versions of themselves.

    Parameters
    ----------
    val : ndarray or MaskedArray
        Input val
    val2 : ndarray or MaskedArray
        Input val2

    Returns
    -------
    mask, val, val2: ndarray or None
        Mask: (None or bool ndarray), val, val2: ndarray
    """
    def get_as_filled_ndarray(mask, val):
        """
        Fill the given MaskedArray ``val`` from the first non-masked
        element in the array.  This ensures that upstream Time initialization
        will succeed.

        Note that nothing happens if there are no masked elements.
        """
        fill_value = None

        if np.any(val.mask):
            # Final mask is the logical-or of inputs
            mask = mask | val.mask

            # First unmasked element.  If all elements are masked then
            # use fill_value=None from above which will use val.fill_value.
            # As long as the user has set this appropriately then all will
            # be fine.
            val_unmasked = val.compressed()  # 1-d ndarray of unmasked values
            if len(val_unmasked) > 0:
                fill_value = val_unmasked[0]

        # Fill the input ``val``.  If fill_value is None then this just returns
        # an ndarray view of val (no copy).
        val = val.filled(fill_value)

        return mask, val

    mask = False
    if isinstance(val, np.ma.MaskedArray):
        mask, val = get_as_filled_ndarray(mask, val)
    if isinstance(val2, np.ma.MaskedArray):
        mask, val2 = get_as_filled_ndarray(mask, val2)

    return mask, val, val2


class OperandTypeError(TypeError):
    def __init__(self, left, right, op=None):
        op_string = '' if op is None else f' for {op}'
        super().__init__(
            "Unsupported operand type(s){}: "
            "'{}' and '{}'".format(op_string,
                                   left.__class__.__name__,
                                   right.__class__.__name__))


def update_leap_seconds(files=None):
    """If the current ERFA leap second table is out of date, try to update it.

    Uses `astropy.utils.iers.LeapSeconds.auto_open` to try to find an
    up-to-date table.  See that routine for the definition of "out of date".

    In order to make it safe to call this any time, all exceptions are turned
    into warnings,

    Parameters
    ----------
    files : list of path, optional
        List of files/URLs to attempt to open.  By default, uses defined by
        `astropy.utils.iers.LeapSeconds.auto_open`, which includes the table
        used by ERFA itself, so if that is up to date, nothing will happen.

    Returns
    -------
    n_update : int
        Number of items updated.

    """
    try:
        from astropy.utils import iers

        table = iers.LeapSeconds.auto_open(files)
        return erfa.leap_seconds.update(table)

    except Exception as exc:
        warn("leap-second auto-update failed due to the following "
             f"exception: {exc!r}", AstropyWarning)
        return 0
