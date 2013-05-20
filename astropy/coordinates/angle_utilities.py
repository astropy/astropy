# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains utility functions that are for internal use in
astropy.coordinates.angles. Mainly they are conversions from one format
of data to another.
"""
from __future__ import unicode_literals

import math
import os
from warnings import warn

from .errors import *
from ..utils import format_exception
from .. import units as u


class _AngleParser(object):
    """
    Parses the various angle formats including:

       * 01:02:30.43 degrees
       * 1 2 0 hours
       * 1°2′3″
       * 1d2m3s
       * -1h2m3s

    This class should not be used directly.  Use `parse_angle`
    instead.
    """
    def __init__(self):
        if '_parser' not in _AngleParser.__dict__:
            _AngleParser._parser, _AngleParser._lexer = self._make_parser()

    @classmethod
    def _get_simple_unit_names(cls):
        simple_units = set(u.radian.find_equivalent_units())
        simple_units.remove(u.deg)
        simple_units.remove(u.hourangle)
        simple_unit_names = set()
        for unit in simple_units:
            simple_unit_names.update(unit.names)
        return list(simple_unit_names)

    @classmethod
    def _make_parser(cls):
        from ..extern.ply import lex, yacc

        # List of token names.
        tokens = (
            'SIGN',
            'UINT',
            'UFLOAT',
            'COLON',
            'DEGREE',
            'HOUR',
            'MINUTE',
            'SECOND',
            'SIMPLE_UNIT'
        )

        # NOTE THE ORDERING OF THESE RULES IS IMPORTANT!!
        # Regular expression rules for simple tokens
        def t_UFLOAT(t):
            r'((\d+\.\d*)|(\.\d+))([eE][+-]?\d+)?'
            t.value = float(t.value)
            return t

        def t_UINT(t):
            r'\d+'
            t.value = int(t.value)
            return t

        def t_SIGN(t):
            r'[+-]'
            t.value = float(t.value + '1')
            return t

        t_COLON = ':'
        t_DEGREE = r'd(eg(ree(s)?)?)?|°'
        t_HOUR = r'hour(s)?|h(r)?|ʰ'
        t_MINUTE = r'm(in(ute(s)?)?)?|′|\''
        t_SECOND = r's(ec(ond(s)?)?)?|″|\"'
        t_SIMPLE_UNIT = '|'.join(
            '({0})'.format(x) for x in cls._get_simple_unit_names())

        # A string containing ignored characters (spaces)
        t_ignore  = ' '

        # Error handling rule
        def t_error(t):
            raise ValueError(
                "Invalid character at col {0}".format(t.lexpos))

        # Build the lexer
        try:
            from . import angle_lextab
            lexer = lex.lex(optimize=True, lextab=angle_lextab)
        except ImportError:
            lexer = lex.lex(optimize=True, lextab='angle_lextab',
                            outputdir=os.path.dirname(__file__))

        def p_angle(p):
            '''
            angle : hms
                  | dms
                  | simple
            '''
            p[0] = p[1]

        def p_sign(p):
            '''
            sign : SIGN
                 |
            '''
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1.0

        def p_ufloat(p):
            '''
            ufloat : UFLOAT
                   | UINT
            '''
            p[0] = float(p[1])

        def p_colon(p):
            '''
            colon : sign UINT COLON UINT
                  | sign UINT COLON UINT COLON ufloat
            '''
            if len(p) == 5:
                p[0] = (p[1] * p[2], p[4], 0.0)
            elif len(p) == 7:
                p[0] = (p[1] * p[2], p[4], p[6])

        def p_spaced(p):
            '''
            spaced : sign UINT UINT
                   | sign UINT UINT ufloat
            '''
            if len(p) == 4:
                p[0] = (p[1] * p[2], p[3], 0.0)
            elif len(p) == 5:
                p[0] = (p[1] * p[2], p[3], p[4])

        def p_generic(p):
            '''
            generic : colon
                    | spaced
                    | sign UFLOAT
                    | sign UINT
            '''
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = p[1] * p[2]

        def p_hms(p):
            '''
            hms : sign UINT HOUR UINT
                | sign UINT HOUR UINT MINUTE
                | sign UINT HOUR UINT MINUTE ufloat
                | sign UINT HOUR UINT MINUTE ufloat SECOND
                | generic HOUR
            '''
            if len(p) == 3:
                p[0] = (p[1], u.hourangle)
            elif len(p) in (5, 6):
                p[0] = ((p[1] * p[2], p[4], 0.0), u.hourangle)
            elif len(p) in (7, 8):
                p[0] = ((p[1] * p[2], p[4], p[6]), u.hourangle)

        def p_dms(p):
            '''
            dms : sign UINT DEGREE UINT
                | sign UINT DEGREE UINT MINUTE
                | sign UINT DEGREE UINT MINUTE ufloat
                | sign UINT DEGREE UINT MINUTE ufloat SECOND
                | generic DEGREE
            '''
            if len(p) == 3:
                p[0] = (p[1], u.degree)
            elif len(p) in (5, 6):
                p[0] = ((p[1] * p[2], p[4], 0.0), u.degree)
            elif len(p) in (7, 8):
                p[0] = ((p[1] * p[2], p[4], p[6]), u.degree)

        def p_simple(p):
            '''
            simple : generic
                   | generic SIMPLE_UNIT
            '''
            if len(p) == 2:
                p[0] = (p[1], None)
            else:
                p[0] = (p[1], p[2])

        def p_error(p):
            raise ValueError

        try:
            from . import angle_parsetab
            parser = yacc.yacc(debug=False, tabmodule=angle_parsetab,
                               write_tables=False)
        except ImportError:
            parser = yacc.yacc(debug=False, tabmodule='angle_parsetab',
                               outputdir=os.path.dirname(__file__))

        return parser, lexer

    def parse(self, angle, unit, debug=False):
        try:
            found_angle, found_unit = self._parser.parse(
                angle, lexer=self._lexer, debug=debug)
        except ValueError as e:
            if str(e):
                raise ValueError("{0} in angle {1!r}".format(
                    str(e), angle))
            else:
                raise ValueError(
                    "Syntax error parsing angle {0!r}".format(angle))

        if unit is not None:
            unit = u.Unit(unit)
            if (found_unit is not None and
                found_unit is not unit):
                raise u.UnitsException(
                    "Unit in string ({0}) does not match requested unit "
                    "({1})".format(found_unit, unit))
        else:
            if found_unit is None:
                raise u.UnitsException("No unit specified")
            else:
                unit = found_unit

        return found_angle, unit


def _check_hour_range(hrs):
    """
    Checks that the given value is in the range (-24, 24).
    """
    if math.fabs(hrs) == 24.:
        warn(IllegalHourWarning(hrs, 'Treating as 24 hr'))
    elif not -24. < hrs < 24.:
        raise IllegalHourError(hrs)


def _check_minute_range(min):
    """
    Checks that the given value is in the range [0,60].  If the value
    is equal to 60, then a warning is raised.
    """
    if min == 60.:
        warn(IllegalMinuteWarning(min, 'Treating as 0 min, +1 hr/deg'))
    elif not 0. <= min < 60.:
        # "Error: minutes not in range [0,60) ({0}).".format(min))
        raise IllegalMinuteError(min)


def _check_second_range(sec):
    """
    Checks that the given value is in the range [0,60].  If the value
    is equal to 60, then a warning is raised.
    """
    if sec == 60.:
        warn(IllegalSecondWarning(sec, 'Treating as 0 sec, +1 min'))
    elif not 0. <= sec < 60.:
        # "Error: seconds not in range [0,60) ({0}).".format(sec))
        raise IllegalSecondError(sec)


def check_hms_ranges(h, m, s):
    """
    Checks that the given hour, minute and second are all within
    reasonable range.
    """
    _check_hour_range(h)
    _check_minute_range(m)
    _check_second_range(s)
    return None


def parse_angle(angle, unit=None, debug=False):
    """
    Parses an input string value into an angle value.

    Parameters
    ----------
    angle : str
        A string representing the angle.  May be in one of the following forms:

            * 01:02:30.43 degrees
            * 1 2 0 hours
            * 1°2′3″
            * 1d2m3s
            * -1h2m3s

    unit : `~astropy.units.UnitBase` instance, optional
        The unit used to interpret the string.  If `unit` is not
        provided, the unit must be explicitly represented in the
        string, either at the end or as number separators.

    debug : bool, optional
        If `True`, print debugging information from the parser.

    Returns
    -------
    value, unit : tuple
        `value` is the value as a floating point number or three-part
        tuple, and `unit` is a `Unit` instance which is either the
        unit passed in or the one explicitly mentioned in the input
        string.
    """
    return _AngleParser().parse(angle, unit, debug=debug)


def degrees_to_dms(d):
    """
    Convert a floating-point degree value into a ``(degree, arcminute,
    arcsecond)`` tuple.
    """
    sign = math.copysign(1.0, d)

    (df, d) = math.modf(abs(d))  # (degree fraction, degree)
    (mf, m) = math.modf(df * 60.)  # (minute fraction, minute)
    s = mf * 60.

    _check_minute_range(m)
    _check_second_range(s)

    return (float(sign * d), int(sign * m), sign * s)


def dms_to_degrees(d, m, s):
    """
    Convert degrees, arcminute, arcsecond to a float degrees value.
    """

    _check_minute_range(m)
    _check_second_range(s)

    # determine sign
    sign = math.copysign(1.0, d)

    # TODO: This will fail if d or m have values after the decimal
    # place

    try:
        d = int(abs(d))
        m = int(abs(m))
        s = float(abs(s))
    except ValueError:
        raise ValueError(format_exception(
            "{func}: dms values ({1[0]},{2[1]},{3[2]}) could not be "
            "converted to numbers.", d, m, s))

    return sign * (d + m / 60. + s / 3600.)


def hms_to_hours(h, m, s):
    """
    Convert hour, minute, second to a float hour value.
    """

    check_hms_ranges(h, m, s)

    # determine sign
    sign = math.copysign(1.0, h)

    # TODO: This will fail if d or m have values after the decimal
    # place

    try:
        h = int(abs(h))
        m = int(abs(m))
        s = float(abs(s))
    except ValueError:
        raise ValueError(format_exception(
            "{func}: HMS values ({1[0]},{2[1]},{3[2]}) could not be "
            "converted to numbers.", h, m, s))

    return sign * (h + m / 60. + s / 3600.)


def hms_to_degrees(h, m, s):
    """
    Convert hour, minute, second to a float degrees value.
    """

    return hms_to_hours(h, m, s) * 15.


def hms_to_radians(h, m, s):
    """
    Convert hour, minute, second to a float radians value.
    """

    return u.degree.to(u.radian, hms_to_degrees(h, m, s))


def hms_to_dms(h, m, s):
    """
    Convert degrees, arcminutes, arcseconds to an ``(hour, minute, second)``
    tuple.
    """

    return degrees_to_dms(hms_to_degrees(h, m, s))


def hours_to_decimal(h):
    """
    Convert any parseable hour value into a float value.
    """
    from . import angles
    return angles.Angle(h, unit=u.hourangle).hour


def hours_to_radians(h):
    """
    Convert an angle in Hours to Radians.
    """

    return u.hourangle.to(u.radian, h)


def hours_to_hms(h):
    """
    Convert an floating-point hour value into an ``(hour, minute,
    second)`` tuple.
    """

    sign = math.copysign(1.0, h)

    (hf, h) = math.modf(abs(h))  # (degree fraction, degree)
    (mf, m) = math.modf(hf * 60.0)  # (minute fraction, minute)
    s = mf * 60.0

    check_hms_ranges(h, m, s)  # throws exception if out of range

    return (float(sign * h), int(sign * m), sign * s)


def radians_to_degrees(r):
    """
    Convert an angle in Radians to Degrees.
    """
    return u.radian.to(u.degree, r)


def radians_to_hours(r):
    """
    Convert an angle in Radians to Hours.
    """
    return u.radian.to(u.hourangle, r)


def radians_to_hms(r):
    """
    Convert an angle in Radians to an ``(hour, minute, second)`` tuple.
    """

    hours = radians_to_hours(r)
    return hours_to_hms(hours)


def radians_to_dms(r):
    """
    Convert an angle in Radians to an ``(degree, arcminute,
    arcsecond)`` tuple.
    """

    degrees = u.radian.to(u.degree, r)
    return degrees_to_dms(degrees)


def hours_to_string(h, precision=5, pad=False, sep=('h', 'm', 's')):
    """
    Takes a decimal hour value and returns a string formatted as hms with
    separator specified by the 'sep' parameter.

    TODO: More detailed description here!
    """

    if pad:
        if h < 0:
            pad = 3
        else:
            pad = 2
    else:
        pad = 0

    if not isinstance(sep, tuple):
        # Note: This will convert 'hms' to ('h', 'm', 's'); a potentially nice
        # shortcut
        sep = tuple(sep)

    if len(sep) == 1:
        sep = sep + (sep[0], '')
    elif len(sep) == 2:
        sep = sep + ('',)
    elif len(sep) != 3:
        raise ValueError(
            "Invalid separator specification for converting angle to string.")

    literal = ('{0:0{pad}.0f}{sep[0]}{1:02d}{sep[1]}{2:0{width}.{precision}f}'
               '{sep[2]}')
    h, m, s = hours_to_hms(h)
    return literal.format(h, abs(m), abs(s), sep=sep, pad=pad,
                          width=(precision + 3), precision=precision)


def degrees_to_string(d, precision=5, pad=False, sep=':'):
    """
    Takes a decimal hour value and returns a string formatted as dms with
    separator specified by the 'sep' parameter.
    """

    if pad:
        if d < 0:
            pad = 3
        else:
            pad = 2
    else:
        pad = 0

    if not isinstance(sep, tuple):
        sep = tuple(sep)

    if len(sep) == 1:
        sep = sep + (sep[0], '')
    elif len(sep) == 2:
        sep = sep + ('',)
    elif len(sep) != 3:
        raise ValueError(
            "Invalid separator specification for converting angle to string.")

    literal = ('{0:0{pad}.0f}{sep[0]}{1:02d}{sep[1]}{2:0{width}.{precision}f}'
               '{sep[2]}')
    d, m, s = degrees_to_dms(d)
    return literal.format(d, abs(m), abs(s), sep=sep, pad=pad,
                          width=(precision + 3), precision=precision)


#<----------Spherical angular distances------------->
def small_angle_sphere_dist(lon1, lat1, lon2, lat2):
    """
    Euclidean angular distance "on a sphere" - only valid on sphere in the
    small-angle approximation.

    .. warning::
        Do not use this unless you know small-angle is a valid approximation
        for your problem and performance is a major conern.  In general this
        is very wrong.

    Inputs must be in radians.
    """

    from math import cos

    dlat = lat2 - lat1
    dlon = (lon2 - lon1) * cos((lat1 + lat2) / 2.)

    return (dlat ** 2 + dlon ** 2) ** 0.5


def simple_sphere_dist(lon1, lat1, lon2, lat2):
    """
    Simple formula for angular distance on a sphere: numerically unstable
    for small distances.

    Inputs must be in radians.
    """

    # FIXME: array: use numpy functions
    from math import acos, sin, cos

    cdlon = cos(lon2 - lon1)
    return acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(-lat2) * cdlon)


def haversine_sphere_dist(lon1, lat1, lon2, lat2):
    """
    Haversine formula for angular distance on a sphere: more stable at poles

    Inputs must be in radians.
    """

    # FIXME: array: use numpy functions
    from math import asin, sin, cos

    sdlat = sin((lat2 - lat1) / 2)
    sdlon = sin((lon2 - lon1) / 2)
    coslats = cos(lat1) * cos(lat2)

    return 2 * asin((sdlat ** 2 + coslats * sdlon ** 2) ** 0.5)


def haversine_atan_sphere_dist(lon1, lat1, lon2, lat2):
    """
    Haversine formula for angular distance on a sphere: more stable at poles.
    This version uses arctan instead of arcsin and thus does better with sign
    conventions.

    Inputs must be in radians.
    """

    # FIXME: array: use numpy functions
    from math import atan2, sin, cos

    sdlat = sin((lat2 - lat1) / 2)
    sdlon = sin((lon2 - lon1) / 2)
    coslats = cos(lat1) * cos(lat2)

    numerator = sdlat ** 2 + coslats * sdlon ** 2

    return 2 * atan2(numerator ** 0.5, (1 - numerator) ** 0.5)


def vincenty_sphere_dist(lon1, lat1, lon2, lat2):
    """
    Vincenty formula for angular distance on a sphere: stable at poles and
    antipodes but more complex/computationally expensive.

    Note that this is the only version actually used in the `AngularSeparation`
    classes, so the other `*_spher_dist` functions are only for possible
    future internal use.

    Inputs must be in radians.
    """
    #FIXME: array: use numpy functions
    from math import atan2, sin, cos

    sdlon = sin(lon2 - lon1)
    cdlon = cos(lon2 - lon1)
    slat1 = sin(lat1)
    slat2 = sin(lat2)
    clat1 = cos(lat1)
    clat2 = cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return atan2((num1 ** 2 + num2 ** 2) ** 0.5, denominator)
