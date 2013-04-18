# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains utility functions that are for internal use in
astropy.coordinates.angles. Mainly they are conversions from one format
of data to another.
"""
from __future__ import unicode_literals

import re
import math
import inspect # NB: get the function name with: inspect.stack()[0][3]
from warnings import warn

from .errors import *
from ..utils import format_exception


def _check_hour_range(hrs):
    ''' Checks that the given value is in the range (-24,24). '''
    if math.fabs(hrs) == 24.:
        warn(IllegalHourWarning(hrs))
    elif not -24. < hrs < 24.:
        raise IllegalHourError(hrs)


def _check_minute_range(min):
    ''' Checks that the given value is in the range [0,60].
    If the value is equal to 60, then a warning is raised. '''
    if min == 60.:
        warn(IllegalMinuteWarning(min))
    elif not 0. <= min < 60.:
        # "Error: minutes not in range [0,60) ({0}).".format(min))
        raise IllegalMinuteError(min)


def _check_second_range(sec):
    ''' Checks that the given value is in the range [0,60].
    If the value is equal to 60, then a warning is raised. '''
    if sec == 60.:
        warn(IllegalSecondWarning(sec))
    elif not 0. <= sec < 60.:
        # "Error: seconds not in range [0,60) ({0}).".format(sec))
        raise IllegalSecondError(sec)


def check_hms_ranges(h, m, s):
    _check_hour_range(h)
    _check_minute_range(m)
    _check_second_range(s)
    return None

# these regexes are used in parse_degrees
# accept these as (one or more repeated) delimiters: :, whitespace, /
_dms_div_regex_str = '[:|/|\t|\-|\sDdMmSs]{1,2}'
# Look for a pattern where d,m,s is specified
_dms_regex = re.compile('^([+-]{0,1}\d{1,3})' + _dms_div_regex_str +
                        '(\d{1,2})' + _dms_div_regex_str +
                        '(\d{1,2}[\.0-9]*)' + '[Ss]{0,1}' + '$')
# look for a pattern where only d,m is specified
_dm_regex = re.compile('^([+-]{0,1}\d{1,3})' + _dms_div_regex_str +
                       '(\d{1,2})' + '[Mm]{0,1}' + '$')


def parse_degrees(degrees, output_dms=False):
    """
    Parses an input "degrees" value into decimal degrees or a
    degree,arcminute,arcsecond tuple.

    Convert degrees given in any parseable format (float, string, or Angle)
    into degrees, arcminutes, and arcseconds components or decimal degrees.

    Parameters
    ----------
    degrees : float, int, str
        If a string, accepts values in these formats:
            * [+|-]DD:MM:SS.sss (string), e.g. +12:52:32.423 or -12:52:32.423
            * DD.dddd (float, string), e.g. 12.542326
            * DD MM SS.sss (string, array), e.g. +12 52 32.423
        Whitespace may be spaces and/or tabs.
    output_dms : bool
        If True, returns a tuple of (degree, arcminute, arcsecond)

    Returns
    -------
    deg : float or tuple
         Returns degrees in decimal form unless the keyword "output_dms" is
         True, in which case a tuple (d, m, s).
    """

    from .angles import Angle

    # either a string or a float
    x = degrees

    if isinstance(x, float) or isinstance(x, int):
        parsed_degrees = float(x)

    elif isinstance(x, basestring):
        x = x.strip()

        string_parsed = False

        # See if the string is just a float or int value.
        try:
            parsed_degrees = float(x)
            string_parsed = True
        except ValueError:
            pass

        if not string_parsed:
            try:
                elems = _dms_regex.search(x).groups()
                # We need to convert elems[0] to float to preserve sign for -0.0
                parsed_degrees = dms_to_degrees(
                    float(elems[0]), int(elems[1]), float(elems[2]))
                string_parsed = True
            except AttributeError:
                # regular expression did not match - try again below
                # make sure to let things like IllegalMinuteError, etc. through
                pass

        if not string_parsed:
            try:
                elems = _dm_regex.search(x).groups()
                # We need to convert elems[0] to float to preserve sign for -0.0
                parsed_degrees = dms_to_degrees(
                    float(elems[0]), int(elems[1]), 0.0)
                string_parsed = True
            except AttributeError:
                # regular expression did not match - try again below
                # make sure to let things like IllegalMinuteError, etc. through
                pass

        if not string_parsed:
            # look for a '°' symbol
            for unitStr in ["degrees", "degree", "deg", "d", "°"]:
                x = x.replace(unitStr, '')
                try:
                    parsed_degrees = float(x)
                    string_parsed = True
                except ValueError:
                    pass

        if not string_parsed:
            raise ValueError(format_exception(
                "{func}: Invalid input string! ('{0}')", x))

    elif isinstance(x, Angle):
        parsed_degrees = x.degrees

    elif isinstance(x, tuple):
        parsed_degrees = dms_to_degrees(*x)

    else:
        raise ValueError(format_exception(
            "{func}: could not parse value of {0}.", type(x)))

    return degrees_to_dms(parsed_degrees) if output_dms else parsed_degrees


# these regexes are used in parse_hours
# accept these as (one or more repeated) delimiters: :, whitespace, /
_hms_div_regex_str = '[:|/|\t|\-|\sHhMmSs]{1,2}'
# Look for a pattern where h,m,s is specified
_hms_regex = re.compile('^([+-]{0,1}\d{1,2})' + _hms_div_regex_str +
                        '(\d{1,2})' + _hms_div_regex_str +
                        '(\d{1,2}[\.0-9]*)' + '[Ss]{0,1}' + '$')
# look for a pattern where only h,m is specified
_hm_regex = re.compile('^([+-]{0,1}\d{1,2})' + _hms_div_regex_str +
                       '(\d{1,2})' + '[Mm]{0,1}' + '$')


def parse_hours(hours, output_hms=False):
    """
    Returns an hour value (as a decimal or HMS tuple) from the integer, float,
    or string provided.

    Convert hours given in any parseable format (float, string, tuple, list, or
    Angle) into hour, minute, and seconds components or decimal hours.

    Parameters
    ----------
    hours : float, str, int
        If a string, accepts values in these formats:
            * HH:MM:SS.sss (string), e.g. 12:52:32.423
            * HH.dddd (float, string), e.g. 12.542326
            * HH MM SS.sss (string, array), e.g. 12 52 32.423
        Surrounding whitespace in a string value is allowed.
    output_hms : bool
        If True, returns a tuple of (hour, minute, second)

    Returns
    -------
    hrs : float or tuple
         Returns degrees in hours form unless the keyword "output_dms" is
         True, in which case a tuple (h, m, s).
    """
    import datetime as py_datetime

    from .angles import Angle

    # either a string or a float
    x = hours

    if isinstance(x, float) or isinstance(x, int):
        parsed_hours = x
        parsed_hms = hours_to_hms(parsed_hours)

    elif isinstance(x, basestring):
        x = x.strip()

        try:
            parsed_hours = float(x)
            parsed_hms = hours_to_hms(parsed_hours)
        except ValueError:

            string_parsed = False

            try:
                elems = _hms_regex.search(x).groups()
                string_parsed = True
            except:
                pass  # try again below

            if string_parsed:
                h, m, s = float(elems[0]), int(elems[1]), float(elems[2])
                parsed_hours = hms_to_hours(h, m, s)
                parsed_hms = (h, m, s)

            else:
                try:
                    elems = _hm_regex.search(x).groups()
                    string_parsed = True
                except:
                    raise ValueError(format_exception(
                        "{func}: Invalid input string, can't parse to "
                        "HMS. ({0})", x))
                h, m, s = float(elems[0]), int(elems[1]), 0.0
                parsed_hours = hms_to_hours(h, m, s)
                parsed_hms = (h, m, s)

    elif isinstance(x, Angle):
        parsed_hours = x.hours
        parsed_hms = hours_to_hms(parsed_hours)

    #TODO: make a decision as to whether or not to allow datetime objects
    #elif isinstance(x, py_datetime.datetime):
    #    parsed_hours = datetimeToDecimalTime(x)
    #    parsed_hms = hours_to_hms(parsed_hours)

    elif isinstance(x, tuple):
        if len(x) == 3:
            parsed_hours = hms_to_hours(*x)
            parsed_hms = x
        else:
            raise ValueError(format_exception(
                "{filename}:{func}: Incorrect number of values given, expected "
                "(h,m,s), got: {0}", x))

    elif isinstance(x, list):
        if len(x) == 3:
            try:
                h = float(x[0])
                m = float(x[1])
                s = float(x[2])
            except ValueError:
                raise ValueError(format_exception(
                    "{filename}:{func}: Array values ([h,m,s] expected) "
                    "could not be coerced into floats. {0}", x))

            parsed_hours = hms_to_hours(h, m, s)
            parsed_hms = (h, m, s)
            if output_hms:
                return (h, m, s)
            else:
                return hms_to_hours(h, m, s)

        else:
            raise ValueError(format_exception(
                "{filename}:{func}: Array given must contain exactly three "
                "elements ([h,m,s]), provided: {0}", x))

    else:
        raise ValueError(
            "parse_hours: could not parse value of type "
            "{0}.".format(type(x).__name__))

    if output_hms:
        return parsed_hms
    else:
        return parsed_hours


def parse_radians(radians):
    """
    Parses an input "radians" value into a float number.

    Convert radians given in any parseable format (float or Angle) into float
    radians.

    ..Note::
        This function is mostly for consistency with the other "parse"
        functions, like parse_hours and parse_degrees.

    Parameters
    ----------
    radians : float, int, Angle
        The input angle.
    """

    from .angles import Angle

    x = radians

    if type(x) in [float, int]:
        return float(x)
    elif isinstance(x, Angle):
        return x.radians
    else:
        raise ValueError(format_exception(
            "{func}: could not parse value of type {0}.", type(x).__name__))


def degrees_to_dms(d):
    """
    Convert any parseable degree value (see: parse_degrees) into a
    ``(degree, arcminute, arcsecond)`` tuple.
    """
    sign = math.copysign(1.0, d)

    (df, d) = math.modf(abs(d))  # (degree fraction, degree)
    (mf, m) = math.modf(df * 60.)  # (minute fraction, minute)
    s = mf * 60.

    _check_minute_range(m)
    _check_second_range(s)

    return (float(sign * d), int(sign * m), sign * s)


def dms_to_degrees(d, m, s):
    """ Convert degrees, arcminute, arcsecond to a float degrees value. """

    _check_minute_range(m)
    _check_second_range(s)

    # determine sign
    sign = math.copysign(1.0, d)

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
    """ Convert hour, minute, second to a float hour value. """

    check_hms_ranges(h, m, s)

    # determine sign
    sign = math.copysign(1.0, h)

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
    """ Convert hour, minute, second to a float degrees value. """

    return hms_to_hours(h, m, s) * 15.


def hms_to_radians(h, m, s):
    """ Convert hour, minute, second to a float radians value. """

    return math.radians(hms_to_degrees(h, m, s))


def hms_to_dms(h, m, s):
    """
    Convert degrees, arcminutes, arcseconds to an ``(hour, minute, second)``
    tuple.
    """

    return degrees_to_dms(hms_to_degrees(h, m, s))


def hours_to_decimal(h):
    """
    Convert any parseable hour value (see: parse_hours) into a float value.
    """

    return parse_hours(h, output_hms=False)


def hours_to_radians(h):
    """ Convert an angle in Hours to Radians. """

    return math.radians(h * 15.)


def hours_to_hms(h):
    """
    Convert any parseable hour value (see: parse_hours) into an
    ``(hour, minute, second)`` tuple.
    """

    sign = math.copysign(1.0, h)

    (hf, h) = math.modf(abs(h))  # (degree fraction, degree)
    (mf, m) = math.modf(hf * 60.0)  # (minute fraction, minute)
    s = mf * 60.0

    check_hms_ranges(h, m, s)  # throws exception if out of range

    return (float(sign * h), int(sign * m), sign * s)


def radians_to_degrees(r):
    """ Convert an angle in Radians to Degrees """

    try:
        r = float(r)
    except ValueError:
        raise ValueError(format_exception(
            "{func}: degree value ({1[0]}) could not be converted to a "
            "float.", r))

    return math.degrees(r)


def radians_to_hours(r):
    """ Convert an angle in Radians to Hours """

    return radians_to_degrees(r) / 15.


def radians_to_hms(r):
    """ Convert an angle in Radians to an hour,minute,second tuple """

    hours = radians_to_hours(r)
    return hours_to_hms(hours)


def radians_to_dms(r):
    """ Convert an angle in Radians to an degree,arcminute,arcsecond tuple """

    degrees = math.degrees(r)
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
