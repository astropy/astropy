# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This module includes files automatically generated from ply (these end in
# _lextab.py and _parsetab.py). To generate these files, remove them from this
# folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/coordinates
#
# You can then commit the changes to the re-generated _lextab.py and
# _parsetab.py files.

"""
This module contains utility functions that are for internal use in
astropy.coordinates.angles. Mainly they are conversions from one format
of data to another.
"""

import os
import threading
from warnings import warn

import numpy as np

from .errors import (IllegalHourWarning, IllegalHourError,
                     IllegalMinuteWarning, IllegalMinuteError,
                     IllegalSecondWarning, IllegalSecondError)
from astropy.utils import format_exception
from astropy import units as u


TAB_HEADER = """# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/coordinates
#
# You can then commit the changes to this file.

"""


class _AngleParser:
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
    # For safe multi-threaded operation all class (but not instance)
    # members that carry state should be thread-local. They are stored
    # in the following class member
    _thread_local = threading.local()

    def __init__(self):
        # TODO: in principle, the parser should be invalidated if we change unit
        # system (from CDS to FITS, say).  Might want to keep a link to the
        # unit_registry used, and regenerate the parser/lexer if it changes.
        # Alternatively, perhaps one should not worry at all and just pre-
        # generate the parser for each release (as done for unit formats).
        # For some discussion of this problem, see
        # https://github.com/astropy/astropy/issues/5350#issuecomment-248770151
        if '_parser' not in _AngleParser._thread_local.__dict__:
            (_AngleParser._thread_local._parser,
             _AngleParser._thread_local._lexer) = self._make_parser()

    @classmethod
    def _get_simple_unit_names(cls):
        simple_units = set(
            u.radian.find_equivalent_units(include_prefix_units=True))
        simple_unit_names = set()
        # We filter out degree and hourangle, since those are treated
        # separately.
        for unit in simple_units:
            if unit != u.deg and unit != u.hourangle:
                simple_unit_names.update(unit.names)
        return sorted(simple_unit_names)

    @classmethod
    def _make_parser(cls):
        from astropy.extern.ply import lex, yacc

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
            r'((\d+\.\d*)|(\.\d+))([eE][+-−]?\d+)?'
            # The above includes Unicode "MINUS SIGN" \u2212.  It is
            # important to include the hyphen last, or the regex will
            # treat this as a range.
            t.value = float(t.value.replace('−', '-'))
            return t

        def t_UINT(t):
            r'\d+'
            t.value = int(t.value)
            return t

        def t_SIGN(t):
            r'[+−-]'
            # The above include Unicode "MINUS SIGN" \u2212.  It is
            # important to include the hyphen last, or the regex will
            # treat this as a range.
            if t.value == '+':
                t.value = 1.0
            else:
                t.value = -1.0
            return t

        def t_SIMPLE_UNIT(t):
            t.value = u.Unit(t.value)
            return t
        t_SIMPLE_UNIT.__doc__ = '|'.join(
            f'(?:{x})' for x in cls._get_simple_unit_names())

        t_COLON = ':'
        t_DEGREE = r'd(eg(ree(s)?)?)?|°'
        t_HOUR = r'hour(s)?|h(r)?|ʰ'
        t_MINUTE = r'm(in(ute(s)?)?)?|′|\'|ᵐ'
        t_SECOND = r's(ec(ond(s)?)?)?|″|\"|ˢ'

        # A string containing ignored characters (spaces)
        t_ignore = ' '

        # Error handling rule
        def t_error(t):
            raise ValueError(
                f"Invalid character at col {t.lexpos}")

        lexer_exists = os.path.exists(os.path.join(os.path.dirname(__file__),
                                      'angle_lextab.py'))

        # Build the lexer
        lexer = lex.lex(optimize=True, lextab='angle_lextab',
                        outputdir=os.path.dirname(__file__))

        if not lexer_exists:
            cls._add_tab_header('angle_lextab')

        def p_angle(p):
            '''
            angle : hms
                  | dms
                  | arcsecond
                  | arcminute
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
            colon : sign UINT COLON ufloat
                  | sign UINT COLON UINT COLON ufloat
            '''
            if len(p) == 5:
                p[0] = (p[1] * p[2], p[4])
            elif len(p) == 7:
                p[0] = (p[1] * p[2], p[4], p[6])

        def p_spaced(p):
            '''
            spaced : sign UINT ufloat
                   | sign UINT UINT ufloat
            '''
            if len(p) == 4:
                p[0] = (p[1] * p[2], p[3])
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
            hms : sign UINT HOUR
                | sign UINT HOUR ufloat
                | sign UINT HOUR UINT MINUTE
                | sign UINT HOUR UFLOAT MINUTE
                | sign UINT HOUR UINT MINUTE ufloat
                | sign UINT HOUR UINT MINUTE ufloat SECOND
                | generic HOUR
            '''
            if len(p) == 3:
                p[0] = (p[1], u.hourangle)
            elif len(p) == 4:
                p[0] = (p[1] * p[2], u.hourangle)
            elif len(p) in (5, 6):
                p[0] = ((p[1] * p[2], p[4]), u.hourangle)
            elif len(p) in (7, 8):
                p[0] = ((p[1] * p[2], p[4], p[6]), u.hourangle)

        def p_dms(p):
            '''
            dms : sign UINT DEGREE
                | sign UINT DEGREE ufloat
                | sign UINT DEGREE UINT MINUTE
                | sign UINT DEGREE UFLOAT MINUTE
                | sign UINT DEGREE UINT MINUTE ufloat
                | sign UINT DEGREE UINT MINUTE ufloat SECOND
                | generic DEGREE
            '''
            if len(p) == 3:
                p[0] = (p[1], u.degree)
            elif len(p) == 4:
                p[0] = (p[1] * p[2], u.degree)
            elif len(p) in (5, 6):
                p[0] = ((p[1] * p[2], p[4]), u.degree)
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

        def p_arcsecond(p):
            '''
            arcsecond : generic SECOND
            '''
            p[0] = (p[1], u.arcsecond)

        def p_arcminute(p):
            '''
            arcminute : generic MINUTE
            '''
            p[0] = (p[1], u.arcminute)

        def p_error(p):
            raise ValueError

        parser_exists = os.path.exists(os.path.join(os.path.dirname(__file__),
                                       'angle_parsetab.py'))

        parser = yacc.yacc(debug=False, tabmodule='angle_parsetab',
                           outputdir=os.path.dirname(__file__),
                           write_tables=True)

        if not parser_exists:
            cls._add_tab_header('angle_parsetab')

        return parser, lexer

    @classmethod
    def _add_tab_header(cls, name):

        lextab_file = os.path.join(os.path.dirname(__file__), name + '.py')

        with open(lextab_file, 'r') as f:
            contents = f.read()

        with open(lextab_file, 'w') as f:
            f.write(TAB_HEADER)
            f.write(contents)

    def parse(self, angle, unit, debug=False):
        try:
            found_angle, found_unit = self._thread_local._parser.parse(
                angle, lexer=self._thread_local._lexer, debug=debug)
        except ValueError as e:
            if str(e):
                raise ValueError("{} in angle {!r}".format(
                    str(e), angle))
            else:
                raise ValueError(
                    f"Syntax error parsing angle {angle!r}")

        if unit is None and found_unit is None:
            raise u.UnitsError("No unit specified")

        return found_angle, found_unit


def _check_hour_range(hrs):
    """
    Checks that the given value is in the range (-24, 24).
    """
    if np.any(np.abs(hrs) == 24.):
        warn(IllegalHourWarning(hrs, 'Treating as 24 hr'))
    elif np.any(hrs < -24.) or np.any(hrs > 24.):
        raise IllegalHourError(hrs)


def _check_minute_range(m):
    """
    Checks that the given value is in the range [0,60].  If the value
    is equal to 60, then a warning is raised.
    """
    if np.any(m == 60.):
        warn(IllegalMinuteWarning(m, 'Treating as 0 min, +1 hr/deg'))
    elif np.any(m < -60.) or np.any(m > 60.):
        # "Error: minutes not in range [-60,60) ({0}).".format(min))
        raise IllegalMinuteError(m)


def _check_second_range(sec):
    """
    Checks that the given value is in the range [0,60].  If the value
    is equal to 60, then a warning is raised.
    """
    if np.any(sec == 60.):
        warn(IllegalSecondWarning(sec, 'Treating as 0 sec, +1 min'))
    elif sec is None:
        pass
    elif np.any(sec < -60.) or np.any(sec > 60.):
        # "Error: seconds not in range [-60,60) ({0}).".format(sec))
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
        The unit used to interpret the string.  If ``unit`` is not
        provided, the unit must be explicitly represented in the
        string, either at the end or as number separators.

    debug : bool, optional
        If `True`, print debugging information from the parser.

    Returns
    -------
    value, unit : tuple
        ``value`` is the value as a floating point number or three-part
        tuple, and ``unit`` is a `Unit` instance which is either the
        unit passed in or the one explicitly mentioned in the input
        string.
    """
    return _AngleParser().parse(angle, unit, debug=debug)


def degrees_to_dms(d):
    """
    Convert a floating-point degree value into a ``(degree, arcminute,
    arcsecond)`` tuple.
    """
    sign = np.copysign(1.0, d)

    (df, d) = np.modf(np.abs(d))  # (degree fraction, degree)
    (mf, m) = np.modf(df * 60.)  # (minute fraction, minute)
    s = mf * 60.

    return np.floor(sign * d), sign * np.floor(m), sign * s


def dms_to_degrees(d, m, s=None):
    """
    Convert degrees, arcminute, arcsecond to a float degrees value.
    """

    _check_minute_range(m)
    _check_second_range(s)

    # determine sign
    sign = np.copysign(1.0, d)

    try:
        d = np.floor(np.abs(d))
        if s is None:
            m = np.abs(m)
            s = 0
        else:
            m = np.floor(np.abs(m))
            s = np.abs(s)
    except ValueError:
        raise ValueError(format_exception(
            "{func}: dms values ({1[0]},{2[1]},{3[2]}) could not be "
            "converted to numbers.", d, m, s))

    return sign * (d + m / 60. + s / 3600.)


def hms_to_hours(h, m, s=None):
    """
    Convert hour, minute, second to a float hour value.
    """

    check_hms_ranges(h, m, s)

    # determine sign
    sign = np.copysign(1.0, h)

    try:
        h = np.floor(np.abs(h))
        if s is None:
            m = np.abs(m)
            s = 0
        else:
            m = np.floor(np.abs(m))
            s = np.abs(s)
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

    sign = np.copysign(1.0, h)

    (hf, h) = np.modf(np.abs(h))  # (degree fraction, degree)
    (mf, m) = np.modf(hf * 60.0)  # (minute fraction, minute)
    s = mf * 60.0

    return (np.floor(sign * h), sign * np.floor(m), sign * s)


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


def sexagesimal_to_string(values, precision=None, pad=False, sep=(':',),
                          fields=3):
    """
    Given an already separated tuple of sexagesimal values, returns
    a string.

    See `hours_to_string` and `degrees_to_string` for a higher-level
    interface to this functionality.
    """

    # Check to see if values[0] is negative, using np.copysign to handle -0
    sign = np.copysign(1.0, values[0])
    # If the coordinates are negative, we need to take the absolute values.
    # We use np.abs because abs(-0) is -0
    # TODO: Is this true? (MHvK, 2018-02-01: not on my system)
    values = [np.abs(value) for value in values]

    if pad:
        if sign == -1:
            pad = 3
        else:
            pad = 2
    else:
        pad = 0

    if not isinstance(sep, tuple):
        sep = tuple(sep)

    if fields < 1 or fields > 3:
        raise ValueError(
            "fields must be 1, 2, or 3")

    if not sep:  # empty string, False, or None, etc.
        sep = ('', '', '')
    elif len(sep) == 1:
        if fields == 3:
            sep = sep + (sep[0], '')
        elif fields == 2:
            sep = sep + ('', '')
        else:
            sep = ('', '', '')
    elif len(sep) == 2:
        sep = sep + ('',)
    elif len(sep) != 3:
        raise ValueError(
            "Invalid separator specification for converting angle to string.")

    # Simplify the expression based on the requested precision.  For
    # example, if the seconds will round up to 60, we should convert
    # it to 0 and carry upwards.  If the field is hidden (by the
    # fields kwarg) we round up around the middle, 30.0.
    if precision is None:
        rounding_thresh = 60.0 - (10.0 ** -4)
    else:
        rounding_thresh = 60.0 - (10.0 ** -precision)

    if fields == 3 and values[2] >= rounding_thresh:
        values[2] = 0.0
        values[1] += 1.0
    elif fields < 3 and values[2] >= 30.0:
        values[1] += 1.0

    if fields >= 2 and values[1] >= 60.0:
        values[1] = 0.0
        values[0] += 1.0
    elif fields < 2 and values[1] >= 30.0:
        values[0] += 1.0

    literal = []
    last_value = ''
    literal.append('{0:0{pad}.0f}{sep[0]}')
    if fields >= 2:
        literal.append('{1:02d}{sep[1]}')
    if fields == 3:
        if precision is None:
            last_value = '{:.4f}'.format(abs(values[2]))
            last_value = last_value.rstrip('0').rstrip('.')
        else:
            last_value = '{0:.{precision}f}'.format(
                abs(values[2]), precision=precision)
        if len(last_value) == 1 or last_value[1] == '.':
            last_value = '0' + last_value
        literal.append('{last_value}{sep[2]}')
    literal = ''.join(literal)
    return literal.format(np.copysign(values[0], sign),
                          int(values[1]), values[2],
                          sep=sep, pad=pad,
                          last_value=last_value)


def hours_to_string(h, precision=5, pad=False, sep=('h', 'm', 's'),
                    fields=3):
    """
    Takes a decimal hour value and returns a string formatted as hms with
    separator specified by the 'sep' parameter.

    ``h`` must be a scalar.
    """
    h, m, s = hours_to_hms(h)
    return sexagesimal_to_string((h, m, s), precision=precision, pad=pad,
                                 sep=sep, fields=fields)


def degrees_to_string(d, precision=5, pad=False, sep=':', fields=3):
    """
    Takes a decimal hour value and returns a string formatted as dms with
    separator specified by the 'sep' parameter.

    ``d`` must be a scalar.
    """
    d, m, s = degrees_to_dms(d)
    return sexagesimal_to_string((d, m, s), precision=precision, pad=pad,
                                 sep=sep, fields=fields)


def angular_separation(lon1, lat1, lon2, lat2):
    """
    Angular separation between two points on a sphere.

    Parameters
    ----------
    lon1, lat1, lon2, lat2 : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the two points. Quantities should be in
        angular units; floats in radians.

    Returns
    -------
    angular separation : `~astropy.units.Quantity` or float
        Type depends on input; `Quantity` in angular units, or float in
        radians.

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1]_,
    which is slightly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    .. [1] https://en.wikipedia.org/wiki/Great-circle_distance
    """

    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.arctan2(np.hypot(num1, num2), denominator)


def position_angle(lon1, lat1, lon2, lat2):
    """
    Position Angle (East of North) between two points on a sphere.

    Parameters
    ----------
    lon1, lat1, lon2, lat2 : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the two points. Quantities should be in
        angular units; floats in radians.

    Returns
    -------
    pa : `~astropy.coordinates.Angle`
        The (positive) position angle of the vector pointing from position 1 to
        position 2.  If any of the angles are arrays, this will contain an array
        following the appropriate `numpy` broadcasting rules.

    """
    from .angles import Angle

    deltalon = lon2 - lon1
    colat = np.cos(lat2)

    x = np.sin(lat2) * np.cos(lat1) - colat * np.sin(lat1) * np.cos(deltalon)
    y = np.sin(deltalon) * colat

    return Angle(np.arctan2(y, x), u.radian).wrap_at(360*u.deg)


def offset_by(lon, lat, posang, distance):
    """
    Point with the given offset from the given point.

    Parameters
    ----------
    lon, lat, posang, distance : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the starting point,
        position angle and distance to the final point.
        Quantities should be in angular units; floats in radians.
        Polar points at lat= +/-90 are treated as limit of +/-(90-epsilon) and same lon.

    Returns
    -------
    lon, lat : `~astropy.coordinates.Angle`
        The position of the final point.  If any of the angles are arrays,
        these will contain arrays following the appropriate `numpy` broadcasting rules.
        0 <= lon < 2pi.

    Notes
    -----
    """
    from .angles import Angle

    # Calculations are done using the spherical trigonometry sine and cosine rules
    # of the triangle A at North Pole,   B at starting point,   C at final point
    # with angles     A (change in lon), B (posang),            C (not used, but negative reciprocal posang)
    # with sides      a (distance),      b (final co-latitude), c (starting colatitude)
    # B, a, c are knowns; A and b are unknowns
    # https://en.wikipedia.org/wiki/Spherical_trigonometry

    cos_a = np.cos(distance)
    sin_a = np.sin(distance)
    cos_c = np.sin(lat)
    sin_c = np.cos(lat)
    cos_B = np.cos(posang)
    sin_B = np.sin(posang)

    # cosine rule: Know two sides: a,c and included angle: B; get unknown side b
    cos_b = cos_c * cos_a + sin_c * sin_a * cos_B
    # sin_b = np.sqrt(1 - cos_b**2)
    # sine rule and cosine rule for A (using both lets arctan2 pick quadrant).
    # multiplying both sin_A and cos_A by x=sin_b * sin_c prevents /0 errors
    # at poles.  Correct for the x=0 multiplication a few lines down.
    # sin_A/sin_a == sin_B/sin_b    # Sine rule
    xsin_A = sin_a * sin_B * sin_c
    # cos_a == cos_b * cos_c + sin_b * sin_c * cos_A  # cosine rule
    xcos_A = cos_a - cos_b * cos_c

    A = Angle(np.arctan2(xsin_A, xcos_A), u.radian)
    # Treat the poles as if they are infinitesimally far from pole but at given lon
    small_sin_c = sin_c < 1e-12
    if small_sin_c.any():
        # For south pole (cos_c = -1), A = posang; for North pole, A=180 deg - posang
        A_pole = (90*u.deg + cos_c*(90*u.deg-Angle(posang, u.radian))).to(u.rad)
        if A.shape:
            # broadcast to ensure the shape is like that of A, which is also
            # affected by the (possible) shapes of lat, posang, and distance.
            small_sin_c = np.broadcast_to(small_sin_c, A.shape)
            A[small_sin_c] = A_pole[small_sin_c]
        else:
            A = A_pole

    outlon = (Angle(lon, u.radian) + A).wrap_at(360.0*u.deg).to(u.deg)
    outlat = Angle(np.arcsin(cos_b), u.radian).to(u.deg)

    return outlon, outlat
