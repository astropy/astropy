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
This module contains formatting functions that are for internal use in
astropy.coordinates.angles. Mainly they are conversions from one format
of data to another.
"""

import threading
from warnings import warn

import numpy as np

from astropy import units as u
from astropy.utils import format_exception, parsing
from astropy.utils.decorators import deprecated

from .errors import (
    IllegalHourError,
    IllegalHourWarning,
    IllegalMinuteError,
    IllegalMinuteWarning,
    IllegalSecondError,
    IllegalSecondWarning,
)


class _AngleParser:
    """
    Parses the various angle formats including:

       * 01:02:30.43 degrees
       * 1 2 0 hours
       * 1°2′3″
       * 1d2m3s
       * -1h2m3s
       * 1°2′3″N

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
        if "_parser" not in _AngleParser._thread_local.__dict__:
            (
                _AngleParser._thread_local._parser,
                _AngleParser._thread_local._lexer,
            ) = self._make_parser()

    @classmethod
    def _get_simple_unit_names(cls):
        simple_units = set(u.radian.find_equivalent_units(include_prefix_units=True))
        simple_unit_names = set()
        # We filter out degree and hourangle, since those are treated
        # separately.
        for unit in simple_units:
            if unit != u.deg and unit != u.hourangle:
                simple_unit_names.update(unit.names)
        return sorted(simple_unit_names)

    @classmethod
    def _make_parser(cls):
        # List of token names.
        tokens = (
            "SIGN",
            "UINT",
            "UFLOAT",
            "COLON",
            "DEGREE",
            "HOUR",
            "MINUTE",
            "SECOND",
            "SIMPLE_UNIT",
            "EASTWEST",
            "NORTHSOUTH",
        )

        # NOTE THE ORDERING OF THESE RULES IS IMPORTANT!!
        # Regular expression rules for simple tokens
        def t_UFLOAT(t):
            r"((\d+\.\d*)|(\.\d+))([eE][+-−]?\d+)?"
            # The above includes Unicode "MINUS SIGN" \u2212.  It is
            # important to include the hyphen last, or the regex will
            # treat this as a range.
            t.value = float(t.value.replace("−", "-"))
            return t

        def t_UINT(t):
            r"\d+"
            t.value = int(t.value)
            return t

        def t_SIGN(t):
            r"[+−-]"
            # The above include Unicode "MINUS SIGN" \u2212.  It is
            # important to include the hyphen last, or the regex will
            # treat this as a range.
            if t.value == "+":
                t.value = 1.0
            else:
                t.value = -1.0
            return t

        def t_EASTWEST(t):
            r"[EW]$"
            t.value = -1.0 if t.value == "W" else 1.0
            return t

        def t_NORTHSOUTH(t):
            r"[NS]$"
            # We cannot use lower-case letters otherwise we'll confuse
            # s[outh] with s[econd]
            t.value = -1.0 if t.value == "S" else 1.0
            return t

        def t_SIMPLE_UNIT(t):
            t.value = u.Unit(t.value)
            return t

        t_SIMPLE_UNIT.__doc__ = "|".join(
            f"(?:{x})" for x in cls._get_simple_unit_names()
        )

        t_COLON = ":"
        t_DEGREE = r"d(eg(ree(s)?)?)?|°"
        t_HOUR = r"hour(s)?|h(r)?|ʰ"
        t_MINUTE = r"m(in(ute(s)?)?)?|′|\'|ᵐ"
        t_SECOND = r"s(ec(ond(s)?)?)?|″|\"|ˢ"

        # A string containing ignored characters (spaces)
        t_ignore = " "

        # Error handling rule
        def t_error(t):
            raise ValueError(f"Invalid character at col {t.lexpos}")

        lexer = parsing.lex(lextab="angle_lextab", package="astropy/coordinates")

        def p_angle(p):
            """
            angle : sign hms eastwest
                  | sign dms dir
                  | sign arcsecond dir
                  | sign arcminute dir
                  | sign simple dir
            """
            sign = p[1] * p[3]
            value, unit = p[2]
            if isinstance(value, tuple):
                p[0] = ((sign * value[0],) + value[1:], unit)
            else:
                p[0] = (sign * value, unit)

        def p_sign(p):
            """
            sign : SIGN
                 |
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1.0

        def p_eastwest(p):
            """
            eastwest : EASTWEST
                     |
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1.0

        def p_dir(p):
            """
            dir : EASTWEST
                | NORTHSOUTH
                |
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1.0

        def p_ufloat(p):
            """
            ufloat : UFLOAT
                   | UINT
            """
            p[0] = p[1]

        def p_colon(p):
            """
            colon : UINT COLON ufloat
                  | UINT COLON UINT COLON ufloat
            """
            if len(p) == 4:
                p[0] = (p[1], p[3])
            elif len(p) == 6:
                p[0] = (p[1], p[3], p[5])

        def p_spaced(p):
            """
            spaced : UINT ufloat
                   | UINT UINT ufloat
            """
            if len(p) == 3:
                p[0] = (p[1], p[2])
            elif len(p) == 4:
                p[0] = (p[1], p[2], p[3])

        def p_generic(p):
            """
            generic : colon
                    | spaced
                    | ufloat
            """
            p[0] = p[1]

        def p_hms(p):
            """
            hms : UINT HOUR
                | UINT HOUR ufloat
                | UINT HOUR UINT MINUTE
                | UINT HOUR UFLOAT MINUTE
                | UINT HOUR UINT MINUTE ufloat
                | UINT HOUR UINT MINUTE ufloat SECOND
                | generic HOUR
            """
            if len(p) == 3:
                p[0] = (p[1], u.hourangle)
            elif len(p) in (4, 5):
                p[0] = ((p[1], p[3]), u.hourangle)
            elif len(p) in (6, 7):
                p[0] = ((p[1], p[3], p[5]), u.hourangle)

        def p_dms(p):
            """
            dms : UINT DEGREE
                | UINT DEGREE ufloat
                | UINT DEGREE UINT MINUTE
                | UINT DEGREE UFLOAT MINUTE
                | UINT DEGREE UINT MINUTE ufloat
                | UINT DEGREE UINT MINUTE ufloat SECOND
                | generic DEGREE
            """
            if len(p) == 3:
                p[0] = (p[1], u.degree)
            elif len(p) in (4, 5):
                p[0] = ((p[1], p[3]), u.degree)
            elif len(p) in (6, 7):
                p[0] = ((p[1], p[3], p[5]), u.degree)

        def p_simple(p):
            """
            simple : generic
                   | generic SIMPLE_UNIT
            """
            if len(p) == 2:
                p[0] = (p[1], None)
            else:
                p[0] = (p[1], p[2])

        def p_arcsecond(p):
            """
            arcsecond : generic SECOND
            """
            p[0] = (p[1], u.arcsecond)

        def p_arcminute(p):
            """
            arcminute : generic MINUTE
            """
            p[0] = (p[1], u.arcminute)

        def p_error(p):
            raise ValueError

        parser = parsing.yacc(tabmodule="angle_parsetab", package="astropy/coordinates")

        return parser, lexer

    def parse(self, angle, unit, debug=False):
        try:
            found_angle, found_unit = self._thread_local._parser.parse(
                angle, lexer=self._thread_local._lexer, debug=debug
            )
        except ValueError as e:
            if str(e):
                raise ValueError(f"{str(e)} in angle {angle!r}") from e
            else:
                raise ValueError(f"Syntax error parsing angle {angle!r}") from e

        if unit is None and found_unit is None:
            raise u.UnitsError("No unit specified")

        return found_angle, found_unit


def _check_hour_range(hrs):
    """
    Checks that the given value is in the range (-24, 24).
    """
    if np.any(np.abs(hrs) == 24.0):
        warn(IllegalHourWarning(hrs, "Treating as 24 hr"))
    elif np.any(hrs < -24.0) or np.any(hrs > 24.0):
        raise IllegalHourError(hrs)


def _check_minute_range(m):
    """
    Checks that the given value is in the range [0,60].  If the value
    is equal to 60, then a warning is raised.
    """
    if np.any(m == 60.0):
        warn(IllegalMinuteWarning(m, "Treating as 0 min, +1 hr/deg"))
    elif np.any(m < -60.0) or np.any(m > 60.0):
        # "Error: minutes not in range [-60,60) ({0}).".format(min))
        raise IllegalMinuteError(m)


def _check_second_range(sec):
    """
    Checks that the given value is in the range [0,60].  If the value
    is equal to 60, then a warning is raised.
    """
    if np.any(sec == 60.0):
        warn(IllegalSecondWarning(sec, "Treating as 0 sec, +1 min"))
    elif sec is None:
        pass
    elif np.any(sec < -60.0) or np.any(sec > 60.0):
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
    (mf, m) = np.modf(df * 60.0)  # (minute fraction, minute)
    s = mf * 60.0

    return np.floor(sign * d), sign * np.floor(m), sign * s


@deprecated(
    since="5.1",
    message=(
        "dms_to_degrees (or creating an Angle with a tuple) has ambiguous "
        "behavior when the degree value is 0. Use {alternative}."
    ),
    alternative=(
        "another way of creating angles instead (e.g. a less "
        "ambiguous string like '-0d1m2.3s')"
    ),
)
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
    except ValueError as err:
        raise ValueError(
            format_exception(
                "{func}: dms values ({1[0]},{2[1]},{3[2]}) could not be "
                "converted to numbers.",
                d,
                m,
                s,
            )
        ) from err

    return sign * (d + m / 60.0 + s / 3600.0)


@deprecated(
    since="5.1",
    message=(
        "hms_to_hours (or creating an Angle with a tuple) has ambiguous "
        "behavior when the hour value is 0. Use {alternative}."
    ),
    alternative=(
        "another way of creating angles instead (e.g. a less "
        "ambiguous string like '-0h1m2.3s')"
    ),
)
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
    except ValueError as err:
        raise ValueError(
            format_exception(
                "{func}: HMS values ({1[0]},{2[1]},{3[2]}) could not be "
                "converted to numbers.",
                h,
                m,
                s,
            )
        ) from err

    return sign * (h + m / 60.0 + s / 3600.0)


def hms_to_degrees(h, m, s):
    """
    Convert hour, minute, second to a float degrees value.
    """
    return hms_to_hours(h, m, s) * 15.0


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


def sexagesimal_to_string(values, precision=None, pad=False, sep=(":",), fields=3):
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
        raise ValueError("fields must be 1, 2, or 3")

    if not sep:  # empty string, False, or None, etc.
        sep = ("", "", "")
    elif len(sep) == 1:
        if fields == 3:
            sep = sep + (sep[0], "")
        elif fields == 2:
            sep = sep + ("", "")
        else:
            sep = ("", "", "")
    elif len(sep) == 2:
        sep = sep + ("",)
    elif len(sep) != 3:
        raise ValueError(
            "Invalid separator specification for converting angle to string."
        )

    # Simplify the expression based on the requested precision.  For
    # example, if the seconds will round up to 60, we should convert
    # it to 0 and carry upwards.  If the field is hidden (by the
    # fields kwarg) we round up around the middle, 30.0.
    if precision is None:
        rounding_thresh = 60.0 - (10.0**-8)
    else:
        rounding_thresh = 60.0 - (10.0**-precision)

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
    last_value = ""
    literal.append("{0:0{pad}.0f}{sep[0]}")
    if fields >= 2:
        literal.append("{1:02d}{sep[1]}")
    if fields == 3:
        if precision is None:
            last_value = f"{abs(values[2]):.8f}"
            last_value = last_value.rstrip("0").rstrip(".")
        else:
            last_value = "{0:.{precision}f}".format(abs(values[2]), precision=precision)
        if len(last_value) == 1 or last_value[1] == ".":
            last_value = "0" + last_value
        literal.append("{last_value}{sep[2]}")
    literal = "".join(literal)
    return literal.format(
        np.copysign(values[0], sign),
        int(values[1]),
        values[2],
        sep=sep,
        pad=pad,
        last_value=last_value,
    )


def hours_to_string(h, precision=5, pad=False, sep=("h", "m", "s"), fields=3):
    """
    Takes a decimal hour value and returns a string formatted as hms with
    separator specified by the 'sep' parameter.

    ``h`` must be a scalar.
    """
    h, m, s = hours_to_hms(h)
    return sexagesimal_to_string(
        (h, m, s), precision=precision, pad=pad, sep=sep, fields=fields
    )


def degrees_to_string(d, precision=5, pad=False, sep=":", fields=3):
    """
    Takes a decimal hour value and returns a string formatted as dms with
    separator specified by the 'sep' parameter.

    ``d`` must be a scalar.
    """
    d, m, s = degrees_to_dms(d)
    return sexagesimal_to_string(
        (d, m, s), precision=precision, pad=pad, sep=sep, fields=fields
    )
