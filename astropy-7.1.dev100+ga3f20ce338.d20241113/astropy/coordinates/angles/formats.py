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
from astropy.utils import parsing

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

        def t_MINUTE(t):
            r"m(in(ute(s)?)?)?|′|\'|ᵐ"
            t.value = u.arcmin
            return t

        def t_SECOND(t):
            r"s(ec(ond(s)?)?)?|″|\"|ˢ"  # codespell:ignore ond
            t.value = u.arcsec
            return t

        t_COLON = ":"
        t_DEGREE = r"d(eg(ree(s)?)?)?|°"
        t_HOUR = r"hour(s)?|h(r)?|ʰ"

        # A string containing ignored characters (spaces)
        t_ignore = " "

        # Error handling rule
        def t_error(t):
            raise ValueError(f"Invalid character at col {t.lexpos}")

        lexer = parsing.lex(lextab="angle_lextab", package="astropy/coordinates/angles")

        def p_angle(p):
            """
            angle : sign hms eastwest
                  | sign dms dir
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
            p[0] = p[1] if len(p) == 2 else 1.0

        def p_eastwest(p):
            """
            eastwest : EASTWEST
                     |
            """
            p[0] = p[1] if len(p) == 2 else 1.0

        def p_dir(p):
            """
            dir : EASTWEST
                | NORTHSOUTH
                |
            """
            p[0] = p[1] if len(p) == 2 else 1.0

        def p_ufloat(p):
            """
            ufloat : UFLOAT
                   | UINT
            """
            p[0] = p[1]

        def p_generic(p):
            """
            generic : ufloat
                    | UINT ufloat
                    | UINT COLON ufloat
                    | UINT UINT ufloat
                    | UINT COLON UINT COLON ufloat
            """
            match p[1:]:
                case [p1]:
                    p[0] = p1
                case [p1, p2] | [p1, ":", p2]:
                    p[0] = (p1, p2)
                case [p1, p2, p3] | [p1, _, p2, _, p3]:
                    p[0] = (p1, p2, p3)

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
                   | generic MINUTE
                   | generic SECOND
                   | generic SIMPLE_UNIT
            """
            p[0] = (p[1], None if len(p) == 2 else p[2])

        def p_error(p):
            raise ValueError

        parser = parsing.yacc(
            tabmodule="angle_parsetab", package="astropy/coordinates/angles"
        )

        return parser, lexer

    def parse(self, angle, unit, debug=False):
        try:
            found_angle, found_unit = self._thread_local._parser.parse(
                angle, lexer=self._thread_local._lexer, debug=debug
            )
        except ValueError as e:
            raise ValueError(
                f"{str(e) or 'syntax error'} parsing angle {angle!r}"
            ) from e

        if unit is None and found_unit is None:
            raise u.UnitsError("No unit specified")

        return found_angle, found_unit


def _check_hour_range(hrs: float) -> None:
    """
    Checks that the given value is in the range [-24,24].  If the value
    is equal to -24 or 24, then a warning is raised.
    """
    if not -24.0 < hrs < 24.0:
        if abs(hrs) != 24.0:
            raise IllegalHourError(hrs)
        warn(IllegalHourWarning(hrs, "Treating as 24 hr"))


def _check_minute_range(m: float) -> None:
    """
    Checks that the given value is in the range [0,60].  If the value
    is equal to 60, then a warning is raised.
    """
    if not 0.0 <= m < 60.0:
        if m != 60.0:
            raise IllegalMinuteError(m)
        warn(IllegalMinuteWarning(m, "Treating as 0 min, +1 hr/deg"))


def _check_second_range(sec: float) -> None:
    """
    Checks that the given value is in the range [0,60].  If the value
    is equal to 60, then a warning is raised.
    """
    if not 0.0 <= sec < 60.0:
        if sec != 60.0:
            raise IllegalSecondError(sec)
        warn(IllegalSecondWarning(sec, "Treating as 0 sec, +1 min"))


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


def _decimal_to_sexagesimal(a, /):
    """
    Convert a floating-point input to a 3 tuple
    - if input is in degrees, the result is (degree, arcminute, arcsecond)
    - if input is in hourangle, the result is (hour, minute, second)
    """
    sign = np.copysign(1.0, a)
    # assuming a in degree, these are (degree fraction, degree)
    (df, d) = np.modf(np.fabs(a))

    # assuming a in degree, these are (arcminute fraction, arcminute)
    (mf, m) = np.modf(df * 60.0)
    s = mf * 60.0

    return np.floor(sign * d), sign * np.floor(m), sign * s


def _decimal_to_sexagesimal_string(
    angle, precision=None, pad=False, sep=(":",), fields=3
):
    """
    Given a floating point angle, convert it to string
    """
    values = _decimal_to_sexagesimal(angle)
    # Check to see if values[0] is negative, using np.copysign to handle -0
    sign = np.copysign(1.0, values[0])
    # If the coordinates are negative, we need to take the absolute values.
    # We use np.abs because abs(-0) is -0
    # TODO: Is this true? (MHvK, 2018-02-01: not on my system)
    values = [np.abs(value) for value in values]

    if pad:
        pad = 3 if sign == -1 else 2
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
    rounding_thresh = 60.0 - (10.0 ** -(8 if precision is None else precision))

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

    literal = f"{np.copysign(values[0], sign):0{pad}.0f}{sep[0]}"
    if fields >= 2:
        literal += f"{int(values[1]):02d}{sep[1]}"
    if fields == 3:
        if precision is None:
            last_value = f"{abs(values[2]):.8f}".rstrip("0").rstrip(".")
        else:
            last_value = f"{abs(values[2]):.{precision}f}"
        if len(last_value) == 1 or last_value[1] == ".":
            last_value = "0" + last_value
        literal += f"{last_value}{sep[2]}"
    return literal
