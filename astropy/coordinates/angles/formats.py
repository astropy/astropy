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
from astropy.utils.compat.numpycompat import NUMPY_LT_2_1

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


def _check_hour_range_array(hrs):
    """Vectorised `_check_hour_range`, reporting at most once for the array."""
    outside = ~((hrs > -24.0) & (hrs < 24.0))
    if not outside.any():
        return
    illegal = outside & (np.abs(hrs) != 24.0)
    if illegal.any():
        raise IllegalHourError(hrs[illegal][0])
    warn(IllegalHourWarning(hrs[outside][0], "Treating as 24 hr"))


def _check_minute_range_array(m):
    """Vectorised `_check_minute_range`, reporting at most once for the array."""
    outside = ~((m >= 0.0) & (m < 60.0))
    if not outside.any():
        return
    illegal = outside & (m != 60.0)
    if illegal.any():
        raise IllegalMinuteError(m[illegal][0])
    warn(IllegalMinuteWarning(m[outside][0], "Treating as 0 min, +1 hr/deg"))


def _check_second_range_array(sec):
    """Vectorised `_check_second_range`, reporting at most once for the array."""
    outside = ~((sec >= 0.0) & (sec < 60.0))
    if not outside.any():
        return
    illegal = outside & (sec != 60.0)
    if illegal.any():
        raise IllegalSecondError(sec[illegal][0])
    warn(IllegalSecondWarning(sec[outside][0], "Treating as 0 sec, +1 min"))


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


# Designators the array parser understands, longest first so that "degrees" is
# replaced before "deg", and "deg" before "d".
_UNIT_DESIGNATORS = (
    "hours", "hour", "hr", "h",
    "degrees", "degree", "deg", "d",
    "minutes", "minute", "min", "m",
    "seconds", "second", "sec", "s",
)  # fmt: skip

# Non-ASCII spellings of a sign or a designator, mapped onto the ASCII ones.
_SYMBOL_REPLACEMENTS = (
    ("\N{MINUS SIGN}", "-"),
    ("\N{DEGREE SIGN}", "d"),
    ("\N{PRIME}", "m"),
    ("\N{DOUBLE PRIME}", "s"),
    ("'", "m"),
    ('"', "s"),
)

# A trailing direction letter is a sign, not a designator.
_DIRECTIONS = (("N", 1.0), ("S", -1.0), ("E", 1.0), ("W", -1.0))

# Below this many elements, parsing one string at a time is quicker: the array
# parser makes a fixed number of passes over the input whatever its length.
# Keeping small inputs on the old path also leaves scalar behaviour untouched.
_MIN_ARRAY_SIZE = 32


def _partition(values, sep):
    """Split each string on the first occurrence of ``sep``.

    Wrapper around `numpy.strings.partition`, which only exists from numpy 2.1.
    """
    if NUMPY_LT_2_1:
        parts = np.char.partition(values, sep)
        return parts[..., 0], parts[..., 1], parts[..., 2]
    return np.strings.partition(values, sep)


def _contains(values, sub):
    """Whether any element contains ``sub``.

    Used to skip replacements that would not change anything: scanning for a
    substring is cheaper than rewriting every string in the array.
    """
    return bool((np.strings.find(values, sub) >= 0).any())


def parse_angles(values, unit=None):
    """Parse an array of angle strings without looping over the elements.

    This is a fast path for the common sexagesimal spellings. Anything it does
    not recognise is reported back for `parse_angle` to handle instead, so the
    result is the same either way.

    Parameters
    ----------
    values : `~numpy.ndarray`
        Array of strings, of at least `_MIN_ARRAY_SIZE` elements.
    unit : `~astropy.units.UnitBase` or None
        Unit to assume for strings that do not name one themselves.

    Returns
    -------
    parsed : tuple or None
        ``(angles, unit, unparsed)``, where ``angles`` holds the parsed values
        in ``unit`` and ``unparsed`` flags the elements this parser could not
        handle, which the caller must fill in itself. `None` if the array as a
        whole is not something this parser can take on.
    """
    strings = np.strings.strip(values)

    for symbol, ascii_ in _SYMBOL_REPLACEMENTS:
        if _contains(strings, symbol):
            strings = np.strings.replace(strings, symbol, ascii_)

    # A trailing N/S/E/W flips the sign and is then removed, before the unit
    # designators are touched (S would otherwise look like a seconds marker).
    direction = None
    for letter, sign in _DIRECTIONS:
        if bool(np.strings.endswith(strings, letter).any()):
            if direction is None:
                direction = np.ones(strings.shape)
            found = np.strings.endswith(strings, letter)
            direction = np.where(found, sign, direction)
            strings = np.where(found, np.strings.rstrip(strings, letter), strings)

    # Each string names its own unit, so hours and degrees can be mixed. When
    # they are, the result takes the first element's unit, which is what
    # assembling the individual angles into one array does.
    is_hour = np.strings.find(strings, "h") >= 0
    is_degree = np.strings.find(strings, "d") >= 0
    if is_hour.all():
        found_unit = u.hourangle
    elif is_degree.all():
        found_unit = u.degree
    elif not is_hour.any() and not is_degree.any():
        if unit is None:
            # Nothing names a unit and none was supplied. Hand the whole array
            # back so the per-element parser raises, whatever the reason is.
            return None
        found_unit = unit
    else:
        found_unit = u.hourangle if is_hour[0] else u.degree

    # Spaces are separators here but are ignored by the grammar, so a string
    # mixing them with colons ("06 00:00") is a syntax error rather than
    # something to interpret. Note it now, before both become ":".
    mixed_separators = None
    if _contains(strings, " ") and _contains(strings, ":"):
        mixed_separators = (np.strings.count(strings, " ") > 0) & (
            np.strings.count(strings, ":") > 0
        )

    for designator in _UNIT_DESIGNATORS:
        if _contains(strings, designator):
            strings = np.strings.replace(strings, designator, ":")
    if _contains(strings, " "):
        strings = np.strings.replace(strings, " ", ":")
    strings = np.strings.rstrip(strings, ":")
    while _contains(strings, "::"):
        strings = np.strings.replace(strings, "::", ":")

    # The range checks below depend on how many fields were given, so count
    # them before the padding that follows.
    n_fields = np.strings.count(strings, ":") + 1

    first, _, rest = _partition(strings, ":")
    second, _, rest = _partition(rest, ":")
    third, _, trailing = _partition(rest, ":")

    def _is_number(field):
        digits = np.strings.replace(np.strings.replace(field, ".", ""), "-", "")
        return np.strings.isdigit(digits) | (np.strings.str_len(field) == 0)

    parsed = (
        _is_number(first)
        & _is_number(second)
        & _is_number(third)
        & (np.strings.str_len(trailing) == 0)
        & (np.strings.str_len(first) > 0)
    )
    if mixed_separators is not None:
        parsed &= ~mixed_separators

    angles = np.zeros(strings.shape)
    if not parsed.any():
        return angles, found_unit, ~parsed

    # numpy.loadtxt wants the same number of columns on every line.
    padded = np.where(
        n_fields == 1,
        np.strings.add(strings, ":0:0"),
        np.where(n_fields == 2, np.strings.add(strings, ":0"), strings),
    )
    fields = np.loadtxt(padded[parsed].tolist(), delimiter=":", ndmin=2)
    d, m, s = fields[:, 0], fields[:, 1], fields[:, 2]

    # `parse_angle` only range checks the fields it was actually given.
    given = n_fields[parsed]
    _check_hour_range_array(d[(given > 1) & is_hour[parsed]])
    _check_minute_range_array(m[given > 1])
    _check_second_range_array(s[given > 2])

    values_ = np.abs(d) + m / 60.0 + s / 3600.0
    values_ = np.copysign(values_, d)
    if direction is not None:
        values_ = values_ * direction[parsed]
    # Convert whichever elements did not name the unit the array ended up in.
    other = (is_degree if found_unit is u.hourangle else is_hour)[parsed]
    if other.any():
        other_unit = u.degree if found_unit is u.hourangle else u.hourangle
        converted = u.Quantity(values_, other_unit).to_value(found_unit)
        values_ = np.where(other, converted, values_)

    angles[parsed] = values_
    return angles, found_unit, ~parsed


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
    # Builtin round, not NumPy's, which disagrees with the formatting on ties.
    ndp = 8 if precision is None else precision

    if fields == 3:
        if round(float(values[2]), ndp) >= 60.0:
            values[2] = 0.0
            values[1] += 1.0
    elif fields == 2:
        if values[2] >= 30.0:
            values[1] += 1.0
    # Rounding the seconds into the minutes here too would round twice.
    elif values[1] + values[2] / 60.0 >= 30.0:
        values[0] += 1.0

    if fields >= 2 and values[1] >= 60.0:
        values[1] = 0.0
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
