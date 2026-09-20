# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains formatting functions that are for internal use in
astropy.coordinates.angles. Mainly they are conversions from one format
of data to another.
"""

import functools
import re
from warnings import warn

import numpy as np
from lark import Lark, Token, Transformer, v_args
from lark.exceptions import UnexpectedCharacters, UnexpectedInput

from astropy import units as u
from astropy.utils.parsing import make_parser, token_values

from .errors import (
    IllegalHourError,
    IllegalHourWarning,
    IllegalMinuteError,
    IllegalMinuteWarning,
    IllegalSecondError,
    IllegalSecondWarning,
)

# The regular expression for SIMPLE_UNIT is inserted when the parser is made,
# since it depends on the enabled units.
_GRAMMAR = r"""
angle: sign hms eastwest
     | sign dms dir
     | sign simple dir

sign: SIGN?

eastwest: EASTWEST?

dir: (EASTWEST | NORTHSOUTH)?

?ufloat: UFLOAT
       | UINT

generic: ufloat
       | UINT ufloat
       | UINT ":" ufloat
       | UINT UINT ufloat
       | UINT ":" UINT ":" ufloat

hms: UINT _HOUR
   | UINT _HOUR ufloat
   | UINT _HOUR UINT _MINUTE
   | UINT _HOUR UFLOAT _MINUTE
   | UINT _HOUR UINT _MINUTE ufloat
   | UINT _HOUR UINT _MINUTE ufloat _SECOND
   | generic _HOUR

dms: UINT _DEGREE
   | UINT _DEGREE ufloat
   | UINT _DEGREE UINT _MINUTE
   | UINT _DEGREE UFLOAT _MINUTE
   | UINT _DEGREE UINT _MINUTE ufloat
   | UINT _DEGREE UINT _MINUTE ufloat _SECOND
   | generic _DEGREE

simple: generic
      | generic _MINUTE       -> simple_arcmin
      | generic _SECOND       -> simple_arcsec
      | generic SIMPLE_UNIT   -> simple_unit

// The relative priorities of the terminals matter for the ones that can
// match the same text, since the first matching terminal is used.
// Several terminals include Unicode "MINUS SIGN" −.  It is important
// to include the hyphen last, or the regex will treat this as a range.
UFLOAT.9: /((\d+\.\d*)|(\.\d+))([eE][+-−]?\d+)?/
UINT.8: /\d+/
SIGN.7: /[+−-]/
EASTWEST.6: /[EW]$/
// We cannot use lower-case letters otherwise we'll confuse
// s[outh] with s[econd]
NORTHSOUTH.5: /[NS]$/
SIMPLE_UNIT.4: /__SIMPLE_UNIT__/
_MINUTE.3: /m(in(ute(s)?)?)?|′|'|ᵐ/
_SECOND.2: /s(ec(ond(s)?)?)?|″|"|ˢ/  // codespell:ignore ond
_DEGREE.1: /d(eg(ree(s)?)?)?|°/
_HOUR.1: /hour(s)?|h(r)?|ʰ/
%ignore " "
"""


@v_args(wrapper=token_values)
class _AngleTransformer(Transformer):
    """Turn a parsed angle string into a ``(value, unit)`` tuple.

    Terminal methods are applied as soon as a token is created and set its
    value, rule methods are applied when the rule is reduced.
    """

    def t_UFLOAT(self, t: Token) -> Token:
        t.value = float(t.replace("−", "-"))
        return t

    def t_UINT(self, t: Token) -> Token:
        t.value = int(t)
        return t

    def t_SIGN(self, t: Token) -> Token:
        t.value = 1.0 if t == "+" else -1.0
        return t

    def t_EASTWEST(self, t: Token) -> Token:
        t.value = -1.0 if t == "W" else 1.0
        return t

    def t_NORTHSOUTH(self, t: Token) -> Token:
        t.value = -1.0 if t == "S" else 1.0
        return t

    def t_SIMPLE_UNIT(self, t: Token) -> Token:
        t.value = u.Unit(t.value)
        return t

    def angle(self, sign, value_unit, direction):
        sign = sign * direction
        value, unit = value_unit
        if isinstance(value, tuple):
            return ((sign * value[0],) + value[1:], unit)
        else:
            return (sign * value, unit)

    def sign(self, sign=1.0):
        return sign

    def eastwest(self, direction=1.0):
        return direction

    dir = eastwest

    def generic(self, *values):
        return values[0] if len(values) == 1 else values

    def hms(self, *values):
        return (values[0] if len(values) == 1 else values, u.hourangle)

    def dms(self, *values):
        return (values[0] if len(values) == 1 else values, u.degree)

    def simple(self, value):
        return (value, None)

    def simple_arcmin(self, value):
        return (value, u.arcmin)

    def simple_arcsec(self, value):
        return (value, u.arcsec)

    def simple_unit(self, value, unit):
        return (value, unit)


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
    @functools.cache
    def _make_parser(cls) -> Lark:
        # TODO: in principle, the parser should be invalidated if we change unit
        # system (from CDS to FITS, say).  Might want to keep a link to the
        # unit_registry used, and regenerate the parser if it changes.
        # For some discussion of this problem, see
        # https://github.com/astropy/astropy/issues/5350#issuecomment-248770151
        simple_units = "|".join(
            f"(?:{re.escape(x)})" for x in cls._get_simple_unit_names()
        )
        grammar = _GRAMMAR.replace("__SIMPLE_UNIT__", simple_units)
        return make_parser(grammar, _AngleTransformer(), start="angle")

    def parse(self, angle, unit, debug=False):
        try:
            found_angle, found_unit = self._make_parser().parse(angle)
        except UnexpectedCharacters as e:
            raise ValueError(
                f"Invalid character at col {e.pos_in_stream} parsing angle {angle!r}"
            ) from None
        except UnexpectedInput:
            raise ValueError(f"syntax error parsing angle {angle!r}") from None
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
        No longer has any effect; kept for backwards compatibility.

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
