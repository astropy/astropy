# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Miscellaneous utilities for `astropy.units`.

None of the functions in the module are meant for use outside of the
package.
"""

from __future__ import annotations

import io
import re
from fractions import Fraction
from typing import TYPE_CHECKING, SupportsFloat

import numpy as np
from numpy import finfo

from .errors import UnitScaleError

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping
    from typing import Literal

    from .core import NamedUnit
    from .typing import UnitPower, UnitPowerLike, UnitScale, UnitScaleLike


_float_finfo = finfo(float)
# take float here to ensure comparison with another float is fast
# give a little margin since often multiple calculations happened
_JUST_BELOW_UNITY = float(1.0 - 4.0 * _float_finfo.epsneg)
_JUST_ABOVE_UNITY = float(1.0 + 4.0 * _float_finfo.eps)


def _get_first_sentence(s: str) -> str:
    """
    Get the first sentence from a string and remove any carriage
    returns.
    """
    x = re.match(r".*?\S\.\s", s)
    if x is not None:
        s = x.group(0)
    return s.replace("\n", " ")


def _iter_unit_summary(
    namespace: Mapping[str, object],
) -> Generator[tuple[NamedUnit, str, str, str, Literal["Yes", "No"]], None, None]:
    """
    Generates the ``(unit, doc, represents, aliases, prefixes)``
    tuple used to format the unit summary docs in `generate_unit_summary`.
    """
    from . import core

    # Get all of the units, and keep track of which ones have SI
    # prefixes
    units = []
    has_prefixes = set()
    for key, val in namespace.items():
        # Skip non-unit items
        if not isinstance(val, core.UnitBase):
            continue

        if not isinstance(val, core.NamedUnit):
            raise TypeError(f"{key!r} must be defined with 'def_unit()'")

        # Skip aliases
        if key != val.name:
            continue

        if isinstance(val, core.PrefixUnit):
            # This will return the root unit that is scaled by the prefix
            # attached to it
            has_prefixes.add(val._represents.bases[0].name)
        else:
            units.append(val)

    # Sort alphabetically, case insensitive
    units.sort(key=lambda x: x.name.lower())

    for unit in units:
        doc = _get_first_sentence(unit.__doc__).strip()
        represents = ""
        if isinstance(unit, core.Unit):
            represents = f":math:`{unit._represents.to_string('latex')[1:-1]}`"
        aliases = ", ".join(f"``{x}``" for x in unit.aliases)

        yield (
            unit,
            doc,
            represents,
            aliases,
            "Yes" if unit.name in has_prefixes else "No",
        )


def generate_unit_summary(namespace: Mapping[str, object]) -> str:
    """
    Generates a summary of units from a given namespace.  This is used
    to generate the docstring for the modules that define the actual
    units.

    Parameters
    ----------
    namespace : dict
        A namespace containing units.

    Returns
    -------
    docstring : str
        A docstring containing a summary table of the units.
    """
    docstring = io.StringIO()

    docstring.write(
        """
.. list-table:: Available Units
   :header-rows: 1
   :widths: 10 20 20 20 1

   * - Unit
     - Description
     - Represents
     - Aliases
     - SI Prefixes
"""
    )
    template = """
   * - ``{}``
     - {}
     - {}
     - {}
     - {}
"""
    for unit_summary in _iter_unit_summary(namespace):
        docstring.write(template.format(*unit_summary))

    return docstring.getvalue()


def generate_prefixonly_unit_summary(namespace: Mapping[str, object]) -> str:
    """
    Generates table entries for units in a namespace that are just prefixes
    without the base unit.  Note that this is intended to be used *after*
    `generate_unit_summary` and therefore does not include the table header.

    Parameters
    ----------
    namespace : dict
        A namespace containing units that are prefixes but do *not* have the
        base unit in their namespace.

    Returns
    -------
    docstring : str
        A docstring containing a summary table of the units.
    """
    from .core import PrefixUnit

    faux_namespace = {}
    for unit in namespace.values():
        if isinstance(unit, PrefixUnit):
            base_unit = unit.represents.bases[0]
            faux_namespace[base_unit.name] = base_unit

    docstring = io.StringIO()
    template = """
   * - Prefixes for ``{}``
     - {} prefixes
     - {}
     - {}
     - Only
"""
    for unit_summary in _iter_unit_summary(faux_namespace):
        docstring.write(template.format(*unit_summary))

    return docstring.getvalue()


def is_effectively_unity(value: UnitScaleLike) -> bool | np.bool_:
    # value is *almost* always real, except, e.g., for u.mag**0.5, when
    # it will be complex.  Use try/except to ensure normal case is fast
    try:
        return _JUST_BELOW_UNITY <= value <= _JUST_ABOVE_UNITY
    except TypeError:  # value is complex
        return (
            _JUST_BELOW_UNITY <= value.real <= _JUST_ABOVE_UNITY
            and _JUST_BELOW_UNITY <= value.imag + 1 <= _JUST_ABOVE_UNITY
        )


def sanitize_scale(scale: UnitScaleLike) -> UnitScale:
    if is_effectively_unity(scale):
        return 1.0
    if not scale:
        raise UnitScaleError("cannot create a unit with a scale of 0.")
    if type(scale) is float:  # float is very common, so handle it fast
        return scale
    if isinstance(scale, SupportsFloat):
        return float(scale)

    if abs(scale.real) > abs(scale.imag):
        if is_effectively_unity(scale.imag / scale.real + 1):
            return float(scale.real)
    elif is_effectively_unity(scale.real / scale.imag + 1):
        return complex(0.0, scale.imag)
    return complex(scale)


def maybe_simple_fraction(p: UnitPowerLike, max_denominator: int = 100) -> UnitPower:
    """Fraction very close to x with denominator at most max_denominator.

    The fraction has to be such that fraction/x is unity to within 4 ulp.
    If such a fraction does not exist, returns the float number.

    The algorithm is that of `fractions.Fraction.limit_denominator`, but
    sped up by not creating a fraction to start with.

    If the input is zero, an integer or `fractions.Fraction`, just return it.
    """
    if p.__class__ is int or p.__class__ is Fraction:
        return p
    if p == 0:
        return 0  # p might be some numpy number, but we want a Python int
    n, d = float(p).as_integer_ratio()
    a = n // d
    # Normally, start with 0,1 and 1,0; here we have applied first iteration.
    n0, d0 = 1, 0
    n1, d1 = a, 1
    while d1 <= max_denominator:
        if _JUST_BELOW_UNITY <= n1 / (d1 * p) <= _JUST_ABOVE_UNITY:
            return Fraction(n1, d1)
        n, d = d, n - a * d
        a = n // d
        n0, n1 = n1, n0 + a * n1
        d0, d1 = d1, d0 + a * d1

    return float(p)


def sanitize_power(p: UnitPowerLike) -> UnitPower:
    """Convert the power to a float, an integer, or a Fraction.

    If a fractional power can be represented exactly as a floating point
    number, convert it to a float, to make the math much faster; otherwise,
    retain it as a `fractions.Fraction` object to avoid losing precision.
    Conversely, if the value is indistinguishable from a rational number with a
    low-numbered denominator, convert to a Fraction object.
    If a power can be represented as an integer, use that.

    Parameters
    ----------
    p : float, int, Rational, Fraction
        Power to be converted.
    """
    if p.__class__ is int:
        return p

    denom = getattr(p, "denominator", None)
    if denom is None:
        # This returns either a (simple) Fraction or the same float.
        p = maybe_simple_fraction(p)
        # If still a float, nothing more to be done.
        if isinstance(p, float):
            return p

        # Otherwise, check for simplifications.
        denom = p.denominator

    if denom == 1:
        return int(p.numerator)

    elif (denom & (denom - 1)) == 0:
        # Above is a bit-twiddling hack to see if denom is a power of two.
        # If so, float does not lose precision and will speed things up.
        p = float(p)

    return p


def resolve_fractions(
    a: UnitPowerLike, b: UnitPowerLike
) -> tuple[UnitPowerLike, UnitPowerLike]:
    """
    If either input is a Fraction, convert the other to a Fraction
    (at least if it does not have a ridiculous denominator).
    This ensures that any operation involving a Fraction will use
    rational arithmetic and preserve precision.
    """
    # We short-circuit on the most common cases of int and float, since
    # isinstance(a, Fraction) is very slow for any non-Fraction instances.
    a_is_fraction = (
        a.__class__ is not int and a.__class__ is not float and isinstance(a, Fraction)
    )
    b_is_fraction = (
        b.__class__ is not int and b.__class__ is not float and isinstance(b, Fraction)
    )
    if a_is_fraction and not b_is_fraction:
        b = maybe_simple_fraction(b)
    elif not a_is_fraction and b_is_fraction:
        a = maybe_simple_fraction(a)
    return a, b
