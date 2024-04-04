# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Miscellaneous utilities for `astropy.units`.

None of the functions in the module are meant for use outside of the
package.
"""

import io
import re
from fractions import Fraction

import numpy as np
from numpy import finfo

_float_finfo = finfo(float)
# take float here to ensure comparison with another float is fast
# give a little margin since often multiple calculations happened
_JUST_BELOW_UNITY = float(1.0 - 4.0 * _float_finfo.epsneg)
_JUST_ABOVE_UNITY = float(1.0 + 4.0 * _float_finfo.eps)


def _get_first_sentence(s):
    """
    Get the first sentence from a string and remove any carriage
    returns.
    """
    x = re.match(r".*?\S\.\s", s)
    if x is not None:
        s = x.group(0)
    return s.replace("\n", " ")


def _iter_unit_summary(namespace):
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


def generate_unit_summary(namespace):
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


def generate_prefixonly_unit_summary(namespace):
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
    from . import PrefixUnit

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


def is_effectively_unity(value):
    # value is *almost* always real, except, e.g., for u.mag**0.5, when
    # it will be complex.  Use try/except to ensure normal case is fast
    try:
        return _JUST_BELOW_UNITY <= value <= _JUST_ABOVE_UNITY
    except TypeError:  # value is complex
        return (
            _JUST_BELOW_UNITY <= value.real <= _JUST_ABOVE_UNITY
            and _JUST_BELOW_UNITY <= value.imag + 1 <= _JUST_ABOVE_UNITY
        )


def sanitize_scale(scale):
    if is_effectively_unity(scale):
        return 1.0

    # Maximum speed for regular case where scale is a float.
    if scale.__class__ is float:
        return scale

    # We cannot have numpy scalars, since they don't autoconvert to
    # complex if necessary.  They are also slower.
    if hasattr(scale, "dtype"):
        scale = scale.item()

    # All classes that scale can be (int, float, complex, Fraction)
    # have an "imag" attribute.
    if scale.imag:
        if abs(scale.real) > abs(scale.imag):
            if is_effectively_unity(scale.imag / scale.real + 1):
                return scale.real

        elif is_effectively_unity(scale.real / scale.imag + 1):
            return complex(0.0, scale.imag)

        return scale

    else:
        return scale.real


def maybe_simple_fraction(p, max_denominator=100):
    """Fraction very close to x with denominator at most max_denominator.

    The fraction has to be such that fraction/x is unity to within 4 ulp.
    If such a fraction does not exist, returns the float number.

    The algorithm is that of `fractions.Fraction.limit_denominator`, but
    sped up by not creating a fraction to start with.

    If the input is zero, an integer or `fractions.Fraction`, just return it.
    """
    if p == 0 or p.__class__ is int or p.__class__ is Fraction:
        return p
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

    return p


def validate_power(p):
    """Check that a power can be converted to a floating point value.

    Parameters
    ----------
    p : numerical
        Power to be converted

    Raises
    ------
    ValueError
        If the power is an array in which not all elements are equal.

    Returns
    -------
    p : numerical
        Equals the input unless the input was iterable and all elements
        were the same, in which case it returns the first item.
    """
    if p.__class__ is int or p.__class__ is Fraction:
        return p
    try:
        float(p)
    except Exception:
        p = np.asanyarray(p)
        if ((first := p.flat[0]) == p).all():
            # All the same, now check it is OK.
            float(first)
            return first
        else:
            raise ValueError(
                "Quantities and Units may only be raised to a scalar power"
            ) from None
    else:
        return p


def sanitize_power(p):
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
        p = p.numerator

    elif (denom & (denom - 1)) == 0:
        # Above is a bit-twiddling hack to see if denom is a power of two.
        # If so, float does not lose precision and will speed things up.
        p = float(p)

    return p


def resolve_fractions(a, b):
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


def quantity_asanyarray(a, dtype=None):
    from .quantity import Quantity

    if (
        not isinstance(a, np.ndarray)
        and not np.isscalar(a)
        and any(isinstance(x, Quantity) for x in a)
    ):
        return Quantity(a, dtype=dtype)
    else:
        # skip over some dtype deprecation.
        dtype = np.float64 if dtype is np.inexact else dtype
        return np.asanyarray(a, dtype=dtype)
