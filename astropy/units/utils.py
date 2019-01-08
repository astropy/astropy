# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Miscellaneous utilities for `astropy.units`.

None of the functions in the module are meant for use outside of the
package.
"""


import numbers
import io
import re
from fractions import Fraction

import numpy as np
from numpy import finfo


_float_finfo = finfo(float)
# take float here to ensure comparison with another float is fast
# give a little margin since often multiple calculations happened
_JUST_BELOW_UNITY = float(1.-4.*_float_finfo.epsneg)
_JUST_ABOVE_UNITY = float(1.+4.*_float_finfo.eps)


def _get_first_sentence(s):
    """
    Get the first sentence from a string and remove any carriage
    returns.
    """

    x = re.match(r".*?\S\.\s", s)
    if x is not None:
        s = x.group(0)
    return s.replace('\n', ' ')


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
        represents = ''
        if isinstance(unit, core.Unit):
            represents = ":math:`{0}`".format(
                unit._represents.to_string('latex')[1:-1])
        aliases = ', '.join('``{0}``'.format(x) for x in unit.aliases)

        yield (unit, doc, represents, aliases, 'Yes' if unit.name in has_prefixes else 'No')


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

    docstring.write("""
.. list-table:: Available Units
   :header-rows: 1
   :widths: 10 20 20 20 1

   * - Unit
     - Description
     - Represents
     - Aliases
     - SI Prefixes
""")

    for unit_summary in _iter_unit_summary(namespace):
        docstring.write("""
   * - ``{0}``
     - {1}
     - {2}
     - {3}
     - {4}
""".format(*unit_summary))

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
    for nm, unit in namespace.items():
        if isinstance(unit, PrefixUnit):
            base_unit = unit.represents.bases[0]
            faux_namespace[base_unit.name] = base_unit

    docstring = io.StringIO()

    for unit_summary in _iter_unit_summary(faux_namespace):
        docstring.write("""
   * - Prefixes for ``{0}``
     - {1} prefixes
     - {2}
     - {3}
     - Only
""".format(*unit_summary))

    return docstring.getvalue()


def is_effectively_unity(value):
    # value is *almost* always real, except, e.g., for u.mag**0.5, when
    # it will be complex.  Use try/except to ensure normal case is fast
    try:
        return _JUST_BELOW_UNITY <= value <= _JUST_ABOVE_UNITY
    except TypeError:  # value is complex
        return (_JUST_BELOW_UNITY <= value.real <= _JUST_ABOVE_UNITY and
                _JUST_BELOW_UNITY <= value.imag + 1 <= _JUST_ABOVE_UNITY)


def sanitize_scale(scale):
    if is_effectively_unity(scale):
        return 1.0

    # Maximum speed for regular case where scale is a float.
    if scale.__class__ is float:
        return scale

    # All classes that scale can be (int, float, complex, Fraction)
    # have an "imag" attribute.
    if scale.imag:
        if abs(scale.real) > abs(scale.imag):
            if is_effectively_unity(scale.imag/scale.real + 1):
                return scale.real

        elif is_effectively_unity(scale.real/scale.imag + 1):
            return complex(0., scale.imag)

        return scale

    else:
        return scale.real


def validate_power(p, support_tuples=False):
    """Convert a power to a floating point value, an integer, or a Fraction.

    If a fractional power can be represented exactly as a floating point
    number, convert it to a float, to make the math much faster; otherwise,
    retain it as a `fractions.Fraction` object to avoid losing precision.
    Conversely, if the value is indistinguishable from a rational number with a
    low-numbered denominator, convert to a Fraction object.

    Parameters
    ----------
    p : float, int, Rational, Fraction
        Power to be converted
    """
    denom = getattr(p, 'denominator', None)
    if denom is None:
        try:
            p = float(p)
        except Exception:
            if not np.isscalar(p):
                raise ValueError("Quantities and Units may only be raised "
                                 "to a scalar power")
            else:
                raise

        if (p % 1.0) == 0.0:
            # Denominators of 1 can just be integers.
            p = int(p)
        elif (p * 8.0) % 1.0 == 0.0:
            # Leave alone if the denominator is exactly 2, 4 or 8, since this
            # can be perfectly represented as a float, which means subsequent
            # operations are much faster.
            pass
        else:
            # Convert floats indistinguishable from a rational to Fraction.
            # Here, we do not need to test values that are divisors of a higher
            # number, such as 3, since it is already addressed by 6.
            for i in (10, 9, 7, 6):
                scaled = p * float(i)
                if((scaled + 4. * _float_finfo.eps) % 1.0 <
                   8. * _float_finfo.eps):
                    p = Fraction(int(round(scaled)), i)
                    break

    elif denom == 1:
        p = int(p.numerator)

    elif (denom & (denom - 1)) == 0:
        # Above is a bit-twiddling hack to see if denom is a power of two.
        p = float(p)

    return p


def resolve_fractions(a, b):
    """
    If either input is a Fraction, convert the other to a Fraction.
    This ensures that any operation involving a Fraction will use
    rational arithmetic and preserve precision.
    """
    # We short-circuit on the most common cases of int and float, since
    # isinstance(a, Fraction) is very slow for any non-Fraction instances.
    a_is_fraction = (a.__class__ is not int and a.__class__ is not float and
                     isinstance(a, Fraction))
    b_is_fraction = (b.__class__ is not int and b.__class__ is not float and
                     isinstance(b, Fraction))
    if a_is_fraction and not b_is_fraction:
        b = Fraction(b)
    elif not a_is_fraction and b_is_fraction:
        a = Fraction(a)
    return a, b


def quantity_asanyarray(a, dtype=None):
    from .quantity import Quantity
    if not isinstance(a, np.ndarray) and not np.isscalar(a) and any(isinstance(x, Quantity) for x in a):
        return Quantity(a, dtype=dtype)
    else:
        return np.asanyarray(a, dtype=dtype)
