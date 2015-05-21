# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Miscellaneous utilities for `astropy.units`.

None of the functions in the module are meant for use outside of the
package.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numbers
import io
import re
import warnings

import numpy as np
from numpy import finfo

from ..extern import six
from ..utils.compat.fractions import Fraction
from ..utils.exceptions import AstropyDeprecationWarning

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

    x = re.match(".*?\S\.\s", s)
    if x is not None:
        s = x.group(0)
    return s.replace('\n', ' ')


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

    from . import core

    # Get all of the units, and keep track of which ones have SI
    # prefixes
    units = []
    has_prefixes = set()
    for key, val in list(six.iteritems(namespace)):
        # Skip non-unit items
        if not isinstance(val, core.UnitBase):
            continue

        # Skip aliases
        if key != val.name:
            continue

        if isinstance(val, core.PrefixUnit):
            decomposed = val.decompose()
            if len(decomposed.bases):
                has_prefixes.add(val.decompose().bases[0].name)
            else:
                has_prefixes.add('dimensionless')

        else:
            units.append(val)

    # Sort alphabetically, case insensitive
    units.sort(key=lambda x: x.name.lower())

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

    for unit in units:
        if unit.name in has_prefixes:
            unit_has_prefixes = 'Y'
        else:
            unit_has_prefixes = 'N'
        doc = _get_first_sentence(unit.__doc__).strip()
        represents = ''
        if isinstance(unit, core.Unit):
            represents = ":math:`{0}`".format(
                unit._represents.to_string('latex')[1:-1])
        aliases = ', '.join('``{0}``'.format(x) for x in unit.aliases)
        docstring.write("""
   * - ``{0}``
     - {1}
     - {2}
     - {3}
     - {4}
""".format(unit, doc, represents, aliases, unit_has_prefixes))

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

    if np.iscomplex(scale):  # scale is complex
        if scale == 0.0:
            return 0.0

        if abs(scale.real) > abs(scale.imag):
            if is_effectively_unity(scale.imag/scale.real + 1):
                scale = scale.real
        else:
            if is_effectively_unity(scale.real/scale.imag + 1):
                scale = complex(0., scale.imag)

    return scale


def validate_power(p, support_tuples=False):
    """
    Handles the conversion of a power to a floating point or a
    rational number.

    Parameters
    ----------
    support_tuples : bool, optional
        If `True`, treat 2-tuples as `Fraction` objects.  This
        behavior is deprecated and will be removed in astropy 0.5.
    """
    # For convenience, treat tuples as Fractions
    if support_tuples and isinstance(p, tuple) and len(p) == 2:
        # Deprecated in 0.3.1
        warnings.warn(
            "Using a tuple as a fractional power is deprecated and may be "
            "removed in a future version.  Use Fraction(n, d) instead.",
            AstropyDeprecationWarning)
        p = Fraction(p[0], p[1])

    if isinstance(p, (numbers.Rational, Fraction)):
        # If the fractional power can be represented *exactly* as a
        # floating point number, we convert it to a float, to make the
        # math much faster, otherwise, we retain it as a
        # `fractions.Fraction` object to avoid losing precision.
        denom = p.denominator
        if denom == 1:
            p = int(p.numerator)
        # This is bit-twiddling hack to see if the integer is a
        # power of two
        elif (denom & (denom - 1)) == 0:
            p = float(p)
    else:
        if not np.isscalar(p):
            raise ValueError(
                "Quantities and Units may only be raised to a scalar power")

        p = float(p)

        # If the value is indistinguishable from a rational number
        # with a low-numbered denominator, convert to a Fraction
        # object.  We don't want to convert for denominators that are
        # a power of 2, since those can be perfectly represented, and
        # subsequent operations are much faster if they are retained
        # as floats.  Nor do we need to test values that are divisors
        # of a higher number, such as 3, since it is already addressed
        # by 6.

        # First check for denominator of 1
        if (p % 1.0) == 0.0:
            p = int(p)
        # Leave alone if the denominator is exactly 2, 4 or 8
        elif (p * 8.0) % 1.0 == 0.0:
            pass
        else:
            for i in [10, 9, 7, 6]:
                scaled = p * float(i)
                if((scaled + 4. * _float_finfo.eps) % 1.0 <
                   8. * _float_finfo.eps):
                    p = Fraction(int(round(scaled)), i)
                    break

    return p


def resolve_fractions(a, b):
    """
    If either input is a Fraction, convert the other to a Fraction.
    This ensures that any operation involving a Fraction will use
    rational arithmetic and preserve precision.
    """
    a_is_fraction = isinstance(a, Fraction)
    b_is_fraction = isinstance(b, Fraction)
    if a_is_fraction and not b_is_fraction:
        b = Fraction(b)
    elif not a_is_fraction and b_is_fraction:
        a = Fraction(a)
    return a, b
