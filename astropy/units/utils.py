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

import numpy as np
from numpy import finfo

from ..extern import six
from ..utils.compat.fractions import Fraction


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
    return _JUST_BELOW_UNITY <= value <= _JUST_ABOVE_UNITY


def validate_power(p):
    """
    Handles the conversion of a power to a floating point or a
    rational number.
    """
    # For convenience, treat tuples as Fractions
    if isinstance(p, tuple) and len(p) == 2:
        p = Fraction(p[0], p[1])

    if isinstance(p, numbers.Rational) or isinstance(p, Fraction):
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

        if p % 1.0 == 0.0:
            p = int(p)

    return p
