# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Utilities shared by the different formats.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import warnings
from fractions import Fraction

from ...extern import six
from ...extern.six.moves import zip
from ...utils.misc import did_you_mean


def get_grouped_by_powers(bases, powers):
    """
    Groups the powers and bases in the given
    `~astropy.units.CompositeUnit` into positive powers and
    negative powers for easy display on either side of a solidus.

    Parameters
    ----------
    bases : list of `astropy.units.UnitBase` instances

    powers : list of ints

    Returns
    -------
    positives, negatives : tuple of lists
       Each element in each list is tuple of the form (*base*,
       *power*).  The negatives have the sign of their power reversed
       (i.e. the powers are all positive).
    """
    positive = []
    negative = []
    for base, power in zip(bases, powers):
        if power < 0:
            negative.append((base, -power))
        elif power > 0:
            positive.append((base, power))
        else:
            raise ValueError("Unit with 0 power")
    return positive, negative


def split_mantissa_exponent(v):
    """
    Given a number, split it into its mantissa and base 10 exponent
    parts, each as strings.  If the exponent is too small, it may be
    returned as the empty string.

    The precise rules are based on Python's "general purpose" (`g`)
    formatting.

    Parameters
    ----------
    v : float

    Returns
    -------
    mantissa, exponent : tuple of strings
    """
    x = "{0:.8g}".format(v).split('e')
    if x[0] != '1.' + '0' * (len(x[0]) - 2):
        m = x[0]
    else:
        m = ''

    if len(x) == 2:
        ex = x[1].lstrip("0+")
        if len(ex) > 0 and ex[0] == '-':
            ex = '-' + ex[1:].lstrip('0')
    else:
        ex = ''

    return m, ex


def decompose_to_known_units(unit, func):
    """
    Partially decomposes a unit so it is only composed of units that
    are "known" to a given format.

    Parameters
    ----------
    unit : `~astropy.units.UnitBase` instance

    func : callable
        This function will be called to determine if a given unit is
        "known".  If the unit is not known, this function should raise a
        `ValueError`.

    Returns
    -------
    unit : `~astropy.units.UnitBase` instance
        A flattened unit.
    """
    from .. import core
    if isinstance(unit, core.CompositeUnit):
        new_unit = core.Unit(unit.scale)
        for base, power in zip(unit.bases, unit.powers):
            new_unit = new_unit * decompose_to_known_units(base, func) ** power
        return new_unit
    elif isinstance(unit, core.NamedUnit):
        try:
            func(unit)
        except ValueError:
            if isinstance(unit, core.Unit):
                return decompose_to_known_units(unit._represents, func)
            raise
        return unit


def format_power(power):
    """
    Converts a value for a power (which may be floating point or a
    `fractions.Fraction` object), into a string either looking like
    an integer or a fraction.
    """
    if not isinstance(power, Fraction):
        if power % 1.0 != 0.0:
            frac = Fraction.from_float(power)
            power = frac.limit_denominator(10)
            if power.denominator == 1:
                power = int(power.numerator)
        else:
            power = int(power)
    return six.text_type(power)


def _try_decomposed(unit, format_decomposed):
    represents = getattr(unit, '_represents', None)
    if represents is not None:
        try:
            represents_string = format_decomposed(represents)
        except ValueError:
            pass
        else:
            return represents_string

    decomposed = unit.decompose()
    if decomposed is not unit:
        try:
            decompose_string = format_decomposed(decomposed)
        except ValueError:
            pass
        else:
            return decompose_string

    return None


def did_you_mean_units(s, all_units, deprecated_units, format_decomposed):
    """
    A wrapper around `astropy.utils.misc.did_you_mean` that deals with
    the display of deprecated units.

    Parameters
    ----------
    s : str
        The invalid unit string

    all_units : dict
        A mapping from valid unit names to unit objects.

    deprecated_units : sequence
        The deprecated unit names

    format_decomposed : callable
        A function to turn a decomposed version of the unit into a
        string.  Should return `None` if not possible

    Returns
    -------
    msg : str
        A string message with a list of alternatives, or the empty
        string.
    """
    def fix_deprecated(x):
        if x in deprecated_units:
            results = [x + ' (deprecated)']
            decomposed = _try_decomposed(
                all_units[x], format_decomposed)
            if decomposed is not None:
                results.append(decomposed)
            return results
        return (x,)

    return did_you_mean(s, all_units, fix=fix_deprecated)


def unit_deprecation_warning(s, unit, standard_name, format_decomposed):
    """
    Raises a UnitsWarning about a deprecated unit in a given format.
    Suggests a decomposed alternative if one is available.

    Parameters
    ----------
    s : str
        The deprecated unit name.

    unit : astropy.units.core.UnitBase
        The unit object.

    standard_name : str
        The name of the format for which the unit is deprecated.

    format_decomposed : callable
        A function to turn a decomposed version of the unit into a
        string.  Should return `None` if not possible
    """
    from ..core import UnitsWarning

    message = "The unit '{0}' has been deprecated in the {1} standard.".format(
        s, standard_name)
    decomposed = _try_decomposed(unit, format_decomposed)
    if decomposed is not None:
        message += " Suggested: {0}.".format(decomposed)
    warnings.warn(message, UnitsWarning)
