# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Utilities shared by the different formats.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


from ...extern import six
from ...utils.compat.fractions import Fraction


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
        `~.exceptions.ValueError`.

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
