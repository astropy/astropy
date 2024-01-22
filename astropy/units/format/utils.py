# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Utilities shared by the different formats.
"""


import warnings

from astropy.units.utils import maybe_simple_fraction
from astropy.utils.misc import did_you_mean


def get_grouped_by_powers(bases, powers):
    """
    Groups the powers and bases in the given
    `~astropy.units.CompositeUnit` into positive powers and
    negative powers for easy display on either side of a solidus.

    Parameters
    ----------
    bases : list of `astropy.units.UnitBase` instances

    powers : list of int

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


def split_mantissa_exponent(v, format_spec=".8g"):
    """
    Given a number, split it into its mantissa and base 10 exponent
    parts, each as strings.  If the exponent is too small, it may be
    returned as the empty string.

    Parameters
    ----------
    v : float

    format_spec : str, optional
        Number representation formatting string

    Returns
    -------
    mantissa, exponent : tuple of strings
    """
    x = format(v, format_spec).split("e")

    if len(x) == 2:
        ex = x[1].lstrip("0+")
        if len(ex) > 0 and ex[0] == "-":
            ex = "-" + ex[1:].lstrip("0")
    else:
        ex = ""

    if ex == "" or (x[0] != "1." + "0" * (len(x[0]) - 2)):
        m = x[0]
    else:
        m = ""

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
    from astropy.units import core

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
    else:
        raise TypeError(
            f"unit argument must be a 'NamedUnit' or 'CompositeUnit', not {type(unit)}"
        )


def format_power(power):
    """
    Converts a value for a power (which may be floating point or a
    `fractions.Fraction` object), into a string looking like either
    an integer or a fraction, if the power is close to that.
    """
    if not hasattr(power, "denominator"):
        power = maybe_simple_fraction(power)
        if getattr(power, "denonimator", None) == 1:
            power = power.numerator

    return str(power)


def _try_decomposed(unit, format_decomposed):
    represents = getattr(unit, "_represents", None)
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
            results = [x + " (deprecated)"]
            decomposed = _try_decomposed(all_units[x], format_decomposed)
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
    from astropy.units.core import UnitsWarning

    message = f"The unit '{s}' has been deprecated in the {standard_name} standard."
    decomposed = _try_decomposed(unit, format_decomposed)
    if decomposed is not None:
        message += f" Suggested: {decomposed}."
    warnings.warn(message, UnitsWarning)
