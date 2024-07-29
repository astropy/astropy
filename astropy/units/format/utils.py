# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Utilities shared by the different formats.
"""

from __future__ import annotations

from keyword import iskeyword
from typing import TYPE_CHECKING

from astropy.units.utils import maybe_simple_fraction

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterable, Sequence
    from typing import TypeVar

    from astropy.units import UnitBase
    from astropy.units.typing import Real

    T = TypeVar("T")


def get_grouped_by_powers(
    bases: Sequence[T], powers: Sequence[int]
) -> tuple[list[tuple[T, int]], list[tuple[T, int]]]:
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
    positives, negatives : list of tuple of |Unit| and power
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


def decompose_to_known_units(
    unit: UnitBase, func: Callable[[UnitBase], None]
) -> UnitBase:
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
        for base, power in unit.components():
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


def format_power(power: Real) -> str:
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


def get_non_keyword_units(
    bases: Iterable[str], prefixes: Sequence[str]
) -> Generator[tuple[str, str], None, None]:
    for base in bases:
        for prefix in prefixes:
            if not iskeyword(unit := prefix + base):
                yield unit, base
