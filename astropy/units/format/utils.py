# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Utilities shared by the different formats.
"""

from __future__ import annotations

from keyword import iskeyword
from typing import TYPE_CHECKING

from astropy.units.utils import maybe_simple_fraction

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Sequence

    from astropy.units.typing import Real


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
