# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Utilities shared by the different formats.
"""

from collections.abc import Generator, Iterable, Sequence
from keyword import iskeyword


def get_non_keyword_units(
    bases: Iterable[str], prefixes: Sequence[str]
) -> Generator[tuple[str, str], None, None]:
    for base in bases:
        for prefix in prefixes:
            if not iskeyword(unit := prefix + base):
                yield unit, base
