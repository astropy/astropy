# Licensed under a 3-clause BSD style license - see LICENSE.rst

from enum import StrEnum, auto
from typing import Self


class DeprecatedUnitAction(StrEnum):
    SILENT = auto()
    WARN = auto()
    RAISE = auto()
    CONVERT = auto()

    @classmethod
    def _missing_(cls, value) -> Self | None:
        raise ValueError(
            f"invalid deprecation handling option: {value!r}. Valid options are "
            f"{', '.join(repr(opt.value) for opt in cls)}."
        )
