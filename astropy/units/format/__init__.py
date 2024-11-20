# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A collection of different unit formats.

General usage is by their name in the |Unit| constructor or
in the :meth:`~astropy.units.UnitBase.to_string` method, i.e.,
these classes rarely if ever need to be imported directly.
"""

from __future__ import annotations

import sys
import warnings
from abc import abstractmethod
from typing import TYPE_CHECKING, Protocol, runtime_checkable

from astropy.utils.exceptions import AstropyDeprecationWarning

if TYPE_CHECKING:
    from astropy.units import UnitBase

# This is pretty atrocious, but it will prevent a circular import for those
# formatters that need access to the units.core module An entry for it should
# exist in sys.modules since astropy.units.core imports this module
core = sys.modules["astropy.units.core"]

from .base import Base
from .cds import CDS
from .console import Console
from .fits import FITS
from .generic import Generic
from .latex import Latex, LatexInline
from .ogip import OGIP
from .unicode_format import Unicode
from .vounit import VOUnit

__all__ = [
    "Base",
    "Generic",
    "CDS",
    "Console",
    "FITS",
    "Latex",
    "LatexInline",
    "OGIP",
    "Unicode",
    "VOUnit",
    "get_format",
]


@runtime_checkable
class UnitFormatter(Protocol):
    @classmethod
    @abstractmethod
    def to_string(cls, unit: UnitBase) -> str:
        raise NotImplementedError


@runtime_checkable
class ParsingUnitFormatter(UnitFormatter, Protocol):
    @classmethod
    @abstractmethod
    def parse(cls, s: str) -> UnitBase:
        raise NotImplementedError


def __getattr__(name):
    if name == "Fits":
        warnings.warn(
            AstropyDeprecationWarning(
                'The class "Fits" has been renamed to "FITS" in version 7.0. The old '
                "name is deprecated and may be removed in a future version.\n"
                "        Use FITS instead."
            )
        )
        return FITS
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def _known_formats():
    in_out = []
    out_only = []
    for name, formatter in Base.registry.items():
        if issubclass(formatter, ParsingUnitFormatter):
            in_out.append(name)
        else:
            out_only.append(name)
    return (
        f"Valid formatter names are: {in_out} for input and output, "
        f"and {out_only} for output only."
    )


def get_format(format=None):
    """
    Get a formatter by name.

    Parameters
    ----------
    format : str or `astropy.units.format.Base` subclass
        The name of the format, or the formatter class itself.

    Returns
    -------
    format : `astropy.units.format.Base` subclass
        The requested formatter.
    """
    if format is None:
        return Generic

    if isinstance(format, type) and issubclass(format, UnitFormatter):
        return format
    elif not (isinstance(format, str) or format is None):
        raise TypeError(
            f"Expected a formatter name, not {format!r}.  {_known_formats()}."
        )

    format_lower = format.lower()

    if format_lower in Base.registry:
        return Base.registry[format_lower]

    raise ValueError(f"Unknown format {format!r}.  {_known_formats()}")
