# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A collection of different unit formats.

General usage is by their name in the |Unit| constructor or
in the :meth:`~astropy.units.UnitBase.to_string` method, i.e.,
these classes rarely if ever need to be imported directly.
"""

import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning

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
    "CDS",
    "FITS",
    "OGIP",
    "Base",
    "Console",
    "Generic",
    "Latex",
    "LatexInline",
    "Unicode",
    "VOUnit",
    "get_format",
]


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


def known_formats() -> str:
    return "Valid formatter names are: " + ", ".join(map(repr, Base.registry))


def known_parsers() -> str:
    return "Valid parser names are: " + ", ".join(
        repr(name)
        for name, cls in Base.registry.items()
        if cls.parse.__func__ is not Base.parse.__func__
    )


def get_format(format: str | type[Base] | None = None) -> type[Base]:
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
    if isinstance(format, str):
        try:
            return Base.registry[format.lower()]
        except KeyError:
            raise ValueError(f"Unknown format {format!r}.") from None
    if isinstance(format, type) and issubclass(format, Base):
        return format
    raise TypeError(f"Expected a formatter name, not {format!r}.")
