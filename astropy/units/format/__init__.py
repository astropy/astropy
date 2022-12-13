# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A collection of different unit formats.
"""


# This is pretty atrocious, but it will prevent a circular import for those
# formatters that need access to the units.core module An entry for it should
# exist in sys.modules since astropy.units.core imports this module
import sys

core = sys.modules["astropy.units.core"]

from .base import Base
from .cds import CDS
from .console import Console
from .fits import Fits
from .generic import Generic, Unscaled
from .latex import Latex, LatexInline
from .ogip import OGIP
from .unicode_format import Unicode
from .vounit import VOUnit

__all__ = [
    "Base",
    "Generic",
    "CDS",
    "Console",
    "Fits",
    "Latex",
    "LatexInline",
    "OGIP",
    "Unicode",
    "Unscaled",
    "VOUnit",
    "get_format",
]


def _known_formats():
    inout = [
        name
        for name, cls in Base.registry.items()
        if cls.parse.__func__ is not Base.parse.__func__
    ]
    out_only = [
        name
        for name, cls in Base.registry.items()
        if cls.parse.__func__ is Base.parse.__func__
    ]
    return (
        f"Valid formatter names are: {inout} for input and output, "
        f"and {out_only} for output only."
    )


def get_format(format=None):
    """
    Get a formatter by name.

    Parameters
    ----------
    format : str or `astropy.units.format.Base` instance or subclass
        The name of the format, or the format instance or subclass
        itself.

    Returns
    -------
    format : `astropy.units.format.Base` instance
        The requested formatter.
    """
    if format is None:
        return Generic

    if isinstance(format, type) and issubclass(format, Base):
        return format
    elif not (isinstance(format, str) or format is None):
        raise TypeError(
            f"Formatter must a subclass or instance of a subclass of {Base!r} "
            f"or a string giving the name of the formatter. {_known_formats()}."
        )

    format_lower = format.lower()

    if format_lower in Base.registry:
        return Base.registry[format_lower]

    raise ValueError(f"Unknown format {format!r}.  {_known_formats()}")
