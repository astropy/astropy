# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A collection of different unit formats.
"""


# This is pretty atrocious, but it will prevent a circular import for those
# formatters that need access to the units.core module An entry for it should
# exist in sys.modules since astropy.units.core imports this module
import sys
core = sys.modules['astropy.units.core']

from .base import Base
from .generic import Generic, Unscaled
from .cds import CDS
from .console import Console
from .fits import Fits
from .latex import Latex, LatexInline
from .ogip import OGIP
from .unicode_format import Unicode
from .vounit import VOUnit


__all__ = [
    'Base', 'Generic', 'CDS', 'Console', 'Fits', 'Latex', 'LatexInline',
    'OGIP', 'Unicode', 'Unscaled', 'VOUnit', 'get_format']


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
    if isinstance(format, type) and issubclass(format, Base):
        return format
    elif not (isinstance(format, str) or format is None):
        raise TypeError(
            "Formatter must a subclass or instance of a subclass of {!r} "
            "or a string giving the name of the formatter.  Valid formatter "
            "names are: [{}]".format(Base, ', '.join(Base.registry)))

    if format is None:
        format = 'generic'

    format_lower = format.lower()

    if format_lower in Base.registry:
        return Base.registry[format_lower]

    raise ValueError("Unknown format {!r}.  Valid formatter names are: "
                     "[{}]".format(format, ', '.join(Base.registry)))
