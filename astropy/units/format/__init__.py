# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A collection of different unit formats.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from inspect import isclass

from .base import Base
from .generic import Generic, Unscaled
from .cds import CDS
from .console import Console
from .fits import Fits
from .latex import Latex
from .ogip import OGIP
from .unicode_format import Unicode
from .vounit import VOUnit

__all__ = [
    'Base', 'Generic', 'CDS', 'Console', 'Fits', 'Latex', 'OGIP', 'Unicode',
    'Unscaled', 'VOUnit', 'get_format']

FORMATS = [
    'Generic', 'CDS', 'Console', 'Fits', 'Latex', 'OGIP', 'Unicode',
    'Unscaled', 'VOUnit']


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
        return format()
    elif isinstance(format, Base):
        return format

    if format is None:
        format = 'generic'
    format = format.lower()
    for key in FORMATS:
        val = globals()[key]
        if isclass(val) and (issubclass(val, Base) and key.lower() == format.lower()):
            return val()
    raise ValueError("Unknown format {0!r}".format(format))
