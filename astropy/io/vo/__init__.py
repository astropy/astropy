# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package reads and writes data formats used by the Virtual Observatory (VO)
initiative, particularly the VOTable format.
"""

from __future__ import division, absolute_import

# If we're in the source directory, don't import anything, since that
# requires 2to3 to be run.
from astropy import setup_helpers
if setup_helpers.is_in_build_mode():
    pass
else:
    del setup_helpers
    from .table import parse, parse_single_table
    from .exceptions import (VOWarning, VOTableChangeWarning,
        VOTableSpecWarning, UnimplementedWarning, IOWarning,
        VOTableSpecError)

__all__ = [
    'parse', 'parse_single_table', 'VOWarning', 'VOTableChangeWarning',
    'VOTableSpecWarning', 'UnimplementedWarning', 'IOWarning',
    'VOTableSpecError'
    ]
