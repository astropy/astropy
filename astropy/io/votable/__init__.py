# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package reads and writes data formats used by the Virtual
Observatory (VO) initiative, particularly the VOTable XML format.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

from .table import (
    parse, parse_single_table, validate, from_table, is_votable, writeto)
from .exceptions import (
    VOWarning, VOTableChangeWarning, VOTableSpecWarning, UnimplementedWarning,
    IOWarning, VOTableSpecError)

__all__ = [
    'parse', 'parse_single_table', 'validate', 'from_table',
    'is_votable', 'writeto', 'VOWarning', 'VOTableChangeWarning',
    'VOTableSpecWarning', 'UnimplementedWarning', 'IOWarning',
    'VOTableSpecError']
