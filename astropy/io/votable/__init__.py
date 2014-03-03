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
from ... import config as _config

__all__ = [
    'parse', 'parse_single_table', 'validate', 'from_table',
    'is_votable', 'writeto', 'VOWarning', 'VOTableChangeWarning',
    'VOTableSpecWarning', 'UnimplementedWarning', 'IOWarning',
    'VOTableSpecError', 'conf']


class _Conf(_config.ConfigNamespace):
    pedantic = _config.ConfigItem(
        False,
        'When True, treat fixable violations of the VOTable spec as exceptions.',
        aliases=['astropy.io.votable.table.pedantic'])
conf = _Conf()
