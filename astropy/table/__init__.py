# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six

from .. import config as _config

class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.table`.
    """

    max_lines = _config.ConfigItem(
        25, 'Maximum number of lines for the pretty-printer to use if it '
        'cannot determine the terminal size. Negative numbers mean no limit.',
        aliases=['astropy.table.pprint.max_lines'])
    max_width = _config.ConfigItem(
        80, 'Maximum number of characters for the pretty-printer to use '
        'per line if it cannot determine the terminal size.  Negative numbers '
        'mean no limit.',
        aliases=['astropy.table.pprint.max_width'])
    auto_colname = _config.ConfigItem(
        'col{0}',
        'The template that determines the name of a column if it cannot be '
        'determined. Uses new-style (format method) string formatting.',
        aliases=['astropy.table.column.auto_colname'])
conf = Conf()


from .column import Column, MaskedColumn
from .table import Table, TableColumns, Row, TableFormatter
from .np_utils import TableMergeError
from .operations import join, hstack, vstack

# Import routines that connect readers/writers to astropy.table
from ..io.ascii import connect
from ..io.fits import connect
from ..io.misc import connect
from ..io.votable import connect
