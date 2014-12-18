# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six

from .. import config as _config

class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.table`.
    """

    auto_colname = _config.ConfigItem(
        'col{0}',
        'The template that determines the name of a column if it cannot be '
        'determined. Uses new-style (format method) string formatting.',
        aliases=['astropy.table.column.auto_colname'])
conf = Conf()


from .column import Column, MaskedColumn
from .groups import TableGroups, ColumnGroups
from .table import Table, TableColumns, Row, TableFormatter
from .operations import join, hstack, vstack, unique, TableMergeError

