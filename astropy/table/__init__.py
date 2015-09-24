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
from .table import Table, QTable, TableColumns, Row, TableFormatter, NdarrayMixin
from .operations import join, hstack, vstack, unique, TableMergeError
from .jsviewer import JSViewer

# Import routines that connect readers/writers to astropy.table
from ..io.ascii import connect
from ..io.fits import connect
from ..io.misc import connect
from ..io.votable import connect
from . import jsviewer
