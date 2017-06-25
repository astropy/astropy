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
    default_notebook_table_class = _config.ConfigItem(
        'table-striped table-bordered table-condensed',
        'The table class to be used in Jupyter notebooks when displaying '
        'tables (and not overridden). See <http://getbootstrap.com/css/#tables '
        'for a list of useful bootstrap classes.')
    replace_warnings = _config.ConfigItem(
        ['slice'],
        'List of conditions for issuing a warning when replacing a table '
        "column using setitem, e.g. t['a'] = value.  Allowed options are "
        "'always', 'slice', 'refcount', 'attributes'.",
        'list',
        )
    replace_inplace = _config.ConfigItem(
        False,
        'Always use in-place update of a table column when using setitem, '
        "e.g. t['a'] = value.  This overrides the default behavior of "
        "replacing the column entirely with the new value when possible. "
        "This configuration option will be deprecated and then removed in "
        "subsequent major releases."
        )


conf = Conf()


from .column import Column, MaskedColumn, StringTruncateWarning, ColumnInfo
from .groups import TableGroups, ColumnGroups
from .table import (Table, QTable, TableColumns, Row, TableFormatter,
                    NdarrayMixin, TableReplaceWarning)
from .operations import join, hstack, vstack, unique, TableMergeError
from .bst import BST, FastBST, FastRBT
from .sorted_array import SortedArray
from .serialize import SerializedColumn

# Finally import the formats for the read and write method but delay building
# the documentation until all are loaded. (#5275)
from ..io import registry

with registry.delay_doc_updates(Table):
    # Import routines that connect readers/writers to astropy.table
    from .jsviewer import JSViewer
    from ..io.ascii import connect
    from ..io.fits import connect
    from ..io.misc import connect
    from ..io.votable import connect
