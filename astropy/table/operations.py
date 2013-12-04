"""
High-level table operations:

- join()
- hstack()
- vstack()
"""
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six

from copy import deepcopy
import warnings
import collections

import numpy as np

from ..utils import OrderedDict, metadata

from . import np_utils

__all__ = ['join', 'hstack', 'vstack']


def _merge_col_meta(out, tables, col_name_map, idx_left=0, idx_right=1,
                    metadata_conflicts='warn'):
    """
    Merge column meta data for the ``out`` table.

    This merges column meta, which includes attributes unit, format,
    and description, as well as the actual `meta` atttribute.  It is
    assumed that the ``out`` table was created by merging ``tables``.
    The ``col_name_map`` provides the mapping from col name in ``out``
    back to the original name (which may be different).
    """
    # Set column meta
    attrs = ('unit', 'format', 'description')
    for out_col in six.itervalues(out.columns):
        for idx_table, table in enumerate(tables):
            left_col = out_col
            right_name = col_name_map[out_col.name][idx_table]

            if right_name:
                right_col = table[right_name]
                out_col.meta = metadata.merge(left_col.meta, right_col.meta, metadata_conflicts=metadata_conflicts)
                for attr in attrs:

                    # Pick the metadata item that is not None, or they are both
                    # not None, then if they are equal, there is no conflict,
                    # and if they are different, there is a conflict and we
                    # pick the one on the right (or raise an error).

                    left_attr = getattr(left_col, attr)
                    right_attr = getattr(right_col, attr)

                    if left_attr is None:
                        # This may not seem necessary since merge_attr gets set
                        # to right_attr, but not all objects support != which is
                        # needed for one of the if clauses.
                        merge_attr = right_attr
                    elif right_attr is None:
                        merge_attr = left_attr
                    elif left_attr != right_attr:
                        if metadata_conflicts == 'warn':
                            warnings.warn("In merged column '{0}' the '{1}' attribute does not match "
                                          "({2} != {3}).  Using {3} for merged output"
                                          .format(out_col.name, attr, left_attr, right_attr),
                                          metadata.MergeConflictWarning)
                        elif metadata_conflicts == 'error':
                            raise metadata.MergeConflictError('In merged column {0!r} the {1!r} attribute does not match '
                                          '({2} != {3})'.format(out_col.name, attr, left_attr, right_attr))
                        elif metadata_conflicts != 'silent':
                            raise ValueError('metadata_conflict argument must be one of "silent", "warn", or "error"')
                        merge_attr = right_attr
                    else:  # left_attr == right_attr
                        merge_attr = right_attr

                    setattr(out_col, attr, merge_attr)


def _merge_table_meta(out, tables, metadata_conflicts='warn'):
    out_meta = deepcopy(tables[0].meta)
    for table in tables[1:]:
        out_meta = metadata.merge(out_meta, table.meta, metadata_conflicts=metadata_conflicts)
    out.meta.update(out_meta)


def _get_list_of_tables(tables):
    """
    Check that tables is a Table or sequence of Tables.  Returns the
    corresponding list of Tables.
    """
    from .table import Table, Row

    # Make sure we have a list of things
    if not isinstance(tables, collections.Sequence):
        tables = [tables]

    # Make sure each thing is a Table or Row
    if any(not isinstance(x, (Table, Row)) for x in tables):
        raise TypeError('`tables` arg must be a Table or sequence of Tables or Rows')

    # Convert any Rows to Tables
    tables = [(x if isinstance(x, Table) else Table(x)) for x in tables]

    return tables


def join(left, right, keys=None, join_type='inner',
         uniq_col_name='{col_name}_{table_name}',
         table_names=['1', '2'], metadata_conflicts='warn'):
    """
    Perform a join of the left table with the right table on specified keys.

    Parameters
    ----------
    left : Table object or a value that will initialize a Table object
        Left side table in the join
    right : Table object or a value that will initialize a Table object
        Right side table in the join
    keys : str or list of str
        Name(s) of column(s) used to match rows of left and right tables.
        Default is to use all columns which are common to both tables.
    join_type : str
        Join type ('inner' | 'outer' | 'left' | 'right'), default is 'inner'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2'].
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.
    """
    from .table import Table

    # Try converting inputs to Table as needed
    if not isinstance(left, Table):
        left = Table(left)
    if not isinstance(right, Table):
        right = Table(right)

    col_name_map = OrderedDict()
    out_data = np_utils.join(left._data, right._data, keys, join_type,
                             uniq_col_name, table_names, col_name_map)
    # Create the output (Table or subclass of Table)
    out = Table(out_data)

    # Merge the column and table meta data. Table subclasses might override
    # these methods for custom merge behavior.
    _merge_col_meta(out, [left, right], col_name_map, metadata_conflicts=metadata_conflicts)
    _merge_table_meta(out, [left, right], metadata_conflicts=metadata_conflicts)

    return out


def vstack(tables, join_type='outer', metadata_conflicts='warn'):
    """
    Stack tables vertically (along rows)

    A ``join_type`` of 'exact' means that the tables must all
    have exactly the same column names (though the order can vary).  If
    ``join_type`` is 'inner' then the intersection of common columns will
    be output.  A value of 'outer' (default) means the output will have the union of
    all columns, with table values being masked where no common values are
    available.

    Parameters
    ----------

    tables : Table or list of Table objects
        Table(s) to stack along rows (vertically) with the current table
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'exact'
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.

    Examples
    --------

    To stack two tables along rows do::

      >>> from astropy.table import vstack, Table
      >>> t1 = Table({'a': [1, 2], 'b': [3, 4]}, names=('a', 'b'))
      >>> t2 = Table({'a': [5, 6], 'b': [7, 8]}, names=('a', 'b'))
      >>> print(t1)
       a   b
      --- ---
        1   3
        2   4
      >>> print(t2)
       a   b
      --- ---
        5   7
        6   8
      >>> print(vstack([t1, t2]))
       a   b
      --- ---
        1   3
        2   4
        5   7
        6   8
    """
    from .table import Table

    tables = _get_list_of_tables(tables)  # validates input
    arrays = [table._data for table in tables]
    col_name_map = OrderedDict()

    out_data = np_utils.vstack(arrays, join_type, col_name_map)
    out = Table(out_data)

    # Merge column and table metadata
    _merge_col_meta(out, tables, col_name_map, metadata_conflicts=metadata_conflicts)
    _merge_table_meta(out, tables, metadata_conflicts=metadata_conflicts)

    return out


def hstack(tables, join_type='outer',
           uniq_col_name='{col_name}_{table_name}', table_names=None,
           metadata_conflicts='warn'):
    """
    Stack tables along columns (horizontally)

    A ``join_type`` of 'exact' means that the tables must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be output.  A value of 'outer' (default) means
    the output will have the union of all rows, with table values being
    masked where no common values are available.

    Parameters
    ----------

    tables : List of Table objects
        Tables to stack along columns (horizontally) with the current table
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2', ..].
    col_name_map : empty dict or None
        If passed as a dict then it will be updated in-place with the
        mapping of output to input column names.
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.

    Examples
    --------

    To stack two tables horizontally (along columns) do::

      >>> from astropy.table import Table, hstack
      >>> t1 = Table({'a': [1, 2], 'b': [3, 4]}, names=('a', 'b'))
      >>> t2 = Table({'c': [5, 6], 'd': [7, 8]}, names=('c', 'd'))
      >>> print(t1)
       a   b
      --- ---
        1   3
        2   4
      >>> print(t2)
       c   d
      --- ---
        5   7
        6   8
      >>> print(hstack([t1, t2]))
       a   b   c   d
      --- --- --- ---
        1   3   5   7
        2   4   6   8
    """
    from .table import Table

    tables = _get_list_of_tables(tables)  # validates input
    arrays = [table._data for table in tables]
    col_name_map = OrderedDict()

    out_data = np_utils.hstack(arrays, join_type, uniq_col_name, table_names,
                               col_name_map)
    out = Table(out_data)

    _merge_col_meta(out, tables, col_name_map, metadata_conflicts=metadata_conflicts)
    _merge_table_meta(out, tables, metadata_conflicts=metadata_conflicts)

    return out
