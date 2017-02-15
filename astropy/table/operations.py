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
from ..extern.six.moves import zip, range
from ..units import Quantity

from copy import deepcopy
import warnings
import collections
import itertools
from collections import OrderedDict, Counter

import numpy as np

from ..utils import metadata

from . import _np_utils
from .np_utils import fix_column_name, TableMergeError

__all__ = ['join', 'hstack', 'vstack', 'unique']


def _merge_col_meta(tables, col_name_map, metadata_conflicts='warn'):
    """Merge column meta data from a set of tables.

    This merges column meta, which includes attributes unit, format,
    and description, as well as the actual `meta` attribute.

    Parameters
    ----------
    tables : list of `~astropy.table.Table`
        Tables with columns from which to get meta data.
    col_name_map : dict
        For each output column name key, a list to the input column names
        (`None` if a given table does not have that column).
    metadata_conflicts : {'warn', 'error', 'silent'}
        How to proceed with metadata conflicts. This should be one of:
            * ``'warn'``: pick the last conflicting meta-data value,
              but emit a warning (default);
            * ``'silent'``: silently pick the last conflicting meta-data value;
            * ``'error'``: raise an exception.

    Returns
    -------
    col_attrs_map : dict of dict
        Merged attributes for the output columns.
    """
    # Attributes to merge (other than 'meta' itself)
    attrs = ('meta', 'unit', 'format', 'description')
    # Initialize output attributes dicts
    col_attrs_map = OrderedDict()
    for out_name in col_name_map:
        out_attrs = {}
        for table, right_name in zip(tables, col_name_map[out_name]):
            if right_name is None:
                # input table doesn't have this column.
                continue

            right_col = table[right_name]
            for attr in attrs:
                # Pick the metadata item that is not None, or they are both
                # not None, then if they are equal, there is no conflict,
                # and if they are different, there is a conflict and we
                # pick the one on the right (or raise an error).
                right_attr = getattr(right_col.info, attr, None)
                if not right_attr:
                    continue

                # Use the attribute if we don't have it yet, or if it is a
                # unit and it is equivalent (so there won't be a problem;
                # using here that we always take the right-most occurrence)
                if (attr not in out_attrs or
                    (attr == 'unit' and isinstance(right_col, Quantity))):
                    out_attrs[attr] = right_attr

                elif attr == 'meta':
                    out_attrs[attr] = metadata.merge(out_attrs[attr], right_attr,
                                                     metadata_conflicts=metadata_conflicts)

                elif out_attrs[attr] != right_attr:
                    # we have a real conflict
                    if metadata_conflicts == 'warn':
                        warnings.warn("In merged column '{0}' the '{1}' attribute does not match "
                                      "({2} != {3}).  Using {3} for merged output"
                                      .format(out_name, attr, out_attrs[attr], right_attr),
                                      metadata.MergeConflictWarning)
                    elif metadata_conflicts == 'error':
                        raise metadata.MergeConflictError(
                            'In merged column {0!r} the {1!r} attribute does not match '
                            '({2} != {3})'.format(out_name, attr,
                                                  out_attrs[attr], right_attr))
                    elif metadata_conflicts != 'silent':
                        raise ValueError('metadata_conflicts argument must be one of "silent",'
                                         ' "warn", or "error"')
                    out_attrs[attr] = right_attr

        col_attrs_map[out_name] = out_attrs

    return col_attrs_map


def _merge_table_meta(tables, metadata_conflicts='warn'):
    out_meta = deepcopy(tables[0].meta)
    for table in tables[1:]:
        out_meta = metadata.merge(out_meta, table.meta, metadata_conflicts=metadata_conflicts)
    return out_meta


def _get_list_of_tables(tables):
    """Check that tables is a Table or sequence of Tables.

    Returns the corresponding list of Tables.
    """
    from .table import Table, Row

    # Make sure we have a list of things
    if not isinstance(tables, collections.Sequence):
        tables = [tables]

    # Make sure each thing is a Table or Row
    if any(not isinstance(x, (Table, Row)) for x in tables) or len(tables) == 0:
        raise TypeError('`tables` arg must be a Table or sequence of Tables or Rows')

    # Convert any Rows to Tables
    tables = [(x if isinstance(x, Table) else Table(x)) for x in tables]

    return tables


def _get_out_class(tables):
    """Get deepest subclass from a list of table instances.

    It is assumed that `tables` is a list of at least one element and that they
    are all Table (subclass) instances.  This doesn't handle complicated
    inheritance schemes.
    """
    out_class = tables[0].__class__
    for t in tables[1:]:
        if issubclass(t.__class__, out_class):
            out_class = t.__class__
    return out_class


def join(left, right, keys=None, join_type='inner',
         uniq_col_name='{col_name}_{table_name}',
         table_names=['1', '2'], metadata_conflicts='warn'):
    """Perform a join of the left and right tables on specified keys.

    Parameters
    ----------
    left : `~astropy.table.Table` or a value that will initialize one
        Left side table in the join.
    right : `~astropy.table.Table` or a value that will initialize one
        Right side table in the join.
    keys : str or list of str
        Name(s) of column(s) used to match rows of left and right tables.
        Default is to use all columns which are common to both tables.
    join_type : {'inner', 'outer', 'left', 'right'}, optional
        Default is 'inner'.
    uniq_col_name : str, optional
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str, optional
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2'].
    metadata_conflicts : {'warn', 'error', 'silent'}
        How to proceed with metadata conflicts. This should be one of:
            * ``'warn'``: pick the last conflicting meta-data value,
              but emit a warning (default);
            * ``'silent'``: silently pick the last conflicting meta-data value;
            * ``'error'``: raise an exception.

    Returns
    -------
    joined_table : `~astropy.table.Table` object
        New table containing the result of the join operation.
    """
    from .table import Table

    # Try converting inputs to Table as needed
    if not isinstance(left, Table):
        left = Table(left)
    if not isinstance(right, Table):
        right = Table(right)

    if join_type not in ('inner', 'outer', 'left', 'right'):
        raise ValueError("The 'join_type' argument should be in 'inner', "
                         "'outer', 'left' or 'right' (got '{0}' instead)".
                         format(join_type))

    if keys is None:
        keys = tuple(name for name in left.colnames if name in right.colnames)
        if len(keys) == 0:
            raise TableMergeError('No keys in common between left and right tables')
    elif isinstance(keys, six.string_types):
        keys = (keys,)

    # Check the key columns
    for table, table_label in ((left, 'Left'), (right, 'Right')):
        for name in keys:
            if name not in table.colnames:
                raise TableMergeError('{0} table does not have key column {1!r}'
                                      .format(table_label, name))
            if hasattr(table[name], 'mask') and np.any(table[name].mask):
                raise TableMergeError('{0} key column {1!r} has missing values'
                                      .format(table_label, name))
            if not isinstance(table[name], np.ndarray):
                raise ValueError("non-ndarray column '{0}' not allowed as a key column")

    # Get dict with list of input columns for each output column.
    col_name_map = get_col_name_map([left, right], keys, uniq_col_name, table_names)
    # Get output column attributes by merging input column attributes.
    col_attrs_map = _merge_col_meta([left, right], col_name_map,
                                    metadata_conflicts=metadata_conflicts)
    # Actual join.
    out = _join(left, right, keys, join_type,
                uniq_col_name, table_names, col_name_map, col_attrs_map)

    # Merge the table meta data.
    out_meta = _merge_table_meta([left, right], metadata_conflicts=metadata_conflicts)
    out.meta.update(out_meta)
    return out


def vstack(tables, join_type='outer', metadata_conflicts='warn'):
    """Stack tables vertically (along rows)

    A ``join_type`` of 'exact' means that the tables must all have exactly
    the same column names (though the order can vary).  If ``join_type``
    is 'inner' then the intersection of common columns will be the output.
    A value of 'outer' (default) means the output will have the union of
    all columns, with table values being masked where no common values are
    available.

    Parameters
    ----------
    tables : list of `~astropy.table.Table` objects
        Table(s) to stack along rows (vertically).
    join_type : {'outer', 'exact', 'inner'}, optional
        Join type (see above). Default is 'outer'.
    metadata_conflicts : {'warn', 'error', 'silent'}
        How to proceed with metadata conflicts. This should be one of:
            * ``'warn'``: pick the last conflicting meta-data value,
              but emit a warning (default);
            * ``'silent'``: silently pick the last conflicting meta-data value;
            * ``'error'``: raise an exception.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.

    Notes
    -----
    Both `~astropy.table.Table` and `~astropy.table.QTable` instances can be
    stacked, but for the latter an outer join is not generally possible, since
    this requires output columns with partially masked data, which is not
    supported for `~astropy.units.Quantity` columns.

    The two also differ in how units are treated:

     * For `~astropy.table.Table`, which hold regular `~astropy.table.Column`
       instances, the units are treated as other attributes, i.e., when they
       are not all the same, warnings or errors are given as set by
       ``metadata_conflicts``, but the values in the stacked columns are simply
       those given by the input columns.
     * For `~astropy.table.QTable`, which hold `~astropy.units.Quantity`
       instances, units are considered integral part of the data.  In case of
       conflicts, the output column has the unit of the last conflicting unit,
       and all values are converted to that unit (raising an exception if that
       is not possible).

    Examples
    --------
    To stack two tables along rows do::

      >>> from astropy.table import vstack, Table, QTable
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
      >>> from astropy import units as u
      >>> t3 = QTable({'a': [1., 2.] * u.m})
      >>> t4 = QTable({'a': [3., 4.] * u.km})
      >>> print(vstack([t3, t4]))
        a
        km
      -----
      0.001
      0.002
        3.0
        4.0
      >>> t5 = Table(t3)
      >>> t6 = Table(t4)
      >>> print(vstack([t5, t6]))
       a
       km
      ---
      1.0
      2.0
      3.0
      4.0
    """
    from .column import BaseColumn
    from .table import QTable
    tables = _get_list_of_tables(tables)  # validates input
    if len(tables) == 1:
        return tables[0]  # no point in stacking a single table

    # we can only deal with "standard" columns: Column for Table, and
    # Column or Quantity for QTable.
    if any(table.has_mixin_columns and
           not (isinstance(table, QTable) and
                all(isinstance(col, (BaseColumn, Quantity))
                    for col in six.itervalues(table.columns)))
           for table in tables):
        raise NotImplementedError('vstack not available for tables with mixin columns')

    # Start by assuming an outer match where all names go to output
    names = set(itertools.chain(*[table.colnames for table in tables]))
    # Get dict with list of input columns for each output column.
    col_name_map = get_col_name_map(tables, names)
    if join_type == 'exact':
        # For 'exact', the output must have exactly the same number of columns
        # as each input array.
        for names in six.itervalues(col_name_map):
            if any(x is None for x in names):
                raise TableMergeError('Inconsistent columns in input tables '
                                      "(use 'inner' or 'outer' join_type to "
                                      "allow non-matching columns).")

    elif join_type == 'inner':
        # Keep only those columns where all input arrays have that column.
        col_name_map = OrderedDict((name, in_names) for name, in_names in six.iteritems(col_name_map)
                                   if all(x is not None for x in in_names))
        if len(col_name_map) == 0:
            raise TableMergeError('Input tables have no columns in common')

    elif join_type == 'outer':
        if any(isinstance(table, QTable) and
               any(isinstance(table[name], Quantity) and
                   any(x is None
                       for x in col_name_map[name])
                   for name in table.colnames)
               for table in tables):
            raise NotImplementedError("vstack with join_type 'outer' not "
                                      "available for QTable input with "
                                      "non-matching columns.")

    else:
        raise ValueError("`join_type` must be one of 'inner', 'exact' or "
                         "'outer', not '{0}'".format(join_type))

    # Merge column metadata to create output column attributes.
    col_attrs_map = _merge_col_meta(tables, col_name_map, metadata_conflicts=metadata_conflicts)
    # Create the actual output table.
    out = _vstack(tables, col_name_map, col_attrs_map)
    # Merge input table metadata and use it to set output table metadata.
    out_meta = _merge_table_meta(tables, metadata_conflicts=metadata_conflicts)
    out.meta.update(out_meta)
    return out


def hstack(tables, join_type='outer',
           uniq_col_name='{col_name}_{table_name}', table_names=None,
           metadata_conflicts='warn'):
    """Stack tables along columns (horizontally)

    A ``join_type`` of 'exact' means that the tables must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be the output.  A value of 'outer' (default)
    means the output will have the union of all rows, with table values being
    masked where no common values are available.

    Parameters
    ----------
    tables : List of `~astropy.table.Table` objects
        Tables to stack along columns (horizontally).
    join_type : {'outer', 'exact', 'inner'}, optional
        Join type (see above). Default is 'outer'.
    uniq_col_name : str, optional
        String to generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str, optional
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2', ..].
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'warn'``: pick the last conflicting meta-data value,
              but emit a warning (default);
            * ``'silent'``: silently pick the last conflicting meta-data value;
            * ``'error'``: raise an exception.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.

    Notes
    -----
    Both `~astropy.table.Table` and `~astropy.table.QTable` instances can be
    stacked, but for the latter an outer join is not generally possible, since
    this requires output columns with partially masked data, which is not
    supported for `~astropy.units.Quantity` columns.

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
    # Input validation
    tables = _get_list_of_tables(tables)  # validates input
    if len(tables) == 1:
        return tables[0]  # no point in stacking a single table

    if join_type not in ('inner', 'exact', 'outer'):
        raise ValueError("join_type arg must be 'inner', 'exact' or 'outer', "
                         "not '{0}'".format(join_type))

    if table_names is None:
        table_names = ['{0}'.format(ii + 1) for ii in range(len(tables))]
    elif len(tables) != len(table_names):
        raise ValueError('Number of tables must match number of table_names')

    # Check for different lengths for 'inner' and 'exact'
    if join_type != 'outer' and len(set(len(table) for table in tables)) > 1:
        if join_type == 'exact':
            # For 'exact', all input arrays must have the same length
            raise TableMergeError("Inconsistent number of rows in input tables "
                                  "(use 'inner' or 'outer' join_type to allow "
                                  "non-matching rows)")
        else:
            # For 'inner', we only keep the common rows.
            min_table_len = min(len(table) for table in tables)
            tables = [table[:min_table_len] for table in tables]

    # Get dict with list of input columns for each output column.
    col_name_map = get_col_name_map(tables, [], uniq_col_name, table_names)
    # Unlike vstack/join, no need to merge column metadata, as here every column
    # originates from a single other one,
    out = _hstack(tables, col_name_map)
    # Merge the table meta data.
    out_meta = _merge_table_meta(tables, metadata_conflicts=metadata_conflicts)
    out.meta.update(out_meta)
    return out


def unique(input_table, keys=None, silent=False, keep='first'):
    """
    Returns the unique rows of a table.

    Parameters
    ----------

    input_table : `~astropy.table.Table` object or a value that
        will initialize a `~astropy.table.Table` object
    keys : str or list of str, optional
        Name(s) of column(s) used to create unique rows.
        Default is to use all columns.
    silent : bool, optional
        If `True`, masked value column(s) are silently removed from
        ``keys``. If `False`, an exception is raised when ``keys``
        contains masked value column(s).
        Default is `False`.
    keep : {'first', 'last', 'none'}, optional
        Whether to keep the first or last row for each set of
        duplicates. If 'none', all rows that are duplicate are
        removed, leaving only rows that are already unique in
        the input.
        Default is 'first'.

    Returns
    -------
    unique_table : `~astropy.table.Table` object
        New table containing only the unique rows of ``input_table``.

    Examples
    --------
    >>> from astropy.table import Table
    >>> import numpy as np
    >>> table = Table(data=[[1,2,3,2,3,3],
    ... [2,3,4,5,4,6],
    ... [3,4,5,6,7,8]],
    ... names=['col1', 'col2', 'col3'],
    ... dtype=[np.int32, np.int32, np.int32])
    >>> table
    <Table length=6>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     3     4
        3     4     5
        2     5     6
        3     4     7
        3     6     8
    >>> unique(table, keys='col1')
    <Table length=3>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     3     4
        3     4     5
    >>> unique(table, keys=['col1'], keep='last')
    <Table length=3>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     5     6
        3     6     8
    >>> unique(table, keys=['col1', 'col2'])
    <Table length=5>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     3     4
        2     5     6
        3     4     5
        3     6     8
    >>> unique(table, keys=['col1', 'col2'], keep='none')
    <Table length=4>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     3     4
        2     5     6
        3     6     8
    >>> unique(table, keys=['col1'], keep='none')
    <Table length=1>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3

    """

    if keep not in ('first', 'last', 'none'):
        raise ValueError("'keep' should be one of 'first', 'last', 'none'")

    if isinstance(keys, six.string_types):
        keys = [keys]
    elif keys is None:
        keys = input_table.colnames
    elif len(set(keys)) != len(keys):
        raise ValueError("duplicate key names")

    if input_table.masked:
        for key in keys[:]:
            if np.any(input_table[key].mask):
                if not silent:
                    raise ValueError(
                        "cannot use columns with masked values as keys; "
                        "remove column '{0}' from keys and rerun "
                        "unique()".format(key))
                del keys[keys.index(key)]
        if len(keys) == 0:
            raise ValueError("no column remained in ``keys``; "
                             "unique() cannot work with masked value "
                             "key columns")

    grouped_table = input_table.group_by(keys)
    indices = grouped_table.groups.indices
    if keep == 'first':
        indices = indices[:-1]
    elif keep == 'last':
        indices = indices[1:] - 1
    else:
        indices = indices[:-1][np.diff(indices) == 1]

    return grouped_table[indices]


def get_col_name_map(tables, common_names, uniq_col_name='{col_name}_{table_name}',
                     table_names=None):
    """
    Find the column names mapping when merging the list of tables
    ``tables``.  It is assumed that col names in ``common_names`` are to be
    merged into a single column while the rest will be uniquely represented
    in the output.  The args ``uniq_col_name`` and ``table_names`` specify
    how to rename columns in case of conflicts.

    Returns a dict mapping each output column name to the input(s).  This takes the form
    {outname : (col_name_0, col_name_1, ...), ... }.  For key columns all of input names
    will be present, while for the other non-key columns the value will be (col_name_0,
    None, ..) or (None, col_name_1, ..) etc.
    """

    col_name_map = collections.defaultdict(lambda: [None] * len(tables))
    col_name_list = []

    if table_names is None:
        table_names = [six.text_type(ii + 1) for ii in range(len(tables))]

    for idx, table in enumerate(tables):
        table_name = table_names[idx]
        for name in table.colnames:
            out_name = name

            if name in common_names:
                # If name is in the list of common_names then insert into
                # the column name list, but just once.
                if name not in col_name_list:
                    col_name_list.append(name)
            else:
                # If name is not one of the common column outputs, and it collides
                # with the names in one of the other arrays, then rename
                others = list(tables)
                others.pop(idx)
                if any(name in other.colnames for other in others):
                    out_name = uniq_col_name.format(table_name=table_name, col_name=name)
                col_name_list.append(out_name)

            col_name_map[out_name][idx] = name

    # Check for duplicate output column names
    col_name_count = Counter(col_name_list)
    repeated_names = [name for name, count in six.iteritems(col_name_count) if count > 1]
    if repeated_names:
        raise TableMergeError('Merging column names resulted in duplicates: {0}.  '
                              'Change uniq_col_name or table_names args to fix this.'
                              .format(repeated_names))

    # Convert col_name_map to a regular dict with tuple (immutable) values
    col_name_map = OrderedDict((name, col_name_map[name]) for name in col_name_list)

    return col_name_map


def get_descrs(arrays, col_name_map):
    """
    Find the dtypes descrs resulting from merging the list of arrays' dtypes,
    using the column name mapping ``col_name_map``.

    Return a list of descrs for the output.
    """

    out_descrs = []

    for out_name, in_names in six.iteritems(col_name_map):
        # List of input arrays that contribute to this output column
        in_cols = [arr[name] for arr, name in zip(arrays, in_names) if name is not None]

        # List of names of the columns that contribute to this output column.
        names = [name for name in in_names if name is not None]

        # Output dtype is the superset of all dtypes in in_arrays
        try:
            dtype = common_dtype(in_cols)
        except TableMergeError as tme:
            # Beautify the error message when we are trying to merge columns with incompatible
            # types by including the name of the columns that originated the error.
            raise TableMergeError("The '{0}' columns have incompatible types: {1}"
                                  .format(names[0], tme._incompat_types))

        # Make sure all input shapes are the same
        uniq_shapes = set(col.shape[1:] for col in in_cols)
        if len(uniq_shapes) != 1:
            raise TableMergeError('Key columns {0!r} have different shape'.format(names))
        shape = uniq_shapes.pop()

        out_descrs.append((fix_column_name(out_name), dtype, shape))

    return out_descrs


def common_dtype(cols):
    """
    Use numpy to find the common dtype for a list of columns.

    Only allow columns within the following fundamental numpy data types:
    np.bool_, np.object_, np.number, np.character, np.void
    """
    try:
        return metadata.common_dtype(cols)
    except metadata.MergeConflictError as err:
        tme = TableMergeError('Columns have incompatible types {0}'
                              .format(err._incompat_types))
        tme._incompat_types = err._incompat_types
        raise tme


def _join(left, right, keys=None, join_type='inner',
          uniq_col_name='{col_name}_{table_name}',
          table_names=['1', '2'],
          col_name_map=None, col_attrs_map=None):
    """Perform a join of the left and right Tables on specified keys.

    Parameters
    ----------
    left : Table
        Left side table in the join
    right : Table
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
    col_name_map : dict
        Mapping of output to input column names.

    Returns
    -------
    joined_table : `~astropy.table.Table` object
        New table containing the result of the join operation.
    """
    len_left, len_right = len(left), len(right)

    if len_left == 0 or len_right == 0:
        raise ValueError('input tables for join must both have at least one row')

    # Make an array with just the key columns.  This uses a temporary
    # structured array for efficiency.
    out_descrs = get_descrs([left, right], col_name_map)
    out_keys_dtype = [descr for descr in out_descrs if descr[0] in keys]
    out_keys = np.empty(len_left + len_right, dtype=out_keys_dtype)
    for key in keys:
        out_keys[key][:len_left] = left[key]
        out_keys[key][len_left:] = right[key]
    idx_sort = out_keys.argsort(order=keys)
    out_keys = out_keys[idx_sort]

    # Get all keys
    diffs = np.concatenate(([True], out_keys[1:] != out_keys[:-1], [True]))
    idxs = np.flatnonzero(diffs)

    # Main inner loop in Cython to compute the cartesion product
    # indices for the given join type
    int_join_type = {'inner': 0, 'outer': 1, 'left': 2, 'right': 3}[join_type]
    masked, n_out, left_out, left_mask, right_out, right_mask = \
        _np_utils.join_inner(idxs, idx_sort, len_left, int_join_type)

    # If either of the inputs are masked then the output is masked
    masked = left.masked or right.masked or bool(masked)

    out = _get_out_class([left, right])(masked=masked)

    for (out_name, dtype, shape), col_attrs in zip(out_descrs, six.itervalues(col_attrs_map)):

        left_name, right_name = col_name_map[out_name]
        if left_name and right_name:  # this is a key which comes from left and right
            out[out_name] = out.ColumnClass(length=n_out, name=out_name, dtype=dtype, shape=shape,
                                            **col_attrs)
            # `np.where` would be logical, but is not safe for Quantity.
            # right_mask is True where right is masked, i.e., should not be used.
            out[out_name][~right_mask] = right[right_name].take(right_out)[~right_mask]
            out[out_name][right_mask] = left[left_name].take(left_out)[right_mask]
            continue
        elif left_name:  # out_name came from the left table
            name, array, array_out, array_mask = left_name, left, left_out, left_mask
        elif right_name:
            name, array, array_out, array_mask = right_name, right, right_out, right_mask
        else:
            raise TableMergeError('Unexpected column names (maybe one is ""?)')

        # Finally add the joined column to the output table.
        out[out_name] = array[name][array_out]

        # If the output table is masked then set the output column masking
        # accordingly.  Check for columns that don't support a mask attribute.
        if masked:
            # array_mask is 1-d corresponding to length of output column.  We need
            # make it have the correct shape for broadcasting, i.e. (length, 1, 1, ..).
            # Mixin columns might not have ndim attribute so use len(col.shape).
            array_mask.shape = (out[out_name].shape[0],) + (1,) * (len(out[out_name].shape) - 1)

            if array.masked:
                array_mask = array_mask | array[name].mask[array_out]
            try:
                out[out_name].mask[:] = array_mask
            except ValueError:
                raise ValueError("join requires masking column '{0}' but column"
                                 " type {1} does not support masking"
                                 .format(out_name, out[out_name].__class__.__name__))

    return out


def _vstack(tables, col_name_map=None, col_attrs_map=None):
    """Stack Tables vertically (by rows)

    Parameters
    ----------
    tables : list of Tables
        Tables to stack by rows (vertically)
    col_name_map : dict
        Mapping of output to input column names.
    col_attrs_map : dict
        Attributes for the output columns to be created.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.
    """
    # Trivial case of one input array
    if len(tables) == 1:
        return tables[0]

    # If there are any output columns where one or more input arrays are missing
    # then the output must be masked.  If any input arrays are masked then
    # output is masked.
    masked = any(getattr(table, 'masked', False) for table in tables)
    for names in six.itervalues(col_name_map):
        if any(x is None for x in names):
            masked = True
            break

    lens = [len(table) for table in tables]
    n_rows = sum(lens)
    out = _get_out_class(tables)(masked=masked)
    col_descrs = get_descrs(tables, col_name_map)
    cols = []
    for (name, dtype, shape), out_attrs in zip(col_descrs, col_attrs_map.values()):
        kwargs = dict(name=name, dtype=dtype, shape=shape, length=n_rows)
        if masked:
            kwargs['mask'] = True

        kwargs.update(**out_attrs)
        cols.append(out.ColumnClass(**kwargs))
    out.add_columns(cols, copy=False)

    for out_name, in_names in six.iteritems(col_name_map):
        idx0 = 0
        for name, table in zip(in_names, tables):
            idx1 = idx0 + len(table)
            if name in table.colnames:
                out[out_name][idx0:idx1] = table[name]
            idx0 = idx1

    return out


def _hstack(tables, col_name_map):
    """Stack tables horizontally (by columns)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be the output.  A value of 'outer' means
    the output will have the union of all rows, with array values being
    masked where no common values are available.

    Parameters
    ----------
    tables : List of tables
        Tables to stack by columns (horizontally)
    col_name_map : dict
        Mapping of output to input column names.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.
    """
    table_lens = [len(table) for table in tables]
    # If there are any output rows where one or more input arrays are missing
    # then the output must be masked.  If any input arrays are masked then
    # output is masked.
    masked = (any(getattr(table, 'masked', False) for table in tables) or
              len(set(table_lens)) > 1)

    n_rows = max(table_lens)
    out = _get_out_class(tables)(masked=masked)

    for out_name, in_names in six.iteritems(col_name_map):
        for name, table, table_len in zip(in_names, tables, table_lens):
            if name is None:
                continue

            if n_rows > table_len:
                indices = np.arange(n_rows)
                indices[table_len:] = 0
                out[out_name] = table[name][indices]
                try:
                    out[out_name].mask[table_len:] = True
                except ValueError:
                    raise ValueError("hstack requires masking column '{0}' but column"
                                     " type {1} does not support masking"
                                     .format(out_name, out[out_name].__class__.__name__))
            else:
                out[out_name] = table[name][:n_rows]

    return out
