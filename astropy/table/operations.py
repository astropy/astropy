"""
High-level table operations:

- join()
- setdiff()
- hstack()
- vstack()
- dstack()
"""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy
import collections
import itertools
from collections import OrderedDict, Counter
from collections.abc import Mapping, Sequence

import numpy as np

from astropy.utils import metadata
from .table import Table, QTable, Row, Column, MaskedColumn
from astropy.units import Quantity
from astropy.utils.compat import NUMPY_LT_1_17

from . import _np_utils
from .np_utils import fix_column_name, TableMergeError

__all__ = ['join', 'setdiff', 'hstack', 'vstack', 'unique']


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

    # Make sure we have a list of things
    if not isinstance(tables, Sequence):
        tables = [tables]

    # Make sure there is something to stack
    if len(tables) == 0:
        raise ValueError('no values provided to stack.')

    # Convert inputs (Table, Row, or anything column-like) to Tables.
    # Special case that Quantity converts to a QTable.
    for ii, val in enumerate(tables):
        if isinstance(val, Table):
            pass
        elif isinstance(val, Row):
            tables[ii] = Table(val)
        elif isinstance(val, Quantity):
            tables[ii] = QTable([val])
        else:
            try:
                tables[ii] = Table([val])
            except (ValueError, TypeError):
                raise TypeError('cannot convert {} to table column.'
                                .format(val))

    return tables


def _get_out_class(objs):
    """
    From a list of input objects ``objs`` get merged output object class.

    This is just taken as the deepest subclass. This doesn't handle complicated
    inheritance schemes.
    """
    out_class = objs[0].__class__
    for obj in objs[1:]:
        if issubclass(obj.__class__, out_class):
            out_class = obj.__class__

    if any(not issubclass(out_class, obj.__class__) for obj in objs):
        raise ValueError('unmergeable object classes {}'
                         .format([obj.__class__.__name__ for obj in objs]))

    return out_class


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
        Join type ('inner' | 'outer' | 'left' | 'right' | 'cartesian'), default is 'inner'
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

    Returns
    -------
    joined_table : `~astropy.table.Table` object
        New table containing the result of the join operation.
    """

    # Try converting inputs to Table as needed
    if not isinstance(left, Table):
        left = Table(left)
    if not isinstance(right, Table):
        right = Table(right)

    col_name_map = OrderedDict()
    out = _join(left, right, keys, join_type,
                uniq_col_name, table_names, col_name_map, metadata_conflicts)

    # Merge the column and table meta data. Table subclasses might override
    # these methods for custom merge behavior.
    _merge_table_meta(out, [left, right], metadata_conflicts=metadata_conflicts)

    return out


def setdiff(table1, table2, keys=None):
    """
    Take a set difference of table rows.

    The row set difference will contain all rows in ``table1`` that are not
    present in ``table2``. If the keys parameter is not defined, all columns in
    ``table1`` will be included in the output table.

    Parameters
    ----------
    table1 : `~astropy.table.Table`
        ``table1`` is on the left side of the set difference.
    table2 : `~astropy.table.Table`
        ``table2`` is on the right side of the set difference.
    keys : str or list of str
        Name(s) of column(s) used to match rows of left and right tables.
        Default is to use all columns in ``table1``.

    Returns
    -------
    diff_table : `~astropy.table.Table`
        New table containing the set difference between tables. If the set
        difference is none, an empty table will be returned.

    Examples
    --------
    To get a set difference between two tables::

      >>> from astropy.table import setdiff, Table
      >>> t1 = Table({'a': [1, 4, 9], 'b': ['c', 'd', 'f']}, names=('a', 'b'))
      >>> t2 = Table({'a': [1, 5, 9], 'b': ['c', 'b', 'f']}, names=('a', 'b'))
      >>> print(t1)
       a   b
      --- ---
        1   c
        4   d
        9   f
      >>> print(t2)
       a   b
      --- ---
        1   c
        5   b
        9   f
      >>> print(setdiff(t1, t2))
       a   b
      --- ---
        4   d

      >>> print(setdiff(t2, t1))
       a   b
      --- ---
        5   b
    """
    if keys is None:
        keys = table1.colnames

    #Check that all keys are in table1 and table2
    for tbl, tbl_str in ((table1,'table1'), (table2,'table2')):
        diff_keys = np.setdiff1d(keys, tbl.colnames)
        if len(diff_keys) != 0:
            raise ValueError("The {} columns are missing from {}, cannot take "
                             "a set difference.".format(diff_keys, tbl_str))

    # Make a light internal copy of both tables
    t1 = table1.copy(copy_data=False)
    t1.meta = {}
    t1.keep_columns(keys)
    t1['__index1__'] = np.arange(len(table1))  # Keep track of rows indices

    # Make a light internal copy to avoid touching table2
    t2 = table2.copy(copy_data=False)
    t2.meta = {}
    t2.keep_columns(keys)
    # Dummy column to recover rows after join
    t2['__index2__'] = np.zeros(len(t2), dtype=np.uint8)  # dummy column

    t12 = _join(t1, t2, join_type='left', keys=keys,
                metadata_conflicts='silent')

    # If t12 index2 is masked then that means some rows were in table1 but not table2.
    if hasattr(t12['__index2__'], 'mask'):
        # Define bool mask of table1 rows not in table2
        diff = t12['__index2__'].mask
        # Get the row indices of table1 for those rows
        idx = t12['__index1__'][diff]
        # Select corresponding table1 rows straight from table1 to ensure
        # correct table and column types.
        t12_diff = table1[idx]
    else:
        t12_diff = table1[[]]

    return t12_diff


def dstack(tables, join_type='outer', metadata_conflicts='warn'):
    """
    Stack columns within tables depth-wise

    A ``join_type`` of 'exact' means that the tables must all have exactly
    the same column names (though the order can vary).  If ``join_type``
    is 'inner' then the intersection of common columns will be the output.
    A value of 'outer' (default) means the output will have the union of
    all columns, with table values being masked where no common values are
    available.

    Parameters
    ----------
    tables : Table or list of Table objects
        Table(s) to stack along depth-wise with the current table
        Table columns should have same shape and name for depth-wise stacking
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.

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
      >>> print(dstack([t1, t2]))
      a [2]  b [2]
      ------ ------
      1 .. 5 3 .. 7
      2 .. 6 4 .. 8
    """
    tables = _get_list_of_tables(tables)
    if len(tables) == 1:
        return tables[0]  # no point in stacking a single table

    n_rows = set(len(table) for table in tables)
    if len(n_rows) != 1:
        raise ValueError('Table lengths must all match for dstack')
    n_row = n_rows.pop()

    out = vstack(tables, join_type, metadata_conflicts)

    for name, col in out.columns.items():
        col = out[name]

        # Reshape to so each original column is now in a row.
        # If entries are not 0-dim then those additional shape dims
        # are just carried along.
        # [x x x y y y] => [[x x x],
        #                   [y y y]]
        col.shape = (len(tables), n_row) + col.shape[1:]

        # Transpose the table and row axes to get to
        # [[x, y],
        #  [x, y]
        #  [x, y]]
        axes = np.arange(len(col.shape))
        axes[:2] = [1, 0]

        # This temporarily makes `out` be corrupted (columns of different
        # length) but it all works out in the end.
        out.columns.__setitem__(name, col.transpose(axes), validated=True)

    return out


def vstack(tables, join_type='outer', metadata_conflicts='warn'):
    """
    Stack tables vertically (along rows)

    A ``join_type`` of 'exact' means that the tables must all have exactly
    the same column names (though the order can vary).  If ``join_type``
    is 'inner' then the intersection of common columns will be the output.
    A value of 'outer' (default) means the output will have the union of
    all columns, with table values being masked where no common values are
    available.

    Parameters
    ----------
    tables : Table or list of Table objects
        Table(s) to stack along rows (vertically) with the current table
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.

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
    tables = _get_list_of_tables(tables)  # validates input
    if len(tables) == 1:
        return tables[0]  # no point in stacking a single table
    col_name_map = OrderedDict()

    out = _vstack(tables, join_type, col_name_map, metadata_conflicts)

    # Merge table metadata
    _merge_table_meta(out, tables, metadata_conflicts=metadata_conflicts)

    return out


def hstack(tables, join_type='outer',
           uniq_col_name='{col_name}_{table_name}', table_names=None,
           metadata_conflicts='warn'):
    """
    Stack tables along columns (horizontally)

    A ``join_type`` of 'exact' means that the tables must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be the output.  A value of 'outer' (default)
    means the output will have the union of all rows, with table values being
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
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.

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
    tables = _get_list_of_tables(tables)  # validates input
    if len(tables) == 1:
        return tables[0]  # no point in stacking a single table
    col_name_map = OrderedDict()

    out = _hstack(tables, join_type, uniq_col_name, table_names,
                  col_name_map)

    _merge_table_meta(out, tables, metadata_conflicts=metadata_conflicts)

    return out


def unique(input_table, keys=None, silent=False, keep='first'):
    """
    Returns the unique rows of a table.

    Parameters
    ----------

    input_table : `~astropy.table.Table` object or a value that
        will initialize a `~astropy.table.Table` object
    keys : str or list of str
        Name(s) of column(s) used to create unique rows.
        Default is to use all columns.
    keep : one of 'first', 'last' or 'none'
        Whether to keep the first or last row for each set of
        duplicates. If 'none', all rows that are duplicate are
        removed, leaving only rows that are already unique in
        the input.
        Default is 'first'.
    silent : bool
        If `True`, masked value column(s) are silently removed from
        ``keys``. If `False`, an exception is raised when ``keys``
        contains masked value column(s).
        Default is `False`.

    Returns
    -------
    unique_table : `~astropy.table.Table` object
        New table containing only the unique rows of ``input_table``.

    Examples
    --------
    >>> from astropy.table import unique, Table
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

    if isinstance(keys, str):
        keys = [keys]
    if keys is None:
        keys = input_table.colnames
    else:
        if len(set(keys)) != len(keys):
            raise ValueError("duplicate key names")

    # Check for columns with masked values
    nkeys = 0
    for key in keys[:]:
        col = input_table[key]
        if hasattr(col, 'mask') and np.any(col.mask):
            if not silent:
                raise ValueError(
                    "cannot use columns with masked values as keys; "
                    "remove column '{}' from keys and rerun "
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


def get_col_name_map(arrays, common_names, uniq_col_name='{col_name}_{table_name}',
                     table_names=None):
    """
    Find the column names mapping when merging the list of tables
    ``arrays``.  It is assumed that col names in ``common_names`` are to be
    merged into a single column while the rest will be uniquely represented
    in the output.  The args ``uniq_col_name`` and ``table_names`` specify
    how to rename columns in case of conflicts.

    Returns a dict mapping each output column name to the input(s).  This takes the form
    {outname : (col_name_0, col_name_1, ...), ... }.  For key columns all of input names
    will be present, while for the other non-key columns the value will be (col_name_0,
    None, ..) or (None, col_name_1, ..) etc.
    """

    col_name_map = collections.defaultdict(lambda: [None] * len(arrays))
    col_name_list = []

    if table_names is None:
        table_names = [str(ii + 1) for ii in range(len(arrays))]

    for idx, array in enumerate(arrays):
        table_name = table_names[idx]
        for name in array.colnames:
            out_name = name

            if name in common_names:
                # If name is in the list of common_names then insert into
                # the column name list, but just once.
                if name not in col_name_list:
                    col_name_list.append(name)
            else:
                # If name is not one of the common column outputs, and it collides
                # with the names in one of the other arrays, then rename
                others = list(arrays)
                others.pop(idx)
                if any(name in other.colnames for other in others):
                    out_name = uniq_col_name.format(table_name=table_name, col_name=name)
                col_name_list.append(out_name)

            col_name_map[out_name][idx] = name

    # Check for duplicate output column names
    col_name_count = Counter(col_name_list)
    repeated_names = [name for name, count in col_name_count.items() if count > 1]
    if repeated_names:
        raise TableMergeError('Merging column names resulted in duplicates: {}.  '
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

    for out_name, in_names in col_name_map.items():
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
            raise TableMergeError("The '{}' columns have incompatible types: {}"
                                  .format(names[0], tme._incompat_types))

        # Make sure all input shapes are the same
        uniq_shapes = set(col.shape[1:] for col in in_cols)
        if len(uniq_shapes) != 1:
            raise TableMergeError(f'Key columns {names!r} have different shape')
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
        tme = TableMergeError('Columns have incompatible types {}'
                              .format(err._incompat_types))
        tme._incompat_types = err._incompat_types
        raise tme


def _get_join_sort_idxs(keys, left, right):
    # Go through each of the key columns in order and make columns for
    # a new structured array that represents the lexical ordering of those
    # key columns. This structured array is then argsort'ed. The trick here
    # is that some columns (e.g. Time) may need to be expanded into multiple
    # columns for ordering here.

    ii = 0  # Index for uniquely naming the sort columns
    sort_keys_dtypes = []  # sortable_table dtypes as list of (name, dtype_str, shape) tuples
    sort_keys = []  # sortable_table (structured ndarray) column names
    sort_left = {}  # sortable ndarrays from left table
    sort_right = {}  # sortable ndarray from right table

    for key in keys:
        # get_sortable_arrays() returns a list of ndarrays that can be lexically
        # sorted to represent the order of the column. In most cases this is just
        # a single element of the column itself.
        left_sort_cols = left[key].info.get_sortable_arrays()
        right_sort_cols = right[key].info.get_sortable_arrays()

        if len(left_sort_cols) != len(right_sort_cols):
            # Should never happen because cols are screened beforehand for compatibility
            raise RuntimeError('mismatch in sort cols lengths')

        for left_sort_col, right_sort_col in zip(left_sort_cols, right_sort_cols):
            # Check for consistency of shapes. Mismatch should never happen.
            shape = left_sort_col.shape[1:]
            if shape != right_sort_col.shape[1:]:
                raise RuntimeError('mismatch in shape of left vs. right sort array')

            if shape != ():
                raise ValueError(f'sort key column {key!r} must be 1-d')

            sort_key = str(ii)
            sort_keys.append(sort_key)
            sort_left[sort_key] = left_sort_col
            sort_right[sort_key] = right_sort_col

            # Build up dtypes for the structured array that gets sorted.
            dtype_str = common_dtype([left_sort_col, right_sort_col])
            sort_keys_dtypes.append((sort_key, dtype_str))
            ii += 1

    # Make the empty sortable table and fill it
    len_left = len(left)
    sortable_table = np.empty(len_left + len(right), dtype=sort_keys_dtypes)
    for key in sort_keys:
        sortable_table[key][:len_left] = sort_left[key]
        sortable_table[key][len_left:] = sort_right[key]

    # Finally do the (lexical) argsort and make a new sorted version
    idx_sort = sortable_table.argsort(order=sort_keys)
    sorted_table = sortable_table[idx_sort]

    # Get indexes of unique elements (i.e. the group boundaries)
    diffs = np.concatenate(([True], sorted_table[1:] != sorted_table[:-1], [True]))
    idxs = np.flatnonzero(diffs)

    return idxs, idx_sort


def _join(left, right, keys=None, join_type='inner',
          uniq_col_name='{col_name}_{table_name}',
          table_names=['1', '2'],
          col_name_map=None, metadata_conflicts='warn'):
    """
    Perform a join of the left and right Tables on specified keys.

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
        Join type ('inner' | 'outer' | 'left' | 'right' | 'cartesian'), default is 'inner'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2'].
    col_name_map : empty dict or None
        If passed as a dict then it will be updated in-place with the
        mapping of output to input column names.
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.

    Returns
    -------
    joined_table : `~astropy.table.Table` object
        New table containing the result of the join operation.
    """
    from astropy.time import Time

    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    # Special column name for cartesian join, should never collide with real column
    cartesian_index_name = '__table_cartesian_join_temp_index__'

    if join_type not in ('inner', 'outer', 'left', 'right', 'cartesian'):
        raise ValueError("The 'join_type' argument should be in 'inner', "
                         "'outer', 'left', 'right', or 'cartesian' "
                         "(got '{}' instead)".
                         format(join_type))

    if join_type == 'cartesian':
        if keys:
            raise ValueError('cannot supply keys for a cartesian join')

        # Make light copies of left and right, then add temporary index columns
        # with all the same value so later an outer join turns into a cartesian join.
        left = left.copy(copy_data=False)
        right = right.copy(copy_data=False)
        left[cartesian_index_name] = np.uint8(0)
        right[cartesian_index_name] = np.uint8(0)
        keys = (cartesian_index_name, )

    # If we have a single key, put it in a tuple
    if keys is None:
        keys = tuple(name for name in left.colnames if name in right.colnames)
        if len(keys) == 0:
            raise TableMergeError('No keys in common between left and right tables')
    elif isinstance(keys, str):
        keys = (keys,)

    # Check the key columns
    for arr, arr_label in ((left, 'Left'), (right, 'Right')):
        for name in keys:
            if name not in arr.colnames:
                raise TableMergeError('{} table does not have key column {!r}'
                                      .format(arr_label, name))
            if hasattr(arr[name], 'mask') and np.any(arr[name].mask):
                raise TableMergeError('{} key column {!r} has missing values'
                                      .format(arr_label, name))

    len_left, len_right = len(left), len(right)

    if len_left == 0 or len_right == 0:
        raise ValueError('input tables for join must both have at least one row')

    # Joined array dtype as a list of descr (name, type_str, shape) tuples
    col_name_map = get_col_name_map([left, right], keys, uniq_col_name, table_names)
    out_descrs = get_descrs([left, right], col_name_map)

    try:
        idxs, idx_sort = _get_join_sort_idxs(keys, left, right)
    except NotImplementedError:
        raise TypeError('one or more key columns are not sortable')

    # Main inner loop in Cython to compute the cartesian product
    # indices for the given join type
    int_join_type = {'inner': 0, 'outer': 1, 'left': 2, 'right': 3,
                     'cartesian': 1}[join_type]
    masked, n_out, left_out, left_mask, right_out, right_mask = \
        _np_utils.join_inner(idxs, idx_sort, len_left, int_join_type)

    out = _get_out_class([left, right])()

    for out_name, dtype, shape in out_descrs:
        if out_name == cartesian_index_name:
            continue

        left_name, right_name = col_name_map[out_name]
        if left_name and right_name:  # this is a key which comes from left and right
            cols = [left[left_name], right[right_name]]

            col_cls = _get_out_class(cols)
            if not hasattr(col_cls.info, 'new_like'):
                raise NotImplementedError('join unavailable for mixin column type(s): {}'
                                          .format(col_cls.__name__))

            out[out_name] = col_cls.info.new_like(cols, n_out, metadata_conflicts, out_name)

            if not NUMPY_LT_1_17 or issubclass(col_cls, (Column, Time)):
                out[out_name][:] = np.where(right_mask,
                                            left[left_name].take(left_out),
                                            right[right_name].take(right_out))
            else:
                # np.where does not work for mixin columns (e.g. Quantity) so
                # use a slower workaround.
                non_right_mask = ~right_mask
                if np.any(right_mask):
                    out[out_name][:] = left[left_name].take(left_out)
                if np.any(non_right_mask):
                    out[out_name][non_right_mask] = right[right_name].take(right_out)[non_right_mask]
            continue
        elif left_name:  # out_name came from the left table
            name, array, array_out, array_mask = left_name, left, left_out, left_mask
        elif right_name:
            name, array, array_out, array_mask = right_name, right, right_out, right_mask
        else:
            raise TableMergeError('Unexpected column names (maybe one is ""?)')

        # Select the correct elements from the original table
        col = array[name][array_out]

        # If the output column is masked then set the output column masking
        # accordingly.  Check for columns that don't support a mask attribute.
        if masked and np.any(array_mask):
            # If col is a Column but not MaskedColumn then upgrade at this point
            # because masking is required.
            if isinstance(col, Column) and not isinstance(col, MaskedColumn):
                col = out.MaskedColumn(col, copy=False)

            # array_mask is 1-d corresponding to length of output column.  We need
            # make it have the correct shape for broadcasting, i.e. (length, 1, 1, ..).
            # Mixin columns might not have ndim attribute so use len(col.shape).
            array_mask.shape = (col.shape[0],) + (1,) * (len(col.shape) - 1)

            # Now broadcast to the correct final shape
            array_mask = np.broadcast_to(array_mask, col.shape)

            try:
                col[array_mask] = col.info.mask_val
            except Exception:  # Not clear how different classes will fail here
                raise NotImplementedError(
                    "join requires masking column '{}' but column"
                    " type {} does not support masking"
                    .format(out_name, col.__class__.__name__))

        # Set the output table column to the new joined column
        out[out_name] = col

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, Mapping):
        _col_name_map.update(col_name_map)

    return out


def _vstack(arrays, join_type='outer', col_name_map=None, metadata_conflicts='warn'):
    """
    Stack Tables vertically (by rows)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same column names (though the order can vary).  If
    ``join_type`` is 'inner' then the intersection of common columns will
    be the output.  A value of 'outer' means the output will have the union of
    all columns, with array values being masked where no common values are
    available.

    Parameters
    ----------
    arrays : list of Tables
        Tables to stack by rows (vertically)
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    col_name_map : empty dict or None
        If passed as a dict then it will be updated in-place with the
        mapping of output to input column names.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.
    """
    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    # Input validation
    if join_type not in ('inner', 'exact', 'outer'):
        raise ValueError("`join_type` arg must be one of 'inner', 'exact' or 'outer'")

    # Trivial case of one input array
    if len(arrays) == 1:
        return arrays[0]

    # Start by assuming an outer match where all names go to output
    names = set(itertools.chain(*[arr.colnames for arr in arrays]))
    col_name_map = get_col_name_map(arrays, names)

    # If require_match is True then the output must have exactly the same
    # number of columns as each input array
    if join_type == 'exact':
        for names in col_name_map.values():
            if any(x is None for x in names):
                raise TableMergeError('Inconsistent columns in input arrays '
                                      "(use 'inner' or 'outer' join_type to "
                                      "allow non-matching columns)")
        join_type = 'outer'

    # For an inner join, keep only columns where all input arrays have that column
    if join_type == 'inner':
        col_name_map = OrderedDict((name, in_names) for name, in_names in col_name_map.items()
                                   if all(x is not None for x in in_names))
        if len(col_name_map) == 0:
            raise TableMergeError('Input arrays have no columns in common')

    lens = [len(arr) for arr in arrays]
    n_rows = sum(lens)
    out = _get_out_class(arrays)()

    for out_name, in_names in col_name_map.items():
        # List of input arrays that contribute to this output column
        cols = [arr[name] for arr, name in zip(arrays, in_names) if name is not None]

        col_cls = _get_out_class(cols)
        if not hasattr(col_cls.info, 'new_like'):
            raise NotImplementedError('vstack unavailable for mixin column type(s): {}'
                                      .format(col_cls.__name__))
        try:
            col = col_cls.info.new_like(cols, n_rows, metadata_conflicts, out_name)
        except metadata.MergeConflictError as err:
            # Beautify the error message when we are trying to merge columns with incompatible
            # types by including the name of the columns that originated the error.
            raise TableMergeError("The '{}' columns have incompatible types: {}"
                                  .format(out_name, err._incompat_types))

        idx0 = 0
        for name, array in zip(in_names, arrays):
            idx1 = idx0 + len(array)
            if name in array.colnames:
                col[idx0:idx1] = array[name]
            else:
                # If col is a Column but not MaskedColumn then upgrade at this point
                # because masking is required.
                if isinstance(col, Column) and not isinstance(col, MaskedColumn):
                    col = out.MaskedColumn(col, copy=False)

                try:
                    col[idx0:idx1] = col.info.mask_val
                except Exception:
                    raise NotImplementedError(
                        "vstack requires masking column '{}' but column"
                        " type {} does not support masking"
                        .format(out_name, col.__class__.__name__))
            idx0 = idx1

        out[out_name] = col

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, Mapping):
        _col_name_map.update(col_name_map)

    return out


def _hstack(arrays, join_type='outer', uniq_col_name='{col_name}_{table_name}',
           table_names=None, col_name_map=None):
    """
    Stack tables horizontally (by columns)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be the output.  A value of 'outer' means
    the output will have the union of all rows, with array values being
    masked where no common values are available.

    Parameters
    ----------
    arrays : List of tables
        Tables to stack by columns (horizontally)
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2', ..].

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.
    """

    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    # Input validation
    if join_type not in ('inner', 'exact', 'outer'):
        raise ValueError("join_type arg must be either 'inner', 'exact' or 'outer'")

    if table_names is None:
        table_names = ['{}'.format(ii + 1) for ii in range(len(arrays))]
    if len(arrays) != len(table_names):
        raise ValueError('Number of arrays must match number of table_names')

    # Trivial case of one input arrays
    if len(arrays) == 1:
        return arrays[0]

    col_name_map = get_col_name_map(arrays, [], uniq_col_name, table_names)

    # If require_match is True then all input arrays must have the same length
    arr_lens = [len(arr) for arr in arrays]
    if join_type == 'exact':
        if len(set(arr_lens)) > 1:
            raise TableMergeError("Inconsistent number of rows in input arrays "
                                  "(use 'inner' or 'outer' join_type to allow "
                                  "non-matching rows)")
        join_type = 'outer'

    # For an inner join, keep only the common rows
    if join_type == 'inner':
        min_arr_len = min(arr_lens)
        if len(set(arr_lens)) > 1:
            arrays = [arr[:min_arr_len] for arr in arrays]
        arr_lens = [min_arr_len for arr in arrays]

    # If there are any output rows where one or more input arrays are missing
    # then the output must be masked.  If any input arrays are masked then
    # output is masked.

    n_rows = max(arr_lens)
    out = _get_out_class(arrays)()

    for out_name, in_names in col_name_map.items():
        for name, array, arr_len in zip(in_names, arrays, arr_lens):
            if name is None:
                continue

            if n_rows > arr_len:
                indices = np.arange(n_rows)
                indices[arr_len:] = 0
                col = array[name][indices]

                # If col is a Column but not MaskedColumn then upgrade at this point
                # because masking is required.
                if isinstance(col, Column) and not isinstance(col, MaskedColumn):
                    col = out.MaskedColumn(col, copy=False)

                try:
                    col[arr_len:] = col.info.mask_val
                except Exception:
                    raise NotImplementedError(
                        "hstack requires masking column '{}' but column"
                        " type {} does not support masking"
                        .format(out_name, col.__class__.__name__))
            else:
                col = array[name][:n_rows]

            out[out_name] = col

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, Mapping):
        _col_name_map.update(col_name_map)

    return out
