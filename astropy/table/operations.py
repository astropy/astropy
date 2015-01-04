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
from ..extern.six.moves import zip

from copy import deepcopy
import warnings
import collections
import itertools

import numpy as np
from numpy import ma

from ..utils import OrderedDict, metadata

from . import _np_utils
from .np_utils import fix_column_name, TableMergeError

__all__ = ['join', 'hstack', 'vstack', 'unique']


def _merge_col_meta(out, tables, col_name_map, idx_left=0, idx_right=1,
                    metadata_conflicts='warn'):
    """
    Merge column meta data for the ``out`` table.

    This merges column meta, which includes attributes unit, format,
    and description, as well as the actual `meta` attribute.  It is
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
                out_col.meta = metadata.merge(left_col.meta, right_col.meta,
                                              metadata_conflicts=metadata_conflicts)
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
                            raise metadata.MergeConflictError(
                                'In merged column {0!r} the {1!r} attribute does not match '
                                '({2} != {3})'.format(out_col.name, attr, left_attr, right_attr))
                        elif metadata_conflicts != 'silent':
                            raise ValueError('metadata_conflict argument must be one of "silent",'
                                             ' "warn", or "error"')
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

    col_name_map = OrderedDict()
    out = _join(left, right, keys, join_type,
                uniq_col_name, table_names, col_name_map)

    # Merge the column and table meta data. Table subclasses might override
    # these methods for custom merge behavior.
    _merge_col_meta(out, [left, right], col_name_map, metadata_conflicts=metadata_conflicts)
    _merge_table_meta(out, [left, right], metadata_conflicts=metadata_conflicts)

    return out


def vstack(tables, join_type='outer', metadata_conflicts='warn'):
    """
    Stack tables vertically (along rows)

    A ``join_type`` of 'exact' means that the tables must all have exactly
    the same column names (though the order can vary).  If ``join_type``
    is 'inner' then the intersection of common columns will be output.
    A value of 'outer' (default) means the output will have the union of
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
    col_name_map = OrderedDict()

    out = _vstack(tables, join_type, col_name_map)

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
    col_name_map = OrderedDict()

    out = _hstack(tables, join_type, uniq_col_name, table_names,
                  col_name_map)

    _merge_col_meta(out, tables, col_name_map, metadata_conflicts=metadata_conflicts)
    _merge_table_meta(out, tables, metadata_conflicts=metadata_conflicts)

    return out


def unique(input_table, keys=None, silent=False):
    """
    Returns the unique rows of a table.

    Parameters
    ----------

    input_table : `~astropy.table.Table` object or a value that
    will initialize a `~astropy.table.Table` object
        Input table.
    keys : str or list of str
        Name(s) of column(s) used to unique rows.
        Default is to use all columns.
    silent : boolean
        If `True` masked value column(s) are silently removed from
        ``keys``. If `False` an exception is raised when ``keys`` contains
        masked value column(s).
        Default is `False`.

    Returns
    -------
    unique_table : `~astropy.table.Table` object
        Table containing only the unique rays of ``input_table``.

    """

    if keys is None:
        keys = input_table.colnames

    if input_table.masked:
        if isinstance(keys, six.string_types):
            keys = [keys, ]
        for i, key in enumerate(keys):
            if np.any(input_table[key].mask):
                if not silent:
                    raise ValueError("Cannot unique masked value key columns, "
                                     "remove column '{0}' from keys and rerun "
                                     "unique.".format(key))
                del keys[i]
        if len(keys) == 0:
            raise ValueError("No column remained in ``keys``, unique cannot "
                             "work with masked value key columns.")

    grouped_table = input_table.group_by(keys)
    unique_table = grouped_table[grouped_table.groups.indices[:-1]]

    return unique_table


def _counter(iterable):
    """
    Count instances of each unique value in ``iterable``.  Returns a dict
    with the counts.  Would use collections.Counter but this isn't available in 2.6.
    """
    counts = collections.defaultdict(int)
    for val in iterable:
        counts[val] += 1
    return counts


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
        table_names = [six.text_type(ii + 1) for ii in range(len(arrays))]

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
    col_name_count = _counter(col_name_list)
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
            raise TableMergeError('Key columns {0!r} have different shape'.format(name))
        shape = uniq_shapes.pop()

        out_descrs.append((fix_column_name(out_name), dtype, shape))

    return out_descrs


def common_dtype(cols):
    """
    Use numpy to find the common dtype for a list of columns.

    Only allow columns within the following fundamental numpy data types:
    np.bool_, np.object_, np.number, np.character, np.void
    """
    np_types = (np.bool_, np.object_, np.number, np.character, np.void)
    uniq_types = set(tuple(issubclass(col.dtype.type, np_type) for np_type in np_types)
                     for col in cols)
    if len(uniq_types) > 1:
        # Embed into the exception the actual list of incompatible types.
        incompat_types = [col.dtype.name for col in cols]
        tme = TableMergeError('Columns have incompatible types {0}'
                              .format(incompat_types))
        tme._incompat_types = incompat_types
        raise tme

    arrs = [np.empty(1, dtype=col.dtype) for col in cols]

    # For string-type arrays need to explicitly fill in non-zero
    # values or the final arr_common = .. step is unpredictable.
    for arr in arrs:
        if arr.dtype.kind in ('S', 'U'):
            arr[0] = '0' * arr.itemsize

    arr_common = np.array([arr[0] for arr in arrs])
    return arr_common.dtype.str


def _join(left, right, keys=None, join_type='inner',
         uniq_col_name='{col_name}_{table_name}',
         table_names=['1', '2'],
         col_name_map=None):
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
        Join type ('inner' | 'outer' | 'left' | 'right'), default is 'inner'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2'].
    col_name_map : empty dict or None
        If passed as a dict then it will be updated in-place with the
        mapping of output to input column names.

    Returns
    -------
    joined_table : `~astropy.table.Table` object
        New table containing the result of the join operation.
    """
    from .table import Table

    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    if join_type not in ('inner', 'outer', 'left', 'right'):
        raise ValueError("The 'join_type' argument should be in 'inner', "
                         "'outer', 'left' or 'right' (got '{0}' instead)".
                         format(join_type))

    # If we have a single key, put it in a tuple
    if keys is None:
        keys = tuple(name for name in left.colnames if name in right.colnames)
        if len(keys) == 0:
            raise TableMergeError('No keys in common between left and right tables')
    elif isinstance(keys, six.string_types):
        keys = (keys,)

    # Check the key columns
    for arr, arr_label in ((left, 'Left'), (right, 'Right')):
        for name in keys:
            if name not in arr.colnames:
                raise TableMergeError('{0} table does not have key column {1!r}'
                                      .format(arr_label, name))
            if hasattr(arr[name], 'mask') and np.any(arr[name].mask):
                raise TableMergeError('{0} key column {1!r} has missing values'
                                      .format(arr_label, name))

    len_left, len_right = len(left), len(right)

    if len_left == 0 or len_right == 0:
        raise ValueError('input tables for join must both have at least one row')


    # Joined array dtype as a list of descr (name, type_str, shape) tuples
    col_name_map = get_col_name_map([left, right], keys, uniq_col_name, table_names)
    out_descrs = get_descrs([left, right], col_name_map)

    # Make an array with just the key columns.  This uses a temporary
    # structured array for efficiency.
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
    if left.masked or right.masked:
        masked = True
    masked = bool(masked)

    out = Table(masked=masked)

    for out_name, dtype, shape in out_descrs:
        out[out_name] = out.ColumnClass(length=n_out, name=out_name, dtype=dtype, shape=shape)

        left_name, right_name = col_name_map[out_name]
        if left_name and right_name:  # this is a key which comes from left and right
            out[out_name] = np.where(right_mask,
                                     left[left_name].take(left_out),
                                     right[right_name].take(right_out))
            continue
        elif left_name:  # out_name came from the left table
            name, array, array_out, array_mask = left_name, left, left_out, left_mask
        elif right_name:
            name, array, array_out, array_mask = right_name, right, right_out, right_mask
        else:
            raise TableMergeError('Unexpected column names (maybe one is ""?)')

        out[out_name] = array[name].take(array_out, axis=0)
        if masked:
            if array.masked:
                array_mask = array_mask | array[name].mask.take(array_out)
            out[out_name].mask = array_mask

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, collections.Mapping):
        _col_name_map.update(col_name_map)

    return out


def _vstack(arrays, join_type='inner', col_name_map=None):
    """
    Stack Tables vertically (by rows)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same column names (though the order can vary).  If
    ``join_type`` is 'inner' then the intersection of common columns will
    be output.  A value of 'outer' means the output will have the union of
    all columns, with array values being masked where no common values are
    available.

    Parameters
    ----------
    arrays : list of Tables
        Tables to stack by rows (vertically)
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'exact'
    col_name_map : empty dict or None
        If passed as a dict then it will be updated in-place with the
        mapping of output to input column names.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.
    """
    from .table import Table

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
        for names in six.itervalues(col_name_map):
            if any(x is None for x in names):
                raise TableMergeError('Inconsistent columns in input arrays '
                                      "(use 'inner' or 'outer' join_type to "
                                      "allow non-matching columns)")
        join_type = 'outer'

    # For an inner join, keep only columns where all input arrays have that column
    if join_type == 'inner':
        col_name_map = OrderedDict((name, in_names) for name, in_names in six.iteritems(col_name_map)
                                   if all(x is not None for x in in_names))
        if len(col_name_map) == 0:
            raise TableMergeError('Input arrays have no columns in common')

    # If there are any output columns where one or more input arrays are missing
    # then the output must be masked.  If any input arrays are masked then
    # output is masked.
    masked = any(getattr(arr, 'masked', False) for arr in arrays)
    for names in six.itervalues(col_name_map):
        if any(x is None for x in names):
            masked = True
            break

    lens = [len(arr) for arr in arrays]
    n_rows = sum(lens)
    out = Table(masked=masked)
    out_descrs = get_descrs(arrays, col_name_map)
    for out_descr in out_descrs:
        name = out_descr[0]
        dtype = out_descr[1:]
        if masked:
            out[name] = ma.array(data=np.zeros(n_rows, dtype),
                                 mask=np.ones(n_rows, ma.make_mask_descr(dtype)))
        else:
            out[name] = np.empty(n_rows, dtype=dtype)

    for out_name, in_names in six.iteritems(col_name_map):
        idx0 = 0
        for name, array in zip(in_names, arrays):
            idx1 = idx0 + len(array)
            if name in array.colnames:
                out[out_name][idx0:idx1] = array[name]
            idx0 = idx1

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, collections.Mapping):
        _col_name_map.update(col_name_map)

    return out


def _hstack(arrays, join_type='exact', uniq_col_name='{col_name}_{table_name}',
           table_names=None, col_name_map=None):
    """
    Stack tables horizontally (by columns)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be output.  A value of 'outer' means
    the output will have the union of all rows, with array values being
    masked where no common values are available.

    Parameters
    ----------
    arrays : List of tables
        Tables to stack by columns (horizontally)
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'exact'
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
    from .table import Table

    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    # Input validation
    if join_type not in ('inner', 'exact', 'outer'):
        raise ValueError("join_type arg must be either 'inner', 'exact' or 'outer'")

    if table_names is None:
        table_names = ['{0}'.format(ii + 1) for ii in range(len(arrays))]
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
    masked = any(getattr(arr, 'masked', False) for arr in arrays) or len(set(arr_lens)) > 1

    n_rows = max(arr_lens)
    out = Table(masked=masked)
    out_descrs = get_descrs(arrays, col_name_map)

    for out_descr in out_descrs:
        name = out_descr[0]
        dtype = out_descr[1:]
        if masked:
            # Adapted from ma.all_masked() code.  Here the array is filled with
            # zeros instead of empty.  This avoids the bug reported here:
            # https://github.com/numpy/numpy/issues/3276
            out[name] = ma.array(data=np.zeros(n_rows, dtype),
                                 mask=np.ones(n_rows, ma.make_mask_descr(dtype)))
        else:
            out[name] = np.empty(n_rows, dtype=dtype)

    for out_name, in_names in six.iteritems(col_name_map):
        for name, array, arr_len in zip(in_names, arrays, arr_lens):
            if name is not None:
                out[out_name][:arr_len] = array[name]

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, collections.Mapping):
        _col_name_map.update(col_name_map)

    return out
