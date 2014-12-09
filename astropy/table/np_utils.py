"""
High-level operations for numpy structured arrays.

join():  Perform a database join of two numpy ndarrays.
hstack(): Horizontally stack a list of numpy ndarrays.
vstack(): Vertically stack a list of numpy ndarrays.

Some code and inspriration taken from numpy.lib.recfunctions.join_by().
Redistribution license restrictions apply.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six
from ..extern.six.moves import zip as izip
from ..utils.decorators import deprecated

from itertools import chain
import collections

import numpy as np
import numpy.ma as ma

from . import _np_utils
from ..utils import OrderedDict

__all__ = ['join', 'hstack', 'vstack', 'TableMergeError']

DEPRECATION_MESSAGE = ('The %(func)s %(obj_type)s is deprecated and may '
                       'be removed in a future version. '
                       'Contact the Astropy developers if you need '
                       'continued support for this function.')

class TableMergeError(ValueError):
    pass


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
    Find the column names mapping when merging the list of structured ndarrays
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
        for name in array.dtype.names:
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
                if any(name in other.dtype.names for other in others):
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
        in_cols = [arr[name] for arr, name in izip(arrays, in_names) if name is not None]

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
    Use numpy to find the common dtype for a list of structured ndarray columns.

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


@deprecated('1.0', message=DEPRECATION_MESSAGE)
def join(left, right, keys=None, join_type='inner',
         uniq_col_name='{col_name}_{table_name}',
         table_names=['1', '2'],
         col_name_map=None):
    """
    Perform a join of the left and right numpy structured array on specified keys.

    Parameters
    ----------
    left : structured array
        Left side table in the join
    right : structured array
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
    """
    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    if join_type not in ('inner', 'outer', 'left', 'right'):
        raise ValueError("The 'join_type' argument should be in 'inner', "
                         "'outer', 'left' or 'right' (got '{0}' instead)".
                         format(join_type))

    # If we have a single key, put it in a tuple
    if keys is None:
        keys = tuple(name for name in left.dtype.names if name in right.dtype.names)
        if len(keys) == 0:
            raise TableMergeError('No keys in common between left and right tables')
    elif isinstance(keys, six.string_types):
        keys = (keys,)

    # Check the key columns
    for arr, arr_label in ((left, 'Left'), (right, 'Right')):
        for name in keys:
            if name not in arr.dtype.names:
                raise TableMergeError('{0} table does not have key column {1!r}'
                                      .format(arr_label, name))
            if hasattr(arr[name], 'mask') and np.any(arr[name].mask):
                raise TableMergeError('{0} key column {1!r} has missing values'
                                      .format(arr_label, name))

    # Make sure we work with ravelled arrays
    left = left.ravel()
    right = right.ravel()
    len_left, len_right = len(left), len(right)
    left_names, right_names = left.dtype.names, right.dtype.names

    # Joined array dtype as a list of descr (name, type_str, shape) tuples
    col_name_map = get_col_name_map([left, right], keys, uniq_col_name, table_names)
    out_descrs = get_descrs([left, right], col_name_map)

    # Make an array with just the key columns
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
    if any(isinstance(array, ma.MaskedArray) for array in (left, right)):
        masked = True

    if masked:
        out = ma.empty(n_out, dtype=out_descrs)
    else:
        out = np.empty(n_out, dtype=out_descrs)

    # If either input array was zero length then stub a new version
    # with one row.  In this case the corresponding left_out or right_out
    # will contain all zeros with mask set to true.  This allows the
    # take(*_out) method calls to work as expected.
    if len(left) == 0:
        left = left.__class__(1, dtype=left.dtype)
    if len(right) == 0:
        right = right.__class__(1, dtype=right.dtype)

    for out_name, left_right_names in six.iteritems(col_name_map):
        left_name, right_name = left_right_names

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
            if isinstance(array, ma.MaskedArray):
                array_mask = array_mask | array[name].mask.take(array_out)
            out[out_name].mask = array_mask

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, collections.Mapping):
        _col_name_map.update(col_name_map)

    return out


def _check_for_sequence_of_structured_arrays(arrays):
    err = '`arrays` arg must be a sequence (e.g. list) of structured arrays'
    if not isinstance(arrays, collections.Sequence):
        raise TypeError(err)
    for array in arrays:
        # Must be structured array
        if not isinstance(array, np.ndarray) or array.dtype.names is None:
            raise TypeError(err)
    if len(arrays) == 0:
        raise ValueError('`arrays` arg must include at least one array')


@deprecated('1.0', message=DEPRECATION_MESSAGE)
def vstack(arrays, join_type='inner', col_name_map=None):
    """
    Stack structured arrays vertically (by rows)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same column names (though the order can vary).  If
    ``join_type`` is 'inner' then the intersection of common columns will
    be output.  A value of 'outer' means the output will have the union of
    all columns, with array values being masked where no common values are
    available.

    Parameters
    ----------

    arrays : list of structured arrays
        Structured array(s) to stack by rows (vertically)
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'exact'
    col_name_map : empty dict or None
        If passed as a dict then it will be updated in-place with the
        mapping of output to input column names.

    Examples
    --------

    To stack two structured arrays by rows do::

      >>> from astropy.table import np_utils
      >>> t1 = np.array([(1, 2),
      ...                (3, 4)], dtype=[(str('a'), 'i4'), (str('b'), 'i4')])
      >>> t2 = np.array([(5, 6),
      ...                (7, 8)], dtype=[(str('a'), 'i4'), (str('b'), 'i4')])
      >>> np_utils.vstack([t1, t2])
      array([(1, 2),
             (3, 4),
             (5, 6),
             (7, 8)],
            dtype=[('a', '<i4'), ('b', '<i4')])
    """
    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    # Input validation
    if join_type not in ('inner', 'exact', 'outer'):
        raise ValueError("`join_type` arg must be one of 'inner', 'exact' or 'outer'")

    _check_for_sequence_of_structured_arrays(arrays)

    # Trivial case of one input array
    if len(arrays) == 1:
        return arrays[0]

    # Start by assuming an outer match where all names go to output
    names = set(chain(*[arr.dtype.names for arr in arrays]))
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
    masked = any(isinstance(arr, ma.MaskedArray) for arr in arrays)
    for names in six.itervalues(col_name_map):
        if any(x is None for x in names):
            masked = True
            break

    lens = [len(arr) for arr in arrays]
    n_rows = sum(lens)
    out_descrs = get_descrs(arrays, col_name_map)
    if masked:
        # Make a masked array with all values initially masked.  Note
        # that setting an array value automatically unmasks it.
        # See comment in hstack for heritage of this code.
        out = ma.masked_array(np.zeros(n_rows, out_descrs),
                              mask=np.ones(n_rows, ma.make_mask_descr(out_descrs)))
    else:
        out = np.empty(n_rows, dtype=out_descrs)

    for out_name, in_names in six.iteritems(col_name_map):
        idx0 = 0
        for name, array in izip(in_names, arrays):
            idx1 = idx0 + len(array)
            if name in array.dtype.names:
                out[out_name][idx0:idx1] = array[name]
            idx0 = idx1

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, collections.Mapping):
        _col_name_map.update(col_name_map)

    return out


@deprecated('1.0', message=DEPRECATION_MESSAGE)
def hstack(arrays, join_type='exact', uniq_col_name='{col_name}_{table_name}',
           table_names=None, col_name_map=None):
    """
    Stack structured arrays by horizontally (by columns)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be output.  A value of 'outer' means
    the output will have the union of all rows, with array values being
    masked where no common values are available.

    Parameters
    ----------

    arrays : List of structured array objects
        Structured arrays to stack by columns (horizontally)
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'exact'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2', ..].

    Examples
    --------

    To stack two arrays horizontally (by columns) do::

      >>> from astropy.table import np_utils
      >>> t1 = np.array([(1, 2),
      ...                (3, 4)], dtype=[(str('a'), 'i4'), (str('b'), 'i4')])
      >>> t2 = np.array([(5, 6),
      ...                (7, 8)], dtype=[(str('c'), 'i4'), (str('d'), 'i4')])
      >>> np_utils.hstack([t1, t2])
      array([(1, 2, 5, 6),
             (3, 4, 7, 8)],
            dtype=[('a', '<i4'), ('b', '<i4'), ('c', '<i4'), ('d', '<i4')])
    """
    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    # Input validation
    if join_type not in ('inner', 'exact', 'outer'):
        raise ValueError("join_type arg must be either 'inner', 'exact' or 'outer'")
    _check_for_sequence_of_structured_arrays(arrays)

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

    # For an inner join, keep only columns where all input arrays have that column
    if join_type == 'inner':
        min_arr_len = min(arr_lens)
        arrays = [arr[:min_arr_len] for arr in arrays]
        arr_lens = [min_arr_len for arr in arrays]

    # If there are any output rows where one or more input arrays are missing
    # then the output must be masked.  If any input arrays are masked then
    # output is masked.
    masked = (any(isinstance(arr, ma.MaskedArray) for arr in arrays) or
              len(set(arr_lens)) > 1)

    n_rows = max(arr_lens)
    out_descrs = get_descrs(arrays, col_name_map)
    if masked:
        # Adapted from ma.all_masked() code.  Here the array is filled with
        # zeros instead of empty.  This avoids the bug reported here:
        # https://github.com/numpy/numpy/issues/3276
        out = ma.masked_array(np.zeros(n_rows, out_descrs),
                              mask=np.ones(n_rows, ma.make_mask_descr(out_descrs)))
    else:
        out = np.empty(n_rows, dtype=out_descrs)

    for out_name, in_names in six.iteritems(col_name_map):
        for name, array, arr_len in izip(in_names, arrays, arr_lens):
            if name is not None:
                out[out_name][:arr_len] = array[name]

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, collections.Mapping):
        _col_name_map.update(col_name_map)

    return out


@deprecated('1.0', message=DEPRECATION_MESSAGE)
def get_groups(table, keys):
    """
    Get groups for numpy structured array on specified keys.

    Parameters
    ----------
    table : structured array
        Table to group
    keys : str or list of str
        Name(s) of column(s) used to match rows of table.
    """
    if isinstance(keys, six.string_types):
        keys = (keys,)

    # Check the key columns
    for name in keys:
        if name not in table.dtype.names:
            raise TableMergeError('Table does not have key column {1!r}'
                                  .format(name))
        if hasattr(table[name], 'mask') and np.any(table[name].mask):
            raise TableMergeError('{0} key column {1!r} has missing values'
                                  .format(name))

    # Make sure we work with ravelled arrays
    table = table.ravel()
    len_table = len(table)

    # oined array dtype as a list of descr (name, type_str, shape) tuples
    col_name_map = get_col_name_map([table], keys)
    out_descrs = get_descrs([table], col_name_map)

    # Make an array with just the key columns
    out_keys_dtype = [descr for descr in out_descrs if descr[0] in keys]
    out_keys = np.empty(len_table, dtype=out_keys_dtype)
    for key in keys:
        out_keys[key] = table[key]
    idx_sort = out_keys.argsort(order=keys)
    out_keys = out_keys[idx_sort]

    # Get all keys
    diffs = np.concatenate(([True], out_keys[1:] != out_keys[:-1], [True]))
    idxs = np.flatnonzero(diffs)

    return idxs, out_keys


def fix_column_name(val):
    """
    Fixes column names so that they are compatible with Numpy on
    Python 2.  Raises a ValueError exception if the column name
    contains Unicode characters, which can not reasonably be used as a
    column name.
    """
    if val is not None:
        try:
            val = str(val)
        except UnicodeEncodeError:
            if not six.PY3:
                raise ValueError(
                    "Column names must not contain Unicode characters "
                    "on Python 2")
            raise

    return val


def recarray_fromrecords(rec_list):
    """
    Partial replacement for `~numpy.core.records.fromrecords` which includes
    a workaround for the bug with unicode arrays described at:
    https://github.com/astropy/astropy/issues/3052

    This should not serve as a full replacement for the original function;
    this only does enough to fulfill the needs of the table module.
    """

    # Note: This is just copying what Numpy does for converting arbitrary rows
    # to column arrays in the recarray module; it could be there is a better
    # way
    nfields = len(rec_list[0])
    obj = np.array(rec_list, dtype=object)
    array_list = [np.array(obj[..., i].tolist()) for i in range(nfields)]
    formats = []
    for obj in array_list:
        formats.append(obj.dtype.str)
    formats = ','.join(formats)
    return np.rec.fromarrays(array_list, formats=formats)
