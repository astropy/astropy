"""
Perform a database join of two numpy ndarrays.

Some code taken from numpy.lib.recfunctions.join_by().
Redistribution license restrictions apply.
"""

import collections

import numpy as np
import numpy.ma as ma
import _np_utils

def common_dtype(col0, col1):
    """
    Use numpy to find the common dtype for two structured ndarray columns.
    """
    arr0 = np.empty(1, dtype=col0.dtype)
    arr1 = np.empty(1, dtype=col1.dtype)
    for arr in (arr0, arr1):
        if arr.dtype.kind in ('S', 'U'):
            arr[0] = '0' * arr.itemsize
    arr_common = np.array([arr0[0], arr1[0]])
    return arr_common.dtype.str


def get_merge_descrs(left, right, keys, uniq_col_name='{col_name}_{table_name}',
                     table_names=['1', '2']):
    """
    Find the dtypes descrs resulting from merging the left and right dtypes,
    assuming the names in `keys` are common in the output.

    Return a list of descrs for the output and a dict mapping
    each output column name to the input(s).  This takes the form
    { outname : (left_col_name, right_col_name), ... }.
    For key columns both the left and right input names will be
    present, while for the other non-key columns the value will
    be either (left_col_name, None) or (None, right_col_name).
    """

    out_descrs = []
    col_name_map = collections.defaultdict(lambda: list([None, None]))

    # TODO: should probably refactor this into one routine that determines
    # col_name_map and a second that gets the merged descrs.

    for array, other in ((left, right),
                         (right, left)):
        idx = (0 if array is left else 1)
        table_name = table_names[idx]
        for descr in array.dtype.descr:
            name = descr[0]
            shape = array[name].shape[1:]
            out_descr = [name, descr[1], shape]

            if name in keys:
                out_name = name  # note: in future there may be diff keys for left / right
                if array is right:
                    col_name_map[out_name][1] = name
                    # Skip over keys for the right array, already was done for the left
                    continue
                else:
                    col_name_map[out_name][0] = name
                if shape != other[name].shape[1:]:
                    raise ValueError('Key columns {0!r} have different shape'.format(name))
                out_descr[1] = common_dtype(array[name], other[name])
            elif name in other.dtype.names:
                out_descr[0] = uniq_col_name.format(table_name=table_name, col_name=name)

            col_name_map[out_descr[0]][idx] = name
            out_descrs.append(tuple(out_descr))

    # Convert col_name_map to a regular dict with tuple (immutable) values
    col_name_map = dict((key, tuple(val)) for key, val in col_name_map.items())

    return out_descrs, col_name_map


def join(left, right, keys=None, join_type='inner',
         uniq_col_name='{col_name}_{table_name}',
         table_names=['1', '2'],
         col_name_map=None):
    """
    Perform a join of the left and right numpy structured array on specified keys.

    Parameters
    ----------
    left : structured ndarray
        Left side table in the join
    right : structured ndarray
        Right side table in the join
    keys : str or list of str
        Column(s) used to match rows of left and right tables.  Default
        is to use all columns which are common to both tables.
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

    # This function improves on np.lib.recfunctions.join_by():
    #   key values do not have to be unique (ok with cartesion join)
    #   key columns do not have to be in the same order
    #   key columns can have non-trivial shape

    if join_type not in ('inner', 'outer', 'left', 'right'):
        raise ValueError("The 'join_type' argument should be in 'inner', "
                         "'outer', 'left' or 'right' (got '{0}' instead)".
                         format(join_type))

    # If we have a single key, put it in a tuple
    if keys is None:
        keys = tuple(name for name in left.dtype.names if name in right.dtype.names)
        if len(keys) == 0:
            raise ValueError('No keys in common between left and right tables')
    elif isinstance(keys, basestring):
        keys = (keys,)

    # Check the keys
    for name in keys:
        if name not in left.dtype.names:
            raise ValueError('left does not have key field %s' % name)
        if name not in right.dtype.names:
            raise ValueError('right does not have key field %s' % name)

    # Make sure we work with ravelled arrays
    left = left.ravel()
    right = right.ravel()
    len_left, len_right = len(left), len(right)
    left_names, right_names = left.dtype.names, right.dtype.names

    # Joined array dtype as a list of descr (name, type_str, shape) tuples
    out_descrs, _col_name_map = get_merge_descrs(left, right, keys, uniq_col_name, table_names)
    # If col_name_map supplied as a dict input, then update.
    if isinstance(col_name_map, dict):
        col_name_map.update(_col_name_map)

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

    for out_name, left_right_names in _col_name_map.items():
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
            raise ValueError('Unexpected column names (maybe one is ""?)')
        out[out_name] = array[name].take(array_out)
        if masked:
            if isinstance(array, ma.MaskedArray):
                array_mask = array_mask | array[name].mask.take(array_out)
            out[out_name].mask = array_mask

    return out
