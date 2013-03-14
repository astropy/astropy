"""
Perform a database join of two numpy ndarrays.

Some code taken from numpy.lib.recfunctions.join_by().
Redistribution license restrictions apply.
"""

import numpy as np
import numpy.ma as ma
import cyjoin


# np.lib.recfunctions.join_by() limitations:
#   key columns must be in the same order
#   key columns cannot have non-trivial shape

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


def get_merge_descrs(left, right, keys, uniq_col_name, table_names):
    """
    Find the dtypes descrs resulting from merging the left and right dtypes,
    assuming the names in `keys` are common in the output.
    """

    out_descrs = []
    for array, other, table_name in ((left, right, table_names[0]),
                                   (right, left, table_names[1])):
        for descr in array.dtype.descr:
            name = descr[0]
            shape = array[name].shape[1:]
            out_descr = [name, descr[1], shape]

            if name in keys:
                if array is right:
                    # Skip over keys for the right array, already done with left
                    continue
                if shape != other[name].shape[1:]:
                    raise ValueError('Key columns {0!r} have different shape'.format(name))
                out_descr[1] = common_dtype(array[name], other[name])
            elif name in other.dtype.names:
                out_descr[0] = uniq_col_name.format(table_name=table_name, col_name=name)

            out_descrs.append(tuple(out_descr))

    return out_descrs


def join(left, right, keys=None, join_type='inner',
         uniq_col_name='{col_name}_{table_name}',
         table_names=['1', '2']):
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

    """

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
    out_descrs = get_merge_descrs(left, right, keys, uniq_col_name, table_names)

    # Make an array with just the key columns
    aux_dtype = [descr for descr in out_descrs if descr[0] in keys]
    aux = np.empty(len_left + len_right, dtype=aux_dtype)
    for key in keys:
        aux[key][:len_left] = left[key]
        aux[key][len_left:] = right[key]
    idx_sort = aux.argsort(order=keys)
    aux = aux[idx_sort]

    # Get all keys
    diffs = np.concatenate(([True], aux[1:] != aux[:-1], [True]))
    idxs = np.flatnonzero(diffs)

    # Main inner loop in Cython to compute the cartesion product
    # indices for the given join type
    int_join_type = {'inner': 0, 'outer': 1, 'left': 2, 'right': 3}[join_type]
    masked, n_out, left_out, left_mask, right_out, right_mask = \
        cyjoin.join_inner(idxs, idx_sort, len_left, int_join_type)

    if masked:
        out = ma.empty(n_out, dtype=out_descrs)
    else:
        out = np.empty(n_out, dtype=out_descrs)

    out_names = out.dtype.names

    for array, array_out, array_mask, table_name in [
            (left, left_out, left_mask, table_names[0]),
            (right, right_out, right_mask, table_names[1])]:
        names = (name for name in array.dtype.names if name not in keys)
        for name in names:
            # TODO: return the mapping of input to output col names when generating
            # dtype to avoid repeating the logic here.
            if name not in out_names:
                out_name = uniq_col_name.format(table_name=table_name, col_name=name)
            else:
                out_name = name
            out[out_name] = array[name].take(array_out)
            if masked and name not in keys:
                out[out_name].mask = array_mask

    for name in keys:
        out[name] = np.where(right_mask, left[name].take(left_out), right[name].take(right_out))

    return out
