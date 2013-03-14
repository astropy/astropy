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


def get_merge_descrs(left, right, keys, left_fix, right_fix):
    """
    Find the dtypes descrs resulting from merging the left and right dtypes,
    assuming the names in `keys` are common in the output.
    """

    out_descrs = []
    for array, other, name_fix in ((left, right, left_fix),
                                   (right, left, right_fix)):
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
                out_descr[0] = name + name_fix

            out_descrs.append(tuple(out_descr))

    return out_descrs


def join(left, right, keys=None, jointype='inner', left_fix='_l', right_fix='_r'):
    # Check jointype
    if jointype not in ('inner', 'outer', 'left', 'right'):
        raise ValueError("The 'jointype' argument should be in 'inner', "
                         "'outer', 'left' or 'right' (got '{0}' instead)".
                         format(jointype))

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
    out_descrs = get_merge_descrs(left, right, keys, left_fix, right_fix)

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
    int_jointype = {'inner': 0, 'outer': 1, 'left': 2, 'right': 3}[jointype]
    masked, n_out, left_out, left_mask, right_out, right_mask = \
        cyjoin.join_inner(idxs, idx_sort, len_left, int_jointype)

    if masked:
        out = ma.empty(n_out, dtype=out_descrs)
    else:
        out = np.empty(n_out, dtype=out_descrs)

    out_names = out.dtype.names

    for array, array_out, array_mask, name_fix in [
            (left, left_out, left_mask, left_fix),
            (right, right_out, right_mask, right_fix)]:
        names = (name for name in array.dtype.names if name not in keys)
        for name in names:
            out_name = (name if name in out_names else name + name_fix)
            out[out_name] = array[name].take(array_out)
            if masked and name not in keys:
                out[out_name].mask = array_mask

    for name in keys:
        out[name] = np.where(right_mask, left[name].take(left_out), right[name].take(right_out))

    return out
