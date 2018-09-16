"""
High-level operations for numpy structured arrays.

Some code and inspiration taken from numpy.lib.recfunctions.join_by().
Redistribution license restrictions apply.
"""

from itertools import chain
import collections
from collections import OrderedDict, Counter
from collections.abc import Sequence

import numpy as np
import numpy.ma as ma

from . import _np_utils

__all__ = ['TableMergeError']


class TableMergeError(ValueError):
    pass


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
        table_names = [str(ii + 1) for ii in range(len(arrays))]

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
    col_name_count = Counter(col_name_list)
    repeated_names = [name for name, count in col_name_count.items() if count > 1]
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


def _check_for_sequence_of_structured_arrays(arrays):
    err = '`arrays` arg must be a sequence (e.g. list) of structured arrays'
    if not isinstance(arrays, Sequence):
        raise TypeError(err)
    for array in arrays:
        # Must be structured array
        if not isinstance(array, np.ndarray) or array.dtype.names is None:
            raise TypeError(err)
    if len(arrays) == 0:
        raise ValueError('`arrays` arg must include at least one array')


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
