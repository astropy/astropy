import numpy as np
import warnings

from itertools import izip


def group_by(table, keys):
    """
    Get groups for numpy structured array on specified keys.

    Parameters
    ----------
    table : structured array
        Table to group
    keys : str, list of str, `Table`, or Numpy array
        Grouping key specifier

    Returns
    -------
    idxs, idx_sort : numpy arrays
    """
    from .table import Table

    # Pre-convert string to tuple of strings, or Table to the underlying structured array
    if isinstance(keys, basestring):
        keys = (keys,)
    elif isinstance(keys, Table):
        keys = keys._data

    if isinstance(keys, (list, tuple)):
        for name in keys:
            if name not in table.colnames:
                raise ValueError('Table does not have key column {0!r}'.format(name))
            if table.masked and np.any(table[name].mask):
                raise ValueError('Missing values in key column {0!r} are not allowed'.format(name))

        keys = tuple(keys)
        table_keys = table[keys]._data

    elif isinstance(keys, np.ndarray):
        table_keys = keys
        if len(table_keys) != len(table):
            raise ValueError('Input keys array length {0} does not match table length {1}'
                             .format(len(table_keys), len(table)))
        keys = ()

    else:
        raise TypeError('Keys input must be string, list, tuple or numpy array, but got {0}'
                        .format(type(keys)))

    idx_sort = table_keys.argsort()
    table_keys = table_keys[idx_sort]

    # Get all keys
    diffs = np.concatenate(([True], table_keys[1:] != table_keys[:-1], [True]))
    indices = np.flatnonzero(diffs)

    # Note the use of copy=False because a copy is already made with table[idx_sort]
    out = table.__class__(table[idx_sort], copy=False)
    out.groups = TableGroups(out, indices=indices, group_keys=keys)

    return out


class BaseGroups(object):
    """
    A class to represent groups within a table of heterogeneous data.

      - ``group_keys``: list of key column names
      - ``group_indices``: index values corresponding to group boundaries
      - ``aggregate()``: method to create new table by aggregating within groups
    """
    def __init__(self, parent, indices=None, group_keys=None):
        self.parent = parent  # parent Table or Column
        self._indices = (np.array([0, len(self.parent)]) if indices is None else indices)
        self._group_keys = group_keys or ()

    def values(self):
        i0s = self.indices[:-1]
        i1s = self.indices[1:]
        for i0, i1 in izip(i0s, i1s):
            yield self.parent[i0:i1]

    def keys(self):
        i0s = self.indices[:-1]
        i1s = self.indices[1:]
        for i0, i1 in izip(i0s, i1s):
            key_vals = tuple(self.parent[key][i0] for key in self.parent.group_keys)
            if len(key_vals) == 0:
                key_vals = None
            elif len(key_vals) == 1:
                key_vals = key_vals[0]
            yield key_vals

    @property
    def indices(self):
        return self._indices

    @property
    def group_keys(self):
        return self._group_keys

    def __getitem__(self, item):
        if isinstance(item, int):
            i0, i1 = self.indices[item], self.indices[item + 1]
            return self.parent[i0:i1]
        elif isinstance(item, slice):
            raise NotImplementedError()


class TableGroups(BaseGroups):
    def aggregate(self, func):
        """
        Aggregate each group in the Table into a single row by applying the reduction
        function ``func`` to group values in each column.

        Parameters
        ----------
        func : function
            Function that reduces an array of values to a single value

        Returns
        -------
        out : Table
            New table with the aggregated rows.
        """
        idxs0, idxs1 = self.indices[:-1], self.indices[1:]
        out_cols = []

        for col in self.parent.columns.values():
            # For key columns just pick off first in each group since they are identical
            if col.name in self.group_keys:
                vals = col.take(idxs0)
            else:
                try:
                    vals = np.array([func(col[i0: i1]) for i0, i1 in izip(idxs0, idxs1)])
                except Exception:
                    warnings.warn("Cannot aggregate column '{0}'"
                                  .format(col.name))
                    continue

            out_cols.append((col, vals))

        out_cols = [col.__class__(data=vals, name=col.name, description=col.description,
                                  unit=col.unit, format=col.format, meta=col.meta)
                    for col, vals in out_cols]

        return self.parent.__class__(out_cols, meta=self.parent.meta)
