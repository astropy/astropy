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
    out._groups = TableGroups(out, indices=indices, group_keys=keys)

    return out


class BaseGroups(object):
    """
    A class to represent groups within a table of heterogeneous data.

      - ``group_keys``: list of key column names
      - ``group_indices``: index values corresponding to group boundaries
      - ``aggregate()``: method to create new table by aggregating within groups
    """
    def values(self):
        i0s, i1s = self.indices[:-1], self.indices[1:]
        for i0, i1 in izip(i0s, i1s):
            yield self.parent[i0:i1]

    def keys(self):
        i0s, i1s = self.indices[:-1], self.indices[1:]
        for i0, i1 in izip(i0s, i1s):
            key_vals = tuple(self.parent[key][i0] for key in self.parent.group_keys)
            if len(key_vals) == 0:
                key_vals = None
            elif len(key_vals) == 1:
                key_vals = key_vals[0]
            yield key_vals

    @property
    def group_keys(self):
        return self._group_keys

    def __getitem__(self, item):
        parent = self.parent_column if isinstance(self, ColumnGroups) else self.parent_table

        if isinstance(item, int):
            i0, i1 = self.indices[item], self.indices[item + 1]
            return parent[i0:i1]
        elif isinstance(item, slice):
            raise NotImplementedError()


class ColumnGroups(BaseGroups):
    def __init__(self, parent_column, indices=None, group_keys=None):
        self.parent_column = parent_column  # parent Column
        self.parent_table = parent_column.parent_table
        self._indices = indices
        self._group_keys = group_keys or ()

    @property
    def indices(self):
        # If the parent column is in a table then use group indices from table
        if self.parent_table:
            return self.parent_table.groups.indices
        else:
            if self._indices is None:
                return np.array([0, len(self.parent_column)])
            else:
                return self._indices

    def aggregate(self, func):
        i0s, i1s = self.indices[:-1], self.indices[1:]
        par_col = self.parent_column
        try:
            vals = np.array([func(par_col[i0: i1]) for i0, i1 in izip(i0s, i1s)])
        except Exception:
            raise TypeError("Cannot aggregate column '{0}'"
                            .format(par_col.name))

        out = par_col.__class__(data=vals, name=par_col.name, description=par_col.description,
                                unit=par_col.unit, format=par_col.format, meta=par_col.meta)
        return out

    def __repr__(self):
        return '<{0} indices={1}>'.format(self.__class__.__name__, self.indices)


class TableGroups(BaseGroups):
    def __init__(self, parent_table, indices=None, group_keys=None):
        self.parent_table = parent_table  # parent Table
        self._indices = indices
        self._group_keys = group_keys or ()

    @property
    def indices(self):
        if self._indices is None:
            return np.array([0, len(self.parent_table)])
        else:
            return self._indices

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
        i0s, i1s = self.indices[:-1], self.indices[1:]
        out_cols = []
        parent_table = self.parent_table

        for col in parent_table.columns.values():
            # For key columns just pick off first in each group since they are identical
            if col.name in self.group_keys:
                # Should just be new_col = col.take(i0s), but there is a bug in
                # MaskedColumn finalize:
                # >>> c = MaskedColumn(data=[1,2], name='a', description='a')
                # >>> print c[1:2].description
                # None
                new_col = col.__class__(data=col.take(i0s), name=col.name,
                                        description=col.description,
                                        unit=col.unit, format=col.format, meta=col.meta)
            else:
                try:
                    new_col = col.groups.aggregate(func)
                except TypeError as err:
                    warnings.warn(str(err))
                    continue

            out_cols.append(new_col)

        return parent_table.__class__(out_cols, meta=parent_table.meta)

    def __repr__(self):
        return '<{0} group_keys={1} indices={2}>'.format(self.__class__.__name__, self.group_keys,
                                                         self.indices)
