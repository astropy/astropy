# Licensed under a 3-clause BSD style license - see LICENSE.rst

import platform
import warnings

import numpy as np

from astropy.utils.exceptions import AstropyUserWarning

from .index import get_index_by_names

__all__ = ["TableGroups", "ColumnGroups"]


def table_group_by(table, keys):
    # index copies are unnecessary and slow down _table_group_by
    with table.index_mode("discard_on_copy"):
        return _table_group_by(table, keys)


def _table_group_by(table, keys):
    """
    Get groups for ``table`` on specified ``keys``.

    Parameters
    ----------
    table : `Table`
        Table to group
    keys : str, list of str, `Table`, or Numpy array
        Grouping key specifier

    Returns
    -------
    grouped_table : Table object with groups attr set accordingly
    """
    from .serialize import represent_mixins_as_columns
    from .table import Table

    # Pre-convert string to tuple of strings, or Table to the underlying structured array
    if isinstance(keys, str):
        keys = (keys,)

    if isinstance(keys, (list, tuple)):
        for name in keys:
            if name not in table.colnames:
                raise ValueError(f"Table does not have key column {name!r}")
            if table.masked and np.any(table[name].mask):
                raise ValueError(
                    f"Missing values in key column {name!r} are not allowed"
                )

        # Make a column slice of the table without copying
        table_keys = table.__class__([table[key] for key in keys], copy=False)

        # If available get a pre-existing index for these columns
        table_index = get_index_by_names(table, keys)
        grouped_by_table_cols = True

    elif isinstance(keys, (np.ndarray, Table)):
        table_keys = keys
        if len(table_keys) != len(table):
            raise ValueError(
                "Input keys array length {} does not match table length {}".format(
                    len(table_keys), len(table)
                )
            )
        table_index = None
        grouped_by_table_cols = False

    else:
        raise TypeError(
            "Keys input must be string, list, tuple, Table or numpy array, but got {}".format(
                type(keys)
            )
        )

    # If there is not already an available index and table_keys is a Table then ensure
    # that all cols (including mixins) are in a form that can sorted with the code below.
    if not table_index and isinstance(table_keys, Table):
        table_keys = represent_mixins_as_columns(table_keys)

    # Get the argsort index `idx_sort`, accounting for particulars
    try:
        # take advantage of index internal sort if possible
        if table_index is not None:
            idx_sort = table_index.sorted_data()
        else:
            idx_sort = table_keys.argsort(kind="mergesort")
        stable_sort = True
    except TypeError:
        # Some versions (likely 1.6 and earlier) of numpy don't support
        # 'mergesort' for all data types.  MacOSX (Darwin) doesn't have a stable
        # sort by default, nor does Windows, while Linux does (or appears to).
        idx_sort = table_keys.argsort()
        stable_sort = platform.system() not in ("Darwin", "Windows")

    # Finally do the actual sort of table_keys values
    table_keys = table_keys[idx_sort]

    # Get all keys
    diffs = np.concatenate(([True], table_keys[1:] != table_keys[:-1], [True]))
    indices = np.flatnonzero(diffs)

    # If the sort is not stable (preserves original table order) then sort idx_sort in
    # place within each group.
    if not stable_sort:
        for i0, i1 in zip(indices[:-1], indices[1:]):
            idx_sort[i0:i1].sort()

    # Make a new table and set the _groups to the appropriate TableGroups object.
    # Take the subset of the original keys at the indices values (group boundaries).
    out = table.__class__(table[idx_sort])
    out_keys = table_keys[indices[:-1]]
    if isinstance(out_keys, Table):
        out_keys.meta["grouped_by_table_cols"] = grouped_by_table_cols
    out._groups = TableGroups(out, indices=indices, keys=out_keys)

    return out


def column_group_by(column, keys):
    """
    Get groups for ``column`` on specified ``keys``.

    Parameters
    ----------
    column : Column object
        Column to group
    keys : Table or Numpy array of same length as col
        Grouping key specifier

    Returns
    -------
    grouped_column : Column object with groups attr set accordingly
    """
    from .serialize import represent_mixins_as_columns
    from .table import Table

    if isinstance(keys, Table):
        keys = represent_mixins_as_columns(keys)
        keys = keys.as_array()

    if not isinstance(keys, np.ndarray):
        raise TypeError(f"Keys input must be numpy array, but got {type(keys)}")

    if len(keys) != len(column):
        raise ValueError(
            "Input keys array length {} does not match column length {}".format(
                len(keys), len(column)
            )
        )

    idx_sort = keys.argsort()
    keys = keys[idx_sort]

    # Get all keys
    diffs = np.concatenate(([True], keys[1:] != keys[:-1], [True]))
    indices = np.flatnonzero(diffs)

    # Make a new column and set the _groups to the appropriate ColumnGroups object.
    # Take the subset of the original keys at the indices values (group boundaries).
    out = column.__class__(column[idx_sort])
    out._groups = ColumnGroups(out, indices=indices, keys=keys[indices[:-1]])

    return out


class BaseGroups:
    """
    A class to represent groups within a table of heterogeneous data.

      - ``keys``: key values corresponding to each group
      - ``indices``: index values in parent table or column corresponding to group boundaries
      - ``aggregate()``: method to create new table by aggregating within groups
    """

    @property
    def parent(self):
        return (
            self.parent_column if isinstance(self, ColumnGroups) else self.parent_table
        )

    def __iter__(self):
        self._iter_index = 0
        return self

    def next(self):
        ii = self._iter_index
        if ii < len(self.indices) - 1:
            i0, i1 = self.indices[ii], self.indices[ii + 1]
            self._iter_index += 1
            return self.parent[i0:i1]
        else:
            raise StopIteration

    __next__ = next

    def __getitem__(self, item):
        parent = self.parent

        if isinstance(item, (int, np.integer)):
            i0, i1 = self.indices[item], self.indices[item + 1]
            out = parent[i0:i1]
            out.groups._keys = parent.groups.keys[item]
        else:
            indices0, indices1 = self.indices[:-1], self.indices[1:]
            try:
                i0s, i1s = indices0[item], indices1[item]
            except Exception as err:
                raise TypeError(
                    "Index item for groups attribute must be a slice, "
                    "numpy mask or int array"
                ) from err
            mask = np.zeros(len(parent), dtype=bool)
            # Is there a way to vectorize this in numpy?
            for i0, i1 in zip(i0s, i1s):
                mask[i0:i1] = True
            out = parent[mask]
            out.groups._keys = parent.groups.keys[item]
            out.groups._indices = np.concatenate([[0], np.cumsum(i1s - i0s)])

        return out

    def __repr__(self):
        return f"<{self.__class__.__name__} indices={self.indices}>"

    def __len__(self):
        return len(self.indices) - 1


class ColumnGroups(BaseGroups):
    def __init__(self, parent_column, indices=None, keys=None):
        self.parent_column = parent_column  # parent Column
        self.parent_table = parent_column.info.parent_table
        self._indices = indices
        self._keys = keys

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

    @property
    def keys(self):
        # If the parent column is in a table then use group indices from table
        if self.parent_table:
            return self.parent_table.groups.keys
        else:
            return self._keys

    def aggregate(self, func):
        from .column import MaskedColumn

        i0s, i1s = self.indices[:-1], self.indices[1:]
        par_col = self.parent_column
        masked = isinstance(par_col, MaskedColumn)
        reduceat = hasattr(func, "reduceat")
        sum_case = func is np.sum
        mean_case = func is np.mean
        try:
            if not masked and (reduceat or sum_case or mean_case):
                if mean_case:
                    vals = np.add.reduceat(par_col, i0s) / np.diff(self.indices)
                else:
                    if sum_case:
                        func = np.add
                    vals = func.reduceat(par_col, i0s)
            else:
                vals = np.array([func(par_col[i0:i1]) for i0, i1 in zip(i0s, i1s)])
            out = par_col.__class__(vals)
        except Exception as err:
            raise TypeError(
                "Cannot aggregate column '{}' with type '{}': {}".format(
                    par_col.info.name, par_col.info.dtype, err
                )
            ) from err

        out_info = out.info
        for attr in ("name", "unit", "format", "description", "meta"):
            try:
                setattr(out_info, attr, getattr(par_col.info, attr))
            except AttributeError:
                pass

        return out

    def filter(self, func):
        """
        Filter groups in the Column based on evaluating function ``func`` on each
        group sub-table.

        The function which is passed to this method must accept one argument:

        - ``column`` : `Column` object

        It must then return either `True` or `False`.  As an example, the following
        will select all column groups with only positive values::

          def all_positive(column):
              if np.any(column < 0):
                  return False
              return True

        Parameters
        ----------
        func : function
            Filter function

        Returns
        -------
        out : Column
            New column with the aggregated rows.
        """
        mask = np.empty(len(self), dtype=bool)
        for i, group_column in enumerate(self):
            mask[i] = func(group_column)

        return self[mask]


class TableGroups(BaseGroups):
    def __init__(self, parent_table, indices=None, keys=None):
        self.parent_table = parent_table  # parent Table
        self._indices = indices
        self._keys = keys

    @property
    def key_colnames(self):
        """
        Return the names of columns in the parent table that were used for grouping.
        """
        # If the table was grouped by key columns *in* the table then treat those columns
        # differently in aggregation.  In this case keys will be a Table with
        # keys.meta['grouped_by_table_cols'] == True.  Keys might not be a Table so we
        # need to handle this.
        grouped_by_table_cols = getattr(self.keys, "meta", {}).get(
            "grouped_by_table_cols", False
        )
        return self.keys.colnames if grouped_by_table_cols else ()

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
        i0s = self.indices[:-1]
        out_cols = []
        parent_table = self.parent_table

        for col in parent_table.columns.values():
            # For key columns just pick off first in each group since they are identical
            if col.info.name in self.key_colnames:
                new_col = col.take(i0s)
            else:
                try:
                    new_col = col.info.groups.aggregate(func)
                except TypeError as err:
                    warnings.warn(str(err), AstropyUserWarning)
                    continue

            out_cols.append(new_col)

        return parent_table.__class__(out_cols, meta=parent_table.meta)

    def filter(self, func):
        """
        Filter groups in the Table based on evaluating function ``func`` on each
        group sub-table.

        The function which is passed to this method must accept two arguments:

        - ``table`` : `Table` object
        - ``key_colnames`` : tuple of column names in ``table`` used as keys for grouping

        It must then return either `True` or `False`.  As an example, the following
        will select all table groups with only positive values in the non-key columns::

          def all_positive(table, key_colnames):
              colnames = [name for name in table.colnames if name not in key_colnames]
              for colname in colnames:
                  if np.any(table[colname] < 0):
                      return False
              return True

        Parameters
        ----------
        func : function
            Filter function

        Returns
        -------
        out : Table
            New table with the aggregated rows.
        """
        mask = np.empty(len(self), dtype=bool)
        key_colnames = self.key_colnames
        for i, group_table in enumerate(self):
            mask[i] = func(group_table, key_colnames)

        return self[mask]

    @property
    def keys(self):
        return self._keys
