# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The Index class can use several implementations as its
engine. Any implementation should implement the following:

__init__(data, row_index) : initialize index based on key/row list pairs
add(key, row) -> None : add (key, row) to existing data
remove(key, data=None) -> boolean : remove data from self[key], or all of
                                    self[key] if data is None
shift_left(row) -> None : decrement row numbers after row
shift_right(row) -> None : increase row numbers >= row
find(key) -> list : list of rows corresponding to key
range(lower, upper, bounds) -> list : rows in self[k] where k is between
                               lower and upper (<= or < based on bounds)
sort() -> None : make row order align with key order
sorted_data() -> list of rows in sorted order (by key)
replace_rows(row_map) -> None : replace row numbers based on slice
items() -> list of tuples of the form (key, data)

Notes
-----
    When a Table is initialized from another Table, indices are
    (deep) copied and their columns are set to the columns of the new Table.

    Column creation:
    Column(c) -> deep copy of indices
    c[[1, 2]] -> deep copy and reordering of indices
    c[1:2] -> reference
    array.view(Column) -> no indices
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from copy import deepcopy
import numpy as np

from ..extern import six
from .bst import BST, FastBST, FastRBT, MinValue, MaxValue
from .sorted_array import SortedArray
from ..time import Time


class QueryError(ValueError):
    '''
    Indicates that a given index cannot handle the supplied query.
    '''
    pass


class Index(object):
    '''
    The Index class makes it possible to maintain indices
    on columns of a Table, so that column values can be queried
    quickly and efficiently. Column values are stored in lexicographic
    sorted order, which allows for binary searching in O(log n).

    Parameters
    ----------
    columns : list or None
        List of columns on which to create an index. If None,
        create an empty index for purposes of deep copying.
    engine : type, instance, or None
        Indexing engine class to use (from among SortedArray, BST,
        FastBST, and FastRBT) or actual engine instance.
        If the supplied argument is None (by default), use SortedArray.
    unique : bool (defaults to False)
        Whether the values of the index must be unique
    '''
    def __new__(cls, *args, **kwargs):
        self = super(Index, cls).__new__(cls)
        self.__init__(*args, **kwargs)
        return SlicedIndex(self, slice(0, 0, None), original=True)

    def __init__(self, columns, engine=None, unique=False):
        from .table import Table, Column

        if engine is not None and not isinstance(engine, type):
            # create from data
            self.engine = engine.__class__
            self.data = engine
            self.columns = columns
            return

        # by default, use SortedArray
        self.engine = engine or SortedArray

        if columns is None: # this creates a special exception for deep copying
            columns = []
            data = []
            row_index = []
        elif len(columns) == 0:
            raise ValueError("Cannot create index without at least one column")
        elif len(columns) == 1:
            col = columns[0]
            row_index = Column(col.argsort())
            data = Table([col[row_index]])
        else:
            num_rows = len(columns[0])

            # replace Time columns with approximate form and remainder
            new_columns = []
            for col in columns:
                if isinstance(col, Time):
                    new_columns.append(col.jd)
                    remainder = col - col.__class__(col.jd, format='jd')
                    new_columns.append(remainder.jd)
                else:
                    new_columns.append(col)

            # sort the table lexicographically and keep row numbers
            table = Table(columns + [np.arange(num_rows)], copy_indices=False)
            sort_columns = new_columns[::-1]
            try:
                lines = table[np.lexsort(sort_columns)]
            except TypeError: # arbitrary mixins might not work with lexsort
                lines = table[table.argsort()]
            data = lines[lines.colnames[:-1]]
            row_index = lines[lines.colnames[-1]]

        self.data = self.engine(data, row_index, unique=unique)
        self.columns = columns

    def __len__(self):
        '''
        Number of rows in index.
        '''
        return len(self.columns[0])

    def replace_col(self, prev_col, new_col):
        '''
        Replace an indexed column with an updated reference.

        Parameters
        ----------
        prev_col : Column
            Column reference to replace
        new_col : Column
            New column reference
        '''
        self.columns[self.col_position(prev_col.info.name)] = new_col

    def reload(self):
        '''
        Recreate the index based on data in self.columns.
        '''
        self.__init__(self.columns, engine=self.engine)

    def col_position(self, col_name):
        '''
        Return the position of col_name in self.columns.

        Parameters
        ----------
        col_name : str
            Name of column to look up
        '''
        for i, c in enumerate(self.columns):
            if c.info.name == col_name:
                return i
        raise ValueError("Column does not belong to index: {0}".format(col_name))

    def insert_row(self, pos, vals, columns):
        '''
        Insert a new row from the given values.

        Parameters
        ----------
        pos : int
            Position at which to insert row
        vals : list or tuple
            List of values to insert into a new row
        columns : list
            Table column references
        '''
        key = [None] * len(self.columns)
        for i, col in enumerate(columns):
            try:
                key[i] = vals[self.col_position(col.info.name)]
            except ValueError: # not a member of index
                continue
        num_rows = len(self.columns[0])
        if pos < num_rows:
            # shift all rows >= pos to the right
            self.data.shift_right(pos)
        self.data.add(tuple(key), pos)

    def get_row_specifier(self, row_specifier):
        '''
        Return an interable corresponding to the
        input row specifier.

        Parameters
        ----------
        row_specifier : int, list, ndarray, or slice
        '''
        if isinstance(row_specifier, int):
            # single row
            return (row_specifier,)
        elif isinstance(row_specifier, (list, np.ndarray)):
            return row_specifier
        elif isinstance(row_specifier, slice):
            col_len = len(self.columns[0])
            return range(*row_specifier.indices(col_len))
        raise ValueError("Expected int, array of ints, or slice but "
                         "got {0} in remove_rows".format(row_specifier))

    def remove_rows(self, row_specifier):
        '''
        Remove the given rows from the index.

        Parameters
        ----------
        row_specifier : int, list, ndarray, or slice
            Indicates which row(s) to remove
        '''
        rows = []

        # To maintain the correct row order, we loop twice,
        # deleting rows first and then reordering the remaining rows
        for row in self.get_row_specifier(row_specifier):
            self.remove_row(row, reorder=False)
            rows.append(row)
        # second pass - row order is reversed to maintain
        # correct row numbers
        for row in reversed(sorted(rows)):
            self.data.shift_left(row)

    def remove_row(self, row, reorder=True):
        '''
        Remove the given row from the index.

        Parameters
        ----------
        row : int
            Position of row to remove
        reorder : bool
            Whether to reorder indices after removal
        '''
        # for removal, form a key consisting of column values in this row
        if not self.data.remove(tuple([col[row] for col in self.columns]), row):
            raise ValueError("Could not remove row {0} from index".format(row))
        # decrement the row number of all later rows
        if reorder:
            self.data.shift_left(row)

    def find(self, key):
        '''
        Return the row values corresponding to key, in sorted order.

        Parameters
        ----------
        key : tuple
            Values to search for in each column
        '''
        return self.data.find(key)

    def same_prefix(self, key):
        '''
        Return rows whose keys contain the supplied key as a prefix.

        Parameters
        ----------
        key : tuple
            Prefix for which to search
        '''
        return self.same_prefix_range(key, key, (True, True))

    def same_prefix_range(self, lower, upper, bounds=(True, True)):
        '''
        Return rows whose keys have a prefix in the given range.

        Parameters
        ----------
        lower : tuple
            Lower prefix bound
        upper : tuple
            Upper prefix bound
        bounds : tuple (x, y) of bools
            Indicates whether the search should be inclusive or
            exclusive with respect to the endpoints. The first
            argument x corresponds to an inclusive lower bound,
            and the second argument y to an inclusive upper bound.
        '''
        n = len(lower)
        ncols = len(self.columns)
        a = MinValue() if bounds[0] else MaxValue()
        b = MaxValue() if bounds[1] else MinValue()
        # [x, y] search corresponds to [(x, min), (y, max)]
        # (x, y) search corresponds to ((x, max), (x, min))
        lower = lower + tuple((ncols - n) * [a])
        upper = upper + tuple((ncols - n) * [b])
        return self.data.range(lower, upper, bounds)

    def range(self, lower, upper, bounds=(True, True)):
        '''
        Return rows within the given range.

        Parameters
        ----------
        lower : tuple
            Lower prefix bound
        upper : tuple
            Upper prefix bound
        bounds : tuple (x, y) of bools
            Indicates whether the search should be inclusive or
            exclusive with respect to the endpoints. The first
            argument x corresponds to an inclusive lower bound,
            and the second argument y to an inclusive upper bound.
        '''
        return self.data.range(lower, upper, bounds)

    def replace(self, row, col_name, val):
        '''
        Replace the value of a column at a given position.

        Parameters
        ----------
        row : int
            Row number to modify
        col_name : str
            Name of the Column to modify
        val : col.info.dtype
            Value to insert at specified row of col
        '''
        self.remove_row(row, reorder=False)
        key = [c[row] for c in self.columns]
        key[self.col_position(col_name)] = val
        self.data.add(tuple(key), row)

    def replace_rows(self, col_slice):
        '''
        Modify rows in this index to agree with the specified
        slice. For example, given an index
        {'5': 1, '2': 0, '3': 2} on a column ['2', '5', '3'],
        an input col_slice of [2, 0] will result in the relabeling
        {'3': 0, '2': 1} on the sliced column ['3', '2'].

        Parameters
        ----------
        col_slice : list
            Indices to slice
        '''
        row_map = dict((row, i) for i, row in enumerate(col_slice))
        self.data.replace_rows(row_map)

    def sort(self):
        '''
        Make row numbers follow the same sort order as the keys
        of the index.
        '''
        self.data.sort()

    def sorted_data(self):
        '''
        Returns a list of rows in sorted order based on keys;
        essentially acts as an argsort() on columns.
        '''
        return self.data.sorted_data()

    def __getitem__(self, item):
        '''
        Returns a sliced version of this index.

        Parameters
        ----------
        item : slice
            Input slice

        Returns
        -------
        SlicedIndex
            A sliced reference to this index.
        '''
        return SlicedIndex(self, item)

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return str(self)

    def __deepcopy__(self, memo):
        '''
        Return a deep copy of this index.

        Notes
        -----
        The default deep copy must be overridden to perform
        a shallow copy of the index columns, avoiding infinite recursion.

        Parameters
        ----------
        memo : dict
        '''
        num_cols = self.data.num_cols if self.engine == SortedArray else None
        # create an actual Index, not a SlicedIndex
        index = super(Index, Index).__new__(Index)
        index.__init__(None, engine=self.engine)
        index.data = deepcopy(self.data, memo)
        index.columns = self.columns[:] # new list, same columns
        memo[id(self)] = index
        return index


class SlicedIndex(object):
    '''
    This class provides a wrapper around an actual Index object
    to make index slicing function correctly. Since numpy expects
    array slices to provide an actual data view, a SlicedIndex should
    retrieve data directly from the original index and then adapt
    it to the sliced coordinate system as appropriate.

    Parameters
    ----------
    index : Index
        The original Index reference
    index_slice : slice
        The slice to which this SlicedIndex corresponds
    original : bool
        Whether this SlicedIndex represents the original index itself.
        For the most part this is similar to index[:] but certain
        copying operations are avoided, and the slice retains the
        length of the actual index despite modification.
    '''
    def __init__(self, index, index_slice, original=False):
        self.index = index
        self.original = original
        self._frozen = False

        if isinstance(index_slice, tuple):
            self.start, self._stop, self.step = index_slice
        else: # index_slice is an actual slice
            num_rows = len(index.columns[0])
            self.start, self._stop, self.step = index_slice.indices(num_rows)

    @property
    def length(self):
        return 1 + (self.stop - self.start - 1) // self.step

    @property
    def stop(self):
        '''
        The stopping position of the slice, or the end of the
        index if this is an original slice.
        '''
        return len(self.index) if self.original else self._stop

    def __getitem__(self, item):
        '''
        Returns another slice of this Index slice.

        Parameters
        ----------
        item : slice
            Index slice
        '''
        if self.length <= 0:
            # empty slice
            return SlicedIndex(self.index, slice(1, 0))
        start, stop, step = item.indices(self.length)
        new_start = self.orig_coords(start)
        new_stop = self.orig_coords(stop)
        new_step = self.step * step
        return SlicedIndex(self.index, (new_start, new_stop, new_step))

    def sliced_coords(self, rows):
        '''
        Convert the input rows to the sliced coordinate system.

        Parameters
        ----------
        rows : list
            Rows in the original coordinate system

        Returns
        -------
        sliced_rows : list
            Rows in the sliced coordinate system
        '''
        if self.original:
            return rows
        else:
            rows = np.array(rows)
            row0 = rows - self.start
            if self.step != 1:
                correct_mod = np.mod(row0, self.step) == 0
                row0 = row0[correct_mod]
            if self.step > 0:
                ok = (row0 >= 0) & (row0 < self.stop - self.start)
            else:
                ok = (row0 <= 0) & (row0 > self.stop - self.start)
            return row0[ok] // self.step

    def orig_coords(self, row):
        '''
        Convert the input row from sliced coordinates back
        to original coordinates.

        Parameters
        ----------
        row : int
            Row in the sliced coordinate system

        Returns
        -------
        orig_row : int
            Row in the original coordinate system
        '''
        return row if self.original else self.start + row * self.step

    def find(self, key):
        return self.sliced_coords(self.index.find(key))

    def where(self, col_map):
        return self.sliced_coords(self.index.where(col_map))

    def range(self, lower, upper):
        return self.sliced_coords(self.index.range(lower, upper))

    def same_prefix(self, key):
        return self.sliced_coords(self.index.same_prefix(key))

    def sorted_data(self):
        return self.sliced_coords(self.index.sorted_data())

    def replace(self, row, col, val):
        if not self._frozen:
            self.index.replace(self.orig_coords(row), col, val)

    def copy(self):
        if not self.original:
            # replace self.index with a new object reference
            self.index = deepcopy(self.index)
        return self.index

    def insert_row(self, pos, vals, columns):
        if not self._frozen:
            self.copy().insert_row(self.orig_coords(pos), vals,
                                   columns)

    def get_row_specifier(self, row_specifier):
        return [self.orig_coords(x) for x in
                self.index.get_row_specifier(row_specifier)]

    def remove_rows(self, row_specifier):
        if not self._frozen:
            self.copy().remove_rows(row_specifier)

    def replace_rows(self, col_slice):
        if not self._frozen:
            self.index.replace_rows([self.orig_coords(x) for x in col_slice])

    def sort(self):
        if not self._frozen:
            self.copy().sort()

    def __repr__(self):
        if self.original:
            return repr(self.index)
        return 'Index slice {0} of\n{1}'.format(
            (self.start, self.stop, self.step), self.index)

    def __str__(self):
        return repr(self)

    def replace_col(self, prev_col, new_col):
        self.index.replace_col(prev_col, new_col)

    def reload(self):
        self.index.reload()

    def col_position(self, col_name):
        return self.index.col_position(col_name)

    def get_slice(self, col_slice, item):
        '''
        Return a newly created index from the given slice.

        Parameters
        ----------
        col_slice : Column object
            Already existing slice of a single column
        item : list or ndarray
            Slice for retrieval
        '''
        from .table import Table
        if len(self.columns) == 1:
            return Index([col_slice], engine=self.data.__class__)
        t = Table(self.columns, copy_indices=False)
        with t.index_mode('discard_on_copy'):
            new_cols = t[item].columns.values()
        return Index(new_cols, engine=self.data.__class__)

    @property
    def columns(self):
        return self.index.columns

    @property
    def data(self):
        return self.index.data


def get_index(table, table_copy):
    '''
    Inputs a table and some subset of its columns, and
    returns an index corresponding to this subset or None
    if no such index exists.

    Parameters
    ----------
    table : `Table`
        Input table
    table_copy : `Table`
        Subset of the columns in the table argument
    '''
    cols = set(table_copy.columns)
    indices = set()
    for column in cols:
        for index in table[column].info.indices:
            if set([x.info.name for x in index.columns]) == cols:
                return index
    return None


class _IndexModeContext(object):
    '''
    A context manager that allows for special indexing modes, which
    are intended to improve performance. Currently the allowed modes
    are "freeze", in which indices are not modified upon column modification,
    "copy_on_getitem", in which indices are copied upon column slicing,
    and "discard_on_copy", in which indices are discarded upon table
    copying/slicing.
    '''

    _col_subclasses = {}

    def __init__(self, table, mode):
        '''
        Parameters
        ----------
        table : Table
            The table to which the mode should be applied
        mode : str
            Either 'freeze', 'copy_on_getitem', or 'discard_on_copy'.
            In 'discard_on_copy' mode,
            indices are not copied whenever columns or tables are copied.
            In 'freeze' mode, indices are not modified whenever columns are
            modified; at the exit of the context, indices refresh themselves
            based on column values. This mode is intended for scenarios in
            which one intends to make many additions or modifications on an
            indexed column.
            In 'copy_on_getitem' mode, indices are copied when taking column
            slices as well as table slices, so col[i0:i1] will preserve
            indices.
        '''
        self.table = table
        self.mode = mode
        # Used by copy_on_getitem
        self._orig_classes = []
        if mode not in ('freeze', 'discard_on_copy', 'copy_on_getitem'):
            raise ValueError("Expected a mode of either 'freeze', "
                             "'discard_on_copy', or 'copy_on_getitem', got "
                             "'{0}'".format(mode))

    def __enter__(self):
        if self.mode == 'discard_on_copy':
            self.table._copy_indices = False
        elif self.mode == 'copy_on_getitem':
            for col in self.table.columns.values():
                self._orig_classes.append(col.__class__)
                col.__class__ = self._get_copy_on_getitem_shim(col.__class__)
        else:
            for index in self.table.indices:
                index._frozen = True

    def __exit__(self, exc_type, exc_value, traceback):
        if self.mode == 'discard_on_copy':
            self.table._copy_indices = True
        elif self.mode == 'copy_on_getitem':
            for col in reversed(self.table.columns.values()):
                col.__class__ = self._orig_classes.pop()
        else:
            for index in self.table.indices:
                index._frozen = False
                index.reload()

    def _get_copy_on_getitem_shim(self, cls):
        if cls in self._col_subclasses:
            return self._col_subclasses[cls]

        def __getitem__(self, item):
            value = cls.__getitem__(self, item)
            if type(value) is type(self):
                value = self.info.slice_indices(value, item, len(self))

            return value

        clsname = '_{0}WithIndexCopy'.format(cls.__name__)

        new_cls = type(clsname, (cls,), {'__getitem__': __getitem__})

        self._col_subclasses[cls] = new_cls

        return new_cls


class TableIndices(list):
    '''
    A special list of table indices allowing
    for retrieval by column name(s).

    Parameters
    ----------
    lst : list
        List of indices
    '''
    def __init__(self, lst):
        super(TableIndices, self).__init__(lst)

    def __getitem__(self, item):
        '''
        Retrieve an item from the list of indices.

        Parameters
        ----------
        item : int, str, tuple, or list
            Position in list or name(s) of indexed column(s)
        '''
        if isinstance(item, six.string_types):
            item = [item]
        if isinstance(item, (list, tuple)):
            item = list(item)
            for index in self:
                try:
                    for name in item:
                        index.col_position(name)
                    if len(index.columns) == len(item):
                        return index
                except ValueError:
                    pass
            # index search failed
            raise IndexError("No index found for {0}".format(item))

        return super(TableIndices, self).__getitem__(item)


class TableLoc(object):
    '''
    A pseudo-list of Table rows allowing for retrieval
    of rows by indexed column values.

    Parameters
    ----------
    table : Table
        Indexed table to use
    '''
    def __init__(self, table):
        self.table = table
        self.indices = table.indices
        if len(self.indices) == 0:
            raise ValueError("Cannot create TableLoc object with no indices")

    def __getitem__(self, item):
        '''
        Retrieve Table rows by value slice.

        Parameters
        ----------
        item : column element, list, ndarray, slice or tuple
            Can be a value of the table primary index, a list/ndarray
            of such values, or a value slice (both endpoints are included).
            If a tuple is provided, the first element must be
            an index to use instead of the primary key, and the
            second element must be as above.
        '''

        if isinstance(item, tuple):
            key, item = item
        else:
            key = self.table.primary_key

        index = self.indices[key]
        if len(index.columns) > 1:
            raise ValueError("Cannot use .loc on multi-column indices")

        if isinstance(item, slice):
            # None signifies no upper/lower bound
            start = MinValue() if item.start is None else item.start
            stop = MaxValue() if item.stop is None else item.stop
            rows = index.range((start,), (stop,))
        else:
            if not isinstance(item, (list, np.ndarray)): # single element
                item = [item]
            # item should be a list or ndarray of values
            rows = []
            for key in item:
                rows.extend(index.find((key,)))

        if len(rows) == 0: # no matches found
            raise KeyError('No matches found for key {0}'.format(item))
        elif len(rows) == 1: # single row
            return self.table[rows[0]]
        return self.table[rows]


class TableILoc(TableLoc):
    '''
    A variant of TableLoc allowing for row retrieval by
    indexed order rather than data values.

    Parameters
    ----------
    table : Table
        Indexed table to use
    '''
    def __init__(self, table):
        super(TableILoc, self).__init__(table)

    def __getitem__(self, item):
        if isinstance(item, tuple):
            key, item = item
        else:
            key = self.table.primary_key
        index = self.indices[key]
        rows = index.sorted_data()[item]
        table_slice = self.table[rows]

        if len(table_slice) == 0: # no matches found
            raise IndexError('Invalid index for iloc: {0}'.format(item))

        return table_slice
