from copy import deepcopy
import numpy as np

from .bst import BST, RedBlackTree, FastBST, FastRBT
from .sorted_array import SortedArray

'''
The Index class can use several implementations as its
engine. Any implementation should implement the following:

__init__([dict] lines) : initializes based on key/row list pairs
add(key, row) -> None : add (key, row) to existing data
remove(key, data=None) -> boolean : remove data from self[key], or all of
                                    self[key] if data is None
shift_left(row) -> None : decrement row numbers after row
shift_right(row) -> None : increase row numbers >= row
find(key) -> list : list of rows corresponding to key
range(lower, upper, bounds) -> list : rows in self[k] where k is between
                               lower and upper (<= or < based on bounds)
sort() -> list of rows in sorted order (by key)
replace_rows(row_map) -> None : replace row numbers based on slice
items() -> list of tuples of the form (key, data)


Note: when a Table is initialized from another Table, indices are
(deep) copied and their columns are set to the columns of the new Table.

Column creation:
Column(c) -> deep copy of indices
c[[1, 2]] -> deep copy and reordering of indices
c[1:2] -> reference
array.view(Column) -> no indices
'''

class MaxValue(object):
    '''
    Represents an infinite value for purposes
    of tuple comparison.
    '''
    def __gt__(self, other):
        return True

    def __ge__(self, other):
        return True

    def __lt__(self, other):
        return False

    def __le__(self, other):
        return False

    def __repr__(self):
        return "MAX"

    __le__ = __lt__
    __ge__ = __gt__
    __str__ = __repr__


class MinValue(object):
    '''
    The opposite of MaxValue, i.e. a representation of
    negative infinity.
    '''
    def __lt__(self, other):
        return True

    def __le__(self, other):
        return True

    def __gt__(self, other):
        return False

    def __ge__(self, other):
        return False

    def __repr__(self):
        return "MIN"

    __le__ = __lt__
    __ge__ = __gt__
    __str__ = __repr__

class QueryError(ValueError):
    '''
    Indicates that a given index cannot handle the supplied query.
    '''
    pass

class Index:
    '''
    The Index class makes it possible to maintain indices
    on columns of a Table, so that column values can be queried
    quickly and efficiently. Column values are stored in lexicographic
    sorted order, which allows for binary searching in O(log n).
    '''
    def __init__(self, columns, impl=None, col_dtypes=None, data=None):
        from .table import Table
        if data is not None: # create from data
            self.engine = data.__class__
            self.data = data
            self.columns = columns
            return

        # by default, use SortedArray
        self.engine = impl or SortedArray

        if columns is None: # this creates a special exception for deep copying
            columns = []
            lines = []
        elif len(columns) == 0:
            raise ValueError("Cannot create index without at least one column")
        else:
            num_rows = len(columns[0])
            # sort the table lexicographically and keep row numbers
            table = Table(columns + [np.arange(num_rows)])
            lines = table[np.lexsort(columns[::-1])]

        if self.engine == SortedArray:
            self.data = self.engine(lines, col_dtypes=col_dtypes)
        else:
            self.data = self.engine(lines)
        self.columns = columns

    def refresh(self, columns):
        # update index to include correct column references
        self.columns = [columns[x.name] for x in self.columns]

    def col_position(self, col):
        # return the position of col in self.columns
        for i, c in enumerate(self.columns):
            if col.name == c.name:
                return i
        raise ValueError("Column does not belong to index: {0}".format(col))

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
                key[i] = vals[self.col_position(col)]
            except ValueError: # not a member of index
                continue
        # shift all rows >= pos to the right
        self.data.shift_right(pos)
        self.data.add(tuple(key), pos)

    def remove_rows(self, row_specifier):
        # row_specifier must be an int, list of ints, ndarray, or slice
        if isinstance(row_specifier, int):
            # single row
            self.remove_row(row_specifier)
            return
        elif isinstance(row_specifier, (list, np.ndarray)):
            iterable = row_specifier
        elif isinstance(row_specifier, slice):
            col_len = len(self.columns[0])
            iterable = range(*row_specifier.indices(col_len))
        else:
            raise ValueError("Expected int, array of ints, or slice but "
                             "got {0} in remove_rows".format(row_specifier))

        rows = []

        # To maintain the correct row order, we loop twice,
        # deleting rows first and then reordering the remaining rows
        for row in iterable:
            self.remove_row(row, reorder=False)
            rows.append(row)
        # second pass - row order is reversed to maintain
        # correct row numbers
        for row in reversed(sorted(rows)):
            self.data.shift_left(row)

    def remove_row(self, row, reorder=True):
        # for removal, form a key consisting of column values in this row
        if not self.data.remove(tuple([col[row] for col in self.columns]),
                                data=row):
            raise ValueError("Could not remove row {0} from index".format(row))
        # decrement the row number of all later rows
        if reorder:
            self.data.shift_left(row)

    def find(self, key):
        # return the row values corresponding to key, in sorted order
        return self.data.find(key)

    def where(self, col_map):
        # ensure that the keys of col_map form a left prefix of index columns
        # also, a range query can only be on the last of the index columns
        # note: if a range is invalid (upper < lower), there will be no results
        names = [col.name for col in self.columns]
        query_names = col_map.keys()
        if set(names[:len(query_names)]) != set(query_names):
            raise QueryError("Query columns must form a left prefix of "
                             "index columns")
        # query_names is a prefix of index column names
        query_names = names[:len(query_names)]
        for name in query_names[:-1]:
            if isinstance(col_map[name], tuple):
                raise ValueError("Range queries are only valid on the "
                                 "last column of an index")
        base = [col_map[name] for name in query_names[:-1]]
        last_col = query_names[-1]

        if isinstance(col_map[last_col], tuple): # range query
            lower = base + [col_map[last_col][0][0]]
            upper = base + [col_map[last_col][0][1]]
            bounds = col_map[last_col][1]
            # bounds is a tuple of True (<=) or False (<)
            if len(lower) == len(self.columns):
                result = self.data.range(tuple(lower), tuple(upper), bounds)
            else:
                result = self.same_prefix_range(lower, upper, bounds)
        else:
            key = base + [col_map[query_names[-1]]]
            if len(key) == len(self.columns):
                result = self.data.find(tuple(key))
            else:
                result = self.same_prefix(key)
        return sorted(result)

    def range(self, lower, upper):
        # return values between lower and upper
        return self.data.range(lower, upper)

    def same_prefix(self, key):
        # return rows whose keys contain supplied key as a prefix
        return self.same_prefix_range(key, key, (True, True))

    def same_prefix_range(self, lower, upper, bounds):
        '''
        Return rows whose keys have a prefix between lower and upper.
        The parameter bounds is (x, y) where x <-> inclusive lower bound,
        y <-> inclusive upper bound
        '''
        n = len(lower)
        ncols = len(self.columns)
        a = MinValue() if bounds[0] else MaxValue()
        b = MaxValue() if bounds[1] else MinValue()
        # [x, y] search corresponds to [(x, min), (y, max)]
        # (x, y) search corresponds to ((x, max), (x, min))
        lower = tuple(lower + (ncols - n) * [a])
        upper = tuple(upper + (ncols - n) * [b])
        return self.data.range(lower, upper, bounds)

    def replace(self, row, col, val):
        '''
        Replace the value of a column at a given position.

        Parameters
        ----------
        row : int
            Row number to modify
        col : Column
            Column to modify
        val : col.dtype
            Value to insert at specified row of col
        '''
        self.remove_row(row, reorder=False)
        key = [c[row] for c in self.columns]
        key[self.col_position(col)] = val
        self.data.add(tuple(key), row)

    def replace_rows(self, col_slice):
        '''
        Modify rows in this index to agree with the specified
        col_slice. For example, given an index
        {'5': 1, '2': 0, '3': 2} on a column ['2', '5', '3'],
        an input col_slice of [2, 0] will result in the relabeling
        {'3': 0, '2': 1} on the sliced column ['3', '2'].
        '''
        row_map = dict((row, i) for i, row in enumerate(col_slice))
        self.data.replace_rows(row_map)

    def sorted_data(self):
        '''
        Returns a list of rows in sorted order based on keys;
        essentially acts as an argsort() on columns
        '''
        return self.data.sort()

    def __getitem__(self, item):
        # item must be a slice; return sliced copy
        return SlicedIndex(self, item)

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return str(self)

    def __deepcopy__(self, memo):
        # deep copy must be overridden to perform a shallow copy of columns
        num_cols = self.data.num_cols if self.engine == SortedArray else None
        index = Index(None, impl=self.data.__class__, col_dtypes=[x.dtype for
                                                    x in self.columns])
        index.data = deepcopy(self.data, memo)
        index.columns = self.columns[:] # new list, same columns
        memo[id(self)] = index
        return index

class SlicedIndex:
    '''
    This class provides a wrapper around an actual Index object
    to make index slicing function correctly. Since numpy expects
    array slices to provide an actual data view, a SlicedIndex should
    retrieve data directly from the original index and then adapt
    it to the sliced coordinate system as appropriate.
    '''
    def __init__(self, index, index_slice):
        self.index = index
        num_rows = len(index.columns[0])
        if isinstance(index_slice, tuple):
            self.start, self.stop, self.step = index_slice
        else: # index_slice is an actual slice
            self.start, self.stop, self.step = index_slice.indices(num_rows)
        self.length = 1 + (self.stop - self.start - 1) // self.step

    def __getitem__(self, item):
        # item must be a slice
        if self.length <= 0:
            # empty slice
            return SlicedIndex(self.index, slice(1, 0))
        start, stop, step = item.indices(self.length)
        new_start = self.orig_coords(start)
        new_stop = self.orig_coords(stop)
        new_step = self.orig_coords(step)
        return SlicedIndex(self.index, (new_start, new_stop, new_step))

    def sliced_coords(self, rows):
        # Convert the input rows to the sliced coordinate system
        if self.step > 0:
            return [(row - self.start) / self.step for row in rows
                    if self.start <= row < self.stop and
                    (row - self.start) % self.step == 0]
        else:
            return [(row - self.start) / self.step for row in rows
                    if self.stop < row <= self.start and
                    (row - self.start) % self.step == 0]

    def orig_coords(self, row):
        '''
        Convert the input row from sliced coordinates back
        to original coordinates.
        '''
        return self.start + row * self.step

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
        return self.index.replace(self.orig_coords(row), col, val)

    def __repr__(self):
        return 'Index slice {0} of {1}'.format(
            (self.start, self.stop, self.step), self.index)

    ##TODO: adapt other Index methods here


def get_index(table, table_copy):
    '''
    Inputs a table and some subset of its columns, and
    returns an index corresponding to this subset or None
    if no such index exists.
    '''
    cols = set(table_copy.columns)
    indices = set()
    for column in cols:
        for index in table[column].indices:
            if set([x.name for x in index.columns]) == cols:
                return index
    return None

class static_indices:
    # provides a context in which Table indices
    # are not copied
    def __init__(self, table):
        self.table = table

    def __enter__(self):
        self.table._copy_indices = False

    def __exit__(self, exc_type, exc_value, traceback):
        self.table._copy_indices = True
