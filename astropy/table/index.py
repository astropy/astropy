from copy import deepcopy
import numpy as np

from .bst import BST, RedBlackTree, FastBST, FastRBT
from .array import SortedArray

'''
The Index class can use several implementations as its
engine. Any implementation should implement the following:

__init__([dict] lines) : initializes based on key/row list pairs
add(key, row) -> None : add (key, row) to existing data
remove(key, data=None) -> boolean : remove data from self[key], or all of
                                    self[key] if data is None
reorder(row) -> None : decrement row numbers after row
find(key) -> list : list of rows corresponding to key
range(lower, upper) -> list : rows in self[k] where lower<=k<=upper
sort() -> list of rows in sorted order (by key)
replace_rows(row_map) -> None : replace row numbers based on slice


Note: when a Table is initialized from another Table, indices are
(deep) copied and their columns are set to the columns of the new Table.

Column creation:
Column(c) -> deep copy of indices
c[1:2] -> deep copy and reordering of indices
array.view(Column) -> no indices
'''

class MaxValue(object):
    def __gt__(self, other):
        return True

    def __ge__(self, other):
        return True

    def __lt__(self, other):
        return False

    def __le__(self, other):
        return False

class MinValue(object):
    def __lt__(self, other):
        return True

    def __le__(self, other):
        return True

    def __gt__(self, other):
        return False

    def __ge__(self, other):
        return False

class Index:
    def __init__(self, columns, impl=None):
        if impl is None:
            impl = FastRBT

        if columns is None: # this creates a special exception for deep copying
            columns = []
        elif len(columns) == 0:
            raise ValueError("Cannot create index without at least one column")

        # nodes of self.data will be (key val, row index)
        iterable = columns[0] if len(columns) == 1 else zip(*columns)
        lines = {}
        for val, key in enumerate(iterable):
            if not key in lines:
                lines[key] = [val]
            else:
                lines[key].append(val)
        self.data = impl(lines)
        self.columns = columns

    def refresh(self, columns):
        self.columns = [columns[x.name] for x in self.columns]

    def col_position(self, col):
        for i, c in enumerate(self.columns):
            if col.name == c.name:
                return i
        raise ValueError("Column does not belong to index: {0}".format(col))

    def insert_row(self, pos, vals, columns):
        key = [None] * len(self.columns)
        for i, col in enumerate(columns):
            try:
                key[i] = vals[self.col_position(col)]
            except ValueError: # not a member of index
                continue
        self.data.add(tuple(key), pos)

    def remove_rows(self, row_specifier):
        # row_specifier must be an int, list of ints, ndarray, or slice
        if isinstance(row_specifier, int):
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
            self.data.reorder(row)

    def remove_row(self, row, reorder=True):
        if not self.data.remove(tuple([col[row] for col in self.columns]),
                                data=row):
            raise ValueError("Could not remove row {0} from index".format(row))
        # decrement the row number of all later rows
        if reorder:
            self.data.reorder(row)

    def find(self, key):
        return self.data.find(key)

    def where(self, col_map):
        # ensure that the keys of col_map form a left prefix of index columns
        # also, a range query can only be on the last of the index columns
        # note: if a range is invalid (upper < lower), there will be no results
        names = [col.name for col in self.columns]
        query_names = col_map.keys()
        if set(names[:len(query_names)]) != set(query_names):
            raise ValueError("Query columns must form a left prefix of "
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
            lower = base + [col_map[last_col][0]]
            upper = base + [col_map[last_col][1]]
            if len(lower) == len(self.columns):
                return self.data.range(tuple(lower), tuple(upper))
            else:
                return self.same_prefix_range(lower, upper)
        else:
            key = base + [col_map[query_names[-1]]]
            if len(key) == len(self.columns):
                return self.data.find(tuple(key))
            else:
                return self.same_prefix(key)

    def range(self, lower, upper):
        return self.data.range(lower, upper)

    def same_prefix(self, key):
        return self.same_prefix_range(key, key)

    def same_prefix_range(self, lower, upper):
        n = len(lower)
        ncols = len(self.columns)
        lower = tuple(lower + (ncols - n) * [MinValue()])
        upper = tuple(upper + (ncols - n) * [MaxValue()])
        return self.data.range(lower, upper)

    def replace(self, row, col, val):
        self.remove_row(row, reorder=False)
        key = [c[row] for c in self.columns]
        key[self.col_position(col)] = val
        self.data.add(tuple(key), row)

    def replace_rows(self, col_slice):
        row_map = dict((row, i) for i, row in enumerate(col_slice))
        self.data.replace_rows(row_map)

    def sorted_data(self):
        return self.data.sort()

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return str(self)

    def __deepcopy__(self, memo):
        # deep copy must be overridden to perform a shallow copy of columns
        index = Index(None, impl=self.data.__class__)
        index.data = deepcopy(self.data, memo)
        index.columns = self.columns[:] # new list, same columns
        memo[id(self)] = index
        return index

def get_index(table):
    '''
    Returns a sorted table if an index covers the entire table,
    and None otherwise.
    '''
    cols = set(table.columns.values())
    indices = set()
    for column in cols:
        for index in column.indices:
            if set(index.columns) == cols:
                return index
    return None
