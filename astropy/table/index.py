from .bst import BST
from .array import SortedArray

class Index:
    Implementation = SortedArray

    def __init__(self, columns):
        # nodes of self.data will be (key val, row index)
        data = [(c, i) for i, c in enumerate(zip(*columns))]
        self.data = self.Implementation(data)
        self.columns = columns

    def refresh(self, columns):
        self.columns = [columns[x.name] for x in self.columns]

    def col_position(self, col):
        for i, c in enumerate(self.columns):
            if col is c:
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
        if isinstance(row_specifier, int):
            self.remove_row(row_specifier)
        elif isinstance(row_specifier, list): ##TODO: check other iterables
            for row in row_specifier:
                self.remove_row(row)
        else: # must be slice
            max_row = max([row for val, row in row_specifier])
            for row in row_specifier.indices(max_row):
                self.remove_row(row)

    def remove_row(self, row, reorder=True):
        if not self.data.remove(tuple([col[row] for col in self.columns]),
                                data=row):
            raise ValueError("Could not remove row {0} from index".format(row))
        # decrement the row number of all later rows
        if reorder:
            for node in self.data.traverse('inorder'):
                node.data = [x - 1 if x > row else x for x in node.data]

    def find(self, key):
        node = self.data.find(key)
        return [] if node is None else node.data

    def where(self, col_map):
        # ensure that the keys of col_map form a left prefix of index columns
        # also, a range query can only be on the last of the index columns
        # note: if a range is invalid (upper < lower), there will be no results
        names = [col.name for col in self.columns]
        query_names = col_map.keys()
        if set(names[:len(query_names)]) != set(query_names):
            raise ValueError("Query columns must form a left prefix of "
                             "index columns")
        query_names = names[:len(query_names)] # reorder query_names
        for name in query_names[:-1]:
            if isinstance(col_map[name], tuple):
                raise ValueError("Range queries are only valid on the "
                                 "last column of an index")
        base = [col_map[name] for name in query_names[:-1]]
        if isinstance(col_map[query_names[-1]], tuple): # range query
            lower = base + [col_map[query_names[-1]][0]]
            upper = base + [col_map[query_names[-1]][1]]
            lst = [x.data for x in self.data.range(tuple(lower), tuple(upper))]
        else:
            key = base + [col_map[query_names[-1]]]
            if len(key) == len(self.columns):
                lst = [self.find(tuple(key))]
            else:
                lst = [x.data for x in self.data.same_prefix(tuple(key))]
        return [row for l in lst for row in l]

    def range(self, lower, upper):
        return self.data.range(lower, upper)

    def replace(self, row, col, val):
        self.remove_row(row, reorder=False)
        key = [c[row] for c in self.columns]
        key[self.col_position(col)] = val
        self.data.add(tuple(key), row)

    def sorted_data(self):
        lst = [x.data for x in self.data.sort()]
        return [row for l in lst for row in l]

    def nodes(self):
        # for debugging purposes
        return self.data.nodes()

    def __str__(self):
        return str(self.nodes())

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
