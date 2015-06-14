from .bst import BST

class Index:
    Implementation = BST

    def __init__(self, columns):
        # nodes of self.data will be (key val, row index)
        data = [(c, i) for i, c in enumerate(zip(*columns))]
        self.data = self.Implementation(data)
        self.columns = columns

    def insert_row(self, pos, vals, columns):
        key = [vals[columns.index(col)] for col in self.columns]
        self.data.add(tuple(key), pos)

    def remove_rows(self, row_specifier):
        if isinstance(row_specifier, int):
            self.remove_row(row_specifier)
        elif isinstance(row_specifier, list): ##TODO: check other iterables
            for row in row_specifier:
                self.remove_row(row)
        else: # must be slice
            max_row = max([row for val, row in data])
            for row in row_specifier.indices(max_row):
                self.remove_row(row)

    def remove_row(self, row):
        self.data.remove(tuple([col[row] for col in self.columns]), data=row)

    def find(self, *key):
        node = self.data.find(key)
        return [] if node is None else node.data

    def range(self, lower, upper):
        return self.data.range(lower, upper)

    def replace(self, pos, val):
        self.remove_rows(pos)
        self.insert_row(pos, val)

    def sorted_data(self):
        lst = [x.data for x in self.data.traverse('inorder')]
        return [row for l in lst for row in l]

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
