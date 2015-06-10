from bst import BST

class Index:
    Implementation = BST
    def __init__(self, column):
        # nodes of self.data will be (key val, row index)
        data = [(c, i) for i, c in enumerate(column)]
        self.data = self.Implementation(data)

    def insert_row(self, pos, val):
        self.data.add(val, pos)

    def remove_rows(self, row_specifier, col):
        if isinstance(row_specifier, int):
            self.remove_row(row_specifier, col)
        elif isinstance(row_specifier, list): ##TODO: check other iterables
            for row in row_specifier:
                self.remove_row(row, col)
        else: # must be slice
            max_row = max([row for val, row in data])
            for row in row_specifier.indices(max_row):
                self.remove_row(row, col)

    def remove_row(self, row, col):
        self.data.remove(col[row], data=row)

    def find(self, key):
        return self.data.find(key)

    def range(self, lower, upper):
        return self.data.range(lower, upper)

    def replace(self, pos, val, col):
        self.remove_rows(pos, col)
        self.insert_row(pos, val)
