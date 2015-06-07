from bst import BST

class Index:
    Implementation = BST
    def __init__(self, data):
        self.data = self.Implementation(data)
        # nodes of self.data will be (key val, row index)

    def insert_row(self, pos, val):
        self.data.add((val, pos))

    def remove_rows(self, row_specifier):
        if isinstance(row_specifier, int):
            self.data.remove(row_specifier)
        elif isinstance(row_specifier, list): ##TODO: check other iterables
            for row in row_specifier:
                self.data.remove(row)
        else: # must be slice
            max_row = max([row for val, row in data])
            for row in row_specifier.indices(max_row):
                self.data.remove(row)

    def find(self, key):
        return self.data.find(key)

    def range(self, lower, upper):
        return self.data.range(lower, upper)

    def replace(self, pos, val):
        self.remove_rows(pos)
        self.insert_row(pos, val)
