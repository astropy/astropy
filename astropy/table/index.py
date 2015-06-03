DataStructure = lambda x: None # will change

class Index:
    def __init__(self, data):
        self.data = DataStructure(data)
        # nodes of self.data will be (key val, row index)
    def insert_row(self, pos, val):
        pass
    def remove_rows(self, row_specifier):
        pass
    def find(self, key):
        pass
    def range(self, lower, upper):
        pass
    def replace(self, pos, val):
        self.remove_rows(pos)
        self.insert_row(pos, val)
