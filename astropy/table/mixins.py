class ViewColumn(object):
    def __init__(self, col, data=None, name=None):
        from .table import BaseColumn
        if isinstance(col, BaseColumn):
            self.data_name = col.name
            # Make an internal copy of the input column (including metadata etc)
            # but keep the reference to the original data column
            self._col = col.copy(copy_data=False)  # data=col.data ?
        elif isinstance(col, self.__class__):
            self.data_name = col.data_name
            if data is None:
                raise ValueError('Must supply `data` when constucting new ViewColumn '
                                 'from existing ViewColumn object')
            # Since `col` is already a ViewColumn we need to go down one layer and
            # get the existing _col (which is a plain Column or MaskedColumn) and
            # copy it.  In this case we update the data reference.
            self._col = col._col.copy(data=data[col.data_name], copy_data=False)
        self.dependencies = [self.data_name]
        self.name = name
        self.format = col.format
        self.description = col.description
        self.units = col.units

    @property
    def shape(self):
        return self._col.shape

    @property
    def data(self):
        return self._col.data

    def __len__(self):
        return len(self._col)

    # @property
    # def name(self):
    #    return self.col.name

    def __getitem__(self, item):
        return self._col[item]

    def copy(self, data=None, copy_data=False):
        # FIX ME
        return self

    @property
    def format(self):
        return self.__print_format__

    @format.setter
    def format(self, value):
        self.__print_format__ = value

    def __table_replicate__(self, data):
        """
        Replicate the current column but using a new ``data`` ndarray.
        """
        newcol = ViewColumn(self, data, name=self.name)
        return newcol

    def __table_add_column__(self, table, index):
        from .table import TableColumns
        print 'Here in table add column'
        if index is None:
            index = len(table.columns)

        if isinstance(self, ViewColumn):
            columns = TableColumns()

            for i_column, column in enumerate(table.columns.values()):
                if i_column == index:
                    columns[self.name] = self
                columns[column.name] = column
            if index is None or index == len(table.columns):
                columns[self.name] = self
            table.columns = columns
