from .pprint import _pformat_col


class FunctionColumn(object):
    """
    description : str or None
        Full description of column
    units : str or None
        Physical units
    format : str or None or function
        Format string for outputting column values.  This can be an
        "old-style" (``format % value``) or "new-style" (`str.format`)
        format specification string or a function that accepts a single
        value and returns a string.
    meta : dict-like or None
        Meta-data associated with the column
    """
    def __init__(self, func, cols, name=None, description=None, units=None, format=None):
        self._col_dependencies = [col.name for col in cols]
        self.func = func
        self.cols = cols
        self.name = name
        self.__print_format__ = format
        self.description = description
        self.units = units

    @property
    def shape(self):
        return (len(self.parent_table),)  # FIX ME

    @property
    def data(self):
        return self.func(*self.cols)

    def __len__(self):
        return len(self.parent_table)

    def __getitem__(self, item):
        return self.data[item]

    def copy(self, data=None, copy_data=False):
        print 'Calling copy on id={0}'.format(id(self))
        # FIX ME
        self_copy = FunctionColumn(self.func, self.col_names, name=self.name,
                                   description=self.description, units=self.units)
        return self_copy

    def __str__(self):
        lines, n_header = _pformat_col(self)
        return '\n'.join(lines)

    @property
    def format(self):
        return self.__print_format__

    @format.setter
    def format(self, value):
        self.__print_format__ = value

    def __table_replicate__(self, table):
        """
        Replicate the current column but using a new ``table`` (Table object).
        """
        new_col = self.copy()
        new_col.parent_table = table
        return new_col

    def __table_add_column__(self, table, index):
        from .table import TableColumns
        print 'Here in table add column'

        self.parent_table = table
        columns = TableColumns()

        for i_column, column in enumerate(table.columns.values()):
            if i_column == index:
                columns[self.name] = self
            columns[column.name] = column
        else:
            columns[self.name] = self

        table.columns = columns


class ViewColumnOrig(object):
    def __init__(self, col, table=None, name=None):
        from .table import BaseColumn
        if isinstance(col, BaseColumn):
            self.data_name = col.name
            # Make an internal copy of the input column (including metadata etc)
            # but keep the reference to the original data column
            self._col = col.copy(copy_data=False)  # data=col.data ?
        elif isinstance(col, self.__class__):
            self.data_name = col.data_name
            if table is None:
                raise ValueError('Must supply `table` when constucting new ViewColumn '
                                 'from existing ViewColumn object')
            # Since `col` is already a ViewColumn we need to go down one layer and
            # get the existing _col (which is a plain Column or MaskedColumn) and
            # copy it.  In this case we update the data reference.
            self._col = col._col.copy(data=table._data[col.data_name], copy_data=False)
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

    def __table_replicate__(self, table):
        """
        Replicate the current column but using a new ``data`` ndarray.
        """
        newcol = ViewColumn(self, table, name=self.name)
        return newcol

    def __table_add_column__(self, table, index):
        from .table import TableColumns
        print 'Here in table add column'

        if isinstance(self, ViewColumn):
            columns = TableColumns()

            for i_column, column in enumerate(table.columns.values()):
                if i_column == index:
                    columns[self.name] = self
                columns[column.name] = column
            else:
                columns[self.name] = self
            table.columns = columns
