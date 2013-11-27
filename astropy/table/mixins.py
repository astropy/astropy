from .pprint import _pformat_col


class FunctionColumn(object):
    """
    Mixin column that computes a function ``func`` of columns specified by ``col_names``.

    This is a very simple demonstration class which is effectively an ndarray wrapper
    class with a lot of the useful interface API missing.  Maybe inherit from Column?
    """
    def __init__(self, func, col_names, name=None, description=None, units=None, format=None):
        self.func = func
        self.col_names = col_names
        self.name = name
        self.__print_format__ = format
        self.description = description
        self.units = units
        self.parent_table = None

    @property
    def shape(self):
        return self.data.shape if self.data is not None else ()

    @property
    def data(self):
        # self.data is None prior to binding to a parent table
        if self.parent_table is None:
            return None

        # TODO: update _data if other cols have changed since _data was computed.
        if not hasattr(self, '_data'):
            cols = [self.parent_table[col_name] for col_name in self.col_names]
            self._data = self.func(*cols)
        return self._data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.data[item]

    def copy(self, data=None, copy_data=False):
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

    def __table_replicate_column__(self, table, name):
        """
        Replicate the current column but using a new ``table`` (Table object).
        """
        new_col = self.copy()
        new_col.name = name
        new_col.parent_table = table
        return new_col

    def __table_get_columns__(self, name, ColumnClass):
        """
        Get list of columns that represent the internal data for this object assuming the
        mixin column name is ``name``.
        """
        # No internal columns, it relies on already existent (independent) columns
        return []


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
