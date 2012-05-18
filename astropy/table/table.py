import sys
from copy import deepcopy
import collections

import numpy as np

from ..utils import OrderedDict, isiterable
from .structhelper import _drop_fields
from .pprint import _pformat_table, _pformat_col, _more_tabcol
from ..utils.console import color_print

# Python 2 and 3 source compatibility
try:
    unicode
except NameError:
    unicode = basestring = str

AUTO_COLNAME = 'col{0}'


def _auto_names(n_cols):
    return [AUTO_COLNAME.format(i) for i in range(n_cols)]


class TableColumns(OrderedDict):
    """OrderedDict subclass for a set of columns.

    This class enhances item access to provide convenient access to columns
    by name or index, including slice access.  It also handles renaming
    of columns.

    The initialization argument ``cols`` can be any structure that is valid
    for initializing a Python dict.  This includes a dict, list of
    (key, val) tuple pairs, list of [key, val] lists, etc.

    Parameters
    ----------
    cols : dict, list, tuple; optional
        Column objects as data structure that can init dict (see above)
    """

    def __init__(self, cols={}):
        if isinstance(cols, (list, tuple)):
            cols = [(col.name, col) for col in cols]
        super(TableColumns, self).__init__(cols)

    def __getitem__(self, item):
        """Get items from a TableColumns object.
        ::

          tc = TableColumns(cols=[Column('a'), Column('b'), Column('c')])
          tc['a']  # Column('a')
          tc[1] # Column('b')
          tc['a', 'b'] # <TableColumns names=('a', 'b')>
          tc[1:3] # <TableColumns names=('b', 'c')>
        """
        if isinstance(item, basestring):
            return OrderedDict.__getitem__(self, item)
        elif isinstance(item, int):
            return self.values()[item]
        elif isinstance(item, tuple):
            return TableColumns([self[x] for x in item])
        elif isinstance(item, slice):
            return TableColumns([self[x] for x in self.keys()[item]])
        else:
            raise IndexError('Illegal key or index value for TableColumns '
                             'object')

    def __repr__(self):
        names = ("'{0}'".format(x) for x in self.keys())
        return "<TableColumns names=({0})>".format(",".join(names))

    def _rename_column(self, name, new_name):
        if new_name in self:
            raise KeyError("Column {0} already exists".format(new_name))

        mapper = {name: new_name}
        new_names = [mapper.get(name, name) for name in self]
        cols = self.values()
        self.clear()
        self.update(zip(new_names, cols))

    # Define keys and values for Python 2 and 3 source compatibility
    def keys(self):
        return list(OrderedDict.keys(self))

    def values(self):
        return list(OrderedDict.values(self))


class Column(np.ndarray):
    """Define a data column for use in a Table object.

    Parameters
    ----------
    name : str
        Column name and key for reference within Table
    data : list, ndarray or None
        Column data values
    dtype : numpy.dtype compatible value
        Data type for column
    shape : tuple or ()
        Dimensions of a single row element in the column data
    length : int or 0
        Number of row elements in column data
    description : str or None
        Full description of column
    units : str or None
        Physical units
    format : str or None
        Format string for outputting column values.  This can be an
        "old-style" (``format % value``) or "new-style" (`str.format`)
        format specification string.
    meta : dict-like or None
        Meta-data associated with the column

    Examples
    --------
    A Column can be created in two different ways:

    - Provide a ``data`` value and optionally a ``dtype`` value

      Examples::

        col = Column('name', data=[1, 2, 3])         # shape=(3,)
        col = Column('name', data=[[1, 2], [3, 4]])  # shape=(2, 2)
        col = Column('name', data=[1, 2, 3], dtype=float)
        col = Column('name', np.array([1, 2, 3]))
        col = Column('name', ['hello', 'world'])

      The ``dtype`` argument can be any value which is an acceptable
      fixed-size data-type initializer for the numpy.dtype() method.  See
      `<http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html>`_.
      Examples include:

      - Python non-string type (float, int, bool)
      - Numpy non-string type (e.g. np.float32, np.int64, np.bool)
      - Numpy.dtype array-protocol type strings (e.g. 'i4', 'f8', 'S15')

      If no ``dtype`` value is provide then the type is inferred using
      ``np.array(data)``.  When ``data`` is provided then the ``shape``
      and ``length`` arguments are ignored.

    - Provide zero or more of ``dtype``, ``shape``, ``length``

      Examples::

        col = Column('name')
        col = Column('name', dtype=int, length=10, shape=(3,4))

      The default ``dtype`` is ``np.float64`` and the default ``length`` is
      zero.  The ``shape`` argument is the array shape of a single cell in the
      column.  The default ``shape`` is () which means a single value in each
      element.
    """

    def __new__(cls, name, data=None,
                 dtype=None, shape=(), length=0,
                 description=None, units=None, format=None, meta=None):

        if data is None:
            dtype = (np.dtype(dtype).str, shape)
            self_data = np.zeros(length, dtype=dtype)
        elif isinstance(data, Column):
            self_data = np.asarray(data.data, dtype=dtype)
            if description is None:
                description = data.description
            if units is None:
                units = units or data.units
            if format is None:
                format = data.format
            if meta is None:
                meta = deepcopy(data.meta)
        else:
            self_data = np.asarray(data, dtype=dtype)

        self = self_data.view(cls)
        self._name = name
        self.units = units
        self.format = format
        self.description = description
        self.parent_table = None

        self.meta = OrderedDict()
        if meta is not None:
            self.meta.update(meta)

        return self

    def __array_finalize__(self, obj):
        # Obj will be none for direct call to Column() creator
        if obj is None:
            return

        # Self was created from template (e.g. obj[slice] or (obj * 2))
        # or viewcast e.g. obj.view(Column).  In either case we want to
        # init Column attributes for self from obj if possible.
        self.parent_table = None
        for attr in ('name', 'units', 'format', 'description'):
            val = getattr(obj, attr, None)
            setattr(self, attr, val)
        self.meta = deepcopy(getattr(obj, 'meta', {}))

    def _get_name(self):
        return self._name

    def _set_name(self, val):
        if self.parent_table is not None:
            table = self.parent_table
            table.columns._rename_column(self.name, val)
            table._data.dtype.names = table.columns.keys()

        self._name = val

    name = property(_get_name, _set_name)

    @property
    def data(self):
        return self.view(np.ndarray)

    def copy(self, data=None, copy_data=True):
        """Return a copy of the current Column instance.
        """
        if data is None:
            data = self.view(np.ndarray)
            if copy_data:
                data = data.copy()

        return Column(self.name, data, units=self.units, format=self.format,
                      description=self.description, meta=deepcopy(self.meta))

    @property
    def descr(self):
        """Array-interface compliant full description of the column.

        This returns a 3-tuple (name, type, shape) that can always be
        used in a structured array dtype definition.
        """
        return (self.name, self.dtype.str, self.shape[1:])

    def __repr__(self):
        if self.name:
            out = "<Column name={0} units={1} format={2} " \
                "description={3}>\n{4}".format(
                repr(self.name), repr(self.units),
                repr(self.format), repr(self.description), repr(self.data))
        else:
            out = repr(self.data)
        return out

    def attrs_equal(self, col):
        """Compare the column attributes of ``col`` to this object.

        The comparison attributes are: name, units, dtype, format, description,
        and meta.

        Parameters
        ----------
        col: Column
            Comparison column

        Returns
        -------
        equal: boolean
            True if all attributes are equal
        """
        if not isinstance(col, Column):
            raise ValueError('Comparison `col` must be a Column object')

        attrs = ('name', 'units', 'dtype', 'format', 'description', 'meta')
        equal = all(getattr(self, x) == getattr(col, x) for x in attrs)

        return equal

    def pformat(self, max_lines=None, show_name=True, show_units=False):
        """Return a list of formatted string representation of column values.

        If no value of ``max_lines`` is supplied then the height of the screen
        terminal is used to set ``max_lines``.  If the terminal height cannot
        be determined then a default of ``astropy.table.MAX_LINES`` is used.
        If a negative value of ``max_lines`` is supplied then there is no line
        limit applied.

        Parameters
        ----------
        max_lines : int
            Maximum lines of output (header + data rows)

        show_name : bool
            Include column name (default=True)

        show_units : bool
            Include a header row for units (default=False)

        Returns
        -------
        lines : list
            List of lines with header and formatted column values

        """
        lines, n_header = _pformat_col(self, max_lines, show_name, show_units)
        return lines

    def pprint(self, max_lines=None, show_name=True, show_units=False):
        """Print a formatted string representation of column values.

        If no value of ``max_lines`` is supplied then the height of the screen
        terminal is used to set ``max_lines``.  If the terminal height cannot
        be determined then a default of ``astropy.table.MAX_LINES`` is used.
        If a negative value of ``max_lines`` is supplied then there is no line
        limit applied.

        Parameters
        ----------
        max_lines : int
            Maximum number of values in output

        show_name : bool
            Include column name (default=True)

        show_units : bool
            Include a header row for units (default=False)
        """
        lines, n_header = _pformat_col(self, max_lines, show_name, show_units)
        for i, line in enumerate(lines):
            if i < n_header:
                color_print(line, 'red')
            else:
                print line

    def more(self, max_lines=None, show_name=True, show_units=False):
        """Interactively browse column with a paging interface.

        Supported keys::

          f, <space> : forward one page
          b : back one page
          r : refresh same page
          n : next row
          p : previous row
          < : go to beginning
          > : go to end
          q : quit browsing
          h : print this help

        Parameters
        ----------
        max_lines : int
            Maximum number of lines in table output

        show_name : bool
            Include a header row for column names (default=True)

        show_units : bool
            Include a header row for units (default=False)

        """
        _more_tabcol(self, max_lines=max_lines, show_name=show_name,
                     show_units=show_units)

    def __str__(self):
        lines, n_header = _pformat_col(self)
        return '\n'.join(lines)


class Row(object):
    """A class to represent one row of a Table object.

    A Row object is returned when a Table object is indexed with an integer
    or when iterating over a table::

      >>> table = Table([(1, 2), (3, 4)], names=('a', 'b'))
      >>> row = table[1]
      >>> row
      <Row 1 of table
       values=(2, 4)
       dtype=[('a', '<i8'), ('b', '<i8')]>
      >>> row['a']
      2
      >>> row[1]
      4
    """

    def __init__(self, table, index):
        self._table = table
        self._index = index
        self._data = table._data[index]

    def __getitem__(self, item):
        return self.data[item]

    def __setitem__(self, item, val):
        self.data[item] = val

    def __eq__(self, other):
        return self.data == other

    def __ne__(self, other):
        return self.data != other

    def __array__(self, dtype=None):
        """Support converting Row to np.array via np.array(table).

        Coercion to a different dtype via np.array(table, dtype) is not
        supported and will raise a ValueError.
        """
        if dtype is not None:
            raise ValueError('Datatype coercion is not allowed')

        return np.array(self._data)

    def __len__(self):
        return len(self._data)

    @property
    def table(self):
        return self._table

    @property
    def index(self):
        return self._index

    @property
    def data(self):
        return self._data

    @property
    def meta(self):
        return self.table.meta

    @property
    def columns(self):
        return self.table.columns

    @property
    def colnames(self):
        return self.table.colnames

    @property
    def dtype(self):
        return self.data.dtype

    def __repr__(self):
        return "<Row {0} of table\n values={1}\n dtype={2}>".format(
            self.index, self.data, self.dtype)


class Table(object):
    """A class to represent tables of heterogeneous data.

    `Table` provides a class for heterogeneous tabular data, making use of a
    `numpy` structured array internally to store the data values.  A key
    enhancement provided by the `Table` class is the ability to easily modify
    the structure of the table by adding or removing columns, or adding new
    rows of data.  In addition table and column metadata are fully supported.

    `Table` differs from `NDData` by the assumption that the input data
    consists of columns of homogeneous data, where each column has a unique
    identifier and may contain additional metadata such as the data
    units, format, and description.

    Parameters
    ----------
    data : numpy ndarray, dict, list, or Table, optional
        Data to initialize table.
    names : list, optional
        Specify column names
    dtypes : list, optional
        Specify column data types
    meta : dict, optional
        Metadata associated with the table
    copy : boolean, optional
        Copy the input data (default=True).

    """

    def __init__(self, data=None, names=None, dtypes=None, meta=None,
                 copy=True):

        # Set up a placeholder empty table
        self._data = None
        self.columns = TableColumns()
        self.meta = OrderedDict() if meta is None else deepcopy(meta)

        # Must copy if dtypes are changing
        if not copy and dtypes is not None:
            raise ValueError('Cannot specify dtypes when copy=False')

        # Infer the type of the input data and set up the initialization
        # function, number of columns, and potentially the default col names

        default_names = None

        if isinstance(data, (list, tuple)):
            init_func = self._init_from_list
            n_cols = len(data)

        elif isinstance(data, np.ndarray):
            if data.dtype.names:
                init_func = self._init_from_ndarray  # _struct
                n_cols = len(data.dtype.names)
                default_names = data.dtype.names
            else:
                init_func = self._init_from_ndarray  # _homog
                n_cols = data.shape[1]

        elif isinstance(data, dict):
            init_func = self._init_from_dict
            n_cols = len(data.keys())
            default_names = data.keys()

        elif isinstance(data, Table):
            init_func = self._init_from_table
            n_cols = len(data.colnames)
            default_names = data.colnames

        elif data is None:
            if names is None:
                return  # Empty table
            else:
                init_func = self._init_from_list
                n_cols = len(names)
                data = [[]] * n_cols
        else:
            raise ValueError('Data type {0} not allowed to init Table'
                             .format(type(data)))

        # Set up defaults if names and/or dtypes are not specified.
        # A value of None means the actual value will be inferred
        # within the appropriate initialization routine, either from
        # existing specification or auto-generated.

        if names is None:
            names = default_names or [None] * n_cols
        if dtypes is None:
            dtypes = [None] * n_cols
        self._check_names_dtypes(names, dtypes, n_cols)

        # Finally do the real initialization
        init_func(data, names, dtypes, n_cols, copy)

    def __array__(self, dtype=None):
        """Support converting Table to np.array via np.array(table).

        Coercion to a different dtype via np.array(table, dtype) is not
        supported and will raise a ValueError.
        """
        if dtype is not None:
            raise ValueError('Datatype coercion is not allowed')

        # This limitation is because of the following unexpected result that
        # should have made a table copy while changing the column names.
        #
        # >>> d = astropy.table.Table([[1,2],[3,4]])
        # >>> np.array(d, dtype=[('a', 'i8'), ('b', 'i8')])
        # array([(0, 0), (0, 0)],
        #       dtype=[('a', '<i8'), ('b', '<i8')])

        return self._data

    def _check_names_dtypes(self, names, dtypes, n_cols):
        """Make sure that names and dtypes are boths iterable and have
        the same length as data.
        """
        for inp_list, inp_str in ((dtypes, 'dtypes'), (names, 'names')):
            if not isiterable(inp_list):
                raise ValueError('{0} must be a list or None'.format(inp_str))

        if len(names) != n_cols or len(dtypes) != n_cols:
            raise ValueError(
                'Arguments "names" and "dtypes" must match number of columns'
                .format(inp_str))

    def _init_from_list(self, data, names, dtypes, n_cols, copy):
        """Initialize table from a list of columns.  A column can be a
        Column object, np.ndarray, or any other iterable object.
        """
        if not copy:
            raise ValueError('Cannot use copy=False with a list data input')

        cols = []
        def_names = _auto_names(n_cols)
        for col, name, def_name, dtype in zip(data, names, def_names, dtypes):
            if isinstance(col, Column):
                col = Column((name or col.name), col, dtype=dtype)
            elif isinstance(col, np.ndarray) or isiterable(col):
                col = Column((name or def_name), col, dtype=dtype)
            else:
                raise ValueError('Elements in list initialization must be '
                                 'either Column or list-like')
            cols.append(col)

        self._init_from_cols(cols)

    def _init_from_ndarray(self, data, names, dtypes, n_cols, copy):
        """Initialize table from an ndarray structured array"""

        data_names = data.dtype.names or _auto_names(n_cols)
        struct = data.dtype.names is not None
        names = [name or data_names[i] for i, name in enumerate(names)]

        cols = ([data[name] for name in data_names] if struct else
                [data[:, i] for i in range(n_cols)])

        if copy:
            self._init_from_list(cols, names, dtypes, n_cols, copy)
        else:
            dtypes = [(name, col.dtype) for name, col in zip(names, cols)]
            self._data = data.view(dtypes).ravel()
            columns = TableColumns()
            for name in names:
                columns[name] = Column(name, self._data[name])
                columns[name].parent_table = self
            self.columns = columns

    def _init_from_dict(self, data, names, dtypes, n_cols, copy):
        """Initialize table from a dictionary of columns"""

        if not copy:
            raise ValueError('Cannot use copy=False with a dict data input')

        data_list = [data[name] for name in names]
        self._init_from_list(data_list, names, dtypes, n_cols, copy)

    def _init_from_table(self, data, names, dtypes, n_cols, copy):
        """Initialize table from an existing Table object """

        table = data  # data is really a Table, rename for clarity
        data_names = table.colnames
        self.meta = deepcopy(table.meta)
        cols = table.columns.values()

        if copy:
            self._init_from_list(cols, names, dtypes, n_cols, copy)
        else:
            names = [vals[0] or vals[1] for vals in zip(names, data_names)]
            dtypes = [(name, col.dtype) for name, col in zip(names, cols)]
            data = table._data.view(dtypes)

            self._update_table_from_cols(self, data, cols, names)

    def _init_from_cols(self, cols):
        """Initialize table from a list of Column objects"""

        lengths = set(len(col.data) for col in cols)
        if len(lengths) != 1:
            raise ValueError('Inconsistent data column lengths: {0}'
                             .format(lengths))

        names = [col.name for col in cols]
        dtypes = [col.descr for col in cols]
        data = np.empty(lengths.pop(), dtype=dtypes)
        for col in cols:
            data[col.name] = col.data

        self._update_table_from_cols(self, data, cols, names)

    def _new_from_slice(self, slice_):
        """Create a new table as a referenced slice from self."""

        table = Table()
        table.meta = deepcopy(self.meta)
        cols = self.columns.values()
        names = [col.name for col in cols]
        data = self._data[slice_]

        self._update_table_from_cols(table, data, cols, names)

        return table

    @staticmethod
    def _update_table_from_cols(table, data, cols, names):
        """Update the existing ``table`` so that it represents the given
        ``data`` (a structured ndarray) with ``cols`` and ``names``."""

        columns = TableColumns()
        table._data = data

        for name, col in zip(names, cols):
            newcol = col.copy(data=data[name], copy_data=False)
            newcol.name = name
            newcol.parent_table = table
            columns[name] = newcol
        table.columns = columns

    def __repr__(self):
        names = ("'{0}'".format(x) for x in self.colnames)
        s = "<Table rows={0} names=({1})>\n{2}".format(
            self.__len__(),  ','.join(names), repr(self._data))
        return s

    def __str__(self):
        lines, n_header = _pformat_table(self)
        return '\n'.join(lines)

    def pprint(self, max_lines=None, max_width=None, show_name=True,
               show_units=False):
        """Print a formatted string representation of the table.

        If no value of ``max_lines`` is supplied then the height of the screen
        terminal is used to set ``max_lines``.  If the terminal height cannot
        be determined then a default of ``astropy.table.pprint.MAX_LINES`` is
        used.  If a negative value of ``max_lines`` is supplied then there is
        no line limit applied.

        The Same applies for max_width except the default is
        ``astropy.table.pprint.MAX_WIDTH``.

        Parameters
        ----------
        max_lines : int
            Maximum number of lines in table output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_units : bool
            Include a header row for units (default=False)
        """

        lines, n_header = _pformat_table(self, max_lines, max_width, show_name,
                                         show_units)
        for i, line in enumerate(lines):
            if i < n_header:
                color_print(line, 'red')
            else:
                print line

    def pformat(self, max_lines=None, max_width=None, show_name=True,
                show_units=False):
        """Return a list of lines for the formatted string representation of
        the table.

        If no value of ``max_lines`` is supplied then the height of the screen
        terminal is used to set ``max_lines``.  If the terminal height cannot
        be determined then a default of ``astropy.table.pprint.MAX_LINES`` is
        used.  If a negative value of ``max_lines`` is supplied then there is
        no line limit applied.

        The Same applies for max_width except the default is
        ``astropy.table.pprint.MAX_WIDTH``.

        Parameters
        ----------
        max_lines : int or None
            Maximum number of rows to output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_units : bool
            Include a header row for units (default=False)

        Returns
        -------
        lines : list
            Formatted table as a list of strings
        """
        lines, n_header = _pformat_table(self, max_lines, max_width,
                                         show_name, show_units)
        return lines

    def more(self, max_lines=None, max_width=None, show_name=True,
               show_units=False):
        """Interactively browse table with a paging interface.

        Supported keys::

          f, <space> : forward one page
          b : back one page
          r : refresh same page
          n : next row
          p : previous row
          < : go to beginning
          > : go to end
          q : quit browsing
          h : print this help

        Parameters
        ----------
        max_lines : int
            Maximum number of lines in table output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_units : bool
            Include a header row for units (default=False)
        """
        _more_tabcol(self, max_lines, max_width, show_name,
                     show_units)

    def __getitem__(self, item):
        if isinstance(item, basestring):
            return self.columns[item]
        elif isinstance(item, int):
            return Row(self, item)
        elif isinstance(item, tuple):
            if any(x not in set(self.colnames) for x in item):
                raise ValueError('Table column slice must contain only valid '
                                 'column names')
            return Table([self[x] for x in item], meta=deepcopy(self.meta))
        elif (isinstance(item, slice) or isinstance(item, np.ndarray)
              or isinstance(item, list)):
            return self._new_from_slice(item)
        else:
            raise ValueError('Illegal type {0} for table item access'
                             .format(type(item)))

    def __setitem__(self, item, value):
        try:
            self._data[item] = value
        except (ValueError, KeyError, TypeError):
            raise KeyError("Column {0} does not exist".format(item))
        except:
            raise

    def __delitem__(self, item):
        if isinstance(item, basestring):
            self.remove_column(item)
        elif isinstance(item, tuple):
            self.remove_columns(item)

    def __iter__(self):
        self._iter_index = 0
        return self

    def __next__(self):
        """Python 3 iterator"""
        if self._iter_index < len(self._data):
            val = self[self._iter_index]
            self._iter_index += 1
            return val
        else:
            raise StopIteration

    if sys.version_info[0] < 3:  # pragma: py2
        next = __next__

    def field(self, item):
        """Return column[item] for recarray compatibility."""
        return self.columns[item]

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def colnames(self):
        return list(self.columns.keys())

    def keys(self):
        return list(self.columns.keys())

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data)

    def index_column(self, name):
        """
        Return the positional index of column ``name``.

        Parameters
        ----------
        name : str
            column name

        Returns
        -------
        index : int
            Positional index of column ``name``.
        """
        try:
            return self.colnames.index(name)
        except ValueError:
            raise ValueError("Column {0} does not exist".format(name))

    def add_column(self, col, index=None):
        """
        Add a new Column object ``col`` to the table.  If ``index``
        is supplied then insert column before ``index`` position
        in the list of columns, otherwise append column to the end
        of the list.

        Parameters
        ----------
        col : Column
            Column object to add.
        index : int or None
            Insert column before this position or at end (default)
        """
        if index is None:
            index = len(self.columns)
        self.add_columns([col], [index])

    def add_columns(self, cols, indexes=None):
        """
        Add a list of new Column objects ``cols`` to the table.  If a
        corresponding list of ``indexes`` is supplied then insert column before
        each ``index`` position in the *original* list of columns, otherwise
        append columns to the end of the list.

        Parameters
        ----------
        cols : list of Columns
            Column objects to add.
        indexes : list of ints or None
            Insert column before this position or at end (default)
        """
        if indexes is None:
            indexes = [len(self.columns)] * len(cols)
        elif len(indexes) != len(cols):
            raise ValueError('Number of indexes must match number of cols')

        if self._data is None:
            # No existing table data, init from cols
            newcols = cols
        else:
            newcols = list(self.columns.values())
            new_indexes = list(range(len(newcols) + 1))
            for col, index in zip(cols, indexes):
                i = new_indexes.index(index)
                new_indexes.insert(i, None)
                newcols.insert(i, col)

        self._init_from_cols(newcols)

    def remove_column(self, name):
        """
        Remove a column from the table.

        This can also be done with::

          del table[name]

        Parameters
        ----------
        name : str
            Name of column to remove
        """

        self.remove_columns([name])

    def remove_columns(self, names):
        '''
        Remove several columns from the table

        Parameters
        ----------
        names : list
            A list containing the names of the columns to remove
        '''

        for name in names:
            if name not in self.columns:
                raise KeyError("Column {0} does not exist".format(name))

        for name in names:
            self.columns.pop(name)

        self._data = _drop_fields(self._data, names)

    def keep_columns(self, names):
        '''
        Keep only the columns specified (remove the others)

        Parameters
        ----------
        names : list
            A list containing the names of the columns to keep. All other
            columns will be removed.
        '''

        if isinstance(names, basestring):
            names = [names]

        for name in names:
            if name not in self.columns:
                raise KeyError("Column {0} does not exist".format(name))

        remove = list(set(self.keys()) - set(names))

        self.remove_columns(remove)

    def rename_column(self, name, new_name):
        '''
        Rename a column.

        This can also be done directly with by setting the ``name`` attribute
        for a column::

          table[name].name = new_name

        Parameters
        ----------
        name : str
            The current name of the column.
        new_name : str
            The new name for the column
        '''

        if name not in self.keys():
            raise KeyError("Column {0} does not exist".format(name))

        self.columns[name].name = new_name

    def add_row(self, vals=None):
        """Add a new row to the end of the table.

        The ``vals`` argument can be:

        sequence (e.g. tuple or list)
            Column values in the same order as table columns.
        mapping (e.g. dict)
            Keys corresponding to column names.  Missing values will be
            filled with np.zeros for the column dtype.
        None
            All values filled with np.zeros for the column dtype.

        This method requires that the Table object "owns" the underlying array
        data.  In particular one cannot add a row to a Table that was
        initialized with copy=False from an existing array.

        Parameters
        ----------
        vals : tuple, list, dict or None
            Use the specified values in the new row
        """
        newlen = len(self._data) + 1
        self._data.resize((newlen,), refcheck=False)

        if isinstance(vals, collections.Mapping):
            row = self._data[-1]
            for name, val in vals.items():
                try:
                    row[name] = val
                except IndexError:
                    raise ValueError("No column {0} in table".format(name))

        elif isiterable(vals):
            if len(self.columns) != len(vals):
                raise ValueError('Mismatch between number of vals and columns')
            if not isinstance(vals, tuple):
                vals = tuple(vals)
            self._data[-1] = vals

        elif vals is not None:
            raise TypeError('Vals must be an iterable or mapping or None')

        # Add_row() probably corrupted the Column views of self._data.  Rebuild
        # self.columns.  Col.copy() takes an optional data reference that it
        # uses in the copy.
        cols = [c.copy(self._data[c.name]) for c in self.columns.values()]
        self.columns = TableColumns(cols)

    def sort(self, keys):
        '''
        Sort the table according to one or more keys. This operates
        on the existing table and does not return a new table.

        Parameters
        ----------
        keys : str or list of str
            The key(s) to order the table by
        '''
        if type(keys) is not list:
            keys = [keys]
        self._data.sort(order=keys)

    def reverse(self):
        '''
        Reverse the row order of table rows.  The table is reversed
        in place and there are no function arguments.
        '''
        self._data[:] = self._data[::-1].copy()
