# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six
from ..extern.six.moves import zip as izip
from ..extern.six.moves import range as xrange

import warnings
import re

from copy import deepcopy
from distutils import version

import numpy as np
from numpy import ma

from .. import log
from ..io import registry as io_registry
from ..units import Quantity
from ..utils import OrderedDict, isiterable, deprecated
from ..utils.console import color_print
from ..utils.exceptions import AstropyDeprecationWarning
from ..utils.metadata import MetaData
from . import groups
from .pprint import TableFormatter
from .column import BaseColumn, Column, MaskedColumn, _auto_names
from .row import Row
from .np_utils import fix_column_name


# Prior to Numpy 1.6.2, there was a bug (in Numpy) that caused
# sorting of structured arrays containing Unicode columns to
# silently fail.
_NUMPY_VERSION = version.LooseVersion(np.__version__)
_BROKEN_UNICODE_TABLE_SORT = _NUMPY_VERSION < version.LooseVersion('1.6.2')


__doctest_skip__ = ['Table.read', 'Table.write',
                    'Table.convert_bytestring_to_unicode',
                    'Table.convert_unicode_to_bytestring',
                    ]


class TableColumns(OrderedDict):
    """OrderedDict subclass for a set of columns.

    This class enhances item access to provide convenient access to columns
    by name or index, including slice access.  It also handles renaming
    of columns.

    The initialization argument ``cols`` can be a list of ``Column`` objects
    or any structure that is valid for initializing a Python dict.  This
    includes a dict, list of (key, val) tuples or [key, val] lists, etc.

    Parameters
    ----------
    cols : dict, list, tuple; optional
        Column objects as data structure that can init dict (see above)
    """

    def __init__(self, cols={}):
        if isinstance(cols, (list, tuple)):
            # check for Columns in the list
            cols = [((col.name, col) if hasattr(col, 'name') else col)
                    for col in cols]
        super(TableColumns, self).__init__(cols)

    def __getitem__(self, item):
        """Get items from a TableColumns object.
        ::

          tc = TableColumns(cols=[Column(name='a'), Column(name='b'), Column(name='c')])
          tc['a']  # Column('a')
          tc[1] # Column('b')
          tc['a', 'b'] # <TableColumns names=('a', 'b')>
          tc[1:3] # <TableColumns names=('b', 'c')>
        """
        if isinstance(item, six.string_types):
            return OrderedDict.__getitem__(self, item)
        elif isinstance(item, int):
            return self.values()[item]
        elif isinstance(item, tuple):
            return self.__class__([self[x] for x in item])
        elif isinstance(item, slice):
            return self.__class__([self[x] for x in list(self)[item]])
        else:
            raise IndexError('Illegal key or index value for {} object'
                             .format(self.__class__.__name__))

    def __repr__(self):
        names = ("'{0}'".format(x) for x in six.iterkeys(self))
        return "<{1} names=({0})>".format(",".join(names), self.__class__.__name__)

    def _rename_column(self, name, new_name):
        if new_name in self:
            raise KeyError("Column {0} already exists".format(new_name))

        mapper = {name: new_name}
        new_names = [mapper.get(name, name) for name in self]
        cols = list(six.itervalues(self))
        self.clear()
        self.update(list(izip(new_names, cols)))

    # Define keys and values for Python 2 and 3 source compatibility
    def keys(self):
        return list(OrderedDict.keys(self))

    def values(self):
        return list(OrderedDict.values(self))


class Table(object):
    """A class to represent tables of heterogeneous data.

    `Table` provides a class for heterogeneous tabular data, making use of a
    `numpy` structured array internally to store the data values.  A key
    enhancement provided by the `Table` class is the ability to easily modify
    the structure of the table by adding or removing columns, or adding new
    rows of data.  In addition table and column metadata are fully supported.

    `Table` differs from `~astropy.nddata.NDData` by the assumption that the
    input data consists of columns of homogeneous data, where each column
    has a unique identifier and may contain additional metadata such as the
    data unit, format, and description.

    Parameters
    ----------
    data : numpy ndarray, dict, list, or Table, optional
        Data to initialize table.
    masked : bool, optional
        Specify whether the table is masked.
    names : list, optional
        Specify column names
    dtype : list, optional
        Specify column data types
    meta : dict, optional
        Metadata associated with the table.
    copy : bool, optional
        Copy the input data (default=True).
    rows : numpy ndarray, list of lists, optional
        Row-oriented data for table instead of ``data`` argument
    """

    meta = MetaData()

    # Define class attributes for core container objects to allow for subclass
    # customization.
    Row = Row
    Column = Column
    MaskedColumn = MaskedColumn
    TableColumns = TableColumns
    TableFormatter = TableFormatter

    def __init__(self, data=None, masked=None, names=None, dtype=None,
                 meta=None, copy=True, rows=None):

        # Set up a placeholder empty table
        self._data = None
        self._set_masked(masked)
        self.columns = self.TableColumns()
        self.meta = meta
        self.formatter = self.TableFormatter()

        # Must copy if dtype are changing
        if not copy and dtype is not None:
            raise ValueError('Cannot specify dtype when copy=False')

        # Row-oriented input, e.g. list of lists or list of tuples, list of
        # dict, Row instance.  Set data to something that the subsequent code
        # will parse correctly.
        is_list_of_dict = False
        if rows is not None:
            if data is not None:
                raise ValueError('Cannot supply both `data` and `rows` values')
            if all(isinstance(row, dict) for row in rows):
                is_list_of_dict = True  # Avoid doing the all(...) test twice.
                data = rows
            elif isinstance(rows, self.Row):
                data = rows
            else:
                rec_data = np.rec.fromrecords(rows)
                data = [rec_data[name] for name in rec_data.dtype.names]

        # Infer the type of the input data and set up the initialization
        # function, number of columns, and potentially the default col names

        default_names = None

        if isinstance(data, self.Row):
            data = data._table[data._index:data._index + 1]

        if isinstance(data, (list, tuple)):
            init_func = self._init_from_list
            if data and (is_list_of_dict or all(isinstance(row, dict) for row in data)):
                n_cols = len(data[0])
            else:
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
            default_names = list(data)
            n_cols = len(default_names)

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

        # Set up defaults if names and/or dtype are not specified.
        # A value of None means the actual value will be inferred
        # within the appropriate initialization routine, either from
        # existing specification or auto-generated.

        if names is None:
            names = default_names or [None] * n_cols
        if dtype is None:
            dtype = [None] * n_cols

        # Numpy does not support Unicode column names on Python 2, or
        # bytes column names on Python 3, so fix them up now.
        names = [fix_column_name(name) for name in names]

        self._check_names_dtype(names, dtype, n_cols)

        # Finally do the real initialization
        init_func(data, names, dtype, n_cols, copy)

        # Whatever happens above, the masked property should be set to a boolean
        if type(self.masked) != bool:
            raise TypeError("masked property has not been set to True or False")

    def __getstate__(self):
        return (self.columns.values(), self.meta)

    def __setstate__(self, state):
        columns, meta = state
        self.__init__(columns, meta=meta)

    @property
    def mask(self):
        return self._data.mask if self.masked else None

    @mask.setter
    def mask(self, val):
        self._data.mask = val

    @property
    def _mask(self):
        """This is needed due to intricacies in numpy.ma, don't remove it."""
        return self._data.mask

    def filled(self, fill_value=None):
        """Return a copy of self, with masked values filled.

        If input ``fill_value`` supplied then that value is used for all
        masked entries in the table.  Otherwise the individual
        ``fill_value`` defined for each table column is used.

        Parameters
        ----------
        fill_value : str
            If supplied, this ``fill_value`` is used for all masked entries
            in the entire table.

        Returns
        -------
        filled_table : Table
            New table with masked values filled
        """
        if self.masked:
            data = [col.filled(fill_value) for col in six.itervalues(self.columns)]
        else:
            data = self
        return self.__class__(data, meta=deepcopy(self.meta))

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

        return self._data.data if self.masked else self._data

    def _rebuild_table_column_views(self):
        """
        Some table manipulations can corrupt the Column views of self._data.
        This function will cleanly rebuild the columns and self.columns.
        This is a slightly subtle operation, see comments.
        """
        cols = []
        for col in six.itervalues(self.columns):
            # First make a new column based on the name and the original
            # column.  This step is needed because the table manipulation
            # may have changed the table masking so that the original data
            # columns no longer correspond to self.ColumnClass.  This uses
            # data refs, not copies.
            newcol = self.ColumnClass(name=col.name, data=col)

            # Now use the copy() method to copy the column and its metadata,
            # but at the same time set the column data to a view of
            # self._data[col.name].  Somewhat confusingly in this case
            # copy() refers to copying the column attributes, but the data
            # are used by reference.
            newcol = newcol.copy(data=self._data[col.name])

            # Make column aware of the parent table
            newcol.parent_table = self

            cols.append(newcol)

        self.columns = self.TableColumns(cols)

    def _check_names_dtype(self, names, dtype, n_cols):
        """Make sure that names and dtype are boths iterable and have
        the same length as data.
        """
        for inp_list, inp_str in ((dtype, 'dtype'), (names, 'names')):
            if not isiterable(inp_list):
                raise ValueError('{0} must be a list or None'.format(inp_str))

        if len(names) != n_cols or len(dtype) != n_cols:
            raise ValueError(
                'Arguments "names" and "dtype" must match number of columns'
                .format(inp_str))

    def _set_masked_from_cols(self, cols):
        if self.masked is None:
            if any(isinstance(col, (MaskedColumn, ma.MaskedArray)) for col in cols):
                self._set_masked(True)
            else:
                self._set_masked(False)
        elif not self.masked:
            if any(np.any(col.mask) for col in cols if isinstance(col, (MaskedColumn, ma.MaskedArray))):
                self._set_masked(True)

    def _init_from_list(self, data, names, dtype, n_cols, copy):
        """Initialize table from a list of columns.  A column can be a
        Column object, np.ndarray, or any other iterable object.
        """
        if not copy:
            raise ValueError('Cannot use copy=False with a list data input')

        # Set self.masked appropriately, then get class to create column instances.
        self._set_masked_from_cols(data)

        cols = []
        def_names = _auto_names(n_cols)

        if data and all(isinstance(row, dict) for row in data):
            names_from_data = set()
            for row in data:
                names_from_data.update(row)

            cols = {}
            for name in names_from_data:
                cols[name] = []
                for i, row in enumerate(data):
                    try:
                        cols[name].append(row[name])
                    except KeyError:
                        raise ValueError('Row {0} has no value for column {1}'.format(i, name))
            if all(name is None for name in names):
                names = sorted(names_from_data)
            self._init_from_dict(cols, names, dtype, n_cols, copy)
            return

        for col, name, def_name, dtype in zip(data, names, def_names, dtype):
            if isinstance(col, (Column, MaskedColumn)):
                col = self.ColumnClass(name=(name or col.name), data=col, dtype=dtype)
            elif isinstance(col, np.ndarray) or isiterable(col):
                col = self.ColumnClass(name=(name or def_name), data=col, dtype=dtype)
            else:
                raise ValueError('Elements in list initialization must be '
                                 'either Column or list-like')
            cols.append(col)

        self._init_from_cols(cols)

    def _init_from_ndarray(self, data, names, dtype, n_cols, copy):
        """Initialize table from an ndarray structured array"""

        data_names = data.dtype.names or _auto_names(n_cols)
        struct = data.dtype.names is not None
        names = [name or data_names[i] for i, name in enumerate(names)]

        cols = ([data[name] for name in data_names] if struct else
                [data[:, i] for i in range(n_cols)])

        # Set self.masked appropriately, then get class to create column instances.
        self._set_masked_from_cols(cols)

        if copy:
            self._init_from_list(cols, names, dtype, n_cols, copy)
        else:
            dtype = [(name, col.dtype) for name, col in zip(names, cols)]
            self._data = data.view(dtype).ravel()
            columns = self.TableColumns()

            for name in names:
                columns[name] = self.ColumnClass(name=name, data=self._data[name])
                columns[name].parent_table = self
            self.columns = columns

    def _init_from_dict(self, data, names, dtype, n_cols, copy):
        """Initialize table from a dictionary of columns"""

        if not copy:
            raise ValueError('Cannot use copy=False with a dict data input')

        data_list = [data[name] for name in names]
        self._init_from_list(data_list, names, dtype, n_cols, copy)

    def _init_from_table(self, data, names, dtype, n_cols, copy):
        """Initialize table from an existing Table object """

        table = data  # data is really a Table, rename for clarity
        data_names = table.colnames
        self.meta.clear()
        self.meta.update(deepcopy(table.meta))
        cols = list(six.itervalues(table.columns))

        # Set self.masked appropriately from cols
        self._set_masked_from_cols(cols)

        if copy:
            self._init_from_list(cols, names, dtype, n_cols, copy)
        else:
            names = [vals[0] or vals[1] for vals in zip(names, data_names)]
            dtype = [(name, col.dtype) for name, col in zip(names, cols)]
            data = table._data.view(dtype)

            self._update_table_from_cols(self, data, cols, names)

    def _init_from_cols(self, cols):
        """Initialize table from a list of Column objects"""

        lengths = set(len(col.data) for col in cols)
        if len(lengths) != 1:
            raise ValueError('Inconsistent data column lengths: {0}'
                             .format(lengths))

        self._set_masked_from_cols(cols)
        cols = [self.ColumnClass(name=col.name, data=col) for col in cols]

        names = [col.name for col in cols]
        dtype = [col.descr for col in cols]
        empty_init = ma.empty if self.masked else np.empty
        data = empty_init(lengths.pop(), dtype=dtype)
        for col in cols:
            data[col.name] = col.data

        self._update_table_from_cols(self, data, cols, names)

    def _new_from_slice(self, slice_):
        """Create a new table as a referenced slice from self."""

        table = self.__class__(masked=self.masked)
        table.meta.clear()
        table.meta.update(deepcopy(self.meta))
        cols = list(six.itervalues(self.columns))
        names = [col.name for col in cols]
        data = self._data[slice_]

        self._update_table_from_cols(table, data, cols, names)

        return table

    @staticmethod
    def _update_table_from_cols(table, data, cols, names):
        """Update the existing ``table`` so that it represents the given
        ``data`` (a structured ndarray) with ``cols`` and ``names``."""

        columns = table.TableColumns()
        table._data = data

        for name, col in zip(names, cols):
            newcol = col.copy(data=data[name], copy_data=False)
            newcol.name = name
            newcol.parent_table = table
            columns[name] = newcol
        table.columns = columns

    def __repr__(self):
        names = ("'{0}'".format(x) for x in self.colnames)
        if any(col.unit for col in self.columns.values()):
            units = ("{0}".format(
                    col.unit if col.unit is None else '\''+str(col.unit)+'\'')
                    for col in self.columns.values())
            s = "<{3} rows={0} names=({1}) units=({4})>\n{2}".format(
                self.__len__(), ','.join(names), repr(self._data), self.__class__.__name__
                ,','.join(units))
        else:
            s = "<{3} rows={0} names=({1})>\n{2}".format(
                self.__len__(), ','.join(names), repr(self._data), self.__class__.__name__)
        return s

    def __unicode__(self):
        lines, n_header = self.formatter._pformat_table(self)
        return '\n'.join(lines)
    if six.PY3:
        __str__ = __unicode__

    def __bytes__(self):
        return six.text_type(self).encode('utf-8')
    if six.PY2:
        __str__ = __bytes__

    def pprint(self, max_lines=None, max_width=None, show_name=True,
               show_unit=None):
        """Print a formatted string representation of the table.

        If no value of ``max_lines`` is supplied then the height of the
        screen terminal is used to set ``max_lines``.  If the terminal
        height cannot be determined then the default is taken from the
        configuration item ``astropy.conf.max_lines``.  If a negative
        value of ``max_lines`` is supplied then there is no line limit
        applied.

        The same applies for max_width except the configuration item is
        ``astropy.conf.max_width``.

        Parameters
        ----------
        max_lines : int
            Maximum number of lines in table output

        max_width : int or `None`
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        """

        lines, n_header = self.formatter._pformat_table(self, max_lines, max_width, show_name,
                                                        show_unit)
        for i, line in enumerate(lines):
            if i < n_header:
                color_print(line, 'red')
            else:
                print(line)

    def show_in_browser(self,
                        css="table,th,td,tr,tbody {border: 1px solid black; border-collapse: collapse;}",
                        max_lines=5000,
                        jsviewer=False,
                        jskwargs={'use_local_files': True},
                        tableid=None,
                        browser='default'):
        """
        Render the table in HTML and show it in a web browser.

        Parameters
        ----------
        css : string
            A valid CSS string declaring the formatting for the table
        max_lines : int
            Maximum number of rows to export to the table (set low by default
            to avoid memory issues, since the browser view requires duplicating
            the table in memory).  A negative value of ``max_lines`` indicates
            no row limit
        jsviewer : bool
            If `True`, prepends some javascript headers so that the table is
            rendered as a https://datatables.net data table.  This allows
            in-browser searching & sorting.  See `JSViewer
            <http://www.jsviewer.com/>`_
        jskwargs : dict
            Passed to the `JSViewer`_ init.
        tableid : str or `None`
            An html ID tag for the table.  Default is "table{id}", where id is
            the unique integer id of the table object, id(self).
        browser : str
            Any legal browser name, e.g. ``'firefox'``, ``'chrome'``,
            ``'safari'`` (for mac, you may need to use ``'open -a
            "/Applications/Google Chrome.app" %s'`` for Chrome).  If
            ``'default'``, will use the system default browser.
        """

        import os
        import webbrowser
        import tempfile

        # We can't use NamedTemporaryFile here because it gets deleted as
        # soon as it gets garbage collected.

        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, 'table.html')

        with open(path, 'w') as tmp:

            if jsviewer:
                self.write(tmp, format='jsviewer', css=css, max_lines=max_lines,
                           jskwargs=jskwargs, table_id=tableid)
            else:
                self.write(tmp, format='html')

            if browser == 'default':
                webbrowser.open("file://" + path)
            else:
                webbrowser.get(browser).open("file://" + path)

    def pformat(self, max_lines=None, max_width=None, show_name=True,
                show_unit=None, html=False, tableid=None):
        """Return a list of lines for the formatted string representation of
        the table.

        If no value of ``max_lines`` is supplied then the height of the
        screen terminal is used to set ``max_lines``.  If the terminal
        height cannot be determined then the default is taken from the
        configuration item ``astropy.conf.max_lines``.  If a negative
        value of ``max_lines`` is supplied then there is no line limit
        applied.

        The same applies for ``max_width`` except the configuration item  is
        ``astropy.conf.max_width``.

        Parameters
        ----------
        max_lines : int or `None`
            Maximum number of rows to output

        max_width : int or `None`
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        html : bool
            Format the output as an HTML table (default=False)

        tableid : str or `None`
            An ID tag for the table; only used if html is set.  Default is
            "table{id}", where id is the unique integer id of the table object,
            id(self)

        Returns
        -------
        lines : list
            Formatted table as a list of strings

        """
        lines, n_header = self.formatter._pformat_table(self, max_lines, max_width,
                                                        show_name, show_unit, html,
                                                        tableid=tableid)
        return lines

    def more(self, max_lines=None, max_width=None, show_name=True,
             show_unit=None):
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

        max_width : int or `None`
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.
        """
        self.formatter._more_tabcol(self, max_lines, max_width, show_name,
                                    show_unit)

    def _repr_html_(self):
        # Since the user cannot provide input, need a sensible default
        tableid = 'table{id}'.format(id=id(self))
        lines = self.pformat(html=True, tableid=tableid, max_width=-1)
        return ''.join(lines)

    def __getitem__(self, item):
        if isinstance(item, six.string_types):
            return self.columns[item]
        elif isinstance(item, (int, np.integer)):
            return self.Row(self, item)
        elif isinstance(item, (tuple, list)) and all(isinstance(x, six.string_types)
                                                     for x in item):
            bad_names = [x for x in item if x not in self.colnames]
            if bad_names:
                raise ValueError('Slice name(s) {0} not valid column name(s)'
                                 .format(', '.join(bad_names)))
            out = self.__class__([self[x] for x in item], meta=deepcopy(self.meta))
            out._groups = groups.TableGroups(out, indices=self.groups._indices,
                                             keys=self.groups._keys)
            return out
        elif (isinstance(item, slice) or
              isinstance(item, np.ndarray) or
              isinstance(item, list) or
              isinstance(item, tuple) and all(isinstance(x, np.ndarray)
                                              for x in item)):
            # here for the many ways to give a slice; a tuple of ndarray
            # is produced by np.where, as in t[np.where(t['a'] > 2)]
            # For all, a new table is constructed with slice of all columns
            return self._new_from_slice(item)
        else:
            raise ValueError('Illegal type {0} for table item access'
                             .format(type(item)))

    def __setitem__(self, item, value):
        # If the item is a string then it must be the name of a column.
        # If that column doesn't already exist then create it now.
        if isinstance(item, six.string_types) and item not in self.colnames:
            NewColumn = self.MaskedColumn if self.masked else self.Column

            # Make sure value is an ndarray so we can get the dtype
            if not isinstance(value, np.ndarray):
                value = np.asarray(value)

            # Make new column and assign the value.  If the table currently
            # has no rows (len=0) of the value is already a Column then
            # define new column directly from value.  In the latter case
            # this allows for propagation of Column metadata.  Otherwise
            # define a new column with the right length and shape and then
            # set it from value.  This allows for broadcasting, e.g. t['a']
            # = 1.
            if isinstance(value, BaseColumn):
                new_column = value.copy(copy_data=False)
                new_column.name = item
            elif len(self) == 0:
                new_column = NewColumn(name=item, data=value)
            else:
                new_column = NewColumn(name=item, length=len(self), dtype=value.dtype,
                                       shape=value.shape[1:])
                new_column[:] = value

                if isinstance(value, Quantity):
                    new_column.unit = value.unit

            # Now add new column to the table
            self.add_column(new_column)

        elif isinstance(value, Row):
            # Value is another row
            self._data[item] = value.data
        else:
            # Otherwise just delegate to the numpy item setter.
            self._data[item] = value

    def __delitem__(self, item):
        if isinstance(item, six.string_types):
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

    if six.PY2:
        next = __next__

    def field(self, item):
        """Return column[item] for recarray compatibility."""
        return self.columns[item]

    @property
    def masked(self):
        return self._masked

    @masked.setter
    def masked(self, masked):
        raise Exception('Masked attribute is read-only (use t = Table(t, masked=True)'
                        ' to convert to a masked table)')

    def _set_masked(self, masked):
        """
        Set the table masked property.

        Parameters
        ----------
        masked : bool
            State of table masking (`True` or `False`)
        """
        if hasattr(self, '_masked'):
            # The only allowed change is from None to False or True, or False to True
            if self._masked is None and masked in [False, True]:
                self._masked = masked
            elif self._masked is False and masked is True:
                log.info("Upgrading Table to masked Table. Use Table.filled() to convert to unmasked table.")
                self._masked = masked
            elif self._masked is masked:
                raise Exception("Masked attribute is already set to {0}".format(masked))
            else:
                raise Exception("Cannot change masked attribute to {0} once it is set to {1}"
                                .format(masked, self._masked))
        else:
            if masked in [True, False, None]:
                self._masked = masked
            else:
                raise ValueError("masked should be one of True, False, None")
        if self._masked:
            self._column_class = self.MaskedColumn
        else:
            self._column_class = self.Column

    @property
    def ColumnClass(self):
        if self._column_class is None:
            return self.Column
        else:
            return self._column_class

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

    def create_mask(self):
        if isinstance(self._data, ma.MaskedArray):
            raise Exception("data array is already masked")
        else:
            self._data = ma.array(self._data)

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

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Get index of column 'b' of the table::

            >>> t.index_column('b')
            1
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
        index : int or `None`
            Insert column before this position or at end (default)

        Examples
        --------
        Create a table with two columns 'a' and 'b'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> print(t)
             a   b
            --- ---
              1 0.1
              2 0.2
              3 0.3

        Create a third column 'c' and append it to the end of the table::

            >>> col_c = Column(name='c', data=['x', 'y', 'z'])
            >>> t.add_column(col_c)
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Add column 'd' at position 1. Note that the column is inserted
        before the given index::

            >>> col_d = Column(name='d', data=['a', 'b', 'c'])
            >>> t.add_column(col_d, 1)
            >>> print(t)
             a   d   b   c
            --- --- --- ---
              1   a 0.1   x
              2   b 0.2   y
              3   c 0.3   z

        To add several columns use add_columns.
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
        indexes : list of ints or `None`
            Insert column before this position or at end (default)

        Examples
        --------
        Create a table with two columns 'a' and 'b'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> print(t)
             a   b
            --- ---
              1 0.1
              2 0.2
              3 0.3

        Create column 'c' and 'd' and append them to the end of the table::

            >>> col_c = Column(name='c', data=['x', 'y', 'z'])
            >>> col_d = Column(name='d', data=['u', 'v', 'w'])
            >>> t.add_columns([col_c, col_d])
            >>> print(t)
             a   b   c   d
            --- --- --- ---
              1 0.1   x   u
              2 0.2   y   v
              3 0.3   z   w

        Add column 'c' at position 0 and column 'd' at position 1. Note that
        the columns are inserted before the given position::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> col_c = Column(name='c', data=['x', 'y', 'z'])
            >>> col_d = Column(name='d', data=['u', 'v', 'w'])
            >>> t.add_columns([col_c, col_d], [0, 1])
            >>> print(t)
             c   a   d   b
            --- --- --- ---
              x   1   u 0.1
              y   2   v 0.2
              z   3   w 0.3
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

    def remove_row(self, index):
        """
        Remove a row from the table.

        Parameters
        ----------
        index : int
            Index of row to remove

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Remove row 1 from the table::

            >>> t.remove_row(1)
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              3 0.3   z

        To remove several rows at the same time use remove_rows.
        """
        # check the index against the types that work with np.delete
        if not isinstance(index, (six.integer_types, np.integer)):
            raise TypeError("Row index must be an integer")
        self.remove_rows(index)

    def remove_rows(self, row_specifier):
        """
        Remove rows from the table.

        Parameters
        ----------
        row_specifier : slice, int, or array of ints
            Specification for rows to remove

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Remove rows 0 and 2 from the table::

            >>> t.remove_rows([0, 2])
            >>> print(t)
             a   b   c
            --- --- ---
              2 0.2   y


        Note that there are no warnings if the slice operator extends
        outside the data::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> t.remove_rows(slice(10, 20, 1))
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z
        """
        try:
            table = np.delete(self._data, row_specifier, axis=0)
        except (ValueError, IndexError):
            # Numpy <= 1.7 raises ValueError while Numpy >= 1.8 raises IndexError
            raise IndexError('Removing row(s) {0} from table with {1} rows failed'
                             .format(row_specifier, len(self._data)))
        self._data = table

        # after updating the row data, the column views will be out of date
        # and should be updated:
        self._rebuild_table_column_views()

        # Revert groups to default (ungrouped) state
        if hasattr(self, '_groups'):
            del self._groups

    def remove_column(self, name):
        """
        Remove a column from the table.

        This can also be done with::

          del table[name]

        Parameters
        ----------
        name : str
            Name of column to remove

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Remove column 'b' from the table::

            >>> t.remove_column('b')
            >>> print(t)
             a   c
            --- ---
              1   x
              2   y
              3   z

        To remove several columns at the same time use remove_columns.
        """

        self.remove_columns([name])

    def remove_columns(self, names):
        '''
        Remove several columns from the table.

        Parameters
        ----------
        names : list
            A list containing the names of the columns to remove

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...     names=('a', 'b', 'c'))
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Remove columns 'b' and 'c' from the table::

            >>> t.remove_columns(['b', 'c'])
            >>> print(t)
             a
            ---
              1
              2
              3

        Specifying only a single column also works. Remove column 'b' from the table::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...     names=('a', 'b', 'c'))
            >>> t.remove_columns('b')
            >>> print(t)
             a   c
            --- ---
              1   x
              2   y
              3   z

        This gives the same as using remove_column.
        '''
        if isinstance(names, six.string_types):
            names = [names]

        for name in names:
            if name not in self.columns:
                raise KeyError("Column {0} does not exist".format(name))

        for name in names:
            self.columns.pop(name)

        newdtype = [(name, self._data.dtype[name]) for name in self._data.dtype.names
                    if name not in names]
        newdtype = np.dtype(newdtype)

        if newdtype:
            if self.masked:
                table = np.ma.empty(self._data.shape, dtype=newdtype)
            else:
                table = np.empty(self._data.shape, dtype=newdtype)

            for field in newdtype.fields:
                table[field] = self._data[field]
                if self.masked:
                    table[field].fill_value = self._data[field].fill_value
        else:
            table = None

        self._data = table

    def _convert_string_dtype(self, in_kind, out_kind, python3_only):
        """
        Convert string-like columns to/from bytestring and unicode (internal only).

        Parameters
        ----------
        in_kind : str
            Input dtype.kind
        out_kind : str
            Output dtype.kind
        python3_only : bool
            Only do this operation for Python 3
        """
        if python3_only and not six.PY3:
            return

        # If there are no `in_kind` columns then do nothing
        cols = self.columns.values()
        if not any(col.dtype.kind == in_kind for col in cols):
            return

        newcols = []
        for col in cols:
            if col.dtype.kind == in_kind:
                newdtype = re.sub(in_kind, out_kind, col.dtype.str)
                newcol = col.__class__(col, dtype=newdtype)
            else:
                newcol = col
            newcols.append(newcol)

        self._init_from_cols(newcols)

    def convert_bytestring_to_unicode(self, python3_only=False):
        """
        Convert bytestring columns (dtype.kind='S') to unicode (dtype.kind='U') assuming
        ASCII encoding.

        Internally this changes string columns to represent each character in the string
        with a 4-byte UCS-4 equivalent, so it is inefficient for memory but allows Python
        3 scripts to manipulate string arrays with natural syntax.

        The ``python3_only`` parameter is provided as a convenience so that code can
        be written in a Python 2 / 3 compatible way::

          >>> t = Table.read('my_data.fits')
          >>> t.convert_bytestring_to_unicode(python3_only=True)

        Parameters
        ----------
        python3_only : bool
            Only do this operation for Python 3
        """
        self._convert_string_dtype('S', 'U', python3_only)

    def convert_unicode_to_bytestring(self, python3_only=False):
        """
        Convert ASCII-only unicode columns (dtype.kind='U') to bytestring (dtype.kind='S').

        When exporting a unicode string array to a file in Python 3, it may be desirable
        to encode unicode columns as bytestrings.  This routine takes advantage of numpy
        automated conversion which works for strings that are pure ASCII.

        The ``python3_only`` parameter is provided as a convenience so that code can
        be written in a Python 2 / 3 compatible way::

          >>> t.convert_unicode_to_bytestring(python3_only=True)
          >>> t.write('my_data.fits')

        Parameters
        ----------
        python3_only : bool
            Only do this operation for Python 3
        """
        self._convert_string_dtype('U', 'S', python3_only)

    def keep_columns(self, names):
        '''
        Keep only the columns specified (remove the others).

        Parameters
        ----------
        names : list
            A list containing the names of the columns to keep. All other
            columns will be removed.

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3],[0.1, 0.2, 0.3],['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Specifying only a single column name keeps only this column.
        Keep only column 'a' of the table::

            >>> t.keep_columns('a')
            >>> print(t)
             a
            ---
              1
              2
              3

        Specifying a list of column names is keeps is also possible.
        Keep columns 'a' and 'c' of the table::

            >>> t = Table([[1, 2, 3],[0.1, 0.2, 0.3],['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> t.keep_columns(['a', 'c'])
            >>> print(t)
             a   c
            --- ---
              1   x
              2   y
              3   z
        '''

        if isinstance(names, six.string_types):
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

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1,2],[3,4],[5,6]], names=('a','b','c'))
            >>> print(t)
             a   b   c
            --- --- ---
              1   3   5
              2   4   6

        Renaming column 'a' to 'aa'::

            >>> t.rename_column('a' , 'aa')
            >>> print(t)
             aa  b   c
            --- --- ---
              1   3   5
              2   4   6
        '''

        if name not in self.keys():
            raise KeyError("Column {0} does not exist".format(name))

        self.columns[name].name = new_name

    def add_row(self, vals=None, mask=None):
        """Add a new row to the end of the table.

        The ``vals`` argument can be:

        sequence (e.g. tuple or list)
            Column values in the same order as table columns.
        mapping (e.g. dict)
            Keys corresponding to column names.  Missing values will be
            filled with np.zeros for the column dtype.
        `None`
            All values filled with np.zeros for the column dtype.

        This method requires that the Table object "owns" the underlying array
        data.  In particular one cannot add a row to a Table that was
        initialized with copy=False from an existing array.

        The ``mask`` attribute should give (if desired) the mask for the
        values. The type of the mask should match that of the values, i.e. if
        ``vals`` is an iterable, then ``mask`` should also be an iterable
        with the same length, and if ``vals`` is a mapping, then ``mask``
        should be a dictionary.

        Parameters
        ----------
        vals : tuple, list, dict or `None`
            Use the specified values in the new row
        mask : tuple, list, dict or `None`
            Use the specified mask values in the new row

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

           >>> t = Table([[1,2],[4,5],[7,8]], names=('a','b','c'))
           >>> print(t)
            a   b   c
           --- --- ---
             1   4   7
             2   5   8

        Adding a new row with entries '3' in 'a', '6' in 'b' and '9' in 'c'::

           >>> t.add_row([3,6,9])
           >>> print(t)
             a   b   c
             --- --- ---
             1   4   7
             2   5   8
             3   6   9
        """

        def _is_mapping(obj):
            """Minimal checker for mapping (dict-like) interface for obj"""
            attrs = ('__getitem__', '__len__', '__iter__', 'keys', 'values', 'items')
            return all(hasattr(obj, attr) for attr in attrs)

        newlen = len(self._data) + 1

        if vals is None:
            vals = np.zeros(1, dtype=self._data.dtype)[0]

        if mask is not None and not self.masked:
            self._set_masked(True)

        # Create a table with one row to test the operation on
        test_data = (ma.zeros if self.masked else np.zeros)(1, dtype=self._data.dtype)

        if _is_mapping(vals):

            if mask is not None and not _is_mapping(mask):
                raise TypeError("Mismatch between type of vals and mask")

            # Now check that the mask is specified for the same keys as the
            # values, otherwise things get really confusing.
            if mask is not None and set(vals.keys()) != set(mask.keys()):
                raise ValueError('keys in mask should match keys in vals')

            if self.masked:
                # We set the mask to True regardless of whether a mask value
                # is specified or not - that is, any cell where a new row
                # value is not specified should be treated as missing.
                test_data.mask[-1] = (True,) * len(test_data.dtype)

            # First we copy the values
            for name, val in six.iteritems(vals):
                try:
                    test_data[name][-1] = val
                except IndexError:
                    raise ValueError("No column {0} in table".format(name))
                if mask:
                    test_data[name].mask[-1] = mask[name]

        elif isiterable(vals):

            if mask is not None and (not isiterable(mask) or _is_mapping(mask)):
                raise TypeError("Mismatch between type of vals and mask")

            if len(self.columns) != len(vals):
                raise ValueError('Mismatch between number of vals and columns')

            if not isinstance(vals, tuple):
                vals = tuple(vals)

            test_data[-1] = vals

            if mask is not None:

                if len(self.columns) != len(mask):
                    raise ValueError('Mismatch between number of masks and columns')

                if not isinstance(mask, tuple):
                    mask = tuple(mask)

                test_data.mask[-1] = mask

        else:
            raise TypeError('Vals must be an iterable or mapping or None')

        # If no errors have been raised, then the table can be resized
        if self.masked:
            if newlen == 1:
                self._data = ma.empty(1, dtype=self._data.dtype)
            else:
                self._data = ma.resize(self._data, (newlen,))
        else:
            self._data.resize((newlen,), refcheck=False)

        # Assign the new row
        self._data[-1:] = test_data

        self._rebuild_table_column_views()

        # Revert groups to default (ungrouped) state
        if hasattr(self, '_groups'):
            del self._groups

    def argsort(self, keys=None, kind=None):
        """
        Return the indices which would sort the table according to one or
        more key columns.  This simply calls the `numpy.argsort` function on
        the table with the ``order`` parameter set to ``keys``.

        Parameters
        ----------
        keys : str or list of str
            The column name(s) to order the table by
        kind : {'quicksort', 'mergesort', 'heapsort'}, optional
            Sorting algorithm.

        Returns
        -------
        index_array : ndarray, int
            Array of indices that sorts the table by the specified key
            column(s).
        """
        if isinstance(keys, six.string_types):
            keys = [keys]
        kwargs = {}
        if keys:
            kwargs['order'] = keys
        if kind:
            kwargs['kind'] = kind

        data = self._data

        if _BROKEN_UNICODE_TABLE_SORT and keys is not None and any(
                data.dtype[i].kind == 'U' for i in xrange(len(data.dtype))):
            return np.lexsort([data[key] for key in keys[::-1]])
        else:
            return data.argsort(**kwargs)

    def sort(self, keys):
        '''
        Sort the table according to one or more keys. This operates
        on the existing table and does not return a new table.

        Parameters
        ----------
        keys : str or list of str
            The key(s) to order the table by

        Examples
        --------
        Create a table with 3 columns::

            >>> t = Table([['Max', 'Jo', 'John'], ['Miller','Miller','Jackson'],
            ...         [12,15,18]], names=('firstname','name','tel'))
            >>> print(t)
            firstname   name  tel
            --------- ------- ---
                  Max  Miller  12
                   Jo  Miller  15
                 John Jackson  18

        Sorting according to standard sorting rules, first 'name' then 'firstname'::

            >>> t.sort(['name','firstname'])
            >>> print(t)
            firstname   name  tel
            --------- ------- ---
                 John Jackson  18
                   Jo  Miller  15
                  Max  Miller  12
        '''
        if type(keys) is not list:
            keys = [keys]

        data = self._data

        if _BROKEN_UNICODE_TABLE_SORT and any(
                data.dtype[i].kind == 'U' for i in xrange(len(data.dtype))):
            # Use an alternate sort implementation that uses argsort
            ordering = self.argsort(keys=keys)
            data[:] = data[ordering]
        else:
            data.sort(order=keys)

        self._rebuild_table_column_views()

    def reverse(self):
        '''
        Reverse the row order of table rows.  The table is reversed
        in place and there are no function arguments.

        Examples
        --------
        Create a table with three columns::

            >>> t = Table([['Max', 'Jo', 'John'], ['Miller','Miller','Jackson'],
            ...         [12,15,18]], names=('firstname','name','tel'))
            >>> print(t)
            firstname   name  tel
            --------- ------- ---
                  Max  Miller  12
                   Jo  Miller  15
                 John Jackson  18

        Reversing order::

            >>> t.reverse()
            >>> print(t)
            firstname   name  tel
            --------- ------- ---
                 John Jackson  18
                   Jo  Miller  15
                  Max  Miller  12
        '''
        self._data[:] = self._data[::-1].copy()
        self._rebuild_table_column_views()

    @classmethod
    def read(cls, *args, **kwargs):
        """
        Read and parse a data table and return as a Table.

        This function provides the Table interface to the astropy unified I/O
        layer.  This allows easily reading a file in many supported data formats
        using syntax such as::

          >>> from astropy.table import Table
          >>> dat = Table.read('table.dat', format='ascii')
          >>> events = Table.read('events.fits', format='fits')

        The arguments and keywords (other than ``format``) provided to this function are
        passed through to the underlying data reader (e.g. `~astropy.io.ascii.read`).
        """
        return io_registry.read(cls, *args, **kwargs)

    def write(self, *args, **kwargs):
        """
        Write this Table object out in the specified format.

        This function provides the Table interface to the astropy unified I/O
        layer.  This allows easily writing a file in many supported data formats
        using syntax such as::

          >>> from astropy.table import Table
          >>> dat = Table([[1, 2], [3, 4]], names=('a', 'b'))
          >>> dat.write('table.dat', format='ascii')

        The arguments and keywords (other than ``format``) provided to this function are
        passed through to the underlying data reader (e.g. `~astropy.io.ascii.write`).
        """
        io_registry.write(self, *args, **kwargs)

    def copy(self, copy_data=True):
        '''
        Return a copy of the table.


        Parameters
        ----------
        copy_data : bool
            If `True` (the default), copy the underlying data array.
            Otherwise, use the same data array
        '''
        out = self.__class__(self, copy=copy_data)

        # If the current table is grouped then do the same in the copy
        if hasattr(self, '_groups'):
            out._groups = groups.TableGroups(out, indices=self._groups._indices,
                                             keys=self._groups._keys)
        return out

    def __deepcopy__(self, memo=None):
        return self.copy(True)

    def __copy__(self):
        return self.copy(False)

    def __lt__(self, other):
        if six.PY3:
            return super(Table, self).__lt__(other)
        elif six.PY2:
            raise TypeError("unorderable types: Table() < {0}".
                            format(str(type(other))))

    def __gt__(self, other):
        if six.PY3:
            return super(Table, self).__gt__(other)
        elif six.PY2:
            raise TypeError("unorderable types: Table() > {0}".
                            format(str(type(other))))

    def __le__(self, other):
        if six.PY3:
            return super(Table, self).__le__(other)
        elif six.PY2:
            raise TypeError("unorderable types: Table() <= {0}".
                            format(str(type(other))))

    def __ge__(self, other):
        if six.PY3:
            return super(Table, self).__ge__(other)
        else:
            raise TypeError("unorderable types: Table() >= {0}".
                            format(str(type(other))))

    def __eq__(self, other):

        if isinstance(other, Table):
            other = other._data

        if self.masked:
            if isinstance(other, np.ma.MaskedArray):
                result = self._data == other
            else:
                # If mask is True, then by definition the row doesn't match
                # because the other array is not masked.
                false_mask = np.zeros(1, dtype=[(n, bool) for n in self.dtype.names])
                result = (self._data.data == other) & (self.mask == false_mask)
        else:
            if isinstance(other, np.ma.MaskedArray):
                # If mask is True, then by definition the row doesn't match
                # because the other array is not masked.
                false_mask = np.zeros(1, dtype=[(n, bool) for n in other.dtype.names])
                result = (self._data == other.data) & (other.mask == false_mask)
            else:
                result = self._data == other

        return result

    def __ne__(self, other):
        return ~self.__eq__(other)

    @property
    def groups(self):
        if not hasattr(self, '_groups'):
            self._groups = groups.TableGroups(self)
        return self._groups

    def group_by(self, keys):
        """
        Group this table by the specified ``keys``

        This effectively splits the table into groups which correspond to
        unique values of the ``keys`` grouping object.  The output is a new
        `TableGroups` which contains a copy of this table but sorted by row
        according to ``keys``.

        The ``keys`` input to `group_by` can be specified in different ways:

          - String or list of strings corresponding to table column name(s)
          - Numpy array (homogeneous or structured) with same length as this table
          - `Table` with same length as this table

        Parameters
        ----------
        keys : str, list of str, numpy array, or `Table`
            Key grouping object

        Returns
        -------
        out : `Table`
            New table with groups set
        """
        return groups.table_group_by(self, keys)
