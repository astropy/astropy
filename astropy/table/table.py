# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six
from ..extern.six.moves import zip as izip
from ..extern.six.moves import range as xrange

import re

from copy import deepcopy

import numpy as np
from numpy import ma

from .. import log
from ..io import registry as io_registry
from ..units import Quantity
from ..utils import OrderedDict, isiterable, deprecated, minversion
from ..utils.console import color_print
from ..utils.metadata import MetaData
from . import groups
from .pprint import TableFormatter
from .column import (BaseColumn, Column, MaskedColumn, _auto_names, FalseArray,
                     col_getattr, col_setattr, col_copy, _col_update_attrs_from)
from .row import Row
from .np_utils import fix_column_name, recarray_fromrecords


# Prior to Numpy 1.6.2, there was a bug (in Numpy) that caused
# sorting of structured arrays containing Unicode columns to
# silently fail.
_BROKEN_UNICODE_TABLE_SORT = not minversion(np, '1.6.2')


__doctest_skip__ = ['Table.read', 'Table.write',
                    'Table.convert_bytestring_to_unicode',
                    'Table.convert_unicode_to_bytestring',
                    ]


def descr(col):
    """Array-interface compliant full description of a column.

    This returns a 3-tuple (name, type, shape) that can always be
    used in a structured array dtype definition.
    """
    col_dtype_str = col.dtype.str if hasattr(col, 'dtype') else 'O'
    col_shape = col.shape[1:] if hasattr(col, 'shape') else ()
    return (col_getattr(col, 'name'), col_dtype_str, col_shape)


def is_mixin_class(obj):
    """
    Abstraction to determine if ``obj`` should be used as a mixin column
    when input to a table.  This function does not apply to ``Quantity``.

    Parameters
    ----------
    obj : object
        Object being queried for mixin compatibility
    """
    return hasattr(obj, '_astropy_column_attrs')


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
            # `cols` should be a list of two-tuples, but it is allowed to have
            # columns (BaseColumn or mixins) in the list.
            newcols = []
            for col in cols:
                if (isinstance(col, (BaseColumn, Quantity))
                        or is_mixin_class(col)):
                    newcols.append((col_getattr(col, 'name'), col))
                else:
                    newcols.append(col)
            cols = newcols
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

    @property
    @deprecated('0.4', alternative=':attr:`Table.as_array`')
    def _data(self):
        """
        Return a new copy of the table in the form of a structured np.ndarray or
        np.ma.MaskedArray object (as appropriate).

        Prior to version 1.0 of astropy this private property was a modifiable
        view of the table data, but since 1.0 it is a copy.
        """
        return self.as_array()

    def as_array(self):
        """
        Return a new copy of the table in the form of a structured np.ndarray or
        np.ma.MaskedArray object (as appropriate).

        Returns
        -------
        table_array : np.ndarray (unmasked) or np.ma.MaskedArray (masked)
            Copy of table as a numpy structured array
        """
        if len(self.columns) == 0:
            return None

        cols = self.columns.values()
        dtype = [descr(col) for col in cols]
        empty_init = ma.empty if self.masked else np.empty
        data = empty_init(len(self), dtype=dtype)
        for col in cols:
            data[col_getattr(col, 'name')] = col

        return data

    def __init__(self, data=None, masked=None, names=None, dtype=None,
                 meta=None, copy=True, rows=None):

        # Set up a placeholder empty table
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
                rec_data = recarray_fromrecords(rows)
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
        # Dynamic view of available masks
        if self.masked:
            return Table([col.mask for col in self.columns.values()],
                         names=self.colnames, copy=False)
        else:
            return None

    @mask.setter
    def mask(self, val):
        self.mask[:] = val

    @property
    def _mask(self):
        """This is needed so that comparison of a masked Table and a
        MaskedArray works.  The requirement comes from numpy.ma.core
        so don't remove this property."""
        return self.as_array().mask

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

        return self.as_array().data if self.masked else self.as_array()

    def _check_names_dtype(self, names, dtype, n_cols):
        """Make sure that names and dtype are both iterable and have
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

    def _init_from_list_of_dicts(self, data, names, dtype, n_cols, copy):
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

    def _init_from_list(self, data, names, dtype, n_cols, copy):
        """Initialize table from a list of columns.  A column can be a
        Column object, np.ndarray, mixin, or any other iterable object.
        """
        if data and all(isinstance(row, dict) for row in data):
            self._init_from_list_of_dicts(data, names, dtype, n_cols, copy)
            return

        # Set self.masked appropriately, then get class to create column instances.
        self._set_masked_from_cols(data)

        cols = []
        def_names = _auto_names(n_cols)

        for col, name, def_name, dtype in zip(data, names, def_names, dtype):
            if isinstance(col, (Column, MaskedColumn)):
                col = self.ColumnClass(name=(name or col_getattr(col, 'name')),
                                       data=col, dtype=dtype,
                                       copy=copy)
            elif self._is_mixin_column(col):
                # Copy the mixin column attributes if they exist since the copy below
                # may not get this attribute.
                if copy:
                    col = col_copy(col)

                col_setattr(col, 'name', name or col_getattr(col, 'name') or def_name)
            elif isinstance(col, np.ndarray) or isiterable(col):
                col = self.ColumnClass(name=(name or def_name), data=col, dtype=dtype,
                                       copy=copy)
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
            dtype = [(name, col.dtype, col.shape[1:]) for name, col in zip(names, cols)]
            newdata = data.view(dtype).ravel()
            columns = self.TableColumns()

            for name in names:
                columns[name] = self.ColumnClass(name=name, data=newdata[name])
                col_setattr(columns[name], 'parent_table', self)
            self.columns = columns

    def _init_from_dict(self, data, names, dtype, n_cols, copy):
        """Initialize table from a dictionary of columns"""

        # TODO: is this restriction still needed with no ndarray?
        if not copy:
            raise ValueError('Cannot use copy=False with a dict data input')

        data_list = [data[name] for name in names]
        self._init_from_list(data_list, names, dtype, n_cols, copy)

    def _init_from_table(self, data, names, dtype, n_cols, copy):
        """Initialize table from an existing Table object """

        table = data  # data is really a Table, rename for clarity
        self.meta.clear()
        self.meta.update(deepcopy(table.meta))
        cols = list(table.columns.values())

        self._init_from_list(cols, names, dtype, n_cols, copy)

    def _init_from_cols(self, cols):
        """Initialize table from a list of Column objects"""

        lengths = set(len(col) for col in cols)
        if len(lengths) != 1:
            raise ValueError('Inconsistent data column lengths: {0}'
                             .format(lengths))

        # Set the table masking
        self._set_masked_from_cols(cols)

        # Make sure that all Column-based objects have class self.ColumnClass
        newcols = []
        for col in cols:
            # Convert any Columns with units to Quantity for a QTable
            if (isinstance(self, QTable) and isinstance(col, Column)
                    and getattr(col, 'unit', None) is not None):

                qcol = Quantity(col, unit=col.unit, copy=False)
                _col_update_attrs_from(qcol, col, exclude_attrs=['unit', 'dtype', 'parent_table'])

                newcols.append(qcol)
                continue

            if isinstance(col, Column) and not col.__class__ is self.ColumnClass:
                col = self.ColumnClass(col)  # copy attributes and reference data
            newcols.append(col)

        self._update_table_from_cols(self, newcols)

    def _new_from_slice(self, slice_):
        """Create a new table as a referenced slice from self."""

        table = self.__class__(masked=self.masked)
        table.meta.clear()
        table.meta.update(deepcopy(self.meta))
        cols = self.columns.values()
        names = [col_getattr(col, 'name') for col in cols]
        newcols = [col[slice_] for col in cols]

        # Mixin column classes are not responsible for copying column attributes
        # for item/slicing operations.  Do this here in table.
        for name, col, newcol in zip(names, cols, newcols):
            if is_mixin_class(col):
                _col_update_attrs_from(newcol, col, exclude_attrs=['parent_table'])

        self._update_table_from_cols(table, newcols)

        return table

    @staticmethod
    def _update_table_from_cols(table, cols):
        """Update the existing ``table`` so that it represents the given
        ``data`` (a structured ndarray) with ``cols`` and ``names``."""

        colnames = set(col_getattr(col, 'name') for col in cols)
        if None in colnames:
            raise TypeError('Cannot have None for column name')
        if len(colnames) != len(cols):
            raise ValueError('Duplicate column names')

        columns = table.TableColumns((col_getattr(col, 'name'), col) for col in cols)

        for col in cols:
            col_setattr(col, 'parent_table', table)
            if table.masked and not hasattr(col, 'mask'):
                col.mask = FalseArray(col.shape)

        table.columns = columns

    def _base_repr_(self, html=False):
        descr_vals = [self.__class__.__name__]
        for attr, val in (('masked', self.masked),
                          ('length', len(self))):
            descr_vals.append('{0}={1}'.format(attr, repr(val)))

        descr = '<' + ' '.join(descr_vals) + '>\n'

        if html:
            from ..utils.xml.writer import xml_escape
            descr = xml_escape(descr)

        tableid = 'table{id}'.format(id=id(self))
        data_lines, outs = self.formatter._pformat_table(self,
            tableid=tableid, html=html, max_width=(-1 if html else None),
            show_name=True, show_unit=None, show_dtype=True)

        out = descr + '\n'.join(data_lines)
        if six.PY2 and isinstance(out, six.text_type):
            out = out.encode('utf-8')

        return out

    def _repr_html_(self):
        return self._base_repr_(html=True)

    def __repr__(self):
        return self._base_repr_(html=False)

    def __unicode__(self):
        return '\n'.join(self.pformat())
    if six.PY3:
        __str__ = __unicode__

    def __bytes__(self):
        return six.text_type(self).encode('utf-8')
    if six.PY2:
        __str__ = __bytes__

    @property
    def has_mixin_columns(self):
        """
        True if table has any mixin columns (defined as columns that are not Column
        subclasses)
        """
        return any(not isinstance(col, BaseColumn) for col in self.columns.values())

    def _is_mixin_column(self, col, quantity_is_mixin=False):
        """
        Determine if ``col`` meets the protocol for a mixin Table column for
        this table.  By definition a BaseColumn instance is not a mixin.
        """
        if isinstance(col, BaseColumn):
            is_mixin = False
        elif isinstance(col, Quantity):
            is_mixin = quantity_is_mixin
        else:
            is_mixin = is_mixin_class(col)

        return is_mixin

    def pprint(self, max_lines=None, max_width=None, show_name=True,
               show_unit=None, show_dtype=False, align='right'):
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

        show_dtype : bool
            Include a header row for column dtypes (default=True)

        align : str
            Left/right alignment of a column. Default is 'right'.
        """
        lines, outs = self.formatter._pformat_table(self, max_lines, max_width,
                                                    show_name=show_name, show_unit=show_unit,
                                                    show_dtype=show_dtype)
        if outs['show_length']:
            lines.append('Length = {0} rows'.format(len(self)))

        n_header = outs['n_header']

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
                show_unit=None, show_dtype=False, html=False, tableid=None,
                align='right'):
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

        show_dtype : bool
            Include a header row for column dtypes (default=True)

        html : bool
            Format the output as an HTML table (default=False)

        tableid : str or `None`
            An ID tag for the table; only used if html is set.  Default is
            "table{id}", where id is the unique integer id of the table object,
            id(self)
            
        align : str
            Left/right alignment of a column. Default is 'right'.    

        Returns
        -------
        lines : list
            Formatted table as a list of strings

        """

        lines, outs = self.formatter._pformat_table(self, max_lines, max_width,
                                                    show_name=show_name, show_unit=show_unit,
                                                    show_dtype=show_dtype, html=html,
                                                    tableid=tableid)

        if outs['show_length']:
            lines.append('Length = {0} rows'.format(len(self)))

        return lines

    def more(self, max_lines=None, max_width=None, show_name=True,
             show_unit=None, show_dtype=False):
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

        show_dtype : bool
            Include a header row for column dtypes (default=True)
        """
        self.formatter._more_tabcol(self, max_lines, max_width, show_name=show_name,
                                    show_unit=show_unit, show_dtype=show_dtype)

    def __getitem__(self, item):
        if isinstance(item, six.string_types):
            return self.columns[item]
        elif isinstance(item, (int, np.integer)):
            return self.Row(self, item)
        elif (isinstance(item, (tuple, list)) and item and
              all(isinstance(x, six.string_types) for x in item)):
            bad_names = [x for x in item if x not in self.colnames]
            if bad_names:
                raise ValueError('Slice name(s) {0} not valid column name(s)'
                                 .format(', '.join(bad_names)))
            out = self.__class__([self[x] for x in item], meta=deepcopy(self.meta))
            out._groups = groups.TableGroups(out, indices=self.groups._indices,
                                             keys=self.groups._keys)
            return out
        elif ((isinstance(item, np.ndarray) and len(item) == 0) or
              (isinstance(item, (tuple, list)) and not item)):
            # If item is an empty array/list/tuple then return the table with no rows
            return self._new_from_slice([])
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

            # Make sure value has a dtype.  If not make it into a numpy array
            if not hasattr(value, 'dtype') and not self._is_mixin_column(value):
                value = np.asarray(value)

            # Make new column and assign the value.  If the table currently
            # has no rows (len=0) of the value is already a Column then
            # define new column directly from value.  In the latter case
            # this allows for propagation of Column metadata.  Otherwise
            # define a new column with the right length and shape and then
            # set it from value.  This allows for broadcasting, e.g. t['a']
            # = 1.
            name = item
            if isinstance(value, BaseColumn) or self._is_mixin_column(value):
                new_column = col_copy(value)
                col_setattr(new_column, 'name', name)
            elif len(self) == 0:
                new_column = NewColumn(value, name=name)
            else:
                new_column = NewColumn(name=name, length=len(self), dtype=value.dtype,
                                       shape=value.shape[1:])
                new_column[:] = value

                if isinstance(value, Quantity):
                    new_column.unit = value.unit

            # Now add new column to the table
            self.add_columns([new_column], copy=False)

        else:
            n_cols = len(self.columns)

            if isinstance(item, six.string_types):
                # Set an existing column
                self.columns[item][:] = value

            elif isinstance(item, (int, np.integer)):
                # Set the corresponding row assuming value is an iterable.
                if not hasattr(value, '__len__'):
                    raise TypeError('Right side value must be iterable')

                if len(value) != n_cols:
                    raise ValueError('Right side value needs {0} elements (one for each column)'
                                     .format(n_cols))

                for col, val in izip(self.columns.values(), value):
                    col[item] = val

            elif (isinstance(item, slice) or
                  isinstance(item, np.ndarray) or
                  isinstance(item, list) or
                  (isinstance(item, tuple) and  # output from np.where
                   all(isinstance(x, np.ndarray) for x in item))):

                if isinstance(value, Table):
                    vals = (col for col in value.columns.values())

                elif isinstance(value, np.ndarray) and value.dtype.names:
                    vals = (value[name] for name in value.dtype.names)

                elif np.isscalar(value):
                    import itertools
                    vals = itertools.repeat(value, n_cols)

                else:  # Assume this is an iterable that will work
                    if len(value) != n_cols:
                        raise ValueError('Right side value needs {0} elements (one for each column)'
                                         .format(n_cols))
                    vals = value

                for col, val in izip(self.columns.values(), vals):
                    col[item] = val

            else:
                raise ValueError('Illegal type {0} for table item access'
                                 .format(type(item)))

    def __delitem__(self, item):
        if isinstance(item, six.string_types):
            self.remove_column(item)
        elif isinstance(item, tuple):
            self.remove_columns(item)

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
        return np.dtype([descr(col) for col in self.columns.values()])

    @property
    def colnames(self):
        return list(self.columns.keys())

    def keys(self):
        return list(self.columns.keys())

    def __len__(self):
        if len(self.columns) == 0:
            return 0

        lengths = set(len(col) for col in self.columns.values())
        if len(lengths) != 1:
            len_strs = [' {0} : {1}'.format(name, len(col)) for name, col in self.columns.items()]
            raise ValueError('Column length mismatch:\n{0}'.format('\n'.join(len_strs)))

        return lengths.pop()

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

    def add_column(self, col, index=None, rename_duplicate=False):
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
        rename_duplicate : bool
            Uniquify column name if it already exist (default=False)

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

        Add second column named 'b' with rename_duplicate::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> col_b = Column(name='b', data=[1.1, 1.2, 1.3])
            >>> t.add_column(col_b, rename_duplicate=True)
            >>> print(t)
             a   b  b_1
            --- --- ---
              1 0.1 1.1
              2 0.2 1.2
              3 0.3 1.3

        To add several columns use add_columns.
        """
        if index is None:
            index = len(self.columns)
        self.add_columns([col], [index], rename_duplicate=rename_duplicate)

    def add_columns(self, cols, indexes=None, copy=True, rename_duplicate=False):
        """
        Add a list of new Column objects ``cols`` to the table.  If a
        corresponding list of ``indexes`` is supplied then insert column
        before each ``index`` position in the *original* list of columns,
        otherwise append columns to the end of the list.

        Parameters
        ----------
        cols : list of Columns
            Column objects to add.
        indexes : list of ints or `None`
            Insert column before this position or at end (default)
        copy : bool
            Make a copy of the new columns (default=True)
        rename_duplicate : bool
            Uniquify new column names if they duplicate the existing ones
            (default=False)


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

        Add second column 'b' and column 'c' with ``rename_duplicate``::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> col_b = Column(name='b', data=[1.1, 1.2, 1.3])
            >>> col_c = Column(name='c', data=['x', 'y', 'z'])
            >>> t.add_columns([col_b, col_c], rename_duplicate=True)
            >>> print(t)
             a   b  b_1  c
            --- --- --- ---
              1 0.1 1.1  x
              2 0.2 1.2  y
              3 0.3 1.3  z

        """
        if indexes is None:
            indexes = [len(self.columns)] * len(cols)
        elif len(indexes) != len(cols):
            raise ValueError('Number of indexes must match number of cols')

        if copy:
            cols = [col_copy(col) for col in cols]

        if len(self.columns) == 0:
            # No existing table data, init from cols
            newcols = cols
        else:
            newcols = list(self.columns.values())
            new_indexes = list(range(len(newcols) + 1))
            for col, index in zip(cols, indexes):
                i = new_indexes.index(index)
                new_indexes.insert(i, None)
                newcols.insert(i, col)

        if rename_duplicate:
            existing_names = set(self.colnames)
            for col in cols:
                i = 1
                orig_name = col_getattr(col, 'name')
                while col_getattr(col, 'name') in existing_names:
                    # If the column belongs to another table then copy it
                    # before renaming
                    if col_getattr(col, 'parent_table') is not None:
                        col = col_copy(col)
                    new_name = '{0}_{1}'.format(orig_name, i)
                    col_setattr(col, 'name', new_name)
                    i += 1
                existing_names.add(new_name)

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
        keep_mask = np.ones(len(self), dtype=np.bool)
        keep_mask[row_specifier] = False

        columns = self.TableColumns()
        for name, col in self.columns.items():
            newcol = col[keep_mask]
            col_setattr(newcol, 'parent_table', self)
            columns[name] = newcol

        self.columns = columns

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

        TODO: this won't work for mixins

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

        if not isinstance(self.columns[name], BaseColumn):
            raise NotImplementedError('cannot rename a mixin column')

        col_setattr(self.columns[name], 'name', new_name)

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
        self.insert_row(len(self), vals, mask)

    def insert_row(self, index, vals=None, mask=None):
        """Add a new row before the given ``index`` position in the table.

        The ``vals`` argument can be:

        sequence (e.g. tuple or list)
            Column values in the same order as table columns.
        mapping (e.g. dict)
            Keys corresponding to column names.  Missing values will be
            filled with np.zeros for the column dtype.
        `None`
            All values filled with np.zeros for the column dtype.

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
        """
        colnames = self.colnames

        N = len(self)
        if index < -N or index > N:
            raise IndexError("Index {0} is out of bounds for table with length {1}"
                             .format(index, N))
        if index < 0:
            index += N

        def _is_mapping(obj):
            """Minimal checker for mapping (dict-like) interface for obj"""
            attrs = ('__getitem__', '__len__', '__iter__', 'keys', 'values', 'items')
            return all(hasattr(obj, attr) for attr in attrs)

        if mask is not None and not self.masked:
            # Possibly issue upgrade warning and update self.ColumnClass.  This
            # does not change the existing columns.
            self._set_masked(True)

        if _is_mapping(vals) or vals is None:
            # From the vals and/or mask mappings create the corresponding lists
            # that have entries for each table column.
            if mask is not None and not _is_mapping(mask):
                raise TypeError("Mismatch between type of vals and mask")

            # Now check that the mask is specified for the same keys as the
            # values, otherwise things get really confusing.
            if mask is not None and set(vals.keys()) != set(mask.keys()):
                raise ValueError('keys in mask should match keys in vals')

            if vals and any(name not in colnames for name in vals):
                raise ValueError('Keys in vals must all be valid column names')

            vals_list = []
            mask_list = []

            for name in colnames:
                if vals and name in vals:
                    vals_list.append(vals[name])
                    mask_list.append(False if mask is None else mask[name])
                else:
                    col = self[name]
                    if hasattr(col, 'dtype'):
                        # Make a placeholder zero element of the right type which is masked.
                        # This assumes the appropriate insert() method will broadcast a
                        # numpy scalar to the right shape.
                        vals_list.append(np.zeros(shape=(), dtype=col.dtype))

                        # For masked table any unsupplied values are masked by default.
                        mask_list.append(self.masked and vals is not None)
                    else:
                        raise ValueError("Value must be supplied for column '{0}'".format(name))

            vals = vals_list
            mask = mask_list

        if isiterable(vals):
            if mask is not None and (not isiterable(mask) or _is_mapping(mask)):
                raise TypeError("Mismatch between type of vals and mask")

            if len(self.columns) != len(vals):
                raise ValueError('Mismatch between number of vals and columns')

            if mask is not None:
                if len(self.columns) != len(mask):
                    raise ValueError('Mismatch between number of masks and columns')
            else:
                mask = [False] * len(self.columns)

        else:
            raise TypeError('Vals must be an iterable or mapping or None')

        columns = self.TableColumns()
        try:
            # Insert val at index for each column
            for name, col, val, mask_ in izip(colnames, self.columns.values(), vals, mask):
                # If the new row caused a change in self.ColumnClass then
                # Column-based classes need to be converted first.  This is
                # typical for adding a row with mask values to an unmasked table.
                if isinstance(col, Column) and not isinstance(col, self.ColumnClass):
                    col = self.ColumnClass(col, copy=False)

                newcol = col.insert(index, val)
                if not isinstance(newcol, BaseColumn):
                    _col_update_attrs_from(newcol, col)
                    col_setattr(newcol, 'name', name)
                    if self.masked:
                        newcol.mask = FalseArray(newcol.shape)

                if len(newcol) != N + 1:
                    raise ValueError('Incorrect length for column {0} after inserting {1}'
                                     ' (expected {2}, got {3})'
                                     .format(name, val, len(newcol), N + 1))
                col_setattr(newcol, 'parent_table', self)

                # Set mask if needed
                if self.masked:
                    newcol.mask[index] = mask_

                columns[name] = newcol

        except Exception as err:
            raise ValueError("Unable to insert row because of exception in column '{0}':\n{1}"
                             .format(name, err))
        else:
            self.columns = columns

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

        if keys:
            data = self[keys].as_array()
        else:
            data = self.as_array()

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

        indexes = self.argsort(keys)
        for col in self.columns.values():
            col[:] = col.take(indexes, axis=0)

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
        for col in self.columns.values():
            col[:] = col[::-1]

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
            other = other.as_array()

        if self.masked:
            if isinstance(other, np.ma.MaskedArray):
                result = self.as_array() == other
            else:
                # If mask is True, then by definition the row doesn't match
                # because the other array is not masked.
                false_mask = np.zeros(1, dtype=[(n, bool) for n in self.dtype.names])
                result = (self.as_array().data == other) & (self.mask == false_mask)
        else:
            if isinstance(other, np.ma.MaskedArray):
                # If mask is True, then by definition the row doesn't match
                # because the other array is not masked.
                false_mask = np.zeros(1, dtype=[(n, bool) for n in other.dtype.names])
                result = (self.as_array() == other.data) & (other.mask == false_mask)
            else:
                result = self.as_array() == other

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
        if self.has_mixin_columns:
            raise NotImplementedError('group_by not available for tables with mixin columns')

        return groups.table_group_by(self, keys)


class QTable(Table):
    """A class to represent tables of heterogeneous data.

    `QTable` provides a class for heterogeneous tabular data which can be
    easily modified, for instance adding columns or new rows.

    The `QTable` class is identical to `Table` except that columns with an
    associated ``unit`` attribute are converted to `~astropy.units.Quantity`
    objects.

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
    def __init__(self, data=None, masked=None, names=None, dtype=None,
                 meta=None, copy=True, rows=None):
        super(QTable, self).__init__(data, masked, names, dtype, meta, copy, rows)

    def _is_mixin_column(self, col):
        """
        Determine if ``col`` meets the protocol for a mixin Table column for
        this table.  By definition a BaseColumn instance is not a mixin.

        If ``col`` is a string then it refers to a column name in this table.
        """
        return super(QTable, self)._is_mixin_column(col, quantity_is_mixin=True)
