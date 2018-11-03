# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .index import TableIndices, TableLoc, TableILoc, TableLocIndices

import re
import sys
from collections import OrderedDict
from collections.abc import Mapping
import warnings
from copy import deepcopy

import numpy as np
from numpy import ma

from .. import log
from ..io import registry as io_registry
from ..units import Quantity, QuantityInfo
from ..utils import isiterable, ShapedLikeNDArray
from ..utils.console import color_print
from ..utils.metadata import MetaData
from ..utils.data_info import BaseColumnInfo, MixinInfo, ParentDtypeInfo, DataInfo
from ..utils.exceptions import AstropyDeprecationWarning, NoValue

from . import groups
from .pprint import TableFormatter
from .column import (BaseColumn, Column, MaskedColumn, _auto_names, FalseArray,
                     col_copy)
from .row import Row
from .np_utils import fix_column_name, recarray_fromrecords
from .info import TableInfo, serialize_method_as
from .index import Index, _IndexModeContext, get_index
from . import conf


__doctest_skip__ = ['Table.read', 'Table.write',
                    'Table.convert_bytestring_to_unicode',
                    'Table.convert_unicode_to_bytestring',
                    ]


class TableReplaceWarning(UserWarning):
    """
    Warning class for cases when a table column is replaced via the
    Table.__setitem__ syntax e.g. t['a'] = val.

    This does not inherit from AstropyWarning because we want to use
    stacklevel=3 to show the user where the issue occurred in their code.
    """
    pass


def descr(col):
    """Array-interface compliant full description of a column.

    This returns a 3-tuple (name, type, shape) that can always be
    used in a structured array dtype definition.
    """
    col_dtype = 'O' if (col.info.dtype is None) else col.info.dtype
    col_shape = col.shape[1:] if hasattr(col, 'shape') else ()
    return (col.info.name, col_dtype, col_shape)


def has_info_class(obj, cls):
    return hasattr(obj, 'info') and isinstance(obj.info, cls)


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
                if has_info_class(col, BaseColumnInfo):
                    newcols.append((col.info.name, col))
                else:
                    newcols.append(col)
            cols = newcols
        super().__init__(cols)

    def __getitem__(self, item):
        """Get items from a TableColumns object.
        ::

          tc = TableColumns(cols=[Column(name='a'), Column(name='b'), Column(name='c')])
          tc['a']  # Column('a')
          tc[1] # Column('b')
          tc['a', 'b'] # <TableColumns names=('a', 'b')>
          tc[1:3] # <TableColumns names=('b', 'c')>
        """
        if isinstance(item, str):
            return OrderedDict.__getitem__(self, item)
        elif isinstance(item, (int, np.integer)):
            return self.values()[item]
        elif (isinstance(item, np.ndarray) and item.shape == () and item.dtype.kind == 'i'):
            return self.values()[item.item()]
        elif isinstance(item, tuple):
            return self.__class__([self[x] for x in item])
        elif isinstance(item, slice):
            return self.__class__([self[x] for x in list(self)[item]])
        else:
            raise IndexError('Illegal key or index value for {} object'
                             .format(self.__class__.__name__))

    def __setitem__(self, item, value):
        if item in self:
            raise ValueError("Cannot replace column '{0}'.  Use Table.replace_column() instead."
                             .format(item))
        super().__setitem__(item, value)

    def __repr__(self):
        names = ("'{0}'".format(x) for x in self.keys())
        return "<{1} names=({0})>".format(",".join(names), self.__class__.__name__)

    def _rename_column(self, name, new_name):
        if name == new_name:
            return

        if new_name in self:
            raise KeyError("Column {0} already exists".format(new_name))

        mapper = {name: new_name}
        new_names = [mapper.get(name, name) for name in self]
        cols = list(self.values())
        self.clear()
        self.update(list(zip(new_names, cols)))

    # Define keys and values for Python 2 and 3 source compatibility
    def keys(self):
        return list(OrderedDict.keys(self))

    def values(self):
        return list(OrderedDict.values(self))

    def isinstance(self, cls):
        """
        Return a list of columns which are instances of the specified classes.

        Parameters
        ----------
        cls : class or tuple of classes
            Column class (including mixin) or tuple of Column classes.

        Returns
        -------
        col_list : list of Columns
            List of Column objects which are instances of given classes.
        """
        cols = [col for col in self.values() if isinstance(col, cls)]
        return cols

    def not_isinstance(self, cls):
        """
        Return a list of columns which are not instances of the specified classes.

        Parameters
        ----------
        cls : class or tuple of classes
            Column class (including mixin) or tuple of Column classes.

        Returns
        -------
        col_list : list of Columns
            List of Column objects which are not instances of given classes.
        """
        cols = [col for col in self.values() if not isinstance(col, cls)]
        return cols


class Table:
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
    data : numpy ndarray, dict, list, Table, or table-like object, optional
        Data to initialize table.
    masked : bool, optional
        Specify whether the table is masked.
    names : list, optional
        Specify column names.
    dtype : list, optional
        Specify column data types.
    meta : dict, optional
        Metadata associated with the table.
    copy : bool, optional
        Copy the input data. If the input is a Table the ``meta`` is always
        copied regardless of the ``copy`` parameter.
        Default is True.
    rows : numpy ndarray, list of lists, optional
        Row-oriented data for table instead of ``data`` argument.
    copy_indices : bool, optional
        Copy any indices in the input data. Default is True.
    **kwargs : dict, optional
        Additional keyword args when converting table-like object.
    """

    meta = MetaData()

    # Define class attributes for core container objects to allow for subclass
    # customization.
    Row = Row
    Column = Column
    MaskedColumn = MaskedColumn
    TableColumns = TableColumns
    TableFormatter = TableFormatter

    def as_array(self, keep_byteorder=False):
        """
        Return a new copy of the table in the form of a structured np.ndarray or
        np.ma.MaskedArray object (as appropriate).

        Parameters
        ----------
        keep_byteorder : bool, optional
            By default the returned array has all columns in native byte
            order.  However, if this option is `True` this preserves the
            byte order of all columns (if any are non-native).

        Returns
        -------
        table_array : np.ndarray (unmasked) or np.ma.MaskedArray (masked)
            Copy of table as a numpy structured array
        """
        if len(self.columns) == 0:
            return None

        sys_byteorder = ('>', '<')[sys.byteorder == 'little']
        native_order = ('=', sys_byteorder)

        dtype = []

        cols = self.columns.values()

        for col in cols:
            col_descr = descr(col)
            byteorder = col.info.dtype.byteorder

            if not keep_byteorder and byteorder not in native_order:
                new_dt = np.dtype(col_descr[1]).newbyteorder('=')
                col_descr = (col_descr[0], new_dt, col_descr[2])

            dtype.append(col_descr)

        empty_init = ma.empty if self.masked else np.empty
        data = empty_init(len(self), dtype=dtype)
        for col in cols:
            # When assigning from one array into a field of a structured array,
            # Numpy will automatically swap those columns to their destination
            # byte order where applicable
            data[col.info.name] = col

        return data

    def __init__(self, data=None, masked=None, names=None, dtype=None,
                 meta=None, copy=True, rows=None, copy_indices=True,
                 **kwargs):

        # Set up a placeholder empty table
        self._set_masked(masked)
        self.columns = self.TableColumns()
        self.meta = meta
        self.formatter = self.TableFormatter()
        self._copy_indices = True  # copy indices from this Table by default
        self._init_indices = copy_indices  # whether to copy indices in init
        self.primary_key = None

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

        if hasattr(data, '__astropy_table__'):
            # Data object implements the __astropy_table__ interface method.
            # Calling that method returns an appropriate instance of
            # self.__class__ and respects the `copy` arg.  The returned
            # Table object should NOT then be copied (though the meta
            # will be deep-copied anyway).
            data = data.__astropy_table__(self.__class__, copy, **kwargs)
            copy = False
        elif kwargs:
            raise TypeError('__init__() got unexpected keyword argument {!r}'
                            .format(list(kwargs.keys())[0]))

        if (isinstance(data, np.ndarray) and
                data.shape == (0,) and
                not data.dtype.names):
            data = None

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
                if data.shape == ():
                    raise ValueError('Can not initialize a Table with a scalar')
                elif len(data.shape) == 1:
                    data = data[np.newaxis, :]
                n_cols = data.shape[1]

        elif isinstance(data, Mapping):
            init_func = self._init_from_dict
            default_names = list(data)
            n_cols = len(default_names)

        elif isinstance(data, Table):
            init_func = self._init_from_table
            n_cols = len(data.colnames)
            default_names = data.colnames
            # don't copy indices if the input Table is in non-copy mode
            self._init_indices = self._init_indices and data._copy_indices

        elif data is None:
            if names is None:
                if dtype is None:
                    return  # Empty table
                try:
                    # No data nor names but dtype is available.  This must be
                    # valid to initialize a structured array.
                    dtype = np.dtype(dtype)
                    names = dtype.names
                    dtype = [dtype[name] for name in names]
                except Exception:
                    raise ValueError('dtype was specified but could not be '
                                     'parsed for column names')
            # names is guaranteed to be set at this point
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

        # Numpy does not support bytes column names on Python 3, so fix them
        # up now.
        names = [fix_column_name(name) for name in names]

        self._check_names_dtype(names, dtype, n_cols)

        # Finally do the real initialization
        init_func(data, names, dtype, n_cols, copy)

        # Whatever happens above, the masked property should be set to a boolean
        if type(self.masked) is not bool:
            raise TypeError("masked property has not been set to True or False")

    def __getstate__(self):
        columns = OrderedDict((key, col if isinstance(col, BaseColumn) else col_copy(col))
                              for key, col in self.columns.items())
        return (columns, self.meta)

    def __setstate__(self, state):
        columns, meta = state
        self.__init__(columns, meta=meta)

    @property
    def mask(self):
        # Dynamic view of available masks
        if self.masked:
            mask_table = Table([col.mask for col in self.columns.values()],
                               names=self.colnames, copy=False)

            # Set hidden attribute to force inplace setitem so that code like
            # t.mask['a'] = [1, 0, 1] will correctly set the underlying mask.
            # See #5556 for discussion.
            mask_table._setitem_inplace = True
        else:
            mask_table = None

        return mask_table

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
            data = [col.filled(fill_value) for col in self.columns.values()]
        else:
            data = self
        return self.__class__(data, meta=deepcopy(self.meta))

    @property
    def indices(self):
        '''
        Return the indices associated with columns of the table
        as a TableIndices object.
        '''
        lst = []
        for column in self.columns.values():
            for index in column.info.indices:
                if sum([index is x for x in lst]) == 0:  # ensure uniqueness
                    lst.append(index)
        return TableIndices(lst)

    @property
    def loc(self):
        '''
        Return a TableLoc object that can be used for retrieving
        rows by index in a given data range. Note that both loc
        and iloc work only with single-column indices.
        '''
        return TableLoc(self)

    @property
    def loc_indices(self):
        """
        Return a TableLocIndices object that can be used for retrieving
        the row indices corresponding to given table index key value or values.
        """
        return TableLocIndices(self)

    @property
    def iloc(self):
        '''
        Return a TableILoc object that can be used for retrieving
        indexed rows in the order they appear in the index.
        '''
        return TableILoc(self)

    def add_index(self, colnames, engine=None, unique=False):
        '''
        Insert a new index among one or more columns.
        If there are no indices, make this index the
        primary table index.

        Parameters
        ----------
        colnames : str or list
            List of column names (or a single column name) to index
        engine : type or None
            Indexing engine class to use, from among SortedArray, BST,
            FastBST, FastRBT, and SCEngine. If the supplied argument is None
            (by default), use SortedArray.
        unique : bool
            Whether the values of the index must be unique. Default is False.
        '''
        if isinstance(colnames, str):
            colnames = (colnames,)
        columns = self.columns[tuple(colnames)].values()

        # make sure all columns support indexing
        for col in columns:
            if not getattr(col.info, '_supports_indexing', False):
                raise ValueError('Cannot create an index on column "{0}", of '
                                 'type "{1}"'.format(col.info.name, type(col)))

        index = Index(columns, engine=engine, unique=unique)
        if not self.indices:
            self.primary_key = colnames
        for col in columns:
            col.info.indices.append(index)

    def remove_indices(self, colname):
        '''
        Remove all indices involving the given column.
        If the primary index is removed, the new primary
        index will be the most recently added remaining
        index.

        Parameters
        ----------
        colname : str
            Name of column
        '''
        col = self.columns[colname]
        for index in self.indices:
            try:
                index.col_position(col.info.name)
            except ValueError:
                pass
            else:
                for c in index.columns:
                    c.info.indices.remove(index)

    def index_mode(self, mode):
        '''
        Return a context manager for an indexing mode.

        Parameters
        ----------
        mode : str
            Either 'freeze', 'copy_on_getitem', or 'discard_on_copy'.
            In 'discard_on_copy' mode,
            indices are not copied whenever columns or tables are copied.
            In 'freeze' mode, indices are not modified whenever columns are
            modified; at the exit of the context, indices refresh themselves
            based on column values. This mode is intended for scenarios in
            which one intends to make many additions or modifications in an
            indexed column.
            In 'copy_on_getitem' mode, indices are copied when taking column
            slices as well as table slices, so col[i0:i1] will preserve
            indices.
        '''
        return _IndexModeContext(self, mode)

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
            # Structured ndarray gets viewed as a mixin unless already a valid
            # mixin class
            if (isinstance(col, np.ndarray) and len(col.dtype) > 1 and
                    not self._add_as_mixin_column(col)):
                col = col.view(NdarrayMixin)

            if isinstance(col, (Column, MaskedColumn)):
                col = self.ColumnClass(name=(name or col.info.name or def_name),
                                       data=col, dtype=dtype,
                                       copy=copy, copy_indices=self._init_indices)
            elif self._add_as_mixin_column(col):
                # Copy the mixin column attributes if they exist since the copy below
                # may not get this attribute.
                if copy:
                    col = col_copy(col, copy_indices=self._init_indices)

                col.info.name = name or col.info.name or def_name
            elif isinstance(col, np.ndarray) or isiterable(col):
                col = self.ColumnClass(name=(name or def_name), data=col, dtype=dtype,
                                       copy=copy, copy_indices=self._init_indices)
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
                columns[name].info.parent_table = self
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
        self.primary_key = table.primary_key
        cols = list(table.columns.values())

        self._init_from_list(cols, names, dtype, n_cols, copy)

    def _convert_col_for_table(self, col):
        """
        Make sure that all Column objects have correct class for this type of
        Table.  For a base Table this most commonly means setting to
        MaskedColumn if the table is masked.  Table subclasses like QTable
        override this method.
        """
        if col.__class__ is not self.ColumnClass and isinstance(col, Column):
            col = self.ColumnClass(col)  # copy attributes and reference data
        return col

    def _init_from_cols(self, cols):
        """Initialize table from a list of Column or mixin objects"""

        lengths = set(len(col) for col in cols)
        if len(lengths) != 1:
            raise ValueError('Inconsistent data column lengths: {0}'
                             .format(lengths))

        # Set the table masking
        self._set_masked_from_cols(cols)

        # Make sure that all Column-based objects have correct class.  For
        # plain Table this is self.ColumnClass, but for instance QTable will
        # convert columns with units to a Quantity mixin.
        newcols = [self._convert_col_for_table(col) for col in cols]
        self._make_table_from_cols(self, newcols)

        # Deduplicate indices.  It may happen that after pickling or when
        # initing from an existing table that column indices which had been
        # references to a single index object got *copied* into an independent
        # object.  This results in duplicates which will cause downstream problems.
        index_dict = {}
        for col in self.itercols():
            for i, index in enumerate(col.info.indices or []):
                names = tuple(ind_col.info.name for ind_col in index.columns)
                if names in index_dict:
                    col.info.indices[i] = index_dict[names]
                else:
                    index_dict[names] = index

    def _new_from_slice(self, slice_):
        """Create a new table as a referenced slice from self."""

        table = self.__class__(masked=self.masked)
        table.meta.clear()
        table.meta.update(deepcopy(self.meta))
        table.primary_key = self.primary_key
        cols = self.columns.values()

        newcols = []
        for col in cols:
            col.info._copy_indices = self._copy_indices
            newcol = col[slice_]
            if col.info.indices:
                newcol = col.info.slice_indices(newcol, slice_, len(col))
            newcols.append(newcol)
            col.info._copy_indices = True

        self._make_table_from_cols(table, newcols)
        return table

    @staticmethod
    def _make_table_from_cols(table, cols):
        """
        Make ``table`` in-place so that it represents the given list of ``cols``.
        """
        colnames = set(col.info.name for col in cols)
        if None in colnames:
            raise TypeError('Cannot have None for column name')
        if len(colnames) != len(cols):
            raise ValueError('Duplicate column names')

        columns = table.TableColumns((col.info.name, col) for col in cols)

        for col in cols:
            col.info.parent_table = table
            if table.masked and not hasattr(col, 'mask'):
                col.mask = FalseArray(col.shape)

        table.columns = columns

    def itercols(self):
        """
        Iterate over the columns of this table.

        Examples
        --------

        To iterate over the columns of a table::

            >>> t = Table([[1], [2]])
            >>> for col in t.itercols():
            ...     print(col)
            col0
            ----
               1
            col1
            ----
               2

        Using ``itercols()`` is similar to  ``for col in t.columns.values()``
        but is syntactically preferred.
        """
        for colname in self.columns:
            yield self[colname]

    def _base_repr_(self, html=False, descr_vals=None, max_width=None,
                    tableid=None, show_dtype=True, max_lines=None,
                    tableclass=None):
        if descr_vals is None:
            descr_vals = [self.__class__.__name__]
            if self.masked:
                descr_vals.append('masked=True')
            descr_vals.append('length={0}'.format(len(self)))

        descr = ' '.join(descr_vals)
        if html:
            from ..utils.xml.writer import xml_escape
            descr = '<i>{0}</i>\n'.format(xml_escape(descr))
        else:
            descr = '<{0}>\n'.format(descr)

        if tableid is None:
            tableid = 'table{id}'.format(id=id(self))

        data_lines, outs = self.formatter._pformat_table(
            self, tableid=tableid, html=html, max_width=max_width,
            show_name=True, show_unit=None, show_dtype=show_dtype,
            max_lines=max_lines, tableclass=tableclass)

        out = descr + '\n'.join(data_lines)

        return out

    def _repr_html_(self):
        return self._base_repr_(html=True, max_width=-1,
                                tableclass=conf.default_notebook_table_class)

    def __repr__(self):
        return self._base_repr_(html=False, max_width=None)

    def __str__(self):
        return '\n'.join(self.pformat())

    def __bytes__(self):
        return str(self).encode('utf-8')

    @property
    def has_mixin_columns(self):
        """
        True if table has any mixin columns (defined as columns that are not Column
        subclasses).
        """
        return any(has_info_class(col, MixinInfo) for col in self.columns.values())

    def _add_as_mixin_column(self, col):
        """
        Determine if ``col`` should be added to the table directly as
        a mixin column.
        """
        if isinstance(col, BaseColumn):
            return False

        # Is it a mixin but not not Quantity (which gets converted to Column with
        # unit set).
        return has_info_class(col, MixinInfo) and not has_info_class(col, QuantityInfo)

    def pprint(self, max_lines=None, max_width=None, show_name=True,
               show_unit=None, show_dtype=False, align=None):
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
            Maximum number of lines in table output.

        max_width : int or `None`
            Maximum character width of output.

        show_name : bool
            Include a header row for column names. Default is True.

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include a header row for column dtypes. Default is True.

        align : str or list or tuple or `None`
            Left/right alignment of columns. Default is right (None) for all
            columns. Other allowed values are '>', '<', '^', and '0=' for
            right, left, centered, and 0-padded, respectively. A list of
            strings can be provided for alignment of tables with multiple
            columns.
        """
        lines, outs = self.formatter._pformat_table(self, max_lines, max_width,
                                                    show_name=show_name, show_unit=show_unit,
                                                    show_dtype=show_dtype, align=align)
        if outs['show_length']:
            lines.append('Length = {0} rows'.format(len(self)))

        n_header = outs['n_header']

        for i, line in enumerate(lines):
            if i < n_header:
                color_print(line, 'red')
            else:
                print(line)

    def _make_index_row_display_table(self, index_row_name):
        if index_row_name not in self.columns:
            idx_col = self.ColumnClass(name=index_row_name, data=np.arange(len(self)))
            return self.__class__([idx_col] + self.columns.values(),
                                           copy=False)
        else:
            return self

    def show_in_notebook(self, tableid=None, css=None, display_length=50,
                         table_class='astropy-default', show_row_index='idx'):
        """Render the table in HTML and show it in the IPython notebook.

        Parameters
        ----------
        tableid : str or `None`
            An html ID tag for the table.  Default is ``table{id}-XXX``, where
            id is the unique integer id of the table object, id(self), and XXX
            is a random number to avoid conflicts when printing the same table
            multiple times.
        table_class : str or `None`
            A string with a list of HTML classes used to style the table.
            The special default string ('astropy-default') means that the string
            will be retrieved from the configuration item
            ``astropy.table.default_notebook_table_class``. Note that these
            table classes may make use of bootstrap, as this is loaded with the
            notebook.  See `this page <http://getbootstrap.com/css/#tables>`_
            for the list of classes.
        css : string
            A valid CSS string declaring the formatting for the table. Defaults
            to ``astropy.table.jsviewer.DEFAULT_CSS_NB``.
        display_length : int, optional
            Number or rows to show. Defaults to 50.
        show_row_index : str or False
            If this does not evaluate to False, a column with the given name
            will be added to the version of the table that gets displayed.
            This new column shows the index of the row in the table itself,
            even when the displayed table is re-sorted by another column. Note
            that if a column with this name already exists, this option will be
            ignored. Defaults to "idx".

        Notes
        -----
        Currently, unlike `show_in_browser` (with ``jsviewer=True``), this
        method needs to access online javascript code repositories.  This is due
        to modern browsers' limitations on accessing local files.  Hence, if you
        call this method while offline (and don't have a cached version of
        jquery and jquery.dataTables), you will not get the jsviewer features.
        """

        from .jsviewer import JSViewer
        from IPython.display import HTML

        if tableid is None:
            tableid = 'table{0}-{1}'.format(id(self),
                                            np.random.randint(1, 1e6))

        jsv = JSViewer(display_length=display_length)
        if show_row_index:
            display_table = self._make_index_row_display_table(show_row_index)
        else:
            display_table = self
        if table_class == 'astropy-default':
            table_class = conf.default_notebook_table_class
        html = display_table._base_repr_(html=True, max_width=-1, tableid=tableid,
                                         max_lines=-1, show_dtype=False,
                                         tableclass=table_class)

        columns = display_table.columns.values()
        sortable_columns = [i for i, col in enumerate(columns)
                            if col.info.dtype.kind in 'iufc']
        html += jsv.ipynb(tableid, css=css, sort_columns=sortable_columns)
        return HTML(html)

    def show_in_browser(self, max_lines=5000, jsviewer=False,
                        browser='default', jskwargs={'use_local_files': True},
                        tableid=None, table_class="display compact",
                        css=None, show_row_index='idx'):
        """Render the table in HTML and show it in a web browser.

        Parameters
        ----------
        max_lines : int
            Maximum number of rows to export to the table (set low by default
            to avoid memory issues, since the browser view requires duplicating
            the table in memory).  A negative value of ``max_lines`` indicates
            no row limit.
        jsviewer : bool
            If `True`, prepends some javascript headers so that the table is
            rendered as a `DataTables <https://datatables.net>`_ data table.
            This allows in-browser searching & sorting.
        browser : str
            Any legal browser name, e.g. ``'firefox'``, ``'chrome'``,
            ``'safari'`` (for mac, you may need to use ``'open -a
            "/Applications/Google Chrome.app" {}'`` for Chrome).  If
            ``'default'``, will use the system default browser.
        jskwargs : dict
            Passed to the `astropy.table.JSViewer` init. Defaults to
            ``{'use_local_files': True}`` which means that the JavaScript
            libraries will be served from local copies.
        tableid : str or `None`
            An html ID tag for the table.  Default is ``table{id}``, where id
            is the unique integer id of the table object, id(self).
        table_class : str or `None`
            A string with a list of HTML classes used to style the table.
            Default is "display compact", and other possible values can be
            found in https://www.datatables.net/manual/styling/classes
        css : string
            A valid CSS string declaring the formatting for the table. Defaults
            to ``astropy.table.jsviewer.DEFAULT_CSS``.
        show_row_index : str or False
            If this does not evaluate to False, a column with the given name
            will be added to the version of the table that gets displayed.
            This new column shows the index of the row in the table itself,
            even when the displayed table is re-sorted by another column. Note
            that if a column with this name already exists, this option will be
            ignored. Defaults to "idx".
        """

        import os
        import webbrowser
        import tempfile
        from .jsviewer import DEFAULT_CSS
        from urllib.parse import urljoin
        from urllib.request import pathname2url

        if css is None:
            css = DEFAULT_CSS

        # We can't use NamedTemporaryFile here because it gets deleted as
        # soon as it gets garbage collected.
        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, 'table.html')

        with open(path, 'w') as tmp:
            if jsviewer:
                if show_row_index:
                    display_table = self._make_index_row_display_table(show_row_index)
                else:
                    display_table = self
                display_table.write(tmp, format='jsviewer', css=css,
                                    max_lines=max_lines, jskwargs=jskwargs,
                                    table_id=tableid, table_class=table_class)
            else:
                self.write(tmp, format='html')

        try:
            br = webbrowser.get(None if browser == 'default' else browser)
        except webbrowser.Error:
            log.error("Browser '{}' not found.".format(browser))
        else:
            br.open(urljoin('file:', pathname2url(path)))

    def pformat(self, max_lines=None, max_width=None, show_name=True,
                show_unit=None, show_dtype=False, html=False, tableid=None,
                align=None, tableclass=None):
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
            Include a header row for column names. Default is True.

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include a header row for column dtypes. Default is True.

        html : bool
            Format the output as an HTML table. Default is False.

        tableid : str or `None`
            An ID tag for the table; only used if html is set.  Default is
            "table{id}", where id is the unique integer id of the table object,
            id(self)

        align : str or list or tuple or `None`
            Left/right alignment of columns. Default is right (None) for all
            columns. Other allowed values are '>', '<', '^', and '0=' for
            right, left, centered, and 0-padded, respectively. A list of
            strings can be provided for alignment of tables with multiple
            columns.

        tableclass : str or list of str or `None`
            CSS classes for the table; only used if html is set.  Default is
            None.

        Returns
        -------
        lines : list
            Formatted table as a list of strings.

        """

        lines, outs = self.formatter._pformat_table(
            self, max_lines, max_width, show_name=show_name,
            show_unit=show_unit, show_dtype=show_dtype, html=html,
            tableid=tableid, tableclass=tableclass, align=align)

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
            Include a header row for column names. Default is True.

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include a header row for column dtypes. Default is True.
        """
        self.formatter._more_tabcol(self, max_lines, max_width, show_name=show_name,
                                    show_unit=show_unit, show_dtype=show_dtype)

    def __getitem__(self, item):
        if isinstance(item, str):
            return self.columns[item]
        elif isinstance(item, (int, np.integer)):
            return self.Row(self, item)
        elif (isinstance(item, np.ndarray) and item.shape == () and item.dtype.kind == 'i'):
            return self.Row(self, item.item())
        elif self._is_list_or_tuple_of_str(item):
            out = self.__class__([self[x] for x in item],
                                 meta=deepcopy(self.meta),
                                 copy_indices=self._copy_indices)
            out._groups = groups.TableGroups(out, indices=self.groups._indices,
                                             keys=self.groups._keys)
            return out
        elif ((isinstance(item, np.ndarray) and item.size == 0) or
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
        if isinstance(item, str) and item not in self.colnames:
            NewColumn = self.MaskedColumn if self.masked else self.Column
            # If value doesn't have a dtype and won't be added as a mixin then
            # convert to a numpy array.
            if not hasattr(value, 'dtype') and not self._add_as_mixin_column(value):
                value = np.asarray(value)

            # Structured ndarray gets viewed as a mixin (unless already a valid
            # mixin class).
            if (isinstance(value, np.ndarray) and len(value.dtype) > 1 and
                    not self._add_as_mixin_column(value)):
                value = value.view(NdarrayMixin)

            # Make new column and assign the value.  If the table currently
            # has no rows (len=0) of the value is already a Column then
            # define new column directly from value.  In the latter case
            # this allows for propagation of Column metadata.  Otherwise
            # define a new column with the right length and shape and then
            # set it from value.  This allows for broadcasting, e.g. t['a']
            # = 1.
            name = item
            # If this is a column-like object that could be added directly to table
            if isinstance(value, BaseColumn) or self._add_as_mixin_column(value):
                # If we're setting a new column to a scalar, broadcast it.
                # (things will fail in _init_from_cols if this doesn't work)
                if (len(self) > 0 and (getattr(value, 'isscalar', False) or
                                       getattr(value, 'shape', None) == () or
                                       len(value) == 1)):
                    new_shape = (len(self),) + getattr(value, 'shape', ())[1:]
                    if isinstance(value, np.ndarray):
                        value = np.broadcast_to(value, shape=new_shape,
                                                subok=True)
                    elif isinstance(value, ShapedLikeNDArray):
                        value = value._apply(np.broadcast_to, shape=new_shape,
                                             subok=True)

                new_column = col_copy(value)
                new_column.info.name = name

            elif len(self) == 0:
                new_column = NewColumn(value, name=name)
            else:
                new_column = NewColumn(name=name, length=len(self), dtype=value.dtype,
                                       shape=value.shape[1:],
                                       unit=getattr(value, 'unit', None))
                new_column[:] = value

            # Now add new column to the table
            self.add_columns([new_column], copy=False)

        else:
            n_cols = len(self.columns)

            if isinstance(item, str):
                # Set an existing column by first trying to replace, and if
                # this fails do an in-place update.  See definition of mask
                # property for discussion of the _setitem_inplace attribute.
                if (not getattr(self, '_setitem_inplace', False)
                        and not conf.replace_inplace):
                    try:
                        self._replace_column_warnings(item, value)
                        return
                    except Exception:
                        pass
                self.columns[item][:] = value

            elif isinstance(item, (int, np.integer)):
                self._set_row(idx=item, colnames=self.colnames, vals=value)

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

                for col, val in zip(self.columns.values(), vals):
                    col[item] = val

            else:
                raise ValueError('Illegal type {0} for table item access'
                                 .format(type(item)))

    def __delitem__(self, item):
        if isinstance(item, str):
            self.remove_column(item)
        elif isinstance(item, (int, np.integer)):
            self.remove_row(item)
        elif (isinstance(item, (list, tuple, np.ndarray)) and
              all(isinstance(x, str) for x in item)):
            self.remove_columns(item)
        elif (isinstance(item, (list, np.ndarray)) and
              np.asarray(item).dtype.kind == 'i'):
            self.remove_rows(item)
        elif isinstance(item, slice):
            self.remove_rows(item)
        else:
            raise IndexError('illegal key or index value')

    def _ipython_key_completions_(self):
        return self.colnames

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

    @staticmethod
    def _is_list_or_tuple_of_str(names):
        """Check that ``names`` is a tuple or list of strings"""
        return (isinstance(names, (tuple, list)) and names and
                all(isinstance(x, str) for x in names))

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

    def add_column(self, col, index=None, name=None, rename_duplicate=False, copy=True):
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
            Insert column before this position or at end (default).
        name : str
            Column name
        rename_duplicate : bool
            Uniquify column name if it already exist. Default is False.
        copy : bool
            Make a copy of the new column. Default is True.

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

        Add an unnamed column or mixin object in the table using a default name
        or by specifying an explicit name with ``name``. Name can also be overridden::

            >>> t = Table([[1, 2], [0.1, 0.2]], names=('a', 'b'))
            >>> col_c = Column(data=['x', 'y'])
            >>> t.add_column(col_c)
            >>> t.add_column(col_c, name='c')
            >>> col_b = Column(name='b', data=[1.1, 1.2])
            >>> t.add_column(col_b, name='d')
            >>> print(t)
             a   b  col2  c   d
            --- --- ---- --- ---
              1 0.1    x   x 1.1
              2 0.2    y   y 1.2

        To add several columns use add_columns.
        """
        if index is None:
            index = len(self.columns)
        if name is not None:
            name = (name,)

        self.add_columns([col], [index], name, copy=copy, rename_duplicate=rename_duplicate)

    def add_columns(self, cols, indexes=None, names=None, copy=True, rename_duplicate=False):
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
            Insert column before this position or at end (default).
        names : list of str
            Column names
        copy : bool
            Make a copy of the new columns. Default is True.
        rename_duplicate : bool
            Uniquify new column names if they duplicate the existing ones.
            Default is False.


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

        Add unnamed columns or mixin objects in the table using default names
        or by specifying explicit names with ``names``. Names can also be overridden::

            >>> t = Table()
            >>> col_a = Column(data=['x', 'y'])
            >>> col_b = Column(name='b', data=['u', 'v'])
            >>> t.add_columns([col_a, col_b])
            >>> t.add_columns([col_a, col_b], names=['c', 'd'])
            >>> print(t)
            col0  b   c   d
            ---- --- --- ---
               x   u   x   u
               y   v   y   v
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

        if names is None:
            names = (None,) * len(cols)
        elif len(names) != len(cols):
                raise ValueError('Number of names must match number of cols')

        for i, (col, name) in enumerate(zip(cols, names)):
            if name is None:
                if col.info.name is not None:
                    continue
                name = 'col{}'.format(i + len(self.columns))
            if col.info.parent_table is not None:
                col = col_copy(col)
            col.info.name = name

        if rename_duplicate:
            existing_names = set(self.colnames)
            for col in cols:
                i = 1
                orig_name = col.info.name
                if col.info.name in existing_names:
                    # If the column belongs to another table then copy it
                    # before renaming
                    while col.info.name in existing_names:
                        # Iterate until a unique name is found
                        if col.info.parent_table is not None:
                            col = col_copy(col)
                        new_name = '{0}_{1}'.format(orig_name, i)
                        col.info.name = new_name
                        i += 1
                    existing_names.add(new_name)

        self._init_from_cols(newcols)

    def _replace_column_warnings(self, name, col):
        """
        Same as replace_column but issues warnings under various circumstances.
        """
        warns = conf.replace_warnings

        if 'refcount' in warns and name in self.colnames:
            refcount = sys.getrefcount(self[name])

        if name in self.colnames:
            old_col = self[name]

        # This may raise an exception (e.g. t['a'] = 1) in which case none of
        # the downstream code runs.
        self.replace_column(name, col)

        if 'always' in warns:
            warnings.warn("replaced column '{}'".format(name),
                          TableReplaceWarning, stacklevel=3)

        if 'slice' in warns:
            try:
                # Check for ndarray-subclass slice.  An unsliced instance
                # has an ndarray for the base while sliced has the same class
                # as parent.
                if isinstance(old_col.base, old_col.__class__):
                    msg = ("replaced column '{}' which looks like an array slice. "
                           "The new column no longer shares memory with the "
                           "original array.".format(name))
                    warnings.warn(msg, TableReplaceWarning, stacklevel=3)
            except AttributeError:
                pass

        if 'refcount' in warns:
            # Did reference count change?
            new_refcount = sys.getrefcount(self[name])
            if refcount != new_refcount:
                msg = ("replaced column '{}' and the number of references "
                       "to the column changed.".format(name))
                warnings.warn(msg, TableReplaceWarning, stacklevel=3)

        if 'attributes' in warns:
            # Any of the standard column attributes changed?
            changed_attrs = []
            new_col = self[name]
            # Check base DataInfo attributes that any column will have
            for attr in DataInfo.attr_names:
                if getattr(old_col.info, attr) != getattr(new_col.info, attr):
                    changed_attrs.append(attr)

            if changed_attrs:
                msg = ("replaced column '{}' and column attributes {} changed."
                       .format(name, changed_attrs))
                warnings.warn(msg, TableReplaceWarning, stacklevel=3)

    def replace_column(self, name, col):
        """
        Replace column ``name`` with the new ``col`` object.

        Parameters
        ----------
        name : str
            Name of column to replace
        col : column object (list, ndarray, Column, etc)
            New column object to replace the existing column

        Examples
        --------
        Replace column 'a' with a float version of itself::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> float_a = t['a'].astype(float)
            >>> t.replace_column('a', float_a)
        """
        if name not in self.colnames:
            raise ValueError('column name {0} is not in the table'.format(name))

        if self[name].info.indices:
            raise ValueError('cannot replace a table index column')

        t = self.__class__([col], names=[name])
        cols = OrderedDict(self.columns)
        cols[name] = t[name]
        self._init_from_cols(cols.values())

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
        if not isinstance(index, (int, np.integer)):
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
        # Update indices
        for index in self.indices:
            index.remove_rows(row_specifier)

        keep_mask = np.ones(len(self), dtype=bool)
        keep_mask[row_specifier] = False

        columns = self.TableColumns()
        for name, col in self.columns.items():
            newcol = col[keep_mask]
            newcol.info.parent_table = self
            columns[name] = newcol

        self._replace_cols(columns)

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
        if isinstance(names, str):
            names = [names]

        for name in names:
            if name not in self.columns:
                raise KeyError("Column {0} does not exist".format(name))

        for name in names:
            self.columns.pop(name)

    def _convert_string_dtype(self, in_kind, out_kind, encode_decode_func):
        """
        Convert string-like columns to/from bytestring and unicode (internal only).

        Parameters
        ----------
        in_kind : str
            Input dtype.kind
        out_kind : str
            Output dtype.kind
        """

        for col in self.itercols():
            if col.dtype.kind == in_kind:
                try:
                    # This requires ASCII and is faster by a factor of up to ~8, so
                    # try that first.
                    newcol = col.__class__(col, dtype=out_kind)
                except (UnicodeEncodeError, UnicodeDecodeError):
                    newcol = col.__class__(encode_decode_func(col, 'utf-8'))

                    # Quasi-manually copy info attributes.  Unfortunately
                    # DataInfo.__set__ does not do the right thing in this case
                    # so newcol.info = col.info does not get the old info attributes.
                    for attr in col.info.attr_names - col.info._attrs_no_copy - set(['dtype']):
                        value = deepcopy(getattr(col.info, attr))
                        setattr(newcol.info, attr, value)

                self[col.name] = newcol

    def convert_bytestring_to_unicode(self):
        """
        Convert bytestring columns (dtype.kind='S') to unicode (dtype.kind='U')
        using UTF-8 encoding.

        Internally this changes string columns to represent each character
        in the string with a 4-byte UCS-4 equivalent, so it is inefficient
        for memory but allows scripts to manipulate string arrays with
        natural syntax.
        """
        self._convert_string_dtype('S', 'U', np.char.decode)

    def convert_unicode_to_bytestring(self):
        """
        Convert unicode columns (dtype.kind='U') to bytestring (dtype.kind='S')
        using UTF-8 encoding.

        When exporting a unicode string array to a file, it may be desirable
        to encode unicode columns as bytestrings.
        """
        self._convert_string_dtype('U', 'S', np.char.encode)

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

        if isinstance(names, str):
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

        self.columns[name].info.name = new_name

    def _set_row(self, idx, colnames, vals):
        try:
            assert len(vals) == len(colnames)
        except Exception:
            raise ValueError('right hand side must be a sequence of values with '
                             'the same length as the number of selected columns')

        # Keep track of original values before setting each column so that
        # setting row can be transactional.
        orig_vals = []
        cols = self.columns
        try:
            for name, val in zip(colnames, vals):
                orig_vals.append(cols[name][idx])
                cols[name][idx] = val
        except Exception:
            # If anything went wrong first revert the row update then raise
            for name, val in zip(colnames, orig_vals[:-1]):
                cols[name][idx] = val
            raise

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
            for name, col, val, mask_ in zip(colnames, self.columns.values(), vals, mask):
                # If the new row caused a change in self.ColumnClass then
                # Column-based classes need to be converted first.  This is
                # typical for adding a row with mask values to an unmasked table.
                if isinstance(col, Column) and not isinstance(col, self.ColumnClass):
                    col = self.ColumnClass(col, copy=False)

                newcol = col.insert(index, val, axis=0)
                if not isinstance(newcol, BaseColumn):
                    newcol.info.name = name
                    if self.masked:
                        newcol.mask = FalseArray(newcol.shape)

                if len(newcol) != N + 1:
                    raise ValueError('Incorrect length for column {0} after inserting {1}'
                                     ' (expected {2}, got {3})'
                                     .format(name, val, len(newcol), N + 1))
                newcol.info.parent_table = self

                # Set mask if needed
                if self.masked:
                    newcol.mask[index] = mask_

                columns[name] = newcol

            # insert row in indices
            for table_index in self.indices:
                table_index.insert_row(index, vals, self.columns.values())

        except Exception as err:
            raise ValueError("Unable to insert row because of exception in column '{0}':\n{1}"
                             .format(name, err))
        else:
            self._replace_cols(columns)

            # Revert groups to default (ungrouped) state
            if hasattr(self, '_groups'):
                del self._groups

    def _replace_cols(self, columns):
        for col, new_col in zip(self.columns.values(), columns.values()):
            new_col.info.indices = []
            for index in col.info.indices:
                index.columns[index.col_position(col.info.name)] = new_col
                new_col.info.indices.append(index)

        self.columns = columns

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
        if isinstance(keys, str):
            keys = [keys]

        # use index sorted order if possible
        if keys is not None:
            index = get_index(self, self[keys])
            if index is not None:
                return index.sorted_data()

        kwargs = {}
        if keys:
            kwargs['order'] = keys
        if kind:
            kwargs['kind'] = kind

        if keys:
            data = self[keys].as_array()
        else:
            data = self.as_array()

        return data.argsort(**kwargs)

    def sort(self, keys=None):
        '''
        Sort the table according to one or more keys. This operates
        on the existing table and does not return a new table.

        Parameters
        ----------
        keys : str or list of str
            The key(s) to order the table by. If None, use the
            primary index of the Table.

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
        if keys is None:
            if not self.indices:
                raise ValueError("Table sort requires input keys or a table index")
            keys = [x.info.name for x in self.indices[0].columns]

        if isinstance(keys, str):
            keys = [keys]

        indexes = self.argsort(keys)
        sort_index = get_index(self, self[keys])
        if sort_index is not None:
            # avoid inefficient relabelling of sorted index
            prev_frozen = sort_index._frozen
            sort_index._frozen = True

        for col in self.columns.values():
            col[:] = col.take(indexes, axis=0)

        if sort_index is not None:
            # undo index freeze
            sort_index._frozen = prev_frozen
            # now relabel the sort index appropriately
            sort_index.sort()

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
        for index in self.indices:
            index.reverse()

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

        See http://docs.astropy.org/en/stable/io/unified.html for details.

        Parameters
        ----------
        format : str
            File format specifier.
        *args : tuple, optional
            Positional arguments passed through to data reader. If supplied the
            first argument is the input filename.
        **kwargs : dict, optional
            Keyword arguments passed through to data reader.

        Returns
        -------
        out : `Table`
            Table corresponding to file contents

        Notes
        -----
        """
        # The hanging Notes section just above is a section placeholder for
        # import-time processing that collects available formats into an
        # RST table and inserts at the end of the docstring.  DO NOT REMOVE.

        out = io_registry.read(cls, *args, **kwargs)
        # For some readers (e.g., ascii.ecsv), the returned `out` class is not
        # guaranteed to be the same as the desired output `cls`.  If so,
        # try coercing to desired class without copying (io.registry.read
        # would normally do a copy).  The normal case here is swapping
        # Table <=> QTable.
        if cls is not out.__class__:
            try:
                out = cls(out, copy=False)
            except Exception:
                raise TypeError('could not convert reader output to {0} '
                                'class.'.format(cls.__name__))
        return out

    def write(self, *args, **kwargs):
        """Write this Table object out in the specified format.

        This function provides the Table interface to the astropy unified I/O
        layer.  This allows easily writing a file in many supported data formats
        using syntax such as::

          >>> from astropy.table import Table
          >>> dat = Table([[1, 2], [3, 4]], names=('a', 'b'))
          >>> dat.write('table.dat', format='ascii')

        See http://docs.astropy.org/en/stable/io/unified.html for details.

        Parameters
        ----------
        format : str
            File format specifier.
        serialize_method : str, dict, optional
            Serialization method specifier for columns.
        *args : tuple, optional
            Positional arguments passed through to data writer. If supplied the
            first argument is the output filename.
        **kwargs : dict, optional
            Keyword arguments passed through to data writer.

        Notes
        -----
        """
        serialize_method = kwargs.pop('serialize_method', None)
        with serialize_method_as(self, serialize_method):
            io_registry.write(self, *args, **kwargs)

    def copy(self, copy_data=True):
        '''
        Return a copy of the table.

        Parameters
        ----------
        copy_data : bool
            If `True` (the default), copy the underlying data array.
            Otherwise, use the same data array. The ``meta`` is always
            deepcopied regardless of the value for ``copy_data``.
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
        return super().__lt__(other)

    def __gt__(self, other):
        return super().__gt__(other)

    def __le__(self, other):
        return super().__le__(other)

    def __ge__(self, other):
        return super().__ge__(other)

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
        return groups.table_group_by(self, keys)

    def to_pandas(self):
        """
        Return a :class:`pandas.DataFrame` instance

        Returns
        -------
        dataframe : :class:`pandas.DataFrame`
            A pandas :class:`pandas.DataFrame` instance

        Raises
        ------
        ImportError
            If pandas is not installed
        ValueError
            If the Table contains mixin or multi-dimensional columns
        """
        from pandas import DataFrame

        if self.has_mixin_columns:
            raise ValueError("Cannot convert a table with mixin columns to a pandas DataFrame")

        if any(getattr(col, 'ndim', 1) > 1 for col in self.columns.values()):
            raise ValueError("Cannot convert a table with multi-dimensional columns to a pandas DataFrame")

        out = OrderedDict()

        for name, column in self.columns.items():
            if isinstance(column, MaskedColumn) and np.any(column.mask):
                if column.dtype.kind in ['i', 'u']:
                    out[name] = column.astype(float).filled(np.nan)
                    warnings.warn(
                        "converted column '{}' from integer to float".format(
                            name), TableReplaceWarning, stacklevel=3)
                elif column.dtype.kind in ['f', 'c']:
                    out[name] = column.filled(np.nan)
                else:
                    out[name] = column.astype(object).filled(np.nan)
            else:
                out[name] = column

            if out[name].dtype.byteorder not in ('=', '|'):
                out[name] = out[name].byteswap().newbyteorder()

        return DataFrame(out)

    @classmethod
    def from_pandas(cls, dataframe):
        """
        Create a `Table` from a :class:`pandas.DataFrame` instance

        Parameters
        ----------
        dataframe : :class:`pandas.DataFrame`
            The pandas :class:`pandas.DataFrame` instance

        Returns
        -------
        table : `Table`
            A `Table` (or subclass) instance
        """

        out = OrderedDict()

        for name in dataframe.columns:
            column = dataframe[name]
            mask = np.array(column.isnull())
            data = np.array(column)

            if data.dtype.kind == 'O':
                # If all elements of an object array are string-like or np.nan
                # then coerce back to a native numpy str/unicode array.
                string_types = (str, bytes)
                nan = np.nan
                if all(isinstance(x, string_types) or x is nan for x in data):
                    # Force any missing (null) values to b''.  Numpy will
                    # upcast to str/unicode as needed.
                    data[mask] = b''

                    # When the numpy object array is represented as a list then
                    # numpy initializes to the correct string or unicode type.
                    data = np.array([x for x in data])

            if np.any(mask):
                out[name] = MaskedColumn(data=data, name=name, mask=mask)
            else:
                out[name] = Column(data=data, name=name)

        return cls(out)

    info = TableInfo()


class QTable(Table):
    """A class to represent tables of heterogeneous data.

    `QTable` provides a class for heterogeneous tabular data which can be
    easily modified, for instance adding columns or new rows.

    The `QTable` class is identical to `Table` except that columns with an
    associated ``unit`` attribute are converted to `~astropy.units.Quantity`
    objects.

    Parameters
    ----------
    data : numpy ndarray, dict, list, Table, or table-like object, optional
        Data to initialize table.
    masked : bool, optional
        Specify whether the table is masked.
    names : list, optional
        Specify column names.
    dtype : list, optional
        Specify column data types.
    meta : dict, optional
        Metadata associated with the table.
    copy : bool, optional
        Copy the input data. Default is True.
    rows : numpy ndarray, list of lists, optional
        Row-oriented data for table instead of ``data`` argument.
    copy_indices : bool, optional
        Copy any indices in the input data. Default is True.
    **kwargs : dict, optional
        Additional keyword args when converting table-like object.

    """

    def _add_as_mixin_column(self, col):
        """
        Determine if ``col`` should be added to the table directly as
        a mixin column.
        """
        return has_info_class(col, MixinInfo)

    def _convert_col_for_table(self, col):
        if (isinstance(col, Column) and getattr(col, 'unit', None) is not None):
            # We need to turn the column into a quantity, or a subclass
            # identified in the unit (such as u.mag()).
            q_cls = getattr(col.unit, '_quantity_class', Quantity)
            qcol = q_cls(col.data, col.unit, copy=False)
            qcol.info = col.info
            col = qcol
        else:
            col = super()._convert_col_for_table(col)

        return col


class NdarrayMixin(np.ndarray):
    """
    Mixin column class to allow storage of arbitrary numpy
    ndarrays within a Table.  This is a subclass of numpy.ndarray
    and has the same initialization options as ndarray().
    """
    info = ParentDtypeInfo()

    def __new__(cls, obj, *args, **kwargs):
        self = np.array(obj, *args, **kwargs).view(cls)
        if 'info' in getattr(obj, '__dict__', ()):
            self.info = obj.info
        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return

        if callable(super().__array_finalize__):
            super().__array_finalize__(obj)

        # Self was created from template (e.g. obj[slice] or (obj * 2))
        # or viewcast e.g. obj.view(Column).  In either case we want to
        # init Column attributes for self from obj if possible.
        if 'info' in getattr(obj, '__dict__', ()):
            self.info = obj.info

    def __reduce__(self):
        # patch to pickle Quantity objects (ndarray subclasses), see
        # http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html

        object_state = list(super().__reduce__())
        object_state[2] = (object_state[2], self.__dict__)
        return tuple(object_state)

    def __setstate__(self, state):
        # patch to unpickle NdarrayMixin objects (ndarray subclasses), see
        # http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html

        nd_state, own_state = state
        super().__setstate__(nd_state)
        self.__dict__.update(own_state)
