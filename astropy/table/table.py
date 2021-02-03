# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .index import TableIndices, TableLoc, TableILoc, TableLocIndices

import sys
from collections import OrderedDict, defaultdict
from collections.abc import Mapping
import warnings
from copy import deepcopy
import types
import itertools

import numpy as np
from numpy import ma

from astropy import log
from astropy.units import Quantity, QuantityInfo
from astropy.utils import isiterable, ShapedLikeNDArray
from astropy.utils.console import color_print
from astropy.utils.metadata import MetaData, MetaAttribute
from astropy.utils.data_info import BaseColumnInfo, MixinInfo, ParentDtypeInfo, DataInfo
from astropy.utils.decorators import format_doc
from astropy.io.registry import UnifiedReadWriteMethod

from . import groups
from .pprint import TableFormatter
from .column import (BaseColumn, Column, MaskedColumn, _auto_names, FalseArray,
                     col_copy, _convert_sequence_data_to_array)
from .row import Row
from .np_utils import fix_column_name
from .info import TableInfo
from .index import Index, _IndexModeContext, get_index
from .connect import TableRead, TableWrite
from . import conf


_implementation_notes = """
This string has informal notes concerning Table implementation for developers.

Things to remember:

- Table has customizable attributes ColumnClass, Column, MaskedColumn.
  Table.Column is normally just column.Column (same w/ MaskedColumn)
  but in theory they can be different.  Table.ColumnClass is the default
  class used to create new non-mixin columns, and this is a function of
  the Table.masked attribute.  Column creation / manipulation in a Table
  needs to respect these.

- Column objects that get inserted into the Table.columns attribute must
  have the info.parent_table attribute set correctly.  Beware just dropping
  an object into the columns dict since an existing column may
  be part of another Table and have parent_table set to point at that
  table.  Dropping that column into `columns` of this Table will cause
  a problem for the old one so the column object needs to be copied (but
  not necessarily the data).

  Currently replace_column is always making a copy of both object and
  data if parent_table is set.  This could be improved but requires a
  generic way to copy a mixin object but not the data.

- Be aware of column objects that have indices set.

- `cls.ColumnClass` is a property that effectively uses the `masked` attribute
  to choose either `cls.Column` or `cls.MaskedColumn`.
"""

__doctest_skip__ = ['Table.read', 'Table.write', 'Table._read',
                    'Table.convert_bytestring_to_unicode',
                    'Table.convert_unicode_to_bytestring',
                    ]

__doctest_requires__ = {'*pandas': ['pandas>=1.1']}

_pprint_docs = """
    {__doc__}

    Parameters
    ----------
    max_lines : int or `None`
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

_pformat_docs = """
    {__doc__}

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


def _get_names_from_list_of_dict(rows):
    """Return list of column names if ``rows`` is a list of dict that
    defines table data.

    If rows is not a list of dict then return None.
    """
    if rows is None:
        return None

    names = set()
    for row in rows:
        if not isinstance(row, dict):
            return None
        names.update(row)
    return list(names)


# Note to future maintainers: when transitioning this to dict
# be sure to change the OrderedDict ref(s) in Row and in __len__().

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
            return list(self.values())[item]
        elif (isinstance(item, np.ndarray) and item.shape == () and item.dtype.kind == 'i'):
            return list(self.values())[item.item()]
        elif isinstance(item, tuple):
            return self.__class__([self[x] for x in item])
        elif isinstance(item, slice):
            return self.__class__([self[x] for x in list(self)[item]])
        else:
            raise IndexError('Illegal key or index value for {} object'
                             .format(self.__class__.__name__))

    def __setitem__(self, item, value, validated=False):
        """
        Set item in this dict instance, but do not allow directly replacing an
        existing column unless it is already validated (and thus is certain to
        not corrupt the table).

        NOTE: it is easily possible to corrupt a table by directly *adding* a new
        key to the TableColumns attribute of a Table, e.g.
        ``t.columns['jane'] = 'doe'``.

        """
        if item in self and not validated:
            raise ValueError("Cannot replace column '{}'.  Use Table.replace_column() instead."
                             .format(item))
        super().__setitem__(item, value)

    def __repr__(self):
        names = (f"'{x}'" for x in self.keys())
        return "<{1} names=({0})>".format(",".join(names), self.__class__.__name__)

    def _rename_column(self, name, new_name):
        if name == new_name:
            return

        if new_name in self:
            raise KeyError(f"Column {new_name} already exists")

        mapper = {name: new_name}
        new_names = [mapper.get(name, name) for name in self]
        cols = list(self.values())
        self.clear()
        self.update(list(zip(new_names, cols)))

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


class TableReadWrite:
    def __get__(self, instance, owner_cls):
        if instance is None:
            # This is an unbound descriptor on the class
            info = self
            info._parent_cls = owner_cls
        else:
            info = instance.__dict__.get('info')
            if info is None:
                info = instance.__dict__['info'] = self.__class__(bound=True)
            info._parent = instance
        return info


class Table:
    """A class to represent tables of heterogeneous data.

    `~astropy.table.Table` provides a class for heterogeneous tabular data.
    A key enhancement provided by the `~astropy.table.Table` class over
    e.g. a `numpy` structured array is the ability to easily modify the
    structure of the table by adding or removing columns, or adding new
    rows of data.  In addition table and column metadata are fully supported.

    `~astropy.table.Table` differs from `~astropy.nddata.NDData` by the
    assumption that the input data consists of columns of homogeneous data,
    where each column has a unique identifier and may contain additional
    metadata such as the data unit, format, and description.

    See also: https://docs.astropy.org/en/stable/table/

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
    units : list, dict, optional
        List or dict of units to apply to columns.
    descriptions : list, dict, optional
        List or dict of descriptions to apply to columns.
    **kwargs : dict, optional
        Additional keyword args when converting table-like object.
    """

    meta = MetaData(copy=False)

    # Define class attributes for core container objects to allow for subclass
    # customization.
    Row = Row
    Column = Column
    MaskedColumn = MaskedColumn
    TableColumns = TableColumns
    TableFormatter = TableFormatter

    # Unified I/O read and write methods from .connect
    read = UnifiedReadWriteMethod(TableRead)
    write = UnifiedReadWriteMethod(TableWrite)

    def as_array(self, keep_byteorder=False, names=None):
        """
        Return a new copy of the table in the form of a structured np.ndarray or
        np.ma.MaskedArray object (as appropriate).

        Parameters
        ----------
        keep_byteorder : bool, optional
            By default the returned array has all columns in native byte
            order.  However, if this option is `True` this preserves the
            byte order of all columns (if any are non-native).

        names : list, optional:
            List of column names to include for returned structured array.
            Default is to include all table columns.

        Returns
        -------
        table_array : np.ndarray (unmasked) or np.ma.MaskedArray (masked)
            Copy of table as a numpy structured array
        """
        masked = self.masked or self.has_masked_columns or self.has_masked_values
        empty_init = ma.empty if masked else np.empty
        if len(self.columns) == 0:
            return empty_init(0, dtype=None)

        dtype = []

        cols = self.columns.values()

        if names is not None:
            cols = [col for col in cols if col.info.name in names]

        for col in cols:
            col_descr = descr(col)

            if not (col.info.dtype.isnative or keep_byteorder):
                new_dt = np.dtype(col_descr[1]).newbyteorder('=')
                col_descr = (col_descr[0], new_dt, col_descr[2])

            dtype.append(col_descr)

        data = empty_init(len(self), dtype=dtype)
        for col in cols:
            # When assigning from one array into a field of a structured array,
            # Numpy will automatically swap those columns to their destination
            # byte order where applicable
            data[col.info.name] = col

            # For masked out, masked mixin columns need to set output mask attribute.
            if masked and has_info_class(col, MixinInfo) and hasattr(col, 'mask'):
                data[col.info.name].mask = col.mask

        return data

    def __init__(self, data=None, masked=False, names=None, dtype=None,
                 meta=None, copy=True, rows=None, copy_indices=True,
                 units=None, descriptions=None,
                 **kwargs):

        # Set up a placeholder empty table
        self._set_masked(masked)
        self.columns = self.TableColumns()
        self.formatter = self.TableFormatter()
        self._copy_indices = True  # copy indices from this Table by default
        self._init_indices = copy_indices  # whether to copy indices in init
        self.primary_key = None

        # Must copy if dtype are changing
        if not copy and dtype is not None:
            raise ValueError('Cannot specify dtype when copy=False')

        # Specifies list of names found for the case of initializing table with
        # a list of dict. If data are not list of dict then this is None.
        names_from_list_of_dict = None

        # Row-oriented input, e.g. list of lists or list of tuples, list of
        # dict, Row instance.  Set data to something that the subsequent code
        # will parse correctly.
        if rows is not None:
            if data is not None:
                raise ValueError('Cannot supply both `data` and `rows` values')
            if isinstance(rows, types.GeneratorType):
                # Without this then the all(..) test below uses up the generator
                rows = list(rows)

            # Get column names if `rows` is a list of dict, otherwise this is None
            names_from_list_of_dict = _get_names_from_list_of_dict(rows)
            if names_from_list_of_dict:
                data = rows
            elif isinstance(rows, self.Row):
                data = rows
            else:
                data = list(zip(*rows))

        # Infer the type of the input data and set up the initialization
        # function, number of columns, and potentially the default col names

        default_names = None

        # Handle custom (subclass) table attributes that are stored in meta.
        # These are defined as class attributes using the TableAttribute
        # descriptor.  Any such attributes get removed from kwargs here and
        # stored for use after the table is otherwise initialized. Any values
        # provided via kwargs will have precedence over existing values from
        # meta (e.g. from data as a Table or meta via kwargs).
        meta_table_attrs = {}
        if kwargs:
            for attr in list(kwargs):
                descr = getattr(self.__class__, attr, None)
                if isinstance(descr, TableAttribute):
                    meta_table_attrs[attr] = kwargs.pop(attr)

        if hasattr(data, '__astropy_table__'):
            # Data object implements the __astropy_table__ interface method.
            # Calling that method returns an appropriate instance of
            # self.__class__ and respects the `copy` arg.  The returned
            # Table object should NOT then be copied.
            data = data.__astropy_table__(self.__class__, copy, **kwargs)
            copy = False
        elif kwargs:
            raise TypeError('__init__() got unexpected keyword argument {!r}'
                            .format(list(kwargs.keys())[0]))

        if (isinstance(data, np.ndarray)
                and data.shape == (0,)
                and not data.dtype.names):
            data = None

        if isinstance(data, self.Row):
            data = data._table[data._index:data._index + 1]

        if isinstance(data, (list, tuple)):
            # Get column names from `data` if it is a list of dict, otherwise this is None.
            # This might be previously defined if `rows` was supplied as an init arg.
            names_from_list_of_dict = (names_from_list_of_dict
                                       or _get_names_from_list_of_dict(data))
            if names_from_list_of_dict:
                init_func = self._init_from_list_of_dicts
                n_cols = len(names_from_list_of_dict)
            else:
                init_func = self._init_from_list
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
            # If user-input meta is None then use data.meta (if non-trivial)
            if meta is None and data.meta:
                # At this point do NOT deepcopy data.meta as this will happen after
                # table init_func() is called.  But for table input the table meta
                # gets a key copy here if copy=False because later a direct object ref
                # is used.
                meta = data.meta if copy else data.meta.copy()

            # Handle indices on input table. Copy primary key and don't copy indices
            # if the input Table is in non-copy mode.
            self.primary_key = data.primary_key
            self._init_indices = self._init_indices and data._copy_indices

            # Extract default names, n_cols, and then overwrite ``data`` to be the
            # table columns so we can use _init_from_list.
            default_names = data.colnames
            n_cols = len(default_names)
            data = list(data.columns.values())

            init_func = self._init_from_list

        elif data is None:
            if names is None:
                if dtype is None:
                    # Table was initialized as `t = Table()`. Set up for empty
                    # table with names=[], data=[], and n_cols=0.
                    # self._init_from_list() will simply return, giving the
                    # expected empty table.
                    names = []
                else:
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
            raise ValueError('Data type {} not allowed to init Table'
                             .format(type(data)))

        # Set up defaults if names and/or dtype are not specified.
        # A value of None means the actual value will be inferred
        # within the appropriate initialization routine, either from
        # existing specification or auto-generated.

        if dtype is None:
            dtype = [None] * n_cols
        elif isinstance(dtype, np.dtype):
            if default_names is None:
                default_names = dtype.names
            # Convert a numpy dtype input to a list of dtypes for later use.
            dtype = [dtype[name] for name in dtype.names]

        if names is None:
            names = default_names or [None] * n_cols

        # Numpy does not support bytes column names on Python 3, so fix them
        # up now.
        names = [fix_column_name(name) for name in names]

        self._check_names_dtype(names, dtype, n_cols)

        # Finally do the real initialization
        init_func(data, names, dtype, n_cols, copy)

        # Set table meta.  If copy=True then deepcopy meta otherwise use the
        # user-supplied meta directly.
        if meta is not None:
            self.meta = deepcopy(meta) if copy else meta

        # Update meta with TableAttributes supplied as kwargs in Table init.
        # This takes precedence over previously-defined meta.
        if meta_table_attrs:
            for attr, value in meta_table_attrs.items():
                setattr(self, attr, value)

        # Whatever happens above, the masked property should be set to a boolean
        if self.masked not in (None, True, False):
            raise TypeError("masked property must be None, True or False")

        self._set_column_attribute('unit', units)
        self._set_column_attribute('description', descriptions)

    def _set_column_attribute(self, attr, values):
        """Set ``attr`` for columns to ``values``, which can be either a dict (keyed by column
        name) or a dict of name: value pairs.  This is used for handling the ``units`` and
        ``descriptions`` kwargs to ``__init__``.
        """
        if not values:
            return

        if isinstance(values, Row):
            # For a Row object transform to an equivalent dict.
            values = {name: values[name] for name in values.colnames}

        if not isinstance(values, dict):
            # If not a dict map, assume iterable and map to dict if the right length
            if len(values) != len(self.columns):
                raise ValueError(f'sequence of {attr} values must match number of columns')
            values = dict(zip(self.colnames, values))

        for name, value in values.items():
            if name not in self.columns:
                raise ValueError(f'invalid column name {name} for setting {attr} attribute')

            # Special case: ignore unit if it is an empty or blank string
            if attr == 'unit' and isinstance(value, str):
                if value.strip() == '':
                    value = None

            if value not in (np.ma.masked, None):
                setattr(self[name].info, attr, value)

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
        if self.masked or self.has_masked_columns or self.has_masked_values:
            mask_table = Table([getattr(col, 'mask', FalseArray(col.shape))
                                for col in self.itercols()],
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
        """Return copy of self, with masked values filled.

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
        if self.masked or self.has_masked_columns or self.has_masked_values:
            # Get new columns with masked values filled, then create Table with those
            # new cols (copy=False) but deepcopy the meta.
            data = [col.filled(fill_value) if hasattr(col, 'filled') else col
                    for col in self.itercols()]
            return self.__class__(data, meta=deepcopy(self.meta), copy=False)
        else:
            # Return copy of the original object.
            return self.copy()

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
            and SCEngine. If the supplied argument is None
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
                raise ValueError('Cannot create an index on column "{}", of '
                                 'type "{}"'.format(col.info.name, type(col)))

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

        out = self.as_array()
        return out.data if isinstance(out, np.ma.MaskedArray) else out

    def _check_names_dtype(self, names, dtype, n_cols):
        """Make sure that names and dtype are both iterable and have
        the same length as data.
        """
        for inp_list, inp_str in ((dtype, 'dtype'), (names, 'names')):
            if not isiterable(inp_list):
                raise ValueError(f'{inp_str} must be a list or None')

        if len(names) != n_cols or len(dtype) != n_cols:
            raise ValueError(
                'Arguments "names" and "dtype" must match number of columns'
                .format(inp_str))

    def _init_from_list_of_dicts(self, data, names, dtype, n_cols, copy):
        """Initialize table from a list of dictionaries representing rows."""
        # Define placeholder for missing values as a unique object that cannot
        # every occur in user data.
        MISSING = object()

        # Gather column names that exist in the input `data`.
        names_from_data = set()
        for row in data:
            names_from_data.update(row)

        if set(data[0].keys()) == names_from_data:
            names_from_data = list(data[0].keys())
        else:
            names_from_data = sorted(names_from_data)

        # Note: if set(data[0].keys()) != names_from_data, this will give an
        # exception later, so NO need to catch here.

        # Convert list of dict into dict of list (cols), keep track of missing
        # indexes and put in MISSING placeholders in the `cols` lists.
        cols = {}
        missing_indexes = defaultdict(list)
        for name in names_from_data:
            cols[name] = []
            for ii, row in enumerate(data):
                try:
                    val = row[name]
                except KeyError:
                    missing_indexes[name].append(ii)
                    val = MISSING
                cols[name].append(val)

        # Fill the missing entries with first values
        if missing_indexes:
            for name, indexes in missing_indexes.items():
                col = cols[name]
                first_val = next(val for val in col if val is not MISSING)
                for index in indexes:
                    col[index] = first_val

        # prepare initialization
        if all(name is None for name in names):
            names = names_from_data

        self._init_from_dict(cols, names, dtype, n_cols, copy)

        # Mask the missing values if necessary, converting columns to MaskedColumn
        # as needed.
        if missing_indexes:
            for name, indexes in missing_indexes.items():
                col = self[name]
                # Ensure that any Column subclasses with MISSING values can support
                # setting masked values. As of astropy 4.0 the test condition below is
                # always True since _init_from_dict cannot result in mixin columns.
                if isinstance(col, Column) and not isinstance(col, MaskedColumn):
                    self[name] = self.MaskedColumn(col, copy=False)

                # Finally do the masking in a mixin-safe way.
                self[name][indexes] = np.ma.masked
        return

    def _init_from_list(self, data, names, dtype, n_cols, copy):
        """Initialize table from a list of column data.  A column can be a
        Column object, np.ndarray, mixin, or any other iterable object.
        """
        # Special case of initializing an empty table like `t = Table()`. No
        # action required at this point.
        if n_cols == 0:
            return

        cols = []
        default_names = _auto_names(n_cols)

        for col, name, default_name, dtype in zip(data, names, default_names, dtype):
            col = self._convert_data_to_col(col, copy, default_name, dtype, name)

            cols.append(col)

        self._init_from_cols(cols)

    def _convert_data_to_col(self, data, copy=True, default_name=None, dtype=None, name=None):
        """
        Convert any allowed sequence data ``col`` to a column object that can be used
        directly in the self.columns dict.  This could be a Column, MaskedColumn,
        or mixin column.

        The final column name is determined by::

            name or data.info.name or def_name

        If ``data`` has no ``info`` then ``name = name or def_name``.

        The behavior of ``copy`` for Column objects is:
        - copy=True: new class instance with a copy of data and deep copy of meta
        - copy=False: new class instance with same data and a key-only copy of meta

        For mixin columns:
        - copy=True: new class instance with copy of data and deep copy of meta
        - copy=False: original instance (no copy at all)

        Parameters
        ----------
        data : object (column-like sequence)
            Input column data
        copy : bool
            Make a copy
        default_name : str
            Default name
        dtype : np.dtype or None
            Data dtype
        name : str or None
            Column name

        Returns
        -------
        col : Column, MaskedColumn, mixin-column type
            Object that can be used as a column in self
        """
        data_is_mixin = self._is_mixin_for_table(data)
        masked_col_cls = (self.ColumnClass
                          if issubclass(self.ColumnClass, self.MaskedColumn)
                          else self.MaskedColumn)

        try:
            data0_is_mixin = self._is_mixin_for_table(data[0])
        except Exception:
            # Need broad exception, cannot predict what data[0] raises for arbitrary data
            data0_is_mixin = False

        # Structured ndarray gets viewed as a mixin unless already a valid
        # mixin class
        if (not isinstance(data, Column) and not data_is_mixin
                and isinstance(data, np.ndarray) and len(data.dtype) > 1):
            data = data.view(NdarrayMixin)
            data_is_mixin = True

        # Get the final column name using precedence.  Some objects may not
        # have an info attribute.
        if not name:
            if hasattr(data, 'info'):
                name = data.info.name or default_name
            else:
                name = default_name

        if isinstance(data, Column):
            # If self.ColumnClass is a subclass of col, then "upgrade" to ColumnClass,
            # otherwise just use the original class.  The most common case is a
            # table with masked=True and ColumnClass=MaskedColumn.  Then a Column
            # gets upgraded to MaskedColumn, but the converse (pre-4.0) behavior
            # of downgrading from MaskedColumn to Column (for non-masked table)
            # does not happen.
            col_cls = self._get_col_cls_for_table(data)

        elif data_is_mixin:
            # Copy the mixin column attributes if they exist since the copy below
            # may not get this attribute.
            col = col_copy(data, copy_indices=self._init_indices) if copy else data
            col.info.name = name
            return col

        elif data0_is_mixin:
            # Handle case of a sequence of a mixin, e.g. [1*u.m, 2*u.m].
            try:
                col = data[0].__class__(data)
                col.info.name = name
                col.info.indices = []
                return col
            except Exception:
                # If that didn't work for some reason, just turn it into np.array of object
                data = np.array(data, dtype=object)
                col_cls = self.ColumnClass

        elif isinstance(data, np.ma.MaskedArray):
            # Require that col_cls be a subclass of MaskedColumn, remembering
            # that ColumnClass could be a user-defined subclass (though more-likely
            # could be MaskedColumn).
            col_cls = masked_col_cls

        elif data is None:
            # Special case for data passed as the None object (for broadcasting
            # to an object column). Need to turn data into numpy `None` scalar
            # object, otherwise `Column` interprets data=None as no data instead
            # of a object column of `None`.
            data = np.array(None)
            col_cls = self.ColumnClass

        elif not hasattr(data, 'dtype'):
            # `data` is none of the above, convert to numpy array or MaskedArray
            # assuming only that it is a scalar or sequence or N-d nested
            # sequence. This function is relatively intricate and tries to
            # maintain performance for common cases while handling things like
            # list input with embedded np.ma.masked entries. If `data` is a
            # scalar then it gets returned unchanged so the original object gets
            # passed to `Column` later.
            data = _convert_sequence_data_to_array(data, dtype)
            copy = False  # Already made a copy above
            col_cls = masked_col_cls if isinstance(data, np.ma.MaskedArray) else self.ColumnClass

        else:
            col_cls = self.ColumnClass

        try:
            col = col_cls(name=name, data=data, dtype=dtype,
                          copy=copy, copy_indices=self._init_indices)
        except Exception:
            # Broad exception class since we don't know what might go wrong
            raise ValueError('unable to convert data to Column for Table')

        col = self._convert_col_for_table(col)

        return col

    def _init_from_ndarray(self, data, names, dtype, n_cols, copy):
        """Initialize table from an ndarray structured array"""

        data_names = data.dtype.names or _auto_names(n_cols)
        struct = data.dtype.names is not None
        names = [name or data_names[i] for i, name in enumerate(names)]

        cols = ([data[name] for name in data_names] if struct else
                [data[:, i] for i in range(n_cols)])

        self._init_from_list(cols, names, dtype, n_cols, copy)

    def _init_from_dict(self, data, names, dtype, n_cols, copy):
        """Initialize table from a dictionary of columns"""

        data_list = [data[name] for name in names]
        self._init_from_list(data_list, names, dtype, n_cols, copy)

    def _get_col_cls_for_table(self, col):
        """Get the correct column class to use for upgrading any Column-like object.

        For a masked table, ensure any Column-like object is a subclass
        of the table MaskedColumn.

        For unmasked table, ensure any MaskedColumn-like object is a subclass
        of the table MaskedColumn.  If not a MaskedColumn, then ensure that any
        Column-like object is a subclass of the table Column.
        """

        col_cls = col.__class__

        if self.masked:
            if isinstance(col, Column) and not isinstance(col, self.MaskedColumn):
                col_cls = self.MaskedColumn
        else:
            if isinstance(col, MaskedColumn):
                if not isinstance(col, self.MaskedColumn):
                    col_cls = self.MaskedColumn
            elif isinstance(col, Column) and not isinstance(col, self.Column):
                col_cls = self.Column

        return col_cls

    def _convert_col_for_table(self, col):
        """
        Make sure that all Column objects have correct base class for this type of
        Table.  For a base Table this most commonly means setting to
        MaskedColumn if the table is masked.  Table subclasses like QTable
        override this method.
        """
        if isinstance(col, Column) and not isinstance(col, self.ColumnClass):
            col_cls = self._get_col_cls_for_table(col)
            if col_cls is not col.__class__:
                col = col_cls(col, copy=False)

        return col

    def _init_from_cols(self, cols):
        """Initialize table from a list of Column or mixin objects"""

        lengths = set(len(col) for col in cols)
        if len(lengths) > 1:
            raise ValueError('Inconsistent data column lengths: {}'
                             .format(lengths))

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
        if self.meta:
            table.meta = self.meta.copy()  # Shallow copy for slice
        table.primary_key = self.primary_key

        newcols = []
        for col in self.columns.values():
            newcol = col[slice_]

            # Note in line below, use direct attribute access to col.indices for Column
            # instances instead of the generic col.info.indices.  This saves about 4 usec
            # per column.
            if (col if isinstance(col, Column) else col.info).indices:
                # TODO : as far as I can tell the only purpose of setting _copy_indices
                # here is to communicate that to the initial test in `slice_indices`.
                # Why isn't that just sent as an arg to the function?
                col.info._copy_indices = self._copy_indices
                newcol = col.info.slice_indices(newcol, slice_, len(col))

                # Don't understand why this is forcing a value on the original column.
                # Normally col.info does not even have a _copy_indices attribute.  Tests
                # still pass if this line is deleted.  (Each col.info attribute access
                # is expensive).
                col.info._copy_indices = True
            else:
                newcol.info.indices = []

            newcols.append(newcol)

        self._make_table_from_cols(table, newcols, verify=False, names=self.columns.keys())
        return table

    @staticmethod
    def _make_table_from_cols(table, cols, verify=True, names=None):
        """
        Make ``table`` in-place so that it represents the given list of ``cols``.
        """
        if names is None:
            names = [col.info.name for col in cols]

        # Note: we do not test for len(names) == len(cols) if names is not None.  In that
        # case the function is being called by from "trusted" source (e.g. right above here)
        # that is assumed to provide valid inputs.  In that case verify=False.

        if verify:
            if None in names:
                raise TypeError('Cannot have None for column name')
            if len(set(names)) != len(names):
                raise ValueError('Duplicate column names')

        table.columns = table.TableColumns((name, col) for name, col in zip(names, cols))

        for col in cols:
            table._set_col_parent_table_and_mask(col)

    def _set_col_parent_table_and_mask(self, col):
        """
        Set ``col.parent_table = self`` and force ``col`` to have ``mask``
        attribute if the table is masked and ``col.mask`` does not exist.
        """
        # For Column instances it is much faster to do direct attribute access
        # instead of going through .info
        col_info = col if isinstance(col, Column) else col.info
        col_info.parent_table = self

        # Legacy behavior for masked table
        if self.masked and not hasattr(col, 'mask'):
            col.mask = FalseArray(col.shape)

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
            descr_vals.append('length={}'.format(len(self)))

        descr = ' '.join(descr_vals)
        if html:
            from astropy.utils.xml.writer import xml_escape
            descr = '<i>{}</i>\n'.format(xml_escape(descr))
        else:
            descr = f'<{descr}>\n'

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

    @property
    def has_masked_columns(self):
        """True if table has any ``MaskedColumn`` columns.

        This does not check for mixin columns that may have masked values, use the
        ``has_masked_values`` property in that case.

        """
        return any(isinstance(col, MaskedColumn) for col in self.itercols())

    @property
    def has_masked_values(self):
        """True if column in the table has values which are masked.

        This may be relatively slow for large tables as it requires checking the mask
        values of each column.
        """
        for col in self.itercols():
            if hasattr(col, 'mask') and np.any(col.mask):
                return True
        else:
            return False

    def _is_mixin_for_table(self, col):
        """
        Determine if ``col`` should be added to the table directly as
        a mixin column.
        """
        if isinstance(col, BaseColumn):
            return False

        # Is it a mixin but not not Quantity (which gets converted to Column with
        # unit set).
        return has_info_class(col, MixinInfo) and not has_info_class(col, QuantityInfo)

    @format_doc(_pprint_docs)
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

        """
        lines, outs = self.formatter._pformat_table(self, max_lines, max_width,
                                                    show_name=show_name, show_unit=show_unit,
                                                    show_dtype=show_dtype, align=align)
        if outs['show_length']:
            lines.append('Length = {} rows'.format(len(self)))

        n_header = outs['n_header']

        for i, line in enumerate(lines):
            if i < n_header:
                color_print(line, 'red')
            else:
                print(line)

    @format_doc(_pprint_docs)
    def pprint_all(self, max_lines=-1, max_width=-1, show_name=True,
                   show_unit=None, show_dtype=False, align=None):
        """Print a formatted string representation of the entire table.

        This method is the same as `astropy.table.Table.pprint` except that
        the default ``max_lines`` and ``max_width`` are both -1 so that by
        default the entire table is printed instead of restricting to the size
        of the screen terminal.

        """
        return self.pprint(max_lines, max_width, show_name,
                           show_unit, show_dtype, align)

    def _make_index_row_display_table(self, index_row_name):
        if index_row_name not in self.columns:
            idx_col = self.ColumnClass(name=index_row_name, data=np.arange(len(self)))
            return self.__class__([idx_col] + list(self.columns.values()),
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
            notebook.  See `this page <https://getbootstrap.com/css/#tables>`_
            for the list of classes.
        css : str
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
            tableid = 'table{}-{}'.format(id(self),
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
        css : str
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
            log.error(f"Browser '{browser}' not found.")
        else:
            br.open(urljoin('file:', pathname2url(path)))

    @format_doc(_pformat_docs, id="{id}")
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

        """

        lines, outs = self.formatter._pformat_table(
            self, max_lines, max_width, show_name=show_name,
            show_unit=show_unit, show_dtype=show_dtype, html=html,
            tableid=tableid, tableclass=tableclass, align=align)

        if outs['show_length']:
            lines.append('Length = {} rows'.format(len(self)))

        return lines

    @format_doc(_pformat_docs, id="{id}")
    def pformat_all(self, max_lines=-1, max_width=-1, show_name=True,
                    show_unit=None, show_dtype=False, html=False, tableid=None,
                    align=None, tableclass=None):
        """Return a list of lines for the formatted string representation of
        the entire table.

        If no value of ``max_lines`` is supplied then the height of the
        screen terminal is used to set ``max_lines``.  If the terminal
        height cannot be determined then the default is taken from the
        configuration item ``astropy.conf.max_lines``.  If a negative
        value of ``max_lines`` is supplied then there is no line limit
        applied.

        The same applies for ``max_width`` except the configuration item  is
        ``astropy.conf.max_width``.

        """

        return self.pformat(max_lines, max_width, show_name,
                            show_unit, show_dtype, html, tableid,
                            align, tableclass)

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
                                 copy_indices=self._copy_indices)
            out._groups = groups.TableGroups(out, indices=self.groups._indices,
                                             keys=self.groups._keys)
            out.meta = self.meta.copy()  # Shallow copy for meta
            return out
        elif ((isinstance(item, np.ndarray) and item.size == 0)
              or (isinstance(item, (tuple, list)) and not item)):
            # If item is an empty array/list/tuple then return the table with no rows
            return self._new_from_slice([])
        elif (isinstance(item, slice)
              or isinstance(item, np.ndarray)
              or isinstance(item, list)
              or isinstance(item, tuple) and all(isinstance(x, np.ndarray)
                                                 for x in item)):
            # here for the many ways to give a slice; a tuple of ndarray
            # is produced by np.where, as in t[np.where(t['a'] > 2)]
            # For all, a new table is constructed with slice of all columns
            return self._new_from_slice(item)
        else:
            raise ValueError('Illegal type {} for table item access'
                             .format(type(item)))

    def __setitem__(self, item, value):
        # If the item is a string then it must be the name of a column.
        # If that column doesn't already exist then create it now.
        if isinstance(item, str) and item not in self.colnames:
            self.add_column(value, name=item, copy=True)

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

            elif (isinstance(item, slice)
                  or isinstance(item, np.ndarray)
                  or isinstance(item, list)
                  or (isinstance(item, tuple)  # output from np.where
                      and all(isinstance(x, np.ndarray) for x in item))):

                if isinstance(value, Table):
                    vals = (col for col in value.columns.values())

                elif isinstance(value, np.ndarray) and value.dtype.names:
                    vals = (value[name] for name in value.dtype.names)

                elif np.isscalar(value):
                    vals = itertools.repeat(value, n_cols)

                else:  # Assume this is an iterable that will work
                    if len(value) != n_cols:
                        raise ValueError('Right side value needs {} elements (one for each column)'
                                         .format(n_cols))
                    vals = value

                for col, val in zip(self.columns.values(), vals):
                    col[item] = val

            else:
                raise ValueError('Illegal type {} for table item access'
                                 .format(type(item)))

    def __delitem__(self, item):
        if isinstance(item, str):
            self.remove_column(item)
        elif isinstance(item, (int, np.integer)):
            self.remove_row(item)
        elif (isinstance(item, (list, tuple, np.ndarray))
              and all(isinstance(x, str) for x in item)):
            self.remove_columns(item)
        elif (isinstance(item, (list, np.ndarray))
              and np.asarray(item).dtype.kind == 'i'):
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
        if masked in [True, False, None]:
            self._masked = masked
        else:
            raise ValueError("masked should be one of True, False, None")

        self._column_class = self.MaskedColumn if self._masked else self.Column

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
        return (isinstance(names, (tuple, list)) and names
                and all(isinstance(x, str) for x in names))

    def keys(self):
        return list(self.columns.keys())

    def values(self):
        return self.columns.values()

    def items(self):
        return self.columns.items()

    def __len__(self):
        # For performance reasons (esp. in Row) cache the first column name
        # and use that subsequently for the table length.  If might not be
        # available yet or the column might be gone now, in which case
        # try again in the except block.
        try:
            return len(OrderedDict.__getitem__(self.columns, self._first_colname))
        except (AttributeError, KeyError):
            if len(self.columns) == 0:
                return 0

            # Get the first column name
            self._first_colname = next(iter(self.columns))
            return len(self.columns[self._first_colname])

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
            raise ValueError(f"Column {name} does not exist")

    def add_column(self, col, index=None, name=None, rename_duplicate=False, copy=True,
                   default_name=None):
        """
        Add a new column to the table using ``col`` as input.  If ``index``
        is supplied then insert column before ``index`` position
        in the list of columns, otherwise append column to the end
        of the list.

        The ``col`` input can be any data object which is acceptable as a
        `~astropy.table.Table` column object or can be converted.  This includes
        mixin columns and scalar or length=1 objects which get broadcast to match
        the table length.

        To add several columns at once use ``add_columns()`` or simply call
        ``add_column()`` for each one.  There is very little performance difference
        in the two approaches.

        Parameters
        ----------
        col : object
            Data object for the new column
        index : int or `None`
            Insert column before this position or at end (default).
        name : str
            Column name
        rename_duplicate : bool
            Uniquify column name if it already exist. Default is False.
        copy : bool
            Make a copy of the new column. Default is True.
        default_name : str or `None`
            Name to use if both ``name`` and ``col.info.name`` are not available.
            Defaults to ``col{number_of_columns}``.

        Examples
        --------
        Create a table with two columns 'a' and 'b', then create a third column 'c'
        and append it to the end of the table::

            >>> t = Table([[1, 2], [0.1, 0.2]], names=('a', 'b'))
            >>> col_c = Column(name='c', data=['x', 'y'])
            >>> t.add_column(col_c)
            >>> print(t)
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y

        Add column 'd' at position 1. Note that the column is inserted
        before the given index::

            >>> t.add_column(['a', 'b'], name='d', index=1)
            >>> print(t)
             a   d   b   c
            --- --- --- ---
              1   a 0.1   x
              2   b 0.2   y

        Add second column named 'b' with rename_duplicate::

            >>> t = Table([[1, 2], [0.1, 0.2]], names=('a', 'b'))
            >>> t.add_column(1.1, name='b', rename_duplicate=True)
            >>> print(t)
             a   b  b_1
            --- --- ---
              1 0.1 1.1
              2 0.2 1.1

        Add an unnamed column or mixin object in the table using a default name
        or by specifying an explicit name with ``name``. Name can also be overridden::

            >>> t = Table([[1, 2], [0.1, 0.2]], names=('a', 'b'))
            >>> t.add_column(['a', 'b'])
            >>> t.add_column(col_c, name='d')
            >>> print(t)
             a   b  col2  d
            --- --- ---- ---
              1 0.1    a   x
              2 0.2    b   y
        """
        if default_name is None:
            default_name = 'col{}'.format(len(self.columns))

        # Convert col data to acceptable object for insertion into self.columns.
        # Note that along with the lines above and below, this allows broadcasting
        # of scalars to the correct shape for adding to table.
        col = self._convert_data_to_col(col, name=name, copy=copy,
                                        default_name=default_name)

        # Assigning a scalar column to an empty table should result in an
        # exception (see #3811).
        if col.shape == () and len(self) == 0:
            raise TypeError('Empty table cannot have column set to scalar value')
        # Make col data shape correct for scalars.  The second test is to allow
        # broadcasting an N-d element to a column, e.g. t['new'] = [[1, 2]].
        elif (col.shape == () or col.shape[0] == 1) and len(self) > 0:
            new_shape = (len(self),) + getattr(col, 'shape', ())[1:]
            if isinstance(col, np.ndarray):
                col = np.broadcast_to(col, shape=new_shape,
                                      subok=True)
            elif isinstance(col, ShapedLikeNDArray):
                col = col._apply(np.broadcast_to, shape=new_shape,
                                 subok=True)

            # broadcast_to() results in a read-only array.  Apparently it only changes
            # the view to look like the broadcasted array.  So copy.
            col = col_copy(col)

        name = col.info.name

        # Ensure that new column is the right length
        if len(self.columns) > 0 and len(col) != len(self):
            raise ValueError('Inconsistent data column lengths')

        if rename_duplicate:
            orig_name = name
            i = 1
            while name in self.columns:
                # Iterate until a unique name is found
                name = orig_name + '_' + str(i)
                i += 1
            col.info.name = name

        # Set col parent_table weakref and ensure col has mask attribute if table.masked
        self._set_col_parent_table_and_mask(col)

        # Add new column as last column
        self.columns[name] = col

        if index is not None:
            # Move the other cols to the right of the new one
            move_names = self.colnames[index:-1]
            for move_name in move_names:
                self.columns.move_to_end(move_name, last=True)

    def add_columns(self, cols, indexes=None, names=None, copy=True, rename_duplicate=False):
        """
        Add a list of new columns the table using ``cols`` data objects.  If a
        corresponding list of ``indexes`` is supplied then insert column
        before each ``index`` position in the *original* list of columns,
        otherwise append columns to the end of the list.

        The ``cols`` input can include any data objects which are acceptable as
        `~astropy.table.Table` column objects or can be converted.  This includes
        mixin columns and scalar or length=1 objects which get broadcast to match
        the table length.

        From a performance perspective there is little difference between calling
        this method once or looping over the new columns and calling ``add_column()``
        for each column.

        Parameters
        ----------
        cols : list of objects
            List of data objects for the new columns
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
        Create a table with two columns 'a' and 'b', then create columns 'c' and 'd'
        and append them to the end of the table::

            >>> t = Table([[1, 2], [0.1, 0.2]], names=('a', 'b'))
            >>> col_c = Column(name='c', data=['x', 'y'])
            >>> col_d = Column(name='d', data=['u', 'v'])
            >>> t.add_columns([col_c, col_d])
            >>> print(t)
             a   b   c   d
            --- --- --- ---
              1 0.1   x   u
              2 0.2   y   v

        Add column 'c' at position 0 and column 'd' at position 1. Note that
        the columns are inserted before the given position::

            >>> t = Table([[1, 2], [0.1, 0.2]], names=('a', 'b'))
            >>> t.add_columns([['x', 'y'], ['u', 'v']], names=['c', 'd'],
            ...               indexes=[0, 1])
            >>> print(t)
             c   a   d   b
            --- --- --- ---
              x   1   u 0.1
              y   2   v 0.2

        Add second column 'b' and column 'c' with ``rename_duplicate``::

            >>> t = Table([[1, 2], [0.1, 0.2]], names=('a', 'b'))
            >>> t.add_columns([[1.1, 1.2], ['x', 'y']], names=('b', 'c'),
            ...               rename_duplicate=True)
            >>> print(t)
             a   b  b_1  c
            --- --- --- ---
              1 0.1 1.1  x
              2 0.2 1.2  y

        Add unnamed columns or mixin objects in the table using default names
        or by specifying explicit names with ``names``. Names can also be overridden::

            >>> t = Table()
            >>> col_b = Column(name='b', data=['u', 'v'])
            >>> t.add_columns([[1, 2], col_b])
            >>> t.add_columns([[3, 4], col_b], names=['c', 'd'])
            >>> print(t)
            col0  b   c   d
            ---- --- --- ---
               1   u   3   u
               2   v   4   v
        """
        if indexes is None:
            indexes = [len(self.columns)] * len(cols)
        elif len(indexes) != len(cols):
            raise ValueError('Number of indexes must match number of cols')

        if names is None:
            names = (None,) * len(cols)
        elif len(names) != len(cols):
            raise ValueError('Number of names must match number of cols')

        default_names = ['col{}'.format(ii + len(self.columns))
                         for ii in range(len(cols))]

        for ii in reversed(np.argsort(indexes)):
            self.add_column(cols[ii], index=indexes[ii], name=names[ii],
                            default_name=default_names[ii],
                            rename_duplicate=rename_duplicate, copy=copy)

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
            warnings.warn(f"replaced column '{name}'",
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

    def replace_column(self, name, col, copy=True):
        """
        Replace column ``name`` with the new ``col`` object.

        The behavior of ``copy`` for Column objects is:
        - copy=True: new class instance with a copy of data and deep copy of meta
        - copy=False: new class instance with same data and a key-only copy of meta

        For mixin columns:
        - copy=True: new class instance with copy of data and deep copy of meta
        - copy=False: original instance (no copy at all)

        Parameters
        ----------
        name : str
            Name of column to replace
        col : column object (list, ndarray, Column, etc)
            New column object to replace the existing column
        copy : bool
            Make copy of the input ``col``, default=True

        Examples
        --------
        Replace column 'a' with a float version of itself::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> float_a = t['a'].astype(float)
            >>> t.replace_column('a', float_a)
        """
        if name not in self.colnames:
            raise ValueError(f'column name {name} is not in the table')

        if self[name].info.indices:
            raise ValueError('cannot replace a table index column')

        col = self._convert_data_to_col(col, name=name, copy=copy)
        self._set_col_parent_table_and_mask(col)

        # Ensure that new column is the right length, unless it is the only column
        # in which case re-sizing is allowed.
        if len(self.columns) > 1 and len(col) != len(self[name]):
            raise ValueError('length of new column must match table length')

        self.columns.__setitem__(name, col, validated=True)

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

    def iterrows(self, *names):
        """
        Iterate over rows of table returning a tuple of values for each row.

        This method is especially useful when only a subset of columns are needed.

        The ``iterrows`` method can be substantially faster than using the standard
        Table row iteration (e.g. ``for row in tbl:``), since that returns a new
        ``~astropy.table.Row`` object for each row and accessing a column in that
        row (e.g. ``row['col0']``) is slower than tuple access.

        Parameters
        ----------
        names : list
            List of column names (default to all columns if no names provided)

        Returns
        -------
        rows : iterator returning tuples of row values

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table({'a': [1, 2, 3],
            ...            'b': [1.0, 2.5, 3.0],
            ...            'c': ['x', 'y', 'z']})

        To iterate row-wise using column names::

            >>> for a, c in t.iterrows('a', 'c'):
            ...     print(a, c)
            1 x
            2 y
            3 z

        """
        if len(names) == 0:
            names = self.colnames
        else:
            for name in names:
                if name not in self.colnames:
                    raise ValueError(f'{name} is not a valid column name')

        cols = (self[name] for name in names)
        out = zip(*cols)
        return out

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
                raise KeyError(f"Column {name} does not exist")

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
                raise KeyError(f"Column {name} does not exist")

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
            raise KeyError(f"Column {name} does not exist")

        self.columns[name].info.name = new_name

    def rename_columns(self, names, new_names):
        '''
        Rename multiple columns.

        Parameters
        ----------
        names : list, tuple
            A list or tuple of existing column names.
        new_names : list, tuple
            A list or tuple of new column names.

        Examples
        --------
        Create a table with three columns 'a', 'b', 'c'::

            >>> t = Table([[1,2],[3,4],[5,6]], names=('a','b','c'))
            >>> print(t)
              a   b   c
             --- --- ---
              1   3   5
              2   4   6

        Renaming columns 'a' to 'aa' and 'b' to 'bb'::

            >>> names = ('a','b')
            >>> new_names = ('aa','bb')
            >>> t.rename_columns(names, new_names)
            >>> print(t)
             aa  bb   c
            --- --- ---
              1   3   5
              2   4   6
        '''

        if not self._is_list_or_tuple_of_str(names):
            raise TypeError("input 'names' must be a tuple or a list of column names")

        if not self._is_list_or_tuple_of_str(new_names):
            raise TypeError("input 'new_names' must be a tuple or a list of column names")

        if len(names) != len(new_names):
            raise ValueError("input 'names' and 'new_names' list arguments must be the same length")

        for name, new_name in zip(names, new_names):
            self.rename_column(name, new_name)

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
            raise IndexError("Index {} is out of bounds for table with length {}"
                             .format(index, N))
        if index < 0:
            index += N

        def _is_mapping(obj):
            """Minimal checker for mapping (dict-like) interface for obj"""
            attrs = ('__getitem__', '__len__', '__iter__', 'keys', 'values', 'items')
            return all(hasattr(obj, attr) for attr in attrs)

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
                        raise ValueError(f"Value must be supplied for column '{name}'")

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
                # If new val is masked and the existing column does not support masking
                # then upgrade the column to a mask-enabled type: either the table-level
                # default ColumnClass or else MaskedColumn.
                if mask_ and isinstance(col, Column) and not isinstance(col, MaskedColumn):
                    col_cls = (self.ColumnClass
                               if issubclass(self.ColumnClass, self.MaskedColumn)
                               else self.MaskedColumn)
                    col = col_cls(col, copy=False)

                newcol = col.insert(index, val, axis=0)

                if len(newcol) != N + 1:
                    raise ValueError('Incorrect length for column {} after inserting {}'
                                     ' (expected {}, got {})'
                                     .format(name, val, len(newcol), N + 1))
                newcol.info.parent_table = self

                # Set mask if needed and possible
                if mask_:
                    if hasattr(newcol, 'mask'):
                        newcol[index] = np.ma.masked
                    else:
                        raise TypeError("mask was supplied for column '{}' but it does not "
                                        "support masked values".format(col.info.name))

                columns[name] = newcol

            # insert row in indices
            for table_index in self.indices:
                table_index.insert_row(index, vals, self.columns.values())

        except Exception as err:
            raise ValueError("Unable to insert row because of exception in column '{}':\n{}"
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

    def argsort(self, keys=None, kind=None, reverse=False):
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
        reverse : bool
            Sort in reverse order (default=False)

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
            index = get_index(self, names=keys)
            if index is not None:
                idx = np.asarray(index.sorted_data())
                return idx[::-1] if reverse else idx

        kwargs = {}
        if keys:
            # For multiple keys return a structured array which gets sorted,
            # while for a single key return a single ndarray.  Sorting a
            # one-column structured array is slower than ndarray (e.g. a
            # factor of ~6 for a 10 million long random array), and much slower
            # for in principle sortable columns like Time, which get stored as
            # object arrays.
            if len(keys) > 1:
                kwargs['order'] = keys
                data = self.as_array(names=keys)
            else:
                data = self[keys[0]]
        else:
            # No keys provided so sort on all columns.
            data = self.as_array()

        if kind:
            kwargs['kind'] = kind

        # np.argsort will look for a possible .argsort method (e.g., for Time),
        # and if that fails cast to an array and try sorting that way.
        idx = np.argsort(data, **kwargs)

        return idx[::-1] if reverse else idx

    def sort(self, keys=None, reverse=False):
        '''
        Sort the table according to one or more keys. This operates
        on the existing table and does not return a new table.

        Parameters
        ----------
        keys : str or list of str
            The key(s) to order the table by. If None, use the
            primary index of the Table.

        reverse : bool
            Sort in reverse order (default=False)

        Examples
        --------
        Create a table with 3 columns::

            >>> t = Table([['Max', 'Jo', 'John'], ['Miller', 'Miller', 'Jackson'],
            ...            [12, 15, 18]], names=('firstname', 'name', 'tel'))
            >>> print(t)
            firstname   name  tel
            --------- ------- ---
                  Max  Miller  12
                   Jo  Miller  15
                 John Jackson  18

        Sorting according to standard sorting rules, first 'name' then 'firstname'::

            >>> t.sort(['name', 'firstname'])
            >>> print(t)
            firstname   name  tel
            --------- ------- ---
                 John Jackson  18
                   Jo  Miller  15
                  Max  Miller  12

        Sorting according to standard sorting rules, first 'firstname' then 'tel',
        in reverse order::

            >>> t.sort(['firstname', 'tel'], reverse=True)
            >>> print(t)
            firstname   name  tel
            --------- ------- ---
                  Max  Miller  12
                 John Jackson  18
                   Jo  Miller  15
        '''
        if keys is None:
            if not self.indices:
                raise ValueError("Table sort requires input keys or a table index")
            keys = [x.info.name for x in self.indices[0].columns]

        if isinstance(keys, str):
            keys = [keys]

        indexes = self.argsort(keys)

        if reverse:
            indexes = indexes[::-1]

        with self.index_mode('freeze'):
            for name, col in self.columns.items():
                # Make a new sorted column.  This requires that take() also copies
                # relevant info attributes for mixin columns.
                new_col = col.take(indexes, axis=0)

                # First statement in try: will succeed if the column supports an in-place
                # update, and matches the legacy behavior of astropy Table.  However,
                # some mixin classes may not support this, so in that case just drop
                # in the entire new column. See #9553 and #9536 for discussion.
                try:
                    col[:] = new_col
                except Exception:
                    # In-place update failed for some reason, exception class not
                    # predictable for arbitrary mixin.
                    self[col.info.name] = new_col

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
            # First statement in try: will succeed if the column supports an in-place
            # update, and matches the legacy behavior of astropy Table.  However,
            # some mixin classes may not support this, so in that case just drop
            # in the entire new column. See #9836, #9553, and #9536 for discussion.
            new_col = col[::-1]
            try:
                col[:] = new_col
            except Exception:
                # In-place update failed for some reason, exception class not
                # predictable for arbitrary mixin.
                self[col.info.name] = new_col

        for index in self.indices:
            index.reverse()

    def round(self, decimals=0):
        '''
        Round numeric columns in-place to the specified number of decimals.
        Non-numeric columns will be ignored.

        Examples
        --------
        Create three columns with different types:

            >>> t = Table([[1, 4, 5], [-25.55, 12.123, 85],
            ...     ['a', 'b', 'c']], names=('a', 'b', 'c'))
            >>> print(t)
             a    b     c
            --- ------ ---
              1 -25.55   a
              4 12.123   b
              5   85.0   c

        Round them all to 0:

            >>> t.round(0)
            >>> print(t)
             a    b    c
            --- ----- ---
              1 -26.0   a
              4  12.0   b
              5  85.0   c

        Round column 'a' to -1 decimal:

            >>> t.round({'a':-1})
            >>> print(t)
             a    b    c
            --- ----- ---
              0 -26.0   a
              0  12.0   b
              0  85.0   c

        Parameters
        ----------
        decimals: int, dict
            Number of decimals to round the columns to. If a dict is given,
            the columns will be rounded to the number specified as the value.
            If a certain column is not in the dict given, it will remain the
            same.
        '''
        if isinstance(decimals, dict):
            decimal_values = decimals.values()
            column_names = decimals.keys()
        elif isinstance(decimals, int):
            decimal_values = itertools.repeat(decimals)
            column_names = self.colnames
        else:
            raise ValueError("'decimals' argument must be an int or a dict")

        for colname, decimal in zip(column_names, decimal_values):
            col = self.columns[colname]
            if np.issubdtype(col.info.dtype, np.number):
                try:
                    np.around(col, decimals=decimal, out=col)
                except TypeError:
                    # Bug in numpy see https://github.com/numpy/numpy/issues/15438
                    col[()] = np.around(col, decimals=decimal)

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
        return self._rows_equal(other)

    def __ne__(self, other):
        return ~self.__eq__(other)

    def _rows_equal(self, other):
        """
        Row-wise comparison of table with any other object.

        This is actual implementation for __eq__.

        Returns a 1-D boolean numpy array showing result of row-wise comparison.
        This is the same as the ``==`` comparison for tables.

        Parameters
        ----------
        other : Table or DataFrame or ndarray
             An object to compare with table

        Examples
        --------
        Comparing one Table with other::

            >>> t1 = Table([[1,2],[4,5],[7,8]], names=('a','b','c'))
            >>> t2 = Table([[1,2],[4,5],[7,8]], names=('a','b','c'))
            >>> t1._rows_equal(t2)
            array([ True,  True])

        """

        if isinstance(other, Table):
            other = other.as_array()

        if self.has_masked_columns:
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

    def values_equal(self, other):
        """
        Element-wise comparison of table with another table, list, or scalar.

        Returns a ``Table`` with the same columns containing boolean values
        showing result of comparison.

        Parameters
        ----------
        other : Table-like object or list or scalar
             Object to compare with table

        Examples
        --------
        Compare one Table with other::

          >>> t1 = Table([[1, 2], [4, 5], [-7, 8]], names=('a', 'b', 'c'))
          >>> t2 = Table([[1, 2], [-4, 5], [7, 8]], names=('a', 'b', 'c'))
          >>> t1.values_equal(t2)
          <Table length=2>
           a     b     c
          bool  bool  bool
          ---- ----- -----
          True False False
          True  True  True

        """
        if isinstance(other, Table):
            names = other.colnames
        else:
            try:
                other = Table(other, copy=False)
                names = other.colnames
            except Exception:
                # Broadcast other into a dict, so e.g. other = 2 will turn into
                # other = {'a': 2, 'b': 2} and then equality does a
                # column-by-column broadcasting.
                names = self.colnames
                other = {name: other for name in names}

        # Require column names match but do not require same column order
        if set(self.colnames) != set(names):
            raise ValueError('cannot compare tables with different column names')

        eqs = []
        for name in names:
            try:
                np.broadcast(self[name], other[name])  # Check if broadcast-able
                # Catch the numpy FutureWarning related to equality checking,
                # "elementwise comparison failed; returning scalar instead, but
                #  in the future will perform elementwise comparison".  Turn this
                # into an exception since the scalar answer is not what we want.
                with warnings.catch_warnings(record=True) as warns:
                    warnings.simplefilter('always')
                    eq = self[name] == other[name]
                    if (warns and issubclass(warns[-1].category, FutureWarning)
                            and 'elementwise comparison failed' in str(warns[-1].message)):
                        raise FutureWarning(warns[-1].message)
            except Exception as err:
                raise ValueError(f'unable to compare column {name}') from err

            # Be strict about the result from the comparison. E.g. SkyCoord __eq__ is just
            # broken and completely ignores that it should return an array.
            if not (isinstance(eq, np.ndarray)
                    and eq.dtype is np.dtype('bool')
                    and len(eq) == len(self)):
                raise TypeError(f'comparison for column {name} returned {eq} '
                                f'instead of the expected boolean ndarray')

            eqs.append(eq)

        out = Table(eqs, names=names)

        return out

    @property
    def groups(self):
        if not hasattr(self, '_groups'):
            self._groups = groups.TableGroups(self)
        return self._groups

    def group_by(self, keys):
        """
        Group this table by the specified ``keys``

        This effectively splits the table into groups which correspond to unique
        values of the ``keys`` grouping object.  The output is a new
        `~astropy.table.TableGroups` which contains a copy of this table but
        sorted by row according to ``keys``.

        The ``keys`` input to `group_by` can be specified in different ways:

          - String or list of strings corresponding to table column name(s)
          - Numpy array (homogeneous or structured) with same length as this table
          - `~astropy.table.Table` with same length as this table

        Parameters
        ----------
        keys : str, list of str, numpy array, or `~astropy.table.Table`
            Key grouping object

        Returns
        -------
        out : `~astropy.table.Table`
            New table with groups set
        """
        return groups.table_group_by(self, keys)

    def to_pandas(self, index=None, use_nullable_int=True):
        """
        Return a :class:`pandas.DataFrame` instance

        The index of the created DataFrame is controlled by the ``index``
        argument.  For ``index=True`` or the default ``None``, an index will be
        specified for the DataFrame if there is a primary key index on the
        Table *and* if it corresponds to a single column.  If ``index=False``
        then no DataFrame index will be specified.  If ``index`` is the name of
        a column in the table then that will be the DataFrame index.

        In addition to vanilla columns or masked columns, this supports Table
        mixin columns like Quantity, Time, or SkyCoord.  In many cases these
        objects have no analog in pandas and will be converted to a "encoded"
        representation using only Column or MaskedColumn.  The exception is
        Time or TimeDelta columns, which will be converted to the corresponding
        representation in pandas using ``np.datetime64`` or ``np.timedelta64``.
        See the example below.

        Parameters
        ----------
        index : None, bool, str
            Specify DataFrame index mode
        use_nullable_int : bool, default=True
            Convert integer MaskedColumn to pandas nullable integer type.
            If ``use_nullable_int=False`` or the pandas version does not support
            nullable integer types (version < 0.24), then the column is converted
            to float with NaN for missing elements and a warning is issued.

        Returns
        -------
        dataframe : :class:`pandas.DataFrame`
            A pandas :class:`pandas.DataFrame` instance

        Raises
        ------
        ImportError
            If pandas is not installed
        ValueError
            If the Table has multi-dimensional columns

        Examples
        --------
        Here we convert a table with a few mixins to a
        :class:`pandas.DataFrame` instance.

          >>> import pandas as pd
          >>> from astropy.table import QTable
          >>> import astropy.units as u
          >>> from astropy.time import Time, TimeDelta
          >>> from astropy.coordinates import SkyCoord

          >>> q = [1, 2] * u.m
          >>> tm = Time([1998, 2002], format='jyear')
          >>> sc = SkyCoord([5, 6], [7, 8], unit='deg')
          >>> dt = TimeDelta([3, 200] * u.s)

          >>> t = QTable([q, tm, sc, dt], names=['q', 'tm', 'sc', 'dt'])

          >>> df = t.to_pandas(index='tm')
          >>> with pd.option_context('display.max_columns', 20):
          ...     print(df)
                        q  sc.ra  sc.dec              dt
          tm
          1998-01-01  1.0    5.0     7.0 0 days 00:00:03
          2002-01-01  2.0    6.0     8.0 0 days 00:03:20

        """
        from pandas import DataFrame, Series

        if index is not False:
            if index in (None, True):
                # Default is to use the table primary key if available and a single column
                if self.primary_key and len(self.primary_key) == 1:
                    index = self.primary_key[0]
                else:
                    index = False
            else:
                if index not in self.colnames:
                    raise ValueError('index must be None, False, True or a table '
                                     'column name')

        def _encode_mixins(tbl):
            """Encode a Table ``tbl`` that may have mixin columns to a Table with only
            astropy Columns + appropriate meta-data to allow subsequent decoding.
            """
            from . import serialize
            from astropy.time import TimeBase, TimeDelta

            # Convert any Time or TimeDelta columns and pay attention to masking
            time_cols = [col for col in tbl.itercols() if isinstance(col, TimeBase)]
            if time_cols:

                # Make a light copy of table and clear any indices
                new_cols = []
                for col in tbl.itercols():
                    new_col = col_copy(col, copy_indices=False) if col.info.indices else col
                    new_cols.append(new_col)
                tbl = tbl.__class__(new_cols, copy=False)

                for col in time_cols:
                    if isinstance(col, TimeDelta):
                        # Convert to nanoseconds (matches astropy datetime64 support)
                        new_col = (col.sec * 1e9).astype('timedelta64[ns]')
                        nat = np.timedelta64('NaT')
                    else:
                        new_col = col.datetime64.copy()
                        nat = np.datetime64('NaT')
                    if col.masked:
                        new_col[col.mask] = nat
                    tbl[col.info.name] = new_col

            # Convert the table to one with no mixins, only Column objects.
            encode_tbl = serialize.represent_mixins_as_columns(tbl)
            return encode_tbl

        tbl = _encode_mixins(self)

        badcols = [name for name, col in self.columns.items() if len(col.shape) > 1]
        if badcols:
            raise ValueError(
                f'Cannot convert a table with multidimensional columns to a '
                f'pandas DataFrame. Offending columns are: {badcols}\n'
                f'One can filter out such columns using:\n'
                f'names = [name for name in tbl.colnames if len(tbl[name].shape) <= 1]\n'
                f'tbl[names].to_pandas(...)')

        out = OrderedDict()

        for name, column in tbl.columns.items():
            if getattr(column.dtype, 'isnative', True):
                out[name] = column
            else:
                out[name] = column.data.byteswap().newbyteorder('=')

            if isinstance(column, MaskedColumn) and np.any(column.mask):
                if column.dtype.kind in ['i', 'u']:
                    pd_dtype = column.dtype.name
                    if use_nullable_int:
                        # Convert int64 to Int64, uint32 to UInt32, etc for nullable types
                        pd_dtype = pd_dtype.replace('i', 'I').replace('u', 'U')
                    out[name] = Series(out[name], dtype=pd_dtype)

                    # If pandas is older than 0.24 the type may have turned to float
                    if column.dtype.kind != out[name].dtype.kind:
                        warnings.warn(
                            f"converted column '{name}' from {column.dtype} to {out[name].dtype}",
                            TableReplaceWarning, stacklevel=3)
                elif column.dtype.kind not in ['f', 'c']:
                    out[name] = column.astype(object).filled(np.nan)

        kwargs = {'index': out.pop(index)} if index else {}

        return DataFrame(out, **kwargs)

    @classmethod
    def from_pandas(cls, dataframe, index=False, units=None):
        """
        Create a `~astropy.table.Table` from a :class:`pandas.DataFrame` instance

        In addition to converting generic numeric or string columns, this supports
        conversion of pandas Date and Time delta columns to `~astropy.time.Time`
        and `~astropy.time.TimeDelta` columns, respectively.

        Parameters
        ----------
        dataframe : :class:`pandas.DataFrame`
            A pandas :class:`pandas.DataFrame` instance
        index : bool
            Include the index column in the returned table (default=False)
        units: dict
            A dict mapping column names to to a `~astropy.units.Unit`.
            The columns will have the specified unit in the Table.

        Returns
        -------
        table : `~astropy.table.Table`
            A `~astropy.table.Table` (or subclass) instance

        Raises
        ------
        ImportError
            If pandas is not installed

        Examples
        --------
        Here we convert a :class:`pandas.DataFrame` instance
        to a `~astropy.table.QTable`.

          >>> import numpy as np
          >>> import pandas as pd
          >>> from astropy.table import QTable

          >>> time = pd.Series(['1998-01-01', '2002-01-01'], dtype='datetime64[ns]')
          >>> dt = pd.Series(np.array([1, 300], dtype='timedelta64[s]'))
          >>> df = pd.DataFrame({'time': time})
          >>> df['dt'] = dt
          >>> df['x'] = [3., 4.]
          >>> with pd.option_context('display.max_columns', 20):
          ...     print(df)
                  time              dt    x
          0 1998-01-01 0 days 00:00:01  3.0
          1 2002-01-01 0 days 00:05:00  4.0

          >>> QTable.from_pandas(df)
          <QTable length=2>
                    time            dt      x
                   object         object float64
          ----------------------- ------ -------
          1998-01-01T00:00:00.000    1.0     3.0
          2002-01-01T00:00:00.000  300.0     4.0

        """

        out = OrderedDict()

        names = list(dataframe.columns)
        columns = [dataframe[name] for name in names]
        datas = [np.array(column) for column in columns]
        masks = [np.array(column.isnull()) for column in columns]

        if index:
            index_name = dataframe.index.name or 'index'
            while index_name in names:
                index_name = '_' + index_name + '_'
            names.insert(0, index_name)
            columns.insert(0, dataframe.index)
            datas.insert(0, np.array(dataframe.index))
            masks.insert(0, np.zeros(len(dataframe), dtype=bool))

        if units is None:
            units = [None] * len(names)
        else:
            if not isinstance(units, Mapping):
                raise TypeError('Expected a Mapping "column-name" -> "unit"')

            not_found = set(units.keys()) - set(names)
            if not_found:
                warnings.warn('`units` contains additionial columns: {}'.format(
                    not_found
                ))

            units = [units.get(name) for name in names]

        for name, column, data, mask, unit in zip(names, columns, datas, masks, units):

            if column.dtype.kind in ['u', 'i'] and np.any(mask):
                # Special-case support for pandas nullable int
                np_dtype = str(column.dtype).lower()
                data = np.zeros(shape=column.shape, dtype=np_dtype)
                data[~mask] = column[~mask]
                out[name] = MaskedColumn(data=data, name=name, mask=mask, unit=unit, copy=False)
                continue

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

            # Numpy datetime64
            if data.dtype.kind == 'M':
                from astropy.time import Time
                out[name] = Time(data, format='datetime64')
                if np.any(mask):
                    out[name][mask] = np.ma.masked
                out[name].format = 'isot'

            # Numpy timedelta64
            elif data.dtype.kind == 'm':
                from astropy.time import TimeDelta
                data_sec = data.astype('timedelta64[ns]').astype(np.float64) / 1e9
                out[name] = TimeDelta(data_sec, format='sec')
                if np.any(mask):
                    out[name][mask] = np.ma.masked

            else:
                if np.any(mask):
                    out[name] = MaskedColumn(data=data, name=name, mask=mask, unit=unit)
                else:
                    out[name] = Column(data=data, name=name, unit=unit)

        return cls(out)

    info = TableInfo()


class QTable(Table):
    """A class to represent tables of heterogeneous data.

    `~astropy.table.QTable` provides a class for heterogeneous tabular data
    which can be easily modified, for instance adding columns or new rows.

    The `~astropy.table.QTable` class is identical to `~astropy.table.Table`
    except that columns with an associated ``unit`` attribute are converted to
    `~astropy.units.Quantity` objects.

    See also:

    - https://docs.astropy.org/en/stable/table/
    - https://docs.astropy.org/en/stable/table/mixin_columns.html

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

    def _is_mixin_for_table(self, col):
        """
        Determine if ``col`` should be added to the table directly as
        a mixin column.
        """
        return has_info_class(col, MixinInfo)

    def _convert_col_for_table(self, col):
        if isinstance(col, Column) and getattr(col, 'unit', None) is not None:
            # What to do with MaskedColumn with units: leave as MaskedColumn or
            # turn into Quantity and drop mask?  Assuming we have masking support
            # in Quantity someday, let's drop the mask (consistent with legacy
            # behavior) but issue a warning.
            if isinstance(col, MaskedColumn) and np.any(col.mask):
                warnings.warn("dropping mask in Quantity column '{}': "
                              "masked Quantity not supported".format(col.info.name))

            # We need to turn the column into a quantity, or a subclass
            # identified in the unit (such as u.mag()).
            q_cls = getattr(col.unit, '_quantity_class', Quantity)
            qcol = q_cls(col.data, col.unit, copy=False)
            qcol.info = col.info
            qcol.info.indices = col.info.indices
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


class TableAttribute(MetaAttribute):
    """
    Descriptor to define a custom attribute for a Table subclass.

    The value of the ``TableAttribute`` will be stored in a dict named
    ``__attributes__`` that is stored in the table ``meta``.  The attribute
    can be accessed and set in the usual way, and it can be provided when
    creating the object.

    Defining an attribute by this mechanism ensures that it will persist if
    the table is sliced or serialized, for example as a pickle or ECSV file.

    See the `~astropy.utils.metadata.MetaAttribute` documentation for additional
    details.

    Parameters
    ----------
    default : object
        Default value for attribute

    Examples
    --------
      >>> from astropy.table import Table, TableAttribute
      >>> class MyTable(Table):
      ...     identifier = TableAttribute(default=1)
      >>> t = MyTable(identifier=10)
      >>> t.identifier
      10
      >>> t.meta
      OrderedDict([('__attributes__', {'identifier': 10})])
    """
