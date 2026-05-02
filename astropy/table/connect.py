# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import annotations

import itertools
import warnings
from typing import TYPE_CHECKING, Literal, NotRequired, TypedDict

from astropy.io import registry
from astropy.utils.exceptions import AstropyWarning

from .info import serialize_method_as

if TYPE_CHECKING:
    from astropy.table import Table
    from astropy.table.index import SlicedIndex

__all__ = ["TableRead", "TableWrite"]
__doctest_skip__ = ["TableRead", "TableWrite"]


class TableRead(registry.UnifiedReadWrite):
    """Read and parse a data table and return as a Table.

    This function provides the Table interface to the astropy unified I/O
    layer.  This allows easily reading a file in many supported data formats
    using syntax such as::

      >>> from astropy.table import Table
      >>> dat = Table.read('table.dat', format='ascii')
      >>> events = Table.read('events.fits', format='fits')

    Get help on the available readers for ``Table`` using the``help()`` method::

      >>> Table.read.help()  # Get help reading Table and list supported formats
      >>> Table.read.help('fits')  # Get detailed help on Table FITS reader
      >>> Table.read.list_formats()  # Print list of available formats

    See also: https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data reader. If supplied the
        first argument is typically the input filename.
    format : str
        File format specifier.
    units : list, dict, optional
        List or dict of units to apply to columns
    descriptions : list, dict, optional
        List or dict of descriptions to apply to columns
    **kwargs : dict, optional
        Keyword arguments passed through to data reader.

    Returns
    -------
    out : `~astropy.table.Table`
        Table corresponding to file contents

    Notes
    -----
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "read", registry=None)
        # uses default global registry

    def __call__(self, *args, **kwargs):
        cls = self._cls
        units = kwargs.pop("units", None)
        descriptions = kwargs.pop("descriptions", None)

        out = self.registry.read(cls, *args, **kwargs)

        # For some readers (e.g., ascii.ecsv), the returned `out` class is not
        # guaranteed to be the same as the desired output `cls`.  If so,
        # try coercing to desired class without copying (io.registry.read
        # would normally do a copy).  The normal case here is swapping
        # Table <=> QTable.
        if cls is not out.__class__:
            try:
                out = cls(out, copy=False)
            except Exception:
                raise TypeError(
                    f"could not convert reader output to {cls.__name__} class."
                )

        out._set_column_attribute("unit", units)
        out._set_column_attribute("description", descriptions)

        construct_indices(out)

        return out


class TableWrite(registry.UnifiedReadWrite):
    """
    Write this Table object out in the specified format.

    This function provides the Table interface to the astropy unified I/O
    layer.  This allows easily writing a file in many supported data formats
    using syntax such as::

      >>> from astropy.table import Table
      >>> dat = Table([[1, 2], [3, 4]], names=('a', 'b'))
      >>> dat.write('table.dat', format='ascii')

    Get help on the available writers for ``Table`` using the``help()`` method::

      >>> Table.write.help()  # Get help writing Table and list supported formats
      >>> Table.write.help('fits')  # Get detailed help on Table FITS writer
      >>> Table.write.list_formats()  # Print list of available formats

    The ``serialize_method`` argument is explained in the section on
    `Table serialization methods
    <https://docs.astropy.org/en/latest/io/unified.html#table-serialization-methods>`_.

    See also: https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data writer. If supplied the
        first argument is the output filename.
    format : str
        File format specifier.
    serialize_method : str, dict, optional
        Serialization method specifier for columns.
    write_indices : bool, optional
        Write table indices if defined (default is False).
    **kwargs : dict, optional
        Keyword arguments passed through to data writer.

    Notes
    -----
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "write", registry=None)
        # uses default global registry

    def __call__(self, *args, serialize_method=None, write_indices=False, **kwargs):
        tbl = self._instance
        with serialize_method_as(tbl, serialize_method):
            if write_indices and tbl.indices:
                tbl = represent_indices(tbl)
            self.registry.write(tbl, *args, **kwargs)


def construct_indices(tbl: Table) -> None:
    """
    Create indices on a table from column meta ``__indices__`` values.

    This is the inverse of `represent_indices`.

    Parameters
    ----------
    tbl : Table
        Table in which to create indices (in-place).
    """
    # No action if there are no indices defined in table meta. Otherwise get the
    # table indices meta and pop that key off of meta.
    if (indices_meta := tbl.meta.pop("__table_indices__", None)) is None:
        return

    indices: list[SlicedIndex] = []
    row_index_colnames = set()

    primary_key = tuple(indices_meta["primary_key"])
    for index_info in indices_meta["indices"]:
        indices.append(construct_sliced_index(tbl, index_info))
        row_index_colnames.add(index_info["index_colname"])

    for index in indices:
        # Add index to table by adding to each column indices
        for colname in index.id:
            tbl[colname].info.indices.append(index)

    tbl.primary_key = primary_key
    tbl.remove_columns(row_index_colnames)


class IndexInfo(TypedDict):
    colnames: tuple[str, ...]
    index_colname: int
    engine: NotRequired[Literal["SortedArray", "SCEngine"]]
    unique: NotRequired[bool]
    primary: NotRequired[bool]


def construct_sliced_index(tbl: Table, /, index_info: IndexInfo) -> SlicedIndex:
    """
    Construct an index (SlicedIndex) for a table from serialized index information.

    This function reconstructs a SlicedIndex object for a table using the
    provided index metadata, including the engine type, column names, and
    row index column. It is used when deserializing tables with indices
    from formats that store index metadata in column meta information.

    Parameters
    ----------
    tbl : Table
        Table containing the columns and row index data.
    index_info : IndexInfo
        Dictionary containing index metadata, including:
        - 'index_colname': Name of the column with row indices.
        - 'colnames' (tuple[str, ...]): Tuple of column names for the index.
        - 'engine' (str, optional): Name of the index engine class ('SortedArray' or 'SCEngine')
        - 'unique' (bool, optional): Whether the index is unique.
        - 'primary' (bool, optional): Whether index is the primary key index (default=False).

    Returns
    -------
    SlicedIndex
        The reconstructed SlicedIndex object for the table.

    Warns
    -----
    AstropyWarning
        If the engine type is not recognized, a warning is issued and
        no index is created.
    """
    from astropy.table import Table
    from astropy.table.index import ENGINE_CLS_DEFAULT, Index, SlicedIndex
    from astropy.table.soco import SCEngine
    from astropy.table.sorted_array import SortedArray

    match index_info.get("engine", ENGINE_CLS_DEFAULT.__name__):
        case "SortedArray":
            engine_cls = SortedArray
        case "SCEngine":
            engine_cls = SCEngine
        case _:
            warnings.warn(
                f"Unknown index engine {index_info['engine']!r}, "
                "creating index using SortedArray engine",
                AstropyWarning,
                stacklevel=2,
            )
            engine_cls = SortedArray
    row_index_colname = index_info["index_colname"]
    row_index = tbl[row_index_colname]
    colnames = index_info["colnames"]

    # By-reference Table of tbl index columns for input to engine. The current
    # Index.__init__() requires index_columns here to be a Table even if tbl is a
    # QTable. Changing Table to QTable in index.py breaks tests, in particular in
    # TableLoc._get_row_idxs_as_list(): if a slice range has no explicit min or max then
    # a special value from BST is used (MinValue(), MaxValue()), but this does not work
    # with Quantity.searchsorted().
    index_columns = Table([tbl[colname] for colname in colnames], copy=False)

    engine = engine_cls(
        data=index_columns[row_index],  # index columns sorted by index
        row_index=row_index,
        unique=index_info.get("unique", False),
    )
    index = Index(columns=list(index_columns.itercols()), engine=engine)
    sliced_index = SlicedIndex(index, index_slice=slice(None), original=True)

    return sliced_index


def represent_indices(tbl: Table, /) -> Table:
    """
    Return table with indices serialized as new columns and meta data.

    This creates a light copy of the table with columns added for each index that
    contains the row indices in sorted order.  The new columns are given names
    ``__index__0``, ``__index__1``, etc (ensuring these do not conflict with existing
    columns).

    For each column with indices, ``meta['__indices__']`` is create to include
    information about the index. For multi-key indices this results in repeated data in
    the header.

    The key advantage of putting data into column meta is that the FITS writer
    automatically serializes column meta into YAML comments, while table meta is assumed
    to be FITS keyword files.

    See the ``test_indices_serialization_*` tests in ``test_index.py`` for examples.

    Parameters
    ----------
    tbl : Table
        Table with indices to serialize

    Returns
    -------
    Table
        Shallow copy of input table with index columns added and meta updated.
    """
    from astropy.table.index import ENGINE_CLS_DEFAULT

    with tbl.index_mode("discard_on_copy"):
        tbl_out = tbl.copy(copy_data=False)

    indices_info = []
    for index in tbl.indices:
        row_index = index.data.sorted_data()

        # Find unique column name for the index row data.  This must complete in no more
        # than len(tbl_out.colnames) iterations. It is generating a distinct column name
        # in each loop, and tbl_out.colnames has len(tbl_out.colnames) values.
        for ii_index in itertools.count():
            colname = "__index__" + ("" if ii_index == 0 else str(ii_index))
            if colname not in tbl_out.colnames:
                break

        # Make new column for row_index and add meta describing index
        tbl_out[colname] = row_index
        index_info = {
            "index_colname": colname,
            "colnames": list(index.id),
        }
        if type(index.data) is not ENGINE_CLS_DEFAULT:
            index_info["engine"] = type(index.data).__name__
        if index.data.unique:
            index_info["unique"] = True

        indices_info.append(index_info)

        tbl_out.meta["__table_indices__"] = {
            "primary_key": list(tbl.primary_key),
            "indices": indices_info,
        }

    return tbl_out
