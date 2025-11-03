# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

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
    Create indices on a table based on the ``__table_indices__`` meta data.

    This is the inverse of `represent_indices`.

    Parameters
    ----------
    tbl : Table
        Table in which to create indices (in-place).
    """
    indices: list[SlicedIndex] = []

    primary_key = None
    row_index_colnames = set()
    for col in tbl.itercols():
        if not col.info.meta or "__indices__" not in col.info.meta:
            continue

        for index_info in col.info.meta["__indices__"]:
            if index_info.get("primary"):
                primary_key = tuple(index_info["colnames"])
            indices.append(construct_sliced_index(tbl, index_info))
            row_index_colnames.add(index_info["index_colname"])

    # No indices, do nothing
    if not indices:
        return

    if primary_key is None:
        primary_key = indices[0].id

    for index in indices:
        # Add index to table by adding to each column indices and clean up column meta.
        for colname in index.id:
            tbl[colname].info.indices.append(index)
            tbl[colname].info.meta.pop("__indices__", None)

    tbl.primary_key = primary_key
    tbl.remove_columns(row_index_colnames)


def construct_sliced_index(tbl: Table, index_info: dict[str, Any]) -> SlicedIndex:
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
    index_info : dict[str, Any]
        Dictionary containing index metadata, including:
        - 'engine': Name of the index engine class.
        - 'row_index_colname': Name of the column with row indices.
        - 'colnames': Tuple of column names for the index.
        - 'unique': Whether the index is unique.
        - 'primary': Whether index is the primary key index (default=False).

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
    from astropy.table.index import Index, SlicedIndex
    from astropy.table.soco import SCEngine
    from astropy.table.sorted_array import SortedArray

    match index_info["engine"]:
        case "SortedArray":
            engine_cls = SortedArray
        case "SCEngine":
            engine_cls = SCEngine
        case _:
            warnings.warn(
                f'Cannot restore index with engine "{index_info["engine"]}".  '
                "Index not created.",
                AstropyWarning,
            )
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
        unique=index_info["unique"],
    )
    index = Index(columns=list(index_columns.itercols()), engine=engine)
    sliced_index = SlicedIndex(index, index_slice=slice(None), original=True)

    return sliced_index


def represent_indices(tbl: Table) -> Table:
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

    The format is most easily illustrated via the ECSV output in this example::

      >>> from astropy.table import QTable
      >>> import sys
      >>> t = QTable()
      >>> t["a"] = [2, 3, 1]
      >>> t["b"] = [3, 5, 4]
      >>> t.add_index("a")
      >>> t.write(sys.stdout, format="ecsv", write_indices=True)
      # %ECSV 1.0
      # ---
      # datatype:
      # - name: a
      #   datatype: int64
      #   meta: !!omap
      #   - __indices__:
      #     - colnames: !!python/tuple [a]
      #       engine: SortedArray
      #       primary: true
      #       row_index_colname: __index__0
      #       unique: false
      # - {name: b, datatype: int64}
      # - {name: __index__0, datatype: int64}
      # schema: astropy-2.0
      a b __index__0
      2 3 2
      3 5 0
      1 4 1

    Parameters
    ----------
    tbl : Table
        Table with indices to serialize

    Returns
    -------
    Table
        Copy of input table with index columns added and meta updated.
    """
    with tbl.index_mode("discard_on_copy"):
        tbl_out = tbl.copy(copy_data=False)

    indices_info = []
    ii_index = 0
    for index in tbl.indices:
        row_index = index.data.sorted_data()

        # Find unique column name for the index row data
        while True:
            colname = "__index__" + ("" if ii_index == 0 else str(ii_index))
            if colname not in tbl_out.colnames:
                break
            ii_index += 1

        # Make new column for row_index and add meta describing index
        tbl_out[colname] = row_index
        index_info = {
            "index_colname": colname,
            "colnames": list(index.id),
            "engine": index.data.__class__.__name__,
            "unique": index.data.unique,
        }
        if len(tbl.indices) > 1 and index.id == tbl.primary_key:
            index_info["primary"] = True
        indices_info.append(index_info)

    for index_info in indices_info:
        # Store the index information on the first column in the index. This is somewhat
        # arbitrary but eliminates duplication of index_info for multi-column indices.
        col = tbl_out[index_info["colnames"][0]]
        if col.info.meta is None:
            col.info.meta = {}
        col_meta_indices = col.info.meta.setdefault("__indices__", [])
        col_meta_indices.append(index_info)

    return tbl_out
