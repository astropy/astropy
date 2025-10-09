# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from astropy.io import registry
from astropy.utils.exceptions import AstropyWarning

from .info import serialize_method_as

if TYPE_CHECKING:
    from astropy.table import Table

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

        if "__table_indices__" in out.meta:
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

    This is the inverse of `copy_and_serialize_indices`.

    Parameters
    ----------
    tbl : Table
        Table in which to create indices (in-place).
    """
    from astropy.table.index import Index, SlicedIndex
    from astropy.table.soco import SCEngine
    from astropy.table.sorted_array import SortedArray

    indices = []
    indices_meta = tbl.meta["__table_indices__"]

    for index_info in indices_meta["indices"]:
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
        row_index_colname = index_info["row_index_colname"]
        row_index = tbl[row_index_colname]
        colnames = index_info["colnames"]

        # Colnames slice of tbl (cols are by reference, no copy)
        tbl_colnames = tbl[colnames]

        engine = engine_cls(
            data=tbl_colnames[row_index],  # index columns sorted by index
            row_index=row_index,
            unique=index_info["unique"],
        )
        index = Index(columns=list(tbl_colnames.itercols()), engine=engine)
        sliced_index = SlicedIndex(index, index_slice=slice(None), original=True)
        indices.append(sliced_index)

        del tbl[row_index_colname]

    for index in indices:
        # Add index to table by adding to each column indices.
        for colname in index.id:
            tbl[colname].info.indices.append(index)

    tbl.primary_key = indices_meta["primary_key"]
    del tbl.meta["__table_indices__"]


def represent_indices(tbl: Table) -> Table:
    """
    Make a copy of a table with indices serialized as new columns and meta data.

    This creates a light copy of the table with columns added for each index that
    contains the row indices in sorted order.  The new columns are given names
    ``__index__0``, ``__index__1``, etc (ensuring these do not conflict with existing
    columns).  The table ``meta['__table_indices__']`` is also updated to include
    information about the new columns and the index engine.

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
            colname = f"__index__{ii_index}"
            if colname not in tbl_out.colnames:
                break
            ii_index += 1

        # Make new column and add meta
        tbl_out[colname] = row_index
        indices_info.append(
            {
                "row_index_colname": colname,
                "colnames": index.id,
                "engine": index.data.__class__.__name__,
                "unique": index.data.unique,
            }
        )

    tbl_out.meta["__table_indices__"] = {
        "primary_key": tbl.primary_key,
        "indices": indices_info,
    }
    return tbl_out
