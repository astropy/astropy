# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os

import numpy as np

from astropy.io import registry as io_registry
from astropy.table import Table
from astropy.table.column import BaseColumn
from astropy.units import Quantity
from astropy.utils.misc import NOT_OVERWRITING_MSG

from . import from_table, parse
from .tree import TableElement, VOTableFile


def is_votable(origin, filepath, fileobj, *args, **kwargs):
    """
    Reads the header of a file to determine if it is a VOTable file.

    Parameters
    ----------
    origin : str or readable file-like
        Path or file object containing a VOTABLE_ xml file.

    Returns
    -------
    is_votable : bool
        Returns `True` if the given file is a VOTable file.
    """
    from . import is_votable

    if origin == "read":
        if fileobj is not None:
            try:
                result = is_votable(fileobj)
            finally:
                fileobj.seek(0)
            return result
        elif filepath is not None:
            return is_votable(filepath)
        return isinstance(args[0], (VOTableFile, TableElement))

    else:
        return False


def read_table_votable(
    input, table_id=None, use_names_over_ids=False, verify=None, **kwargs
):
    """
    Read a Table object from an VO table file.

    Parameters
    ----------
    input : str or `~astropy.io.votable.tree.VOTableFile` or `~astropy.io.votable.tree.TableElement`
        If a string, the filename to read the table from. If a
        :class:`~astropy.io.votable.tree.VOTableFile` or
        :class:`~astropy.io.votable.tree.TableElement` object, the object to extract
        the table from.

    table_id : str or int, optional
        The table to read in.  If a `str`, it is an ID corresponding
        to the ID of the table in the file (not all VOTable files
        assign IDs to their tables).  If an `int`, it is the index of
        the table in the file, starting at 0.

    use_names_over_ids : bool, optional
        When `True` use the ``name`` attributes of columns as the names
        of columns in the `~astropy.table.Table` instance.  Since names
        are not guaranteed to be unique, this may cause some columns
        to be renamed by appending numbers to the end.  Otherwise
        (default), use the ID attributes as the column names.

    verify : {'ignore', 'warn', 'exception'}, optional
        When ``'exception'``, raise an error when the file violates the spec,
        otherwise either issue a warning (``'warn'``) or silently continue
        (``'ignore'``). Warnings may be controlled using the standard Python
        mechanisms.  See the `warnings` module in the Python standard library
        for more information. When not provided, uses the configuration setting
        ``astropy.io.votable.verify``, which defaults to ``'ignore'``.

    **kwargs
        Additional keyword arguments are passed on to `astropy.io.votable.parse`.
    """
    if not isinstance(input, (VOTableFile, TableElement)):
        input = parse(input, table_id=table_id, verify=verify, **kwargs)

    # Parse all table objects
    table_id_mapping = {}
    tables = []
    if isinstance(input, VOTableFile):
        for table in input.iter_tables():
            if table.ID is not None:
                table_id_mapping[table.ID] = table
            tables.append(table)

        if len(tables) > 1:
            if table_id is None:
                raise ValueError(
                    "Multiple tables found: table id should be set via"
                    " the table_id= argument. The available tables are"
                    f" {', '.join(table_id_mapping)}, or integers less than"
                    f" {len(tables)}."
                )
            elif isinstance(table_id, str):
                if table_id in table_id_mapping:
                    table = table_id_mapping[table_id]
                else:
                    raise ValueError(f"No tables with id={table_id} found")
            elif isinstance(table_id, int):
                if table_id < len(tables):
                    table = tables[table_id]
                else:
                    raise IndexError(
                        f"Table index {table_id} is out of range. {len(tables)} tables"
                        " found"
                    )
        elif len(tables) == 1:
            table = tables[0]
        else:
            raise ValueError("No table found")
    elif isinstance(input, TableElement):
        table = input

    # Convert to an astropy.table.Table object
    return table.to_table(use_names_over_ids=use_names_over_ids)


def write_table_votable(
    input, output, table_id=None, overwrite=False, tabledata_format=None
):
    """
    Write a Table object to an VO table file.

    Parameters
    ----------
    input : Table
        The table to write out.

    output : str
        The filename to write the table to.

    table_id : str, optional
        The table ID to use. If this is not specified, the 'ID' keyword in the
        ``meta`` object of the table will be used.

    overwrite : bool, optional
        Whether to overwrite any existing file without warning.

    tabledata_format : str, optional
        The format of table data to write.  Must be one of ``tabledata``
        (text representation), ``binary`` or ``binary2``.  Default is
        ``tabledata``.  See :ref:`astropy:votable-serialization`.
    """
    # Only those columns which are instances of BaseColumn or Quantity can be written
    unsupported_cols = input.columns.not_isinstance((BaseColumn, Quantity))
    if unsupported_cols:
        unsupported_names = [col.info.name for col in unsupported_cols]
        raise ValueError(
            f"cannot write table with mixin column(s) {unsupported_names} to VOTable"
        )

    # Check if output file already exists
    if isinstance(output, str) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise OSError(NOT_OVERWRITING_MSG.format(output))

    # Create a new VOTable file
    table_file = from_table(input, table_id=table_id)

    # Write out file
    table_file.to_xml(output, tabledata_format=tabledata_format)


io_registry.register_reader("votable", Table, read_table_votable)
io_registry.register_writer("votable", Table, write_table_votable)
io_registry.register_identifier("votable", Table, is_votable)


# VOTable with embedded/linked Parquet file #
def write_table_votable_parquet(input, output, column_metadata, *, overwrite=False):
    """
    This function allows writing a VOTable (XML) with PARQUET
    serialization. This functionality is currently not
    supported by Astropy (with the reason that this method
    requires writing multiple files: a VOTable/XML and
    PARQUET table). This function presents a wrapper, which
    allows to do this. The concept is simple and probably
    can be improved substantially. We first save the PARQUET
    table using Astropy functionality. Then, we create a
    VOTable with binary serialization. The latter is modified
    later to include an external reference to the create
    PARQUET table file.

    Parameters
    ----------
    input : `~astropy.table.Table`
        The table to write out.

    output : str
        The filename to write the table to.

    column_metadata : dict
        Contains the metadata for the columns such as "unit" or
        "ucd" or "utype".
        (Example: {"id": {"unit": "", "ucd": "meta.id", "utype": "none"},
                   "mass": {"unit": "solMass", "ucd": "phys.mass", "utype": "none"}})
    overwrite : bool, optional
        Whether to overwrite any existing file without warning.

    Returns
    -------
    This function creates a VOTable serialized in Parquet.
    Two files are written:
    1. The VOTable (XML file) including the column metadata and a
        ``STREAM`` tag that embeds the PARQUET table.
    2. The PARQUET table itself.

    Both files are stored at the same location. The name of the
    VOTable is ``output``, and the name of the embedded PARQUET
    file is f"{output}.parquet".
    """
    # First save the PARQUET file.
    parquet_filename = f"{output}.parquet"
    path_type = f"file:{'//' if os.path.isabs(parquet_filename) else ''}"

    if os.path.exists(parquet_filename) and not overwrite:
        raise OSError(NOT_OVERWRITING_MSG.format(parquet_filename))
    input.write(parquet_filename, format="parquet", overwrite=overwrite)

    # Second, save table as binary VOT file. We will modify this file
    # later to incorporate the FITS stream. Note that we use here the full
    # table data so we get the datatype and arraysize correct. Later
    # we can maybe make this more efficient and instead write the
    # VOTable file from scratch, especially the FIELDS, which are the
    # most important.
    votablefile = VOTableFile()
    votable = votablefile.from_table(input)

    # Add the fields
    # Maybe there is a smarter way to do this iteratively.
    for field in votable.resources[0].tables[0].fields:
        field.unit = column_metadata[field.name]["unit"]
        field.ucd = column_metadata[field.name]["ucd"]
        field.utype = column_metadata[field.name]["utype"]

    if os.path.exists(output) and not overwrite:
        raise OSError(NOT_OVERWRITING_MSG.format(output))

    votable.to_xml(output, tabledata_format="binary")

    # Now reopen the binary file and replace the binary part with
    # the stream relating to the FITS file. This all is a bit flimsy
    # and needs to be made more bullet-proof.
    with open(output) as f:
        lines = f.readlines()

        # get start and end of <BINARY> tag
        line_start = np.where(["<BINARY>" in line for line in lines])[0][0]
        line_stop = np.where(["</BINARY>" in line for line in lines])[0][0]

        # Add the extension tag
        # We assume here that it is extension #1.
        lines[line_start] = '<PARQUET type="VOTable-remote-file">\n'
        lines[line_start + 1] = f'<STREAM href="{path_type}{parquet_filename}"/>\n'
        lines[line_start + 2] = "</PARQUET>\n"

        # remove last line
        _ = lines.pop(line_stop)

    # write new file
    with open(output, "w") as f:
        f.write("".join(lines))


io_registry.register_writer("votable.parquet", Table, write_table_votable_parquet)
