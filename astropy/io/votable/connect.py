# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os


from . import parse, from_table
from .tree import VOTableFile, Table as VOTable
from astropy.io import registry as io_registry
from astropy.table import Table
from astropy.table.column import BaseColumn
from astropy.units import Quantity
from astropy.utils.misc import NOT_OVERWRITING_MSG


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
    if origin == 'read':
        if fileobj is not None:
            try:
                result = is_votable(fileobj)
            finally:
                fileobj.seek(0)
            return result
        elif filepath is not None:
            return is_votable(filepath)
        elif isinstance(args[0], (VOTableFile, VOTable)):
            return True
        else:
            return False
    else:
        return False


def read_table_votable(input, table_id=None, use_names_over_ids=False,
                       verify=None, **kwargs):
    """
    Read a Table object from an VO table file

    Parameters
    ----------
    input : str or `~astropy.io.votable.tree.VOTableFile` or `~astropy.io.votable.tree.Table`
        If a string, the filename to read the table from. If a
        :class:`~astropy.io.votable.tree.VOTableFile` or
        :class:`~astropy.io.votable.tree.Table` object, the object to extract
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
        Additional keyword arguments are passed on to
        :func:`astropy.io.votable.table.parse`.
    """
    if not isinstance(input, (VOTableFile, VOTable)):
        input = parse(input, table_id=table_id, verify=verify, **kwargs)

    # Parse all table objects
    table_id_mapping = dict()
    tables = []
    if isinstance(input, VOTableFile):
        for table in input.iter_tables():
            if table.ID is not None:
                table_id_mapping[table.ID] = table
            tables.append(table)

        if len(tables) > 1:
            if table_id is None:
                raise ValueError(
                    "Multiple tables found: table id should be set via "
                    "the table_id= argument. The available tables are {}, "
                    'or integers less than {}.'.format(
                        ', '.join(table_id_mapping.keys()), len(tables)))
            elif isinstance(table_id, str):
                if table_id in table_id_mapping:
                    table = table_id_mapping[table_id]
                else:
                    raise ValueError(
                        f"No tables with id={table_id} found")
            elif isinstance(table_id, int):
                if table_id < len(tables):
                    table = tables[table_id]
                else:
                    raise IndexError(
                        "Table index {} is out of range. "
                        "{} tables found".format(
                            table_id, len(tables)))
        elif len(tables) == 1:
            table = tables[0]
        else:
            raise ValueError("No table found")
    elif isinstance(input, VOTable):
        table = input

    # Convert to an astropy.table.Table object
    return table.to_table(use_names_over_ids=use_names_over_ids)


def write_table_votable(input, output, table_id=None, overwrite=False,
                        tabledata_format=None):
    """
    Write a Table object to an VO table file

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
        raise ValueError('cannot write table with mixin column(s) {} to VOTable'
                         .format(unsupported_names))

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


io_registry.register_reader('votable', Table, read_table_votable)
io_registry.register_writer('votable', Table, write_table_votable)
io_registry.register_identifier('votable', Table, is_votable)
