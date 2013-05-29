# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import io
import os

from . import parse, from_table
from .tree import VOTableFile, Table as VOTable
from ...utils import OrderedDict
from .. import registry as io_registry
from ...table import Table


def is_votable(origin, filepath, fileobj, *args, **kwargs):
    """
    Reads the header of a file to determine if it is a VOTable file.

    Parameters
    ----------
    origin : str or readable file-like object
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


def read_table_votable(input, table_id=None, use_names_over_ids=False):
    """
    Read a Table object from an VO table file

    Parameters
    ----------
    input : str or `astropy.io.votable.tree.VOTableFile` or `astropy.io.votable.tree.Table`
        If a string, the filename to read the table from. If a
        :class:`~astropy.io.votable.tree.VOTableFile` or
        :class:`~astropy.io.votable.tree.Table` object, the object to extract
        the table from.

    table_id : str or int, optional
        The table to read in.  If a `str`, it is an ID corresponding
        to the ID of the table in the file (not all VOTable files
        assign IDs to their tables).  If an `int`, it is the index of
        the table in the file, starting at 0.

    use_names_over_ids : boolean, optional
        When `True` use the `name` attributes of columns as the names
        of columns in the `astropy.table.Table` instance.  Since names
        are not guaranteed to be unique, this may cause some columns
        to be renamed by appending numbers to the end.  Otherwise
        (default), use the ID attributes as the column names.
    """
    if not isinstance(input, (VOTableFile, VOTable)):
        input = parse(input, table_id=table_id)

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
                    "the table_id= argument. The available tables are {0}, "
                    'or integers less than {1}.'.format(
                        ', '.join(table_id_mapping.keys()), len(tables)))
            elif isinstance(table_id, basestring):
                if table_id in table_id_mapping:
                    table = table_id_mapping[table_id]
                else:
                    raise ValueError(
                        "No tables with id={0} found".format(table_id))
            elif isinstance(table_id, (int, long)):
                if table_id < len(tables):
                    table = tables[table_id]
                else:
                    raise IndexError(
                        "Table index {0} is out of range. "
                        "{1} tables found".format(
                            table_id, len(tables)))
        elif len(tables) == 1:
            table = tables[0]
        else:
            raise ValueError("No table found")
    elif isinstance(input, VOTable):
        table = input

    # Convert to an astropy.table.Table object
    return table.to_table(use_names_over_ids=use_names_over_ids)


def write_table_votable(input, output, table_id=None, overwrite=False):
    """
    Write a Table object to an VO table file

    Parameters
    ----------
    input : Table
        The table to write out.
    output : str
        The filename to write the table to.
    table_id : str
        The table ID to use. If this is not specified, the 'ID' keyword in the
        ``meta`` object of the table will be used.
    overwrite : bool
        Whether to overwrite any existing file without warning.
    """

    # Check if output file already exists
    if isinstance(output, basestring) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise IOError("File exists: {0}".format(output))

    # Create a new VOTable file
    table_file = from_table(input, table_id=table_id)

    # Write out file
    table_file.to_xml(output)


io_registry.register_reader('votable', Table, read_table_votable)
io_registry.register_writer('votable', Table, write_table_votable)
io_registry.register_identifier('votable', Table, is_votable)
