# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os

from . import parse, from_table
from .tree import VOTableFile, Table
from ...utils import OrderedDict
from ...table import io_registry


def is_votable(origin, args, kwargs):
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
        if isinstance(args[0], basestring) or hasattr(args[0], 'read'):
            return is_votable(args[0])
        elif isinstance(args[0], (VOTableFile, Table)):
            return True
        else:
            return False
    else:
        return False


def read_table_votable(input, table_id=None):
    """
    Read a Table object from an VO table file

    Parameters
    ----------
    input : str or `astropy.io.votable.tree.VOTableFile` or `astropy.io.votable.tree.Table`
        If a string, the filename to read the table from. If a
        :class:`~astropy.io.votable.tree.VOTableFile` or
        :class:`~astropy.io.votable.tree.Table` object, the object to extract
        the table from.
    table_id : str, optional
        The ID of the table to read in.
    """
    if not isinstance(input, (VOTableFile, Table)):
        input = parse(input, table_id=table_id)

    # Parse all table objects
    tables = OrderedDict()
    if isinstance(input, VOTableFile):
        for table in input.iter_tables():
            if table.ID is not None:
                tables[table.ID] = table

        if len(tables) > 1:
            if table_id is None:
                raise ValueError(
                    "Multiple tables found: table id should be set via "
                    "the table_id= argument. The available tables are " +
                    ', '.join(tables.keys()))
            else:
                if table_id in tables:
                    table = tables[table_id]
                else:
                    raise ValueError(
                        "No tables with id={0} found".format(table_id))
        elif len(tables) == 1:
            table = tables[tables.keys()[0]]
        else:
            raise ValueError("No table found")
    elif isinstance(input, Table):
        table = input

    # Convert to an astropy.table.Table object
    return table.to_table()


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


io_registry.register_reader('votable', read_table_votable)
io_registry.register_writer('votable', write_table_votable)
io_registry.register_identifier('votable', is_votable)
