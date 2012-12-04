# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os

from . import parse
from .tree import VOTableFile, Resource
from .tree import Table as VOTable
from ...utils import OrderedDict
from ...table import io_registry


def is_votable(origin, args, kwargs):
    pass


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
    id : str
        The ID of the table to read in
    """

    if isinstance(input, basestring):
        input = parse(input)

    # Parse all table objects
    tables = OrderedDict()
    if isinstance(input, VOTableFile):
        for table in input.iter_tables():
            tables[table.ID] = table

    if len(tables) > 1:
        if table_id is None:
            raise ValueError("Multiple tables found: table id should be set via the id= argument. The available tables are " + ', '.join(tables.keys()))
        else:
            if table_id in tables:
                table = tables[table_id]
            else:
                raise ValueError("No tables with id={0} found".format(table_id))
    elif len(tables) == 1:
        table = tables[tables.keys()[0]]
    else:
        raise ValueError("No table found")

    # Convert to an astropy.table.Table object
    return table.to_table()


def write_table_votable(input, output, table_id=None, compression=False, overwrite=False):
    """
    Write a Table object to an VO table file

    This requires `h5py <http://alfven.org/wp/hdf5-for-python/>`_ to be
    installed.

    Parameters
    ----------
    output : str
        The filename to write the table to.
    compression : bool
        Whether to compress the table inside the HDF5 file.
    table_id : str
        The table ID to use. If this is not specified, the 'ID' keyword in the
        ``meta`` object of the table will be used.
    overwrite : bool
        Whether to overwrite any existing file without warning.
    """

    if os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise IOError("File exists: {0}".format(output))

    # Convert to VOTableFile object
    table_file = VOTableFile()
    resource = Resource()
    table_file.resources.append(resource)
    table = VOTable.from_table(table_file, input)
    table.ID = table_id
    resource.tables.append(table)
    table_file.to_xml(output)


io_registry.register_reader('votable', read_table_votable)
io_registry.register_writer('votable', write_table_votable)
io_registry.register_identifier('votable', is_votable)
