# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
# This file connects ASDF to the astropy.table.Table class

import functools

from astropy.io import registry as io_registry
from astropy.table import Table


def read_table(filename, data_key=None, find_table=None, **kwargs):
    """
    Read a `~astropy.table.Table` object from an ASDF file

    This requires `asdf <https://pypi.org/project/asdf/>`_ to be installed.
    By default, this function will look for a Table object with the key of
    ``data`` in the top-level ASDF tree. The parameters ``data_key`` and
    ``find_key`` can be used to override the default behavior.

    This function is registered as the Table reader for ASDF files with the
    unified I/O interface.

    Parameters
    ----------
    filename : str or :class:`py.lath:local`
        Name of the file to be read
    data_key : str
        Optional top-level key to use for finding the Table in the tree. If not
        provided, uses ``data`` by default. Use of this parameter is not
        compatible with ``find_table``.
    find_table : function
        Optional function to be used for locating the Table in the tree. The
        function takes a single parameter, which is a dictionary representing
        the top of the ASDF tree. The function must return a
        `~astropy.table.Table` instance.

    Returns
    -------
    table : `~astropy.table.Table`
        `~astropy.table.Table` instance
    """
    try:
        import asdf
    except ImportError:
        raise Exception(
            "The asdf module is required to read and write ASDF files")

    if data_key and find_table:
        raise ValueError("Options 'data_key' and 'find_table' are not compatible")

    with asdf.open(filename, **kwargs) as af:
        if find_table:
            return find_table(af.tree)
        else:
            return af[data_key or 'data']


def write_table(table, filename, data_key=None, make_tree=None, **kwargs):
    """
    Write a `~astropy.table.Table` object to an ASDF file.

    This requires `asdf <https://pypi.org/project/asdf/>`_ to be installed.
    By default, this function will write a Table object in the top-level ASDF
    tree using the key of ``data``. The parameters ``data_key`` and
    ``make_tree`` can be used to override the default behavior.

    This function is registered as the Table writer for ASDF files with the
    unified I/O interface.

    Parameters
    ----------
    table : `~astropy.table.Table`
        `~astropy.table.Table` instance to be written
    filename : str or :class:`py.path:local`
        Name of the new ASDF file to be created
    data_key : str
        Optional top-level key in the ASDF tree to use when writing the Table.
        If not provided, uses ``data`` by default. Use of this parameter is not
        compatible with ``make_tree``.
    make_tree : function
        Optional function to be used for creating the ASDF tree. The function
        takes a single parameter, which is the `~astropy.table.Table` instance
        to be written. The function must return a `dict` representing the ASDF
        tree to be created.
    """
    try:
        import asdf
    except ImportError:
        raise Exception(
            "The asdf module is required to read and write ASDF files")

    if data_key and make_tree:
        raise ValueError("Options 'data_key' and 'make_tree' are not compatible")

    if make_tree:
        tree = make_tree(table)
    else:
        tree = {data_key or 'data' : table}

    with asdf.AsdfFile(tree) as af:
        af.write_to(filename, **kwargs)


def asdf_identify(origin, filepath, fileobj, *args, **kwargs):
    try:
        import asdf
    except ImportError:
        return False

    return filepath is not None and filepath.endswith('.asdf')


io_registry.register_reader('asdf', Table, read_table)
io_registry.register_writer('asdf', Table, write_table)
io_registry.register_identifier('asdf', Table, asdf_identify)
