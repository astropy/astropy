# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
# This file connects ASDF to the astropy.table.Table class

import functools

import asdf

from astropy.io import registry as io_registry
from astropy.table import Table


def read_table(filename, data_key=None, find_table=None, **kwargs):

    if data_key and find_table:
        raise ValueError("Options 'data_key' and 'find_table' are not compatible")

    with asdf.open(filename, **kwargs) as af:
        if find_table:
            return find_table(af)
        else:
            return af[data_key or 'data']


def write_table(table, filename, data_key=None, make_tree=None, **kwargs):

    if data_key and make_tree:
        raise ValueError("Options 'data_key' and 'make_tree' are not compatible")

    if make_tree:
        tree = make_tree(table)
    else:
        tree = { data_key or 'data' : table }

    with asdf.AsdfFile(tree) as af:
        af.write_to(filename, **kwargs)


def io_identify(suffix, origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(suffix)


asdf_identify = functools.partial(io_identify, '.asdf')

io_registry.register_reader('asdf', Table, read_table)
io_registry.register_writer('asdf', Table, write_table)
io_registry.register_identifier('asdf', Table, asdf_identify)
