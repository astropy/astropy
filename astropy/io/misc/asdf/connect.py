# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
# This file connects ASDF to the astropy.table.Table class

import functools

import asdf

from astropy.io import registry as io_registry
from astropy.table import Table


def read_table(filename, data_key='data', find_table=None, **kwargs):

    with asdf.open(filename, **kwargs) as af:
        return af[data_key]


def write_table(table, filename, data_key='data', **kwargs):

    tree = { data_key : table }
    with asdf.AsdfFile(tree) as af:
        af.write_to(filename, **kwargs)


def io_identify(suffix, origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(suffix)


asdf_identify = functools.partial(io_identify, '.asdf')

io_registry.register_reader('asdf', Table, read_table)
io_registry.register_writer('asdf', Table, write_table)
io_registry.register_identifier('asdf', Table, asdf_identify)
