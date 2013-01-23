# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects any readers/writers defined in io.misc to the
# astropy.table.Table class

from ...table import io_registry

from .hdf5 import read_table_hdf5, write_table_hdf5, is_hdf5

io_registry.register_reader('hdf5', read_table_hdf5)
io_registry.register_writer('hdf5', write_table_hdf5)
io_registry.register_identifier('hdf5', is_hdf5)
