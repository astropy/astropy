# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects any readers/writers defined in io.misc to the
# astropy.table.Table class

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .hdf5 import read_table_hdf5, write_table_hdf5, is_hdf5
from .html import html_read, html_write, html_identify

from .. import registry as io_registry
from ...table import Table

io_registry.register_reader('hdf5', Table, read_table_hdf5)
io_registry.register_writer('hdf5', Table, write_table_hdf5)
io_registry.register_identifier('hdf5', Table, is_hdf5)
io_registry.register_reader('html', Table, html_read)
io_registry.register_writer('html', Table, html_write)
io_registry.register_identifier('html', Table, html_identify)
