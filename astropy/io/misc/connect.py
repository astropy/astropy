# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects any readers/writers defined in io.misc to the
# astropy.table.Table class

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .hdf5 import HDF5TableIO
