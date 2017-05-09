# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...utils import minversion


__all__ = ['NUMPY_LT_1_9_1', 'NUMPY_LT_1_10', 'NUMPY_LT_1_10_4',
           'NUMPY_LT_1_11', 'NUMPY_LT_1_12']

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_1_9_1 = not minversion('numpy', '1.9.1')
NUMPY_LT_1_10 = not minversion('numpy', '1.10.0')
NUMPY_LT_1_10_4 = not minversion('numpy', '1.10.4')
NUMPY_LT_1_11 = not minversion('numpy', '1.11.0')
NUMPY_LT_1_12 = not minversion('numpy', '1.12')
