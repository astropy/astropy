# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""

import numpy as np
from astropy.utils import minversion

__all__ = ['NUMPY_LT_1_19', 'NUMPY_LT_1_19_1', 'NUMPY_LT_1_20',
           'NUMPY_LT_1_21_1', 'NUMPY_LT_1_22', 'NUMPY_LT_1_22_1',
           'NUMPY_LT_1_23']

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_1_19 = not minversion(np, '1.19')
NUMPY_LT_1_19_1 = not minversion(np, '1.19.1')
NUMPY_LT_1_20 = not minversion(np, '1.20')
NUMPY_LT_1_21_1 = not minversion(np, '1.21.1')
NUMPY_LT_1_22 = not minversion(np, '1.22')
NUMPY_LT_1_22_1 = not minversion(np, '1.22.1')
NUMPY_LT_1_23 = not minversion(np, '1.23dev0')
