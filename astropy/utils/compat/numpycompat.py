# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""
from astropy.utils import minversion


__all__ = ['NUMPY_LT_1_14', 'NUMPY_LT_1_14_1', 'NUMPY_LT_1_14_2',
           'NUMPY_LT_1_16', 'NUMPY_LT_1_17']

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_1_14 = not minversion('numpy', '1.14')
NUMPY_LT_1_14_1 = not minversion('numpy', '1.14.1')
NUMPY_LT_1_14_2 = not minversion('numpy', '1.14.2')
NUMPY_LT_1_16 = not minversion('numpy', '1.16')
NUMPY_LT_1_17 = not minversion('numpy', '1.17')
