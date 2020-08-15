# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""
from astropy.utils import minversion


__all__ = ['NUMPY_LT_1_17', 'NUMPY_LT_1_18', 'NUMPY_LT_1_19',
           'NUMPY_LT_1_20']

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_1_17 = not minversion('numpy', '1.17')
NUMPY_LT_1_18 = not minversion('numpy', '1.18')
NUMPY_LT_1_19 = not minversion('numpy', '1.19')
NUMPY_LT_1_20 = not minversion('numpy', '1.20')
