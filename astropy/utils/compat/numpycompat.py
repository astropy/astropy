# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six
from ...utils import minversion

import numpy as np


__all__ = ['NUMPY_LT_1_6_1', 'NUMPY_LT_1_7', 'NUMPY_LT_1_8', 'NUMPY_LT_1_9',
           'NUMPY_LT_1_9_1']

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_1_6_1 = not minversion(np, '1.6.1')
NUMPY_LE_1_7 = not minversion(np, '1.7.0', inclusive=False)
NUMPY_LT_1_7 = not minversion(np, '1.7.0')
NUMPY_LT_1_8 = not minversion(np, '1.8.0')
NUMPY_LT_1_9 = not minversion(np, '1.9.0')
NUMPY_LT_1_9_1 = not minversion(np, '1.9.1')


def _monkeypatch_unicode_mask_fill_values():
    """
    Numpy <= 1.7.0 on Python 2 does not support Unicode fill values, since
    it assumes that all of the string dtypes are ``S`` and ``V`` (not ``U``).

    This monkey patches the function that validates and corrects a
    fill value to handle this case.
    """
    if NUMPY_LE_1_7 and six.PY2:
        from numpy.ma import core as ma_core
        _check_fill_value_original = ma_core._check_fill_value

        def _check_fill_value(fill_value, ndtype):
            if (not ndtype.fields and
                isinstance(fill_value, six.string_types) and
                ndtype.char in 'SVU'):
                return np.array(fill_value, copy=False, dtype=ndtype)
            return _check_fill_value_original(fill_value, ndtype)

        ma_core._check_fill_value = _check_fill_value


_monkeypatch_unicode_mask_fill_values()
