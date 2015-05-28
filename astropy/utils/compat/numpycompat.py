# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six
from ...utils import minversion


__all__ = ['NUMPY_LT_1_6_1', 'NUMPY_LT_1_7', 'NUMPY_LT_1_8', 'NUMPY_LT_1_9',
           'NUMPY_LT_1_9_1', 'NUMPY_LT_1_10']

# TODO: It might also be nice to have aliases to these named for specific
# features/bugs we're checking for (ex:
# astropy.table.table._BROKEN_UNICODE_TABLE_SORT)
NUMPY_LT_1_6_1 = not minversion('numpy', '1.6.1')
NUMPY_LT_1_7 = not minversion('numpy', '1.7.0')
NUMPY_LT_1_8 = not minversion('numpy', '1.8.0')
NUMPY_LT_1_9 = not minversion('numpy', '1.9.0')
NUMPY_LT_1_9_1 = not minversion('numpy', '1.9.1')
NUMPY_LT_1_10 = not minversion('numpy', '1.9.9999', inclusive=False)


def _monkeypatch_unicode_mask_fill_values():
    """
    Numpy < 1.8.0 on Python 2 does not support Unicode fill values, since
    it assumes that all of the string dtypes are ``S`` and ``V`` (not ``U``).

    This monkey patches the function that validates and corrects a
    fill value to handle this case.
    """
    if NUMPY_LT_1_8 and six.PY2:
        import numpy as np
        from numpy.ma import core as ma_core
        _check_fill_value_original = ma_core._check_fill_value

        def _check_fill_value(fill_value, ndtype):
            if (not ndtype.fields and
                isinstance(fill_value, six.string_types) and
                ndtype.char in 'SVU'):
                return np.array(fill_value, copy=False, dtype=ndtype)
            return _check_fill_value_original(fill_value, ndtype)

        ma_core._check_fill_value = _check_fill_value


def _register_patched_dtype_reduce():
    """
    Numpy < 1.7 has a bug when copying/pickling dtype objects with a
    zero-width void type--i.e. ``np.dtype('V0')``.  Specifically, although
    creating a void type is perfectly valid, it crashes when instantiating
    a dtype using a format string of 'V0', which is what is normally returned
    by dtype.__reduce__() for these dtypes.

    See https://github.com/astropy/astropy/pull/3283#issuecomment-81667461
    """

    if NUMPY_LT_1_7:
        import numpy as np
        import copy_reg

        # Originally this created an alternate constructor that fixed this
        # issue, and returned that constructor from the new reduce_dtype;
        # however that broke pickling since functions can't be pickled, so now
        # we fix the issue directly within the custom __reduce__

        def reduce_dtype(obj):
            info = obj.__reduce__()
            args = info[1]
            if args[0] == 'V0':
                args = ('V',) + args[1:]
                info = (info[0], args) + info[2:]
            return info

        copy_reg.pickle(np.dtype, reduce_dtype)


if not _ASTROPY_SETUP_:
    _monkeypatch_unicode_mask_fill_values()
    _register_patched_dtype_reduce()
