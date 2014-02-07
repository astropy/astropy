# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a collection of monkey patches and workarounds for bugs in
earlier versions of Numpy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six

import numpy as np

# Note: We could have used distutils.version for this comparison,
# but it seems like overkill to import distutils at runtime.
numpy_major, numpy_minor, numpy_rest = np.__version__.split(".", 2)
numpy_major = int(numpy_major)
numpy_minor = int(numpy_minor)


def _monkeypatch_unicode_mask_fill_values():
    """
    Numpy <= 1.7.0 on Python 2 does not support Unicode fill values, since
    it assumes that all of the string dtypes are ``S`` and ``V`` (not ``U``).

    This monkey patches the function that validates and corrects a
    fill value to handle this case.
    """
    if numpy_major == 1 and numpy_minor <= 7 and six.PY2:
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
