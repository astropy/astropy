# -*- coding: utf-8 -*-


# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from ..data_info import dtype_info_name

STRING_TYPE_NAMES = {(True, 'S'): 'bytes',
                     (True, 'U'): 'str'}

DTYPE_TESTS = ((np.array(b'abcd').dtype, STRING_TYPE_NAMES[(True, 'S')] + '4'),
               (np.array(u'abcd').dtype, STRING_TYPE_NAMES[(True, 'U')] + '4'),
               ('S4', STRING_TYPE_NAMES[(True, 'S')] + '4'),
               ('U4', STRING_TYPE_NAMES[(True, 'U')] + '4'),
               (np.void, 'void'),
               (np.int32, 'int32'),
               (bool, 'bool'),
               (float, 'float64'),
               ('<f4', 'float32'),
               ('u8', 'uint64'),
               ('c16', 'complex128'),
               ('object', 'object'))


@pytest.mark.parametrize('input,output', DTYPE_TESTS)
def test_dtype_info_name(input, output):
    """
    Test that dtype_info_name is giving the expected output

    Here the available types::

      'b' boolean
      'i' (signed) integer
      'u' unsigned integer
      'f' floating-point
      'c' complex-floating point
      'O' (Python) objects
      'S', 'a' (byte-)string
      'U' Unicode
      'V' raw data (void)
    """
    assert dtype_info_name(input) == output
