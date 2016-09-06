# -*- coding: utf-8 -*-

# TEST_UNICODE_LITERALS

# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

import numpy as np

from ...extern import six
from ..data_info import dtype_info_name
from ...tests.helper import pytest

STRING_TYPE_NAMES = {(False, 'S'): 'str',  # PY2
                     (False, 'U'): 'unicode',
                     (True, 'S'): 'bytes', # not PY2
                     (True, 'U'): 'str'}

DTYPE_TESTS = ((np.array(b'abcd').dtype, STRING_TYPE_NAMES[(not six.PY2, 'S')] + '4'),
               (np.array(u'abcd').dtype, STRING_TYPE_NAMES[(not six.PY2, 'U')] + '4'),
               ('S4', STRING_TYPE_NAMES[(not six.PY2, 'S')] + '4'),
               ('U4', STRING_TYPE_NAMES[(not six.PY2, 'U')] + '4'),
               (np.void, 'void'),
               (np.int32, 'int32'),
               (np.bool, 'bool'),
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
