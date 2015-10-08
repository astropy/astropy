# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ..nddata import NDData
from ..nduncertainty import StdDevUncertainty, NDUncertainty
from ...tests.helper import pytest, raises


def test_initializing_from_stddevuncertainty():
    u1 = StdDevUncertainty(np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(u1, copy=False)

    assert u1.array is u2.array
