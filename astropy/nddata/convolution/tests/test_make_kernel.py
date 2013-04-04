# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from ....tests.helper import pytest

from ..make_kernel import make_kernel

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

class TestMakeKernel(object):
    """
    Test the make_kernel function
    """

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_airy(self):
        """
        Test kerneltype airy, a.k.a. brickwall
        Checks https://github.com/astropy/astropy/pull/939
        """
        k1 = make_kernel([3, 3], kernelwidth=0.5, kerneltype='airy')
        k2 = make_kernel([3, 3], kernelwidth=0.5, kerneltype='brickwall')
        ref = np.array([[ 0.06375119,  0.12992753,  0.06375119],
                        [ 0.12992753,  0.22528514,  0.12992753],
                        [ 0.06375119,  0.12992753,  0.06375119]])
        assert_allclose(k1, ref, rtol=0, atol=1e-7)
        assert_equal(k1, k2)
