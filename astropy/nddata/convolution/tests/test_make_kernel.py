# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from ....tests.compat import assert_allclose
from ....tests.helper import pytest

from ..make_kernel import make_kernel

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

@pytest.mark.skipif('not HAS_SCIPY')
def test_airy():
    """
    Test kerneltype airy, a.k.a. brickwall
    Checks https://github.com/astropy/astropy/pull/939
    """
    k1 = make_kernel([3, 3], kernelwidth=0.5, kerneltype='airy')
    ref = np.array([[ 0.06375119,  0.12992753,  0.06375119],
                    [ 0.12992753,  0.22528514,  0.12992753],
                    [ 0.06375119,  0.12992753,  0.06375119]])
    assert_allclose(k1, ref, rtol=0, atol=1e-7)
