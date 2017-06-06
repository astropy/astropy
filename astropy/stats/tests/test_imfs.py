
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
#from numpy.testing import assert_equal
#from numpy.testing.utils import assert_allclose
from ...tests.helper import pytest

from .. import imfs

try:
    from scipy.integrate import quad
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.parametrize('pdf',((imfs.Kroupa.pdf,imfs.Salpeter.pdf),))
def test_normalization(pdf):
    abstol = 0.0001
    integral,error = quad(pdf, 0.03, 1000, abstol=abstol)
    assert np.abs(integral-1) < abstol+error

@pytest.mark.parametrize('cdf',((imfs.Kroupa.cdf,imfs.Salpeter.cdf),))
def test_cdf(cdf):
    assert cdf(0.03) == 0
    assert cdf(np.inf) == 1
