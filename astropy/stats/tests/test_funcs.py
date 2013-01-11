# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division

import numpy as np
from ...tests.helper import pytest

from .. import funcs
from ...utils.misc import NumpyRNGContext

try:
    from scipy import stats
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


def test_sigma_clip():
    from numpy.random import randn

    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):  # Amazing, I've got the same combination on my luggage!

        randvar = randn(10000)

        data, mask = funcs.sigma_clip(randvar, 1, 2)
        maskedarr = funcs.sigma_clip(randvar, 1, 2, maout=True)

        assert sum(mask) > 0
        assert data.size < randvar.size
        assert np.all(mask == ~maskedarr.mask)

        #this is actually a silly thing to do, because it uses the standard
        #deviation as the variance, but it tests to make sure these arguments
        #are actually doing something
        data2, mask2 = funcs.sigma_clip(randvar, 1, 2, varfunc=np.std)
        assert not np.all(data == data2)
        assert not np.all(mask == mask2)

        data3, mask3 = funcs.sigma_clip(randvar, 1, 2, cenfunc=np.mean)
        assert not np.all(data == data3)
        assert not np.all(mask == mask3)

        #now just make sure the iters=None method works at all.
        maskedarr = funcs.sigma_clip(randvar, 3, None, maout=True)


@pytest.mark.skipif('not HAS_SCIPY')
def test_compare_to_scipy_sigmaclip():
    from numpy.random import randn
    from numpy.testing import assert_equal

    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):  # Amazing, I've got the same combination on my luggage!

        randvar = randn(10000)

        astropyres = funcs.sigma_clip(randvar, 3, None, np.mean)[0]
        scipyres = stats.sigmaclip(randvar, 3, 3)[0]

        assert_equal(astropyres, scipyres)

