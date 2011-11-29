from __future__ import division

import numpy as np
from .. import alg


def test_sigma_clip():
    from numpy.random import randn

    randvar = randn(10000)

    data, mask = alg.sigma_clip(randvar, 1, 2)
    maskedarr = alg.sigma_clip(randvar, 1, 2, maout=True)

    #in principle it's *possible* that these asserts might fail if the random
    #number generator produces a highly non-gaussian distribution, but it's
    #highly unlikely
    assert sum(mask) > 0
    assert data.size < randvar.size
    assert np.all(mask == ~maskedarr.mask)

    #this is actually a silly thing to do, because it uses the standard
    #deviation as the variance, but it tests to make sure these arguments
    #are actually doing something
    data2, mask2 = alg.sigma_clip(randvar, 1, 2, varfunc=np.std)
    assert not np.all(data == data2)
    assert not np.all(mask == mask2)

    data3, mask3 = alg.sigma_clip(randvar, 1, 2, cenfunc=np.mean)
    assert not np.all(data == data3)
    assert not np.all(mask == mask3)

    #now just make sure the iters=None method works at all.
    maskedarr = alg.sigma_clip(randvar, 3, None, maout=True)
