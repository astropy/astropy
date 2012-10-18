from __future__ import division

import numpy as np
from .. import funcs
from pytest import mark

try:
    from scipy import stats
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


class NumpyRNGContext(object):
    """
    A context manager (for use with the ``with`` statement) that will seed the
    numpy random number generator (RNG) to a specific value, and then restore
    the RNG state back to whatever it was before.

    This is primarily intended for use in the astropy testing suit, but it
    may be useful in ensuring reproducability of monte carlo simulations in a
    science context.

    Parameters
    ----------
    seed : int
        The value to use to seed the numpy RNG

    Examples
    --------
    A tyoical use case might be::

        with NumpyRNGContext(<some seed value you pick>):
            from numpy import random

            randarr = random.randn(100)
            ... run your test using `randarr` ...

        #Any code using numpy.random at this indent level will act just as it
        #would have if it had been before the with statement - e.g. whatever
        #the default seed is.


    """
    def __init__(self, seed):
        self.seed = seed

    def __enter__(self):
        from numpy import random

        self.startstate = random.get_state()
        random.seed(self.seed)

    def __exit__(self, exc_type, exc_value, traceback):
        from numpy import random

        random.set_state(self.startstate)


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


@mark.skipif('not HAS_SCIPY')
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

