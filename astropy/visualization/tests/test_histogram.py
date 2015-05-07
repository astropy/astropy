# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from numpy.testing import assert_allclose

try:
    import matplotlib.pyplot as plt
except ImportError:
    HAS_PLT = False
else:
    HAS_PLT = True

import numpy as np
from .. import hist
from ...stats import histogram


def test_hist_basic(rseed=0):
    rng = np.random.RandomState(rseed)
    x = rng.randn(100)

    for range in [None, (-2, 2)]:
        n1, bins1, patches1 = plt.hist(x, 10, range=range)
        n2, bins2, patches2 = hist(x, 10, range=range)

        assert_allclose(n1, n2)
        assert_allclose(bins1, bins2)


def test_hist_specify_ax(rseed=0):
    rng = np.random.RandomState(rseed)
    x = rng.randn(100)

    fig, ax = plt.subplots(2)
    n1, bins1, patches1 = hist(x, 10, ax=ax[0])
    assert patches1[0].axes is ax[0]

    n2, bins2, patches2 = hist(x, 10, ax=ax[1])
    assert patches2[0].axes is ax[1]


def test_hist_autobin(rseed=0):
    rng = np.random.RandomState(rseed)
    x = rng.randn(100)

    for bintype in [10, np.arange(-3, 3, 10),
                    'knuth', 'scotts', 'freedman', 'blocks', 'adaptive']:
        for range in [None, (-3, 3)]:
            n1, bins1 = histogram(x, bintype, range=range)
            n2, bins2, patches = hist(x, bintype, range=range)
            assert_allclose(n1, n2)
            assert_allclose(bins1, bins2)
