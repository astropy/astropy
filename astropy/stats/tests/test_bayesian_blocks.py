# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ...tests.helper import pytest
from  numpy.testing import assert_allclose
from .. import bayesian_blocks, Events, RegularEvents, PointMeasures


def test_single_change_point(rseed=0):
    rng = np.random.RandomState(rseed)
    x = np.concatenate([rng.rand(100),
                        1 + rng.rand(200)])

    bins = bayesian_blocks(x)

    assert (len(bins) == 3)
    assert_allclose(bins[1], 1, rtol=0.02)


def test_duplicate_events(rseed=0):
    rng = np.random.RandomState(rseed)
    t = rng.rand(100)
    t[80:] = t[:20]

    x = np.ones_like(t)
    x[:20] += 1

    bins1 = bayesian_blocks(t)
    bins2 = bayesian_blocks(t[:80], x[:80])

    assert_allclose(bins1, bins2)


def test_measures_fitness_homoscedastic(rseed=0):
    rng = np.random.RandomState(rseed)
    t = np.linspace(0, 1, 11)
    x = np.exp(-0.5 * (t - 0.5) ** 2 / 0.01 ** 2)
    sigma = 0.05
    x = x + sigma * rng.randn(len(x))

    bins = bayesian_blocks(t, x, sigma, fitness='measures')

    assert_allclose(bins, [0, 0.45, 0.55, 1])


def test_measures_fitness_heteroscedastic():
    rng = np.random.RandomState(1)
    t = np.linspace(0, 1, 11)
    x = np.exp(-0.5 * (t - 0.5) ** 2 / 0.01 ** 2)
    sigma = 0.02 + 0.02 * rng.rand(len(x))
    x = x + sigma * rng.randn(len(x))

    bins = bayesian_blocks(t, x, sigma, fitness='measures')

    assert_allclose(bins, [0, 0.45, 0.55, 1])


def test_regular_events():
    rng = np.random.RandomState(0)
    dt = 0.01
    steps = np.concatenate([np.unique(rng.randint(0, 500, 100)),
                            np.unique(rng.randint(500, 1000, 200))])
    t = dt * steps

    # string fitness
    bins1 = bayesian_blocks(t, fitness='regular_events', dt=dt)
    assert (len(bins1) == 3)
    assert_allclose(bins1[1], 5, rtol=0.05)

    # class name fitness
    bins2 = bayesian_blocks(t, fitness=RegularEvents, dt=dt)
    assert_allclose(bins1, bins2)

    # class instance fitness
    bins3 = bayesian_blocks(t, fitness=RegularEvents(dt=dt))
    assert_allclose(bins1, bins3)


def test_errors():
    rng = np.random.RandomState(0)
    t = rng.rand(100)

    # x must be integer or None for events
    with pytest.raises(ValueError):
        bayesian_blocks(t, fitness='events', x=t)

    # x must be binary for regular events
    with pytest.raises(ValueError):
        bayesian_blocks(t, fitness='regular_events', x=10 * t, dt=1)
    
    # x must be specified for measures
    with pytest.raises(ValueError):
        bayesian_blocks(t, fitness='measures')

    # sigma cannot be specified without x
    with pytest.raises(ValueError):
        bayesian_blocks(t, fitness='events', sigma=0.5)

    # length of x must match length of t
    with pytest.raises(ValueError):
        bayesian_blocks(t, fitness='measures', x=t[:-1])

    # repeated values in t fail when x is specified
    t2 = t.copy()
    t2[1] = t2[0]
    with pytest.raises(ValueError):
        bayesian_blocks(t2, fitness='measures', x=t)

    # sigma must be broadcastable with x
    with pytest.raises(ValueError):
        bayesian_blocks(t, fitness='measures',  x=t, sigma=t[:-1])
    
