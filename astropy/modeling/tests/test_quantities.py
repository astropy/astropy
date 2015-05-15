# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests specifically for models that use units and quantities."""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing import assert_allclose

from ..models import Gaussian1D
from ... import units as u
from ...units import UnitsError, Quantity
from ...tests.helper import pytest


def test_quantities_as_parameters():
    """
    Basic tests for initializing general models (that do not require units)
    with parameters that have units attached.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)
    assert g.amplitude.value == 1.0
    assert g.amplitude.unit is u.J
    assert g.mean.value == 1.0
    assert g.mean.unit is u.m
    assert g.stddev.value == 0.1
    assert g.stddev.unit is u.m


def test_quantity_parameter_arithmetic():
    """
    Test that arithmetic operations with properties that have units return the
    appropriate Quantities.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    assert g.mean + (1 * u.m) == 2 * u.m
    with pytest.raises(UnitsError):
        g.mean + 1
    assert (1 * u.m) + g.mean == 2 * u.m
    with pytest.raises(UnitsError):
        1 + g.mean
    assert g.mean * 2 == (2 * u.m)
    assert 2 * g.mean == (2 * u.m)
    assert g.mean * (2 * u.m) == (2 * (u.m ** 2))
    assert (2 * u.m) * g.mean == (2 * (u.m ** 2))

    assert -g.mean == (-1 * u.m)
    assert abs(-g.mean) == g.mean


def test_quantity_parameter_comparison():
    """
    Basic test of comparison operations on properties with units.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    assert g.mean == 1 * u.m
    assert 1 * u.m == g.mean
    assert g.mean != 1
    assert 1 != g.mean

    assert g.mean < 2 * u.m
    assert 2 * u.m > g.mean
    with pytest.raises(UnitsError):
        g.mean < 2
    with pytest.raises(UnitsError):
        2 > g.mean

    g = Gaussian1D([1, 2] * u.J, [1, 2] * u.m, [0.1, 0.2] * u.m)

    assert g.mean == [1, 2] * u.m
    assert np.all([1, 2] * u.m == g.mean)
    assert g.mean != [1, 2]
    assert np.all([1, 2] != g.mean)
    with pytest.raises(UnitsError):
        g.mean < [3, 4]
    with pytest.raises(UnitsError):
        [3, 4] > g.mean


def test_basic_evaluate_with_quantities():
    """
    Test evaluation of a single model with Quantity parameters, that does
    not explicitly require units.
    """

    g = Gaussian1D(1, 1, 0.1)
    gq = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    assert isinstance(gq(0), Quantity)
    assert gq(0).unit is u.J
    assert g(0) == gq(0).value

    # zero is allowed without explicit units, but other unitless quantities
    # should be an exception
    with pytest.raises(UnitsError):
        gq(1)

    assert gq(1 * u.m).value == g(1)

    # Should get the same numeric result as if we multiplied by 1000
    assert_allclose(gq(0.0005 * u.km).value, g(0.5))
