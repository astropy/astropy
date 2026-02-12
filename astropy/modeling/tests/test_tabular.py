# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from astropy import units as u
from astropy.modeling import models

def test_tabular_uses_quantity():
    """
    Test that uses_quantity reports correctly based on lookup_table units.
    Regression test for #19293.
    """
    x = np.linspace(0, 10, 5)
    y = np.arange(5.0)

    # Unitless lookup_table
    m_unitless = models.Tabular1D(x, y)
    assert not m_unitless.uses_quantity
    assert m_unitless(5.0) == 2.0   
    assert m_unitless(0.0) == 0.0
    assert m_unitless(10.0) == 4.0

    # Quantity lookup_table
    m_qty = models.Tabular1D(x, y * u.m)
    assert m_qty.uses_quantity
    result = m_qty(5.0)
    assert result == 2.0 * u.m
    assert result.unit == u.m

    # Tabular2D unitless
    z = np.outer(y, y)
    m2d_unitless = models.Tabular2D((x, x), z)
    assert not m2d_unitless.uses_quantity
    assert m2d_unitless(5.0, 5.0) == 4.0   

    # Tabular2D with units
    m2d_qty = models.Tabular2D((x, x), z * u.s)
    assert m2d_qty.uses_quantity
    assert m2d_qty(5.0, 5.0).unit == u.s

@pytest.mark.parametrize("method", ["linear", "nearest"])
def test_tabular_evaluation_permissive(method):
    """
    Confirm permissive unit handling (no error on mismatch) is preserved.
    """
    points = np.array([0, 1, 2])
    lookup = np.array([10.0, 20.0, 30.0])

    m = models.Tabular1D(points, lookup, method=method)
    assert not m.uses_quantity

    assert m(1.0 * u.m) == 20.0

    m_qty = models.Tabular1D(points, lookup * u.K, method=method)
    assert m_qty.uses_quantity

    result = m_qty(1.0)
    assert result == 20.0 * u.K
    assert result.unit == u.K