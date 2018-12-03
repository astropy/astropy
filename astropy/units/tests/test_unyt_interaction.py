# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

import numpy as np

from ... import units as u

unyt = pytest.importorskip("unyt")


class TestUnit:
    def test_simple(self):
        assert str(unyt.m == 'm')
        u_m = u.Unit(unyt.m)
        assert u_m is u.m

    def test_dimensionless(self):
        u_dimensionless = u.Unit(unyt.dimensionless)
        assert isinstance(u_dimensionless, u.UnitBase)
        assert u_dimensionless == u.dimensionless_unscaled

    def test_more_complicated(self):
        u_unyt = unyt.m / unyt.s
        u_m_per_s = u.Unit(u_unyt)
        assert isinstance(u_m_per_s, u.UnitBase)
        assert u_m_per_s == u.m / u.s


class TestQuantity:
    @pytest.mark.parametrize('unit', (unyt.m, unyt.dimensionless,
                                      unyt.m/unyt.s))
    @pytest.mark.parametrize('value', (1., np.arange(10.)))
    def test_simple(self, value, unit):
        q_unyt = value * unyt.m
        q = u.Quantity(q_unyt)
        assert type(q) is u.Quantity
        assert q.unit is u.m
        assert np.all(q.value == value)

    def test_conversion(self):
        q_unyt = 1. * u.m
        q = u.Quantity(q_unyt, u.cm)
        assert q.unit == u.cm
        assert q.value == 100.
        q = u.Quantity(q_unyt, unyt.cm)
        assert q.unit == u.cm
        assert q.value == 100.
