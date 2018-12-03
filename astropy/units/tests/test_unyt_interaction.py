# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

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
