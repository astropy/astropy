# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from ... import units as u
from ... import constants as c
from .. import UQuantity


def test_initialisation():
    v1 = UQuantity(5, 2, u.km)
    assert v1.value == 5
    assert v1.unit == u.km
    assert v1.uncertainty == 2
    v2 = UQuantity(5 * u.km, 2)
    assert v2.value == 5
    assert v2.unit == u.km
    assert v2.uncertainty == 2
    v3 = UQuantity(5 * u.km, 2000 * u.m)
    assert v3.value == 5
    assert v3.unit == u.km
    assert v3.uncertainty == 2
    v4 = UQuantity(np.arange(5.), 2., u.km)
    assert np.all(v4.value == np.arange(5.))
    assert v4.unit == u.km
    assert v4.uncertainty == 2
    v5 = UQuantity(np.arange(5.), np.array([1., 2., 1., 2., 1.]), u.km)
    assert np.all(v5.value == np.arange(5.))
    assert v5.unit == u.km
    assert np.all(v5.uncertainty == np.array([1., 2., 1., 2., 1.]))


class TestBasics():
    def setup(self):
        self.v = UQuantity(5., 2., u.km)
        self.a = UQuantity(np.arange(1., 5.), 1., u.s)
        self.b = UQuantity(np.array([1., 2., 3.]), np.array([0.1, 0.2, 0.1]),
                           u.m)

    def test_addition(self):
        c1 = self.v + UQuantity(12, 5, self.v.unit)
        assert c1.value == self.v.value + 12
        assert c1.unit == u.km
        # Uncertainties under addition add in quadrature
        assert np.allclose(c1.uncertainty,
                           np.sqrt(self.v.uncertainty**2 + 5**2))
        # now with different units
        c2 = self.v + UQuantity(12000., 5000., u.m)
        assert c2.value == self.v.value + 12
        assert c2.unit == u.km
        assert np.allclose(c2.uncertainty,
                           np.sqrt(self.v.uncertainty**2 + 5**2))
        # try array
        c3 = self.v + self.b
        assert np.all(c3.nominal_value ==
                      self.v.nominal_value + self.b.nominal_value)
        assert np.allclose(c3.uncertainty,
                           np.sqrt((self.v.uncertainty * self.v.unit)**2 +
                                   (self.b.uncertainty * self.b.unit)**2)
                           .to(c3.unit).value)
        # try adding regular Quantity
        c4 = self.v + 10. * self.v.unit
        assert c4.value == self.v.value + 10.
        assert c4.uncertainty == self.v.uncertainty

    def test_subtraction(self):
        c1 = self.v - UQuantity(12, 5, self.v.unit)
        assert c1.value == self.v.value - 12
        assert c1.unit == u.km
        # Uncertainties under addition add in quadrature
        assert np.allclose(c1.uncertainty,
                           np.sqrt(self.v.uncertainty**2 + 5**2))

    def test_multiplication(self):
        c1 = self.v * self.a
        assert np.all(c1.nominal_value ==
                      self.v.nominal_value * self.a.nominal_value)

        # Fractional uncertainties under multiplication add in quadrature
        assert np.allclose(c1.uncertainty/c1.value,
                           np.sqrt((self.v.uncertainty/self.v.value)**2 +
                                   (self.a.uncertainty/self.a.value)**2))
        # Test multiplication with straight Quantity
        c2 = self.a * (10. * u.s)
        assert np.all(c2.value == self.a.value * 10.)
        assert c2.unit == self.a.unit * u.s
        assert np.all(c2.uncertainty == self.a.uncertainty * 10.)

    def test_division(self):
        c1 = self.v / self.a
        assert np.allclose(c1.value, (self.v.nominal_value /
                                      self.a.nominal_value).to(c1.unit).value)
        # Fractional uncertainties under division add in quadrature
        assert np.allclose(c1.uncertainty/c1.value,
                           np.sqrt((self.v.uncertainty/self.v.value)**2 +
                                   (self.a.uncertainty/self.a.value)**2))


def test_more_complex():
    G = UQuantity(c.G, subok=False)
    m1 = UQuantity(1e15, 1e5, u.kg)
    m2 = UQuantity(100, 10, u.kg)
    r = UQuantity(10000, 500, u.m)
    F = G * (m1 * m2) / r**2
    assert np.allclose(F.value, c.G.si.value * (1e15 * 100) / 10000**2)
    assert F.unit == u.N
    # Uncertainties calculated using partial derivative method
    assert np.allclose(F.uncertainty, np.sqrt(
        (m1.value*m2.value/(r.value**2)*G.uncertainty)**2 +
        (G.value*m2.value/(r.value**2)*m1.uncertainty)**2 +
        (G.value*m1.value/(r.value**2)*m2.uncertainty)**2 +
        (-2*G.value*m1.value*m2.value/(r.value**3)*r.uncertainty)**2))
