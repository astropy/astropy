# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from .. import Variable


def test_initialisation():
    v1 = Variable(5, 2)
    assert v1.nominal_value == 5
    assert v1.uncertainty == 2
    v2 = Variable(5, None)
    assert v2.nominal_value == 5
    assert v2._uncertainty is None
    assert v2.uncertainty == 0
    v4 = Variable(np.arange(5.), 2.)
    assert np.all(v4.nominal_value == np.arange(5.))
    assert v4.uncertainty == 2
    v5 = Variable(np.arange(5.), np.array([1., 2., 1., 2., 1.]))
    assert np.all(v5.nominal_value == np.arange(5.))
    assert np.all(v5.uncertainty == np.array([1., 2., 1., 2., 1.]))


class TestBasics():
    def setup(self):
        self.v = Variable(5., 2.)
        self.a = Variable(np.arange(1., 5.), 1.)
        self.b = Variable(np.array([1., 2., 3.]), np.array([0.1, 0.2, 0.1]))

    def test_addition(self):
        c = self.v + Variable(12, 5)
        assert c.nominal_value == self.v.nominal_value + 12
        # Uncertainties under addition add in quadrature
        assert c.uncertainty == np.sqrt(self.v.uncertainty**2 + 5**2)
        # try array
        c3 = self.v + self.b
        assert np.all(c3.nominal_value ==
                      self.v.nominal_value + self.b.nominal_value)
        assert np.allclose(c3.uncertainty, np.sqrt(self.v.uncertainty**2 +
                                                   self.b.uncertainty**2))
        # Try adding a regular number.
        c4 = self.v + 10.
        assert c4.nominal_value == self.v.nominal_value + 10.
        assert c4.uncertainty == self.v.uncertainty
        # And a variable without an uncertainty.
        c5 = self.v + Variable(10., None)
        assert c5.nominal_value == self.v.nominal_value + 10.
        assert len(c5._uncertainty.derivatives) == 1
        assert c5.uncertainty == self.v.uncertainty

    def test_subtraction(self):
        c = self.v - Variable(12, 5)
        assert c.nominal_value == self.v.nominal_value - 12
        # Uncertainties under addition add in quadrature
        assert c.uncertainty == np.sqrt(self.v.uncertainty**2 + 5**2)

    def test_multiplication(self):
        c = self.v * self.a
        assert np.all(c.nominal_value ==
                      self.v.nominal_value * self.a.nominal_value)

        # Fractional uncertainties under multiplication add in quadrature
        assert np.allclose(c.uncertainty / c.nominal_value, np.sqrt(
            (self.v.uncertainty / self.v.nominal_value)**2 +
            (self.a.uncertainty / self.a.nominal_value)**2))

        # Test multiplication with straight number
        c2 = self.a * 10.
        assert np.all(c2.nominal_value == self.a.nominal_value * 10.)
        assert np.all(c2.uncertainty == self.a.uncertainty * 10.)

    def test_division(self):
        c = self.v / self.a
        assert np.allclose(c.nominal_value, (self.v.nominal_value /
                                             self.a.nominal_value))
        # Fractional uncertainties under division add in quadrature
        assert np.allclose(c.uncertainty/c.nominal_value, np.sqrt(
            (self.v.uncertainty / self.v.nominal_value)**2 +
            (self.a.uncertainty / self.a.nominal_value)**2))

    def test_tracking(self):
        c = Variable(12, 5) + self.v
        assert c.uncertainty == np.sqrt(5**2 + self.v.uncertainty**2)
        c = c - self.v
        assert c.nominal_value == 12
        assert c.uncertainty == 5
        c2 = self.v - self.v
        assert c2.nominal_value == 0
        assert c2.uncertainty == 0
