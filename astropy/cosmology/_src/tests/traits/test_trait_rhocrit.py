import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology._src.traits.rhocrit import CriticalDensity
from astropy.tests.helper import assert_quantity_allclose


class DummyRho(CriticalDensity):
    critical_density0 = 1.0 * u.kg / (u.m**3)

    def efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_rho():
    return DummyRho()


def test_critical_density_returns_quantity(dummy_rho):
    rho = dummy_rho.critical_density(1)
    assert isinstance(rho, u.Quantity)
    assert_quantity_allclose(rho, 1.0 * u.kg / (u.m**3))
