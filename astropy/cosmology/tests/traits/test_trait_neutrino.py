import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology._src.traits.neutrino import NeutrinoComponent
from astropy.tests.helper import assert_quantity_allclose


class DummyNeutrino(NeutrinoComponent):
    Tcmb0 = 2.7255 * u.K
    Ogamma0 = 5e-5

    @property
    def has_massive_nu(self):
        return False

    @property
    def Onu0(self):
        return 0.22710731766 * 3.046 * self.Ogamma0

    def nu_relative_density(self, z):
        return 0.22710731766 * 3.046 * np.ones_like(np.asarray(z))

    def Ogamma(self, z):
        return self.Ogamma0 * (np.asarray(z) + 1.0) ** 4


@pytest.fixture
def dummy_neutrino():
    return DummyNeutrino()


def test_neutrino_onu_and_tnu_basic(dummy_neutrino):
    d = dummy_neutrino
    out = d.Onu(1)
    # scalar-like
    assert np.asarray(out).shape == ()
    # Tnu scales as (1+z)
    assert_quantity_allclose(d.Tnu(1), d.Tnu0 * (1 + 1))
