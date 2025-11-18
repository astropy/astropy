import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.cosmology._src.traits.neutrino import NeutrinoComponent


class _DummyNeutrino(NeutrinoComponent):
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
    return _DummyNeutrino()


def test_neutrino_onu_and_tnu_basic(dummy_neutrino):
    d = dummy_neutrino
    out = d.Onu(1)
    # scalar-like
    assert np.asarray(out).shape == ()
    # Tnu scales as (1+z)
    assert_allclose(d.Tnu(1).to_value(u.K), d.Tnu0.to_value(u.K) * (1 + 1))
