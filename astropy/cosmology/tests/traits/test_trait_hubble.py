import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology._src.traits.hubble import HubbleParameter
from astropy.tests.helper import assert_quantity_allclose


class DummyHubble(HubbleParameter):
    H0 = 70 * u.km / (u.s * u.Mpc)

    def efunc(self, z):
        return np.ones_like(np.asarray(z))

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_hubble():
    return DummyHubble()


def test_hubble_H_and_properties(dummy_hubble):
    h = dummy_hubble
    H1 = h.H(1)
    assert isinstance(H1, u.Quantity)
    assert_quantity_allclose(H1, 70 * u.km / (u.s * u.Mpc))
    # h property
    assert hasattr(h, "h")
    assert isinstance(h.h, (float, np.floating))
