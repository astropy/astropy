import numpy as np
import pytest

from astropy.cosmology._src.traits.photoncomponent import PhotonComponent
from astropy.tests.helper import assert_quantity_allclose

from .helper import is_positional_only


class DummyPhoton(PhotonComponent):
    Ogamma0 = 1e-4

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_photon():
    return DummyPhoton()


def test_photon_signature_and_behavior(dummy_photon):
    assert hasattr(PhotonComponent, "Ogamma")
    assert is_positional_only(PhotonComponent.Ogamma, "z")
    assert_quantity_allclose(dummy_photon.Ogamma(1), 1e-4 * (1 + 1) ** 4)
