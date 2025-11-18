import inspect

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.photoncomponent import PhotonComponent


class _DummyPhoton(PhotonComponent):
    Ogamma0 = 1e-4

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_photon():
    return _DummyPhoton()


def _is_positional_only(func, param_name="z"):
    sig = inspect.signature(func)
    p = sig.parameters.get(param_name)
    return p is not None and p.kind == inspect.Parameter.POSITIONAL_ONLY


def test_photon_signature_and_behavior(dummy_photon):
    assert hasattr(PhotonComponent, "Ogamma")
    assert _is_positional_only(PhotonComponent.Ogamma)
    assert_allclose(dummy_photon.Ogamma(1), 1e-4 * (1 + 1) ** 4)
