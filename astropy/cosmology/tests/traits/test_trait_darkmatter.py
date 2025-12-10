import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.darkmatter import DarkMatterComponent

from .helper import is_positional_only


class DummyDarkMatter(DarkMatterComponent):
    Odm0 = 0.25

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_darkmatter():
    return DummyDarkMatter()


def test_darkmatter_signature_and_behavior(dummy_darkmatter):
    assert hasattr(DarkMatterComponent, "Odm")
    assert is_positional_only(DarkMatterComponent.Odm, "z")
    assert_allclose(dummy_darkmatter.Odm(1), 0.25 * (1 + 1) ** 3)
