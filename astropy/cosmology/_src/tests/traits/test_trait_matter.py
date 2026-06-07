import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.matter import MatterComponent

from .helper import is_positional_only


class DummyMatter(MatterComponent):
    Om0 = 0.3

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


class ZeroMatter(MatterComponent):
    Om0 = 0.0

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_matter():
    return DummyMatter()


@pytest.fixture
def zero_matter():
    return ZeroMatter()


def test_matter_exists_and_signature():
    assert hasattr(MatterComponent, "Om")
    assert is_positional_only(MatterComponent.Om, "z")


def test_matter_scalar_array_quantity_behavior(dummy_matter):
    d = dummy_matter
    # keyword should raise TypeError because z is positional-only
    with pytest.raises(TypeError):
        d.Om(z=1)

    # scalar
    assert_allclose(d.Om(1), 0.3 * (1 + 1) ** 3)

    # array
    zin = np.array([0.0, 1.0])
    out = d.Om(zin)
    assert out.shape == zin.shape
    assert_allclose(out, 0.3 * (1 + zin) ** 3)


def test_matter_zero_case(zero_matter):
    assert_allclose(zero_matter.Om(1), 0.0)
