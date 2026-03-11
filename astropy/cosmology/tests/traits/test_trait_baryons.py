import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.baryons import BaryonComponent

from .helper import is_positional_only


class DummyBaryon(BaryonComponent):
    Ob0 = 0.05

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_baryon():
    return DummyBaryon()


def test_baryon_signature_and_behavior(dummy_baryon):
    assert hasattr(BaryonComponent, "Ob")
    assert is_positional_only(BaryonComponent.Ob, "z")
    assert_allclose(dummy_baryon.Ob(1), 0.05 * (1 + 1) ** 3)
