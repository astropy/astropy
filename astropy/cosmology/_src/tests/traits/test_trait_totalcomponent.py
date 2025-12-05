import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.totalcomponent import TotalComponent


class DummyTotal(TotalComponent):
    @property
    def Otot0(self):
        return 1.0

    def Otot(self, z, /):
        z = np.asarray(z)
        return np.ones_like(z)


@pytest.fixture
def dummy_total():
    return DummyTotal()


def test_totalcomponent_minimal_impl(dummy_total):
    assert_allclose(dummy_total.Otot(1), 1.0)
