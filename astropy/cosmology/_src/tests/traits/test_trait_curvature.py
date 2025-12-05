import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.curvature import CurvatureComponent

from .helper import is_positional_only


class DummyCurvature(CurvatureComponent):
    def __init__(self, ok0):
        self._ok0 = ok0

    @property
    def Ok0(self):
        return self._ok0

    @property
    def is_flat(self):
        return self._ok0 == 0

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_curvature_zero():
    return DummyCurvature(0.0)


@pytest.fixture
def dummy_curvature_nonzero():
    return DummyCurvature(-0.02)


def test_curvature_signature_and_behavior(
    dummy_curvature_zero, dummy_curvature_nonzero
):
    assert hasattr(CurvatureComponent, "Ok")
    assert is_positional_only(CurvatureComponent.Ok, "z")
    assert_allclose(dummy_curvature_zero.Ok(1), 0.0)
    assert_allclose(dummy_curvature_nonzero.Ok(1), -0.02 * (1 + 1) ** 2)
