import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.scale_factor import ScaleFactor


class _DummyScale(ScaleFactor):
    pass


@pytest.fixture
def dummy_scale():
    return _DummyScale()


def test_scale_factor_behavior_and_signature(dummy_scale):
    s = dummy_scale
    assert hasattr(ScaleFactor, "scale_factor")
    assert_allclose(s.scale_factor(1), 1.0 / (1 + 1))
