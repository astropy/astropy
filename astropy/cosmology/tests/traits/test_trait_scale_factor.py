import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.scale_factor import ScaleFactor

from .helper import is_positional_only


class DummyScale(ScaleFactor):
    pass


@pytest.fixture
def dummy_scale():
    return DummyScale()


def test_scale_factor_behavior_and_signature(dummy_scale):
    s = dummy_scale
    assert hasattr(ScaleFactor, "scale_factor")
    # basic value
    assert_allclose(s.scale_factor(1), 1.0 / (1 + 1))

    # scale_factor0 default
    assert s.scale_factor0 == 1 << s.scale_factor0.unit

    # positional-only API
    assert is_positional_only(ScaleFactor.scale_factor, "z")

    # passing as keyword should raise TypeError (positional-only)
    with pytest.raises(TypeError):
        s.scale_factor(z=1)

    # array input returns same-shaped array
    import numpy as np

    arr = np.array([0.0, 1.0, 2.0])
    out = s.scale_factor(arr)
    assert out.shape == arr.shape
