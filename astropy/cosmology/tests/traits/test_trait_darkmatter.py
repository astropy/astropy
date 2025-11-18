import inspect

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.darkmatter import DarkMatterComponent


class _DummyDarkMatter(DarkMatterComponent):
    Odm0 = 0.25

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_darkmatter():
    return _DummyDarkMatter()


def _is_positional_only(func, param_name="z"):
    sig = inspect.signature(func)
    p = sig.parameters.get(param_name)
    return p is not None and p.kind == inspect.Parameter.POSITIONAL_ONLY


def test_darkmatter_signature_and_behavior(dummy_darkmatter):
    assert hasattr(DarkMatterComponent, "Odm")
    assert _is_positional_only(DarkMatterComponent.Odm)
    assert_allclose(dummy_darkmatter.Odm(1), 0.25 * (1 + 1) ** 3)
