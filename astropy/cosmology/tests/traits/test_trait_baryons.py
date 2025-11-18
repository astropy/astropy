import inspect

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.cosmology._src.traits.baryons import BaryonComponent


class _DummyBaryon(BaryonComponent):
    Ob0 = 0.05

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def dummy_baryon():
    return _DummyBaryon()


def _is_positional_only(func, param_name="z"):
    sig = inspect.signature(func)
    p = sig.parameters.get(param_name)
    return p is not None and p.kind == inspect.Parameter.POSITIONAL_ONLY


def test_baryon_signature_and_behavior(dummy_baryon):
    assert hasattr(BaryonComponent, "Ob")
    assert _is_positional_only(BaryonComponent.Ob)
    assert_allclose(dummy_baryon.Ob(1), 0.05 * (1 + 1) ** 3)
