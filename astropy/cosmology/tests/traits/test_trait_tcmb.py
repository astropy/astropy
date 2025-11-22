import pytest
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.cosmology._src.traits.tcmb import TemperatureCMB


class _DummyTcmb(TemperatureCMB):
    Tcmb0 = 2.7 * u.K


@pytest.fixture
def dummy_tcmb():
    return _DummyTcmb()


def test_tcmb_behavior_and_signature(dummy_tcmb):
    t = dummy_tcmb
    assert hasattr(TemperatureCMB, "Tcmb")
    tout = t.Tcmb(1)
    assert tout.unit == u.K
    assert_allclose(tout.to_value(u.K), 2.7 * (1 + 1))
