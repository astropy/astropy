import pytest

import astropy.units as u
from astropy.cosmology._src.traits.tcmb import TemperatureCMB
from astropy.tests.helper import assert_quantity_allclose


class DummyTcmb(TemperatureCMB):
    Tcmb0 = 2.7 * u.K


@pytest.fixture
def dummy_tcmb():
    return DummyTcmb()


def test_tcmb_behavior_and_signature(dummy_tcmb):
    t = dummy_tcmb
    assert hasattr(TemperatureCMB, "Tcmb")
    tout = t.Tcmb(1)
    assert tout.unit == u.K
    assert_quantity_allclose(tout, 2.7 * u.K * (1 + 1))
