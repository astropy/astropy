import numpy as np
import pytest

from astropy.cosmology._src.traits.darkenergy import DarkEnergyComponent

from .helper import is_positional_only


class FailingDarkEnergy(DarkEnergyComponent):
    """A class that fails to implement the required inv_efunc method."""

    Ode0 = 0.7

    def w(self, z, /):
        return -1.0


class ConcreteDarkEnergy(DarkEnergyComponent):
    Ode0 = 0.7

    def w(self, z, /):
        return -1.0

    def inv_efunc(self, z, /):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def concrete_de():
    return ConcreteDarkEnergy()


def test_darkenergy_signature_and_missing_inv_efunc_raises():
    assert hasattr(DarkEnergyComponent, "w")
    assert is_positional_only(DarkEnergyComponent.w, "z")
    # inv_efunc is abstract; FailingDarkEnergy (which lacks it) cannot be instantiated
    with pytest.raises(TypeError, match="abstract"):
        FailingDarkEnergy()


def test_darkenergy_with_inv_efunc(concrete_de):
    pytest.importorskip("scipy")
    # For w = -1 (cosmological constant) the density scale is 1, so Ode==Ode0
    val = concrete_de.Ode(1)
    assert pytest.approx(val) == 0.7
