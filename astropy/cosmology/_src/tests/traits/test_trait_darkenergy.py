import numpy as np
import pytest

from astropy.cosmology._src.traits.darkenergy import DarkEnergyComponent

from .helper import is_positional_only


class MinimalDarkEnergy(DarkEnergyComponent):
    Ode0 = 0.7

    def w(self, z, /):
        return -1.0


class ConcreteDarkEnergy(DarkEnergyComponent):
    Ode0 = 0.7

    def w(self, z, /):
        return -1.0

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def minimal_de():
    return MinimalDarkEnergy()


@pytest.fixture
def concrete_de():
    return ConcreteDarkEnergy()


def test_darkenergy_signature_and_missing_inv_efunc_raises(minimal_de):
    assert hasattr(DarkEnergyComponent, "w")
    assert is_positional_only(DarkEnergyComponent.w, "z")
    # Ode requires inv_efunc; calling should raise NotImplementedError
    with pytest.raises(NotImplementedError):
        minimal_de.Ode(1)


def test_darkenergy_with_inv_efunc(concrete_de):
    pytest.importorskip("scipy")
    # For w = -1 (cosmological constant) the density scale is 1, so Ode==Ode0
    val = concrete_de.Ode(1)
    assert pytest.approx(val) == 0.7
