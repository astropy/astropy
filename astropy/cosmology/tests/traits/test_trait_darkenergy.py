import inspect

import numpy as np
import pytest

from astropy.cosmology._src.traits.darkenergy import DarkEnergyComponent


def _is_positional_only(func, param_name="z"):
    sig = inspect.signature(func)
    p = sig.parameters.get(param_name)
    return p is not None and p.kind == inspect.Parameter.POSITIONAL_ONLY


class _MinimalDarkEnergy(DarkEnergyComponent):
    Ode0 = 0.7

    def w(self, z, /):
        return -1.0


class _ConcreteDarkEnergy(DarkEnergyComponent):
    Ode0 = 0.7

    def w(self, z, /):
        return -1.0

    def inv_efunc(self, z):
        return np.ones_like(np.asarray(z))


@pytest.fixture
def minimal_de():
    return _MinimalDarkEnergy()


@pytest.fixture
def concrete_de():
    return _ConcreteDarkEnergy()


def test_darkenergy_signature_and_missing_inv_efunc_raises(minimal_de):
    assert hasattr(DarkEnergyComponent, "w")
    assert _is_positional_only(DarkEnergyComponent.w)
    # Ode requires inv_efunc; calling should raise NotImplementedError
    with pytest.raises(NotImplementedError):
        minimal_de.Ode(1)


def test_darkenergy_with_inv_efunc(concrete_de):
    # For w = -1 (cosmological constant) the density scale is 1, so Ode==Ode0
    val = concrete_de.Ode(1)
    assert pytest.approx(val) == 0.7
