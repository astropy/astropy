# Licensed under a 3-clause BSD style license - see LICENSE.rst


import pytest

from astropy import units as u
from astropy.units import deprecated

emu = deprecated.emu
GearthRad = deprecated.GearthRad
MjupiterMass = deprecated.MjupiterMass
mjupiterRad = deprecated.MjupiterRad
nearthMass = deprecated.nearthMass


def test_enable():
    with deprecated.enable():
        # `unit in u.Bi.compose()` would use `==` for comparison, but we really
        # do want to check identity, not just equality.
        assert any(unit is emu for unit in u.Bi.compose())


def test_emu():
    assert emu == u.Bi


@pytest.mark.parametrize(
    "unit",
    [emu, GearthRad, MjupiterMass, mjupiterRad, nearthMass],
    ids=lambda x: x.name,
)
def test_deprecated_unit_not_in_main_namespace(unit):
    with pytest.raises(AttributeError):
        getattr(u, unit.name)


@pytest.mark.parametrize(
    "prefixed_unit,base_unit",
    [
        pytest.param(prefixed_unit, base_unit, id=prefixed_unit.name)
        for prefixed_unit, base_unit in [
            (GearthRad, u.earthRad),
            (MjupiterMass, u.jupiterMass),
            (mjupiterRad, u.jupiterRad),
            (nearthMass, u.earthMass),
        ]
    ],
)
def test_deprecated_unit_definition(prefixed_unit, base_unit):
    assert prefixed_unit.represents.bases[0] is base_unit
