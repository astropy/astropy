# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Unit tests for the handling of physical types in `astropy.units`.
"""

import pytest

from astropy import units as u
from astropy.units import physical
from astropy.constants import hbar


unit_physical_type_pairs = [
    (u.m, "length"),
    (u.cm ** 3, "volume"),
    (u.km / u.h, "speed"),
    (u.barn * u.Mpc, "volume"),
    (u.m * u.s ** 8, "unknown"),
    (u.m / u.m, "dimensionless"),
    (hbar.unit, "angular momentum"),
    (u.erg / (u.cm ** 2 * u.s * u.AA), "spectral flux density wav"),
    (u.photon / (u.cm ** 2 * u.s * u.AA), "photon flux density wav"),
    (u.photon / (u.cm ** 2 * u.s * u.Hz), "photon flux density"),
    (u.byte, "data quantity"),
    (u.bit, "data quantity"),
    (u.imperial.mi / u.week, "speed"),
    (u.erg / u.s, "power"),
    (u.C / u.s, "electrical current"),
    (u.C / u.s / u.cm ** 2, "electrical current density"),
    (u.T * u.m ** 2, "magnetic flux"),
    (u.N * u.m, "energy"),
    (u.rad / u.ms, "angular speed"),
    (u.Unit(1), "dimensionless"),
    (u.m ** 2, "area"),
    (u.s, "time"),
    (u.rad, "angle"),
    (u.sr, "solid angle"),
    (u.m / u.s ** 2, "acceleration"),
    (u.Hz, "frequency"),
    (u.g, "mass"),
    (u.mol, "amount of substance"),
    (u.K, "temperature"),
    (u.deg_C, "temperature"),
    (u.imperial.deg_F, "temperature"),
    (u.imperial.deg_R, "temperature"),
    (u.imperial.deg_R / u.m, "temperature_gradient"),
    (u.N, "force"),
    (u.J, "energy"),
    (u.Pa, "pressure"),
    (u.W, "power"),
    (u.kg / u.m ** 3, "mass density"),
    (u.m ** 3 / u.kg, "specific volume"),
    (u.mol / u.m ** 3, "molar concentration"),
    (u.kg * u.m / u.s, "momentum/impulse"),
    (u.kg * u.m ** 2 / u.s, "angular momentum"),
    (u.rad / u.s, "angular speed"),
    (u.rad / u.s ** 2, "angular acceleration"),
    (u.g / (u.m * u.s), "dynamic viscosity"),
    (u.m ** 2 / u.s, "kinematic viscosity"),
    (u.m ** -1, "wavenumber"),
    (u.A, "electrical current"),
    (u.C, "electrical charge"),
    (u.V, "electrical potential"),
    (u.Ohm, "electrical resistance"),
    (u.S, "electrical conductance"),
    (u.F, "electrical capacitance"),
    (u.C * u.m, "electrical dipole moment"),
    (u.A / u.m ** 2, "electrical current density"),
    (u.V / u.m, "electrical field strength"),
    (u.C / u.m ** 2, "electrical flux density"),
    (u.C / u.m ** 3, "electrical charge density"),
    (u.F / u.m, "permittivity"),
    (u.Wb, "magnetic flux"),
    (u.T, "magnetic flux density"),
    (u.A / u.m, "magnetic field strength"),
    (u.H / u.m, "electromagnetic field strength"),
    (u.H, "inductance"),
    (u.cd, "luminous intensity"),
    (u.lm, "luminous flux"),
    (u.lx, "luminous emittance/illuminance"),
    (u.W / u.sr, "radiant intensity"),
    (u.cd / u.m ** 2, "luminance"),
    (u.astrophys.Jy, "spectral flux density"),
    (u.astrophys.R, "photon flux"),
    (u.misc.bit, "data quantity"),
    (u.misc.bit / u.s, "bandwidth"),
    (u.cgs.Franklin, "electrical charge (ESU)"),
    (u.cgs.statampere, "electrical current (ESU)"),
    (u.cgs.Biot, "electrical current (EMU)"),
    (u.cgs.abcoulomb, "electrical charge (EMU)"),
    (u.imperial.btu / (u.s * u.m * u.imperial.deg_F), "thermal conductivity"),
    (u.imperial.cal / u.deg_C, "heat capacity"),
    (u.imperial.cal / u.deg_C / u.g, "specific heat capacity"),
    (u.J * u.m ** -2 * u.s ** -1, "energy flux"),
    (u.W / u.m ** 2, "energy flux"),
    (u.m ** 3 / u.mol, "molar volume"),
    (u.m / u.S, "electrical resistivity"),
    (u.S / u.m, "electrical conductivity"),
    (u.A * u.m ** 2, "magnetic moment"),
    (u.J / u.T, "magnetic moment"),
    (u.yr ** -1 * u.Mpc ** -3, "volumetric rate"),
    (u.m / u.s ** 3, "jerk"),
    (u.m / u.s ** 4, "snap"),
    (u.m / u.s ** 5, "crackle"),
    (u.m / u.s ** 6, "pop"),
    (u.deg_C / u.m, "temperature gradient"),
    (u.imperial.deg_F / u.m, "temperature gradient"),
    (u.imperial.deg_R / u.imperial.ft, "temperature gradient"),
    (u.imperial.Calorie / u.g, "specific energy"),
    (u.mol / u.L / u.s, "reaction rate"),
    (u.imperial.lbf * u.imperial.ft * u.s ** 2, "moment of inertia"),
    (u.mol / u.s, "catalytic activity"),
    (u.imperial.kcal / u.deg_C / u.mol, "molar heat capacity"),
    (u.mol / u.kg, "molality"),
    (u.imperial.inch * u.hr, "absement"),
    (u.imperial.ft ** 3 / u.s, "volumetric flow rate"),
    (u.Hz / u.s, "frequency drift"),
    (u.Pa ** -1, "compressibility"),
    (u.count, "number of counts"),
    (u.electron, "number of electrons"),
    (u.pix, "number of elements on 2D regular grid"),
    (u.vox, "number of elements on 3D regular grid"),
]


@pytest.mark.parametrize("unit, physical_type", unit_physical_type_pairs)
def test_physical_types(unit, physical_type):
    """
    Test that the `physical_type` attribute of `u.Unit` objects provides
    the expected physical type for various units.

    Many of these tests are used to test backwards compatibility.
    """
    if unit.physical_type != physical_type:
        pytest.fail(
            f"{repr(unit)}.physical_type was expected to return "
            f"{repr(physical_type)}, but instead returned "
            f"{unit.physical_type}."
        )


@pytest.mark.parametrize(
    "unit, expected_set",
    [
        (u.m, {"length"}),
        (u.Pa, {"energy density", "pressure", "stress"}),
    ],
)
def test_physical_type_as_set(unit, expected_set):
    """Test making a `physical.PhysicalType` instance into a `set`."""
    resulting_set = set(unit.physical_type)
    assert resulting_set == expected_set


def test_physical_type_id():
    """
    Test that the physical type ID of a `physical.PhysicalType`
    instance matches that of a unit of that physical type.
    """
    time_physical_type = u.s.physical_type
    time_id = time_physical_type._physical_type_id
    hour_id = u.hr._get_physical_type_id()
    assert time_id == hour_id


equivalent_unit_pairs = [
    (u.m, u.m),
    (u.m, u.cm),
    (u.m ** 18, u.pc ** 18),
    (u.N, u.kg * u.m * u.s ** -2),
    (u.barn * u.Mpc, u.cm ** 3),
]


@pytest.mark.parametrize("unit1, unit2", equivalent_unit_pairs)
def test_physical_type_id(unit1, unit2):
    """Test that the physical type IDs of equivalent units match."""
    id1 = physical.get_physical_type(unit1)._physical_type_id
    id2 = physical.get_physical_type(unit2)._physical_type_id
    assert id1 == id2


@pytest.mark.parametrize("unit1, unit2", equivalent_unit_pairs)
def test_physical_type_instance_equality(unit1, unit2):
    """
    Test that `physical.PhysicalType` instances for units of the same
    dimensionality are equal.
    """
    physical_type1 = physical.PhysicalType(unit1, "ptype1")
    physical_type2 = physical.PhysicalType(unit2, "ptype2")
    assert (physical_type1 == physical_type2) is True
    assert (physical_type1 != physical_type2) is False


nonequivalent_unit_pairs = [
    (u.m, u.s),
    (u.m ** 18, u.m ** 19),
    (u.N, u.J),
    (u.barn, u.imperial.deg_F),
]


@pytest.mark.parametrize("unit1, unit2", nonequivalent_unit_pairs)
def test_physical_type_instance_inequality(unit1, unit2):
    """
    Test that `physical.PhysicalType` instances for units with different
    dimensionality are considered unequal
    """
    physical_type1 = physical.PhysicalType(unit1, "ptype1")
    physical_type2 = physical.PhysicalType(unit2, "ptype2")
    assert (physical_type1 != physical_type2) is True
    assert (physical_type1 == physical_type2) is False


def test_physical_type_inequality_other_type():
    """
    Test that a physical type does not turn out to be equal to an object
    of another type.
    """
    assert u.m.physical_type != u.m


def test_physical_type_contains():
    """Test PhysicalType.__contains__."""
    length = u.m.physical_type
    assert "length" in length


def test_physical_type_multiplication():
    """Test PhysicalType.__mul__."""
    length = u.m.physical_type
    time = u.s.physical_type
    speed = (u.m / u.s).physical_type
    assert speed * time == length
    assert speed * u.s == length


def test_physical_type_division():
    """Test PhysicalType.__div__."""
    length = u.m.physical_type
    time = u.s.physical_type
    speed = (u.m / u.s).physical_type
    assert length / time == speed
    assert length / u.s == speed
    assert u.m / time == speed


def test_physical_type_power():
    """Test taking `physical.PhysicalType` instances to different powers."""
    length = u.m.physical_type
    area = u.barn.physical_type
    assert length ** 2 == area
    assert area ** (1 / 2) == length


def test_physical_type_power_typeerror():
    """
    Test that exponentiation with an invalid power raises the
    appropriate exception.
    """
    length = u.m.physical_type
    with pytest.raises(TypeError):
        length ** length


def test_dimensionless_physical_type_equality():
    """Test a dimensionless instance of `physical.PhysicalType`."""
    dimensionless = u.dimensionless_unscaled.physical_type
    assert dimensionless == "dimensionless"


def test_dimensionless_operations():
    dimensionless = u.dimensionless_unscaled.physical_type
    assert dimensionless == dimensionless * dimensionless
    assert dimensionless == dimensionless / dimensionless
    assert dimensionless == dimensionless ** 4


def test_str_for_unique_physical_type():
    """
    Test that a `physical.PhysicalType` instance gets converted to a
    string with `str` correctly.
    """
    length = u.m.physical_type
    speed = (u.m / u.s).physical_type
    assert str(length) == "length"
    assert str(speed) == "speed"


def test_repr_for_unique_physical_type():
    """
    Test that units with only one element in the set of physical types
    gets represented correctly as a string with `repr`.  This behavior
    should mimic a string to match the previous API.
    """
    length = u.m.physical_type
    speed = (u.m / u.s).physical_type
    assert repr(length) == "'length'"
    assert repr(speed) == "'speed'"


def test_str_and_repr_for_multiple_physical_types():
    """
    Test that units with multiple physical types get represented
    correctly as strings with both `str` and `repr`.
    """
    pressure = u.Pa.physical_type
    expected = "{'energy density', 'pressure', 'stress'}"
    assert str(pressure) == expected
    assert repr(pressure) == expected


def test_physical_type_attribute_of_unit():
    """
    Test that a `u.Unit` instance has an attribute named `physical_type`
    which is an instance of `physical.PhysicalType`.
    """
    assert isinstance(u.m.physical_type, physical.PhysicalType)


def test_unknown_unit_physical_type_is_PhysicalType():
    """
    Test that a set of units with an unknown physical type has an
    attribute named `physical_type` that is a `physical.PhysicalType`
    instance.
    """
    unknown_unit = u.s ** 19
    assert isinstance(unknown_unit.physical_type, physical.PhysicalType)


@pytest.mark.parametrize("not_unit", [u.m.physical_type, 1, "1"])
def test_get_physical_type_error(not_unit):
    """
    Test that attempting to get a physical type with
    `physical.get_physical_type` for something that is not a unit raises
    the appropriate exception.
    """
    with pytest.raises(TypeError):
        physical.get_physical_type(not_unit)


@pytest.mark.parametrize(
    "temperature_unit1, temperature_unit2",
    [
        (u.K, u.deg_C),
        (u.K, u.imperial.deg_F),
        (u.deg_C, u.imperial.deg_F),
        (u.K, u.imperial.deg_R),
        (u.imperial.deg_F, u.imperial.deg_R),
    ],
)
def test_physical_type_temperature_units(temperature_unit1, temperature_unit2):
    """
    Because K, °C, & °F have different physical type IDs, test that
    different measurements of temperature are treated as equivalent.
    """
    assert temperature_unit1.physical_type == "temperature"
    assert temperature_unit2.physical_type == "temperature"
    assert temperature_unit1.physical_type == temperature_unit2.physical_type


def test_physical_type_underscore_replacement():
    """
    Test that underscores are treated as spaces in a string
    representation of a physical type.
    """
    pressure = u.Pa.physical_type
    assert pressure == "energy_density"


momentum_representations = [
    "momentum",
    "impulse",
    "momentum/impulse",
    "impulse/momentum",
]


@pytest.mark.parametrize("expected", momentum_representations)
def test_momentum_equality_with_delimiters(expected):
    """
    Test that PhysicalType.__eq__ correctly treats different ways of
    representing momentum and impulse.
    """
    momentum = (u.kg * u.m / u.s).physical_type
    assert momentum == expected


def test_physical_type_iteration():
    """Test iterating through different physical type names."""
    pressure = u.Pa.physical_type
    physical_type_names = [physical_type_name for physical_type_name in pressure]
    assert physical_type_names == ["energy density", "pressure", "stress"]


def test_redefining_physical_type_warning():
    """
    Test that trying to redefine a physical type with
    `physical.def_physical_type` issues the appropriate warning.
    """
    weird_unit = u.s ** 42
    physical.def_physical_type(weird_unit, "original physical type")
    with pytest.raises(ValueError):
        physical.def_physical_type(weird_unit, "new physical type")


def test_redundant_physical_type():
    """
    Test that a physical type name already in use cannot be assigned for
    another unit (excluding `"unknown"`).
    """
    with pytest.raises(ValueError):
        physical.def_physical_type(u.m ** 23, "length")
