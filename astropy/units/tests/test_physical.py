# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Unit tests for the handling of physical types in `astropy.units`.
"""

import pickle

import pytest

from astropy import units as u
from astropy.constants import hbar
from astropy.units import physical

unit_physical_type_pairs = [
    (u.m, "length"),
    (u.cm**3, "volume"),
    (u.km / u.h, "speed"),
    (u.barn * u.Mpc, "volume"),
    (u.m * u.s**8, "unknown"),
    (u.m / u.m, "dimensionless"),
    (hbar.unit, "angular momentum"),
    (u.erg / (u.cm**2 * u.s * u.AA), "spectral flux density wav"),
    (u.photon / (u.cm**2 * u.s * u.AA), "photon flux density wav"),
    (u.photon / (u.cm**2 * u.s * u.Hz), "photon flux density"),
    (u.Jy / u.sr, "surface brightness"),
    (u.J * u.m**-3 * u.s**-1 * u.sr**-1, "surface brightness wav"),
    (u.photon / u.Hz / u.cm**2 / u.s / u.sr, "photon surface brightness"),
    (u.photon / u.AA / u.cm**2 / u.s / u.sr, "photon surface brightness wav"),
    (u.byte, "data quantity"),
    (u.bit, "data quantity"),
    (u.imperial.mi / u.week, "speed"),
    (u.erg / u.s, "power"),
    (u.C / u.s, "electrical current"),
    (u.C / u.s / u.cm**2, "electrical current density"),
    (u.T * u.m**2, "magnetic flux"),
    (u.N * u.m, "energy"),
    (u.rad / u.ms, "angular speed"),
    (u.Unit(1), "dimensionless"),
    (u.m**2, "area"),
    (u.s, "time"),
    (u.rad, "angle"),
    (u.sr, "solid angle"),
    (u.m / u.s**2, "acceleration"),
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
    (u.kg / u.m**3, "mass density"),
    (u.m**3 / u.kg, "specific volume"),
    (u.mol / u.m**3, "molar concentration"),
    (u.kg * u.m / u.s, "momentum/impulse"),
    (u.kg * u.m**2 / u.s, "angular momentum"),
    (u.rad / u.s, "angular speed"),
    (u.rad / u.s**2, "angular acceleration"),
    (u.g / (u.m * u.s), "dynamic viscosity"),
    (u.m**2 / u.s, "kinematic viscosity"),
    (u.m**-1, "wavenumber"),
    (u.A, "electrical current"),
    (u.C, "electrical charge"),
    (u.V, "electrical potential"),
    (u.Ohm, "electrical resistance"),
    (u.S, "electrical conductance"),
    (u.F, "electrical capacitance"),
    (u.C * u.m, "electrical dipole moment"),
    (u.A / u.m**2, "electrical current density"),
    (u.V / u.m, "electrical field strength"),
    (u.C / u.m**2, "electrical flux density"),
    (u.C / u.m**3, "electrical charge density"),
    (u.F / u.m, "permittivity"),
    (u.Wb, "magnetic flux"),
    (u.Wb**2, "magnetic helicity"),
    (u.T, "magnetic flux density"),
    (u.A / u.m, "magnetic field strength"),
    (u.H / u.m, "electromagnetic field strength"),
    (u.H, "inductance"),
    (u.cd, "luminous intensity"),
    (u.lm, "luminous flux"),
    (u.lx, "luminous emittance/illuminance"),
    (u.W / u.sr, "radiant intensity"),
    (u.cd / u.m**2, "luminance"),
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
    (u.J * u.m**-2 * u.s**-1, "energy flux"),
    (u.W / u.m**2, "energy flux"),
    (u.m**3 / u.mol, "molar volume"),
    (u.m / u.S, "electrical resistivity"),
    (u.S / u.m, "electrical conductivity"),
    (u.A * u.m**2, "magnetic moment"),
    (u.J / u.T, "magnetic moment"),
    (u.yr**-1 * u.Mpc**-3, "volumetric rate"),
    (u.m / u.s**3, "jerk"),
    (u.m / u.s**4, "snap"),
    (u.m / u.s**5, "crackle"),
    (u.m / u.s**6, "pop"),
    (u.deg_C / u.m, "temperature gradient"),
    (u.imperial.deg_F / u.m, "temperature gradient"),
    (u.imperial.deg_R / u.imperial.ft, "temperature gradient"),
    (u.imperial.Calorie / u.g, "specific energy"),
    (u.mol / u.L / u.s, "reaction rate"),
    (u.imperial.lbf * u.imperial.ft * u.s**2, "moment of inertia"),
    (u.mol / u.s, "catalytic activity"),
    (u.imperial.kcal / u.deg_C / u.mol, "molar heat capacity"),
    (u.mol / u.kg, "molality"),
    (u.imperial.inch * u.hr, "absement"),
    (u.imperial.ft**3 / u.s, "volumetric flow rate"),
    (u.Hz / u.s, "frequency drift"),
    (u.Pa**-1, "compressibility"),
    (u.dimensionless_unscaled, "dimensionless"),
]


@pytest.mark.parametrize("unit, physical_type", unit_physical_type_pairs)
def test_physical_type_names(unit, physical_type):
    """
    Test that the `physical_type` attribute of `u.Unit` objects provides
    the expected physical type for various units.

    Many of these tests are used to test backwards compatibility.
    """
    assert unit.physical_type == physical_type, (
        f"{unit!r}.physical_type was expected to return "
        f"{physical_type!r}, but instead returned {unit.physical_type!r}."
    )


length = u.m.physical_type
time = u.s.physical_type
speed = (u.m / u.s).physical_type
area = (u.m**2).physical_type
wavenumber = (u.m**-1).physical_type
dimensionless = u.dimensionless_unscaled.physical_type
pressure = u.Pa.physical_type
momentum = (u.kg * u.m / u.s).physical_type


@pytest.mark.parametrize(
    "physical_type_representation, physical_type_name",
    [
        (1.0, "dimensionless"),
        (u.m, "length"),
        ("work", "work"),
        (5 * u.m, "length"),
        (length, length),
        (u.Pa, "energy_density"),  # attribute-accessible name
        ("energy_density", "energy_density"),  # attribute-accessible name
    ],
)
def test_getting_physical_type(physical_type_representation, physical_type_name):
    """Test different ways of getting a physical type."""
    physical_type = physical.get_physical_type(physical_type_representation)
    assert isinstance(physical_type, physical.PhysicalType)
    assert physical_type == physical_type_name


@pytest.mark.parametrize(
    "argument, exception",
    [
        ("unknown", ValueError),
        ("not a name of a physical type", ValueError),
        ({"this set cannot be made into a Quantity"}, TypeError),
    ],
)
def test_getting_physical_type_exceptions(argument, exception):
    """
    Test that `get_physical_type` raises appropriate exceptions when
    provided with invalid arguments.
    """
    with pytest.raises(exception):
        physical.get_physical_type(argument)


def test_physical_type_cannot_become_quantity():
    """
    Test that `PhysicalType` instances cannot be cast into `Quantity`
    objects.  A failure in this test could be related to failures
    in subsequent tests.
    """
    with pytest.raises(TypeError):
        u.Quantity(u.m.physical_type, u.m)


# left term, right term, operator, expected value
operation_parameters = [
    (length, length, "__eq__", True),
    (length, area, "__eq__", False),
    (length, "length", "__eq__", True),
    ("length", length, "__eq__", NotImplemented),
    (dimensionless, dimensionless, "__eq__", True),
    (momentum, "momentum/impulse", "__eq__", True),  # test delimiters in names
    (pressure, "energy_density", "__eq__", True),  # test underscores in names
    ((u.m**8).physical_type, "unknown", "__eq__", True),
    ((u.m**8).physical_type, (u.m**9).physical_type, "__eq__", False),
    (length, length, "__ne__", False),
    (speed, time, "__ne__", True),
    (pressure, dimensionless, "__ne__", True),
    (length, u.m, "__eq__", NotImplemented),
    (length, length, "__mul__", area),
    (speed, time, "__mul__", length),
    (speed, time, "__rmul__", length),
    (length, time, "__truediv__", speed),
    (area, length, "__truediv__", length),
    (length, area, "__rtruediv__", length),
    (dimensionless, dimensionless, "__mul__", dimensionless),
    (dimensionless, dimensionless, "__truediv__", dimensionless),
    (length, 2, "__pow__", area),
    (area, 0.5, "__pow__", length),
    (dimensionless, 4, "__pow__", dimensionless),
    (u.m, length, "__mul__", NotImplemented),
    (3.2, length, "__mul__", NotImplemented),
    (u.m, time, "__truediv__", NotImplemented),
    (3.2, length, "__truediv__", NotImplemented),
    (length, u.m, "__mul__", area),
    (length, u.m, "__rmul__", area),
    (speed, u.s, "__mul__", length),
    (length, 1, "__mul__", length),
    (length, 1, "__rmul__", length),
    (length, u.s, "__truediv__", speed),
    (area, 1, "__truediv__", area),
    (time, u.m, "__rtruediv__", speed),
    (length, 1.0, "__rtruediv__", wavenumber),
    (length, 2, "__pow__", area),
    (length, 32, "__mul__", NotImplemented),
    (length, 0, "__rmul__", NotImplemented),
    (length, 3.2, "__truediv__", NotImplemented),
    (length, -1, "__rtruediv__", NotImplemented),
    (length, "length", "__mul__", area),
    (length, "length", "__rmul__", area),
    (area, "length", "__truediv__", length),
    (length, "area", "__rtruediv__", length),
]


@pytest.mark.parametrize("left, right, operator, expected", operation_parameters)
def test_physical_type_operations(left, right, operator, expected):
    """
    Test that `PhysicalType` dunder methods that require another
    argument behave as intended.
    """
    assert getattr(left, operator)(right) == expected


unit_with_physical_type_set = [
    (u.m, {"length"}),
    (u.kg * u.m / u.s, {"impulse", "momentum"}),
    (u.Pa, {"energy density", "pressure", "stress"}),
]


@pytest.mark.parametrize("unit, expected_set", unit_with_physical_type_set)
def test_physical_type_as_set(unit, expected_set):
    """Test making a `physical.PhysicalType` instance into a `set`."""
    resulting_set = set(unit.physical_type)
    assert resulting_set == expected_set


def test_physical_type_iteration():
    """Test iterating through different physical type names."""
    physical_type_names = list(pressure)
    assert physical_type_names == ["energy density", "pressure", "stress"]


def test_physical_type_in():
    """
    Test that `in` works as expected for `PhysicalType` objects with one
    or multiple names.
    """
    assert "length" in length
    assert "pressure" in pressure


equivalent_unit_pairs = [
    (u.m, u.m),
    (u.m, u.cm),
    (u.N, u.kg * u.m * u.s**-2),
    (u.barn * u.Mpc, u.cm**3),
    (u.K, u.deg_C),
    (u.K, u.imperial.deg_R),
    (u.K, u.imperial.deg_F),
    (u.deg_C, u.imperial.deg_F),
    (u.m**18, u.pc**18),
]


@pytest.mark.parametrize("unit1, unit2", equivalent_unit_pairs)
def test_physical_type_instance_equality(unit1, unit2):
    """
    Test that `physical.PhysicalType` instances for units of the same
    dimensionality are equal.
    """
    assert (unit1.physical_type == unit2.physical_type) is True
    assert (unit1.physical_type != unit2.physical_type) is False


@pytest.mark.parametrize("unit1, unit2", equivalent_unit_pairs)
def test_get_physical_type_equivalent_pairs(unit1, unit2):
    """
    Test that `get_physical_type` retrieves the same `PhysicalType`
    instances for equivalent physical types, except for unknown types
    which are not cataloged.
    """
    physical_type1 = physical.get_physical_type(unit1)
    physical_type2 = physical.get_physical_type(unit2)
    assert physical_type1 == physical_type2
    if physical_type1 != "unknown":
        assert physical_type1 is physical_type2


nonequivalent_unit_pairs = [
    (u.m, u.s),
    (u.m**18, u.m**19),
    (u.N, u.J),
    (u.barn, u.imperial.deg_F),
]


@pytest.mark.parametrize("unit1, unit2", nonequivalent_unit_pairs)
def test_physical_type_instance_inequality(unit1, unit2):
    """
    Test that `physical.PhysicalType` instances for units with different
    dimensionality are considered unequal.
    """
    physical_type1 = physical.PhysicalType(unit1, "ptype1")
    physical_type2 = physical.PhysicalType(unit2, "ptype2")
    assert (physical_type1 != physical_type2) is True
    assert (physical_type1 == physical_type2) is False


physical_type_with_expected_str = [
    (length, "length"),
    (speed, "speed/velocity"),
    (pressure, "energy density/pressure/stress"),
    (u.deg_C.physical_type, "temperature"),
    ((u.J / u.K / u.kg).physical_type, "specific entropy/specific heat capacity"),
]

physical_type_with_expected_repr = [
    (length, "PhysicalType('length')"),
    (speed, "PhysicalType({'speed', 'velocity'})"),
    (pressure, "PhysicalType({'energy density', 'pressure', 'stress'})"),
    (u.deg_C.physical_type, "PhysicalType('temperature')"),
    (
        (u.J / u.K / u.kg).physical_type,
        "PhysicalType({'specific entropy', 'specific heat capacity'})",
    ),
]


@pytest.mark.parametrize("physical_type, expected_str", physical_type_with_expected_str)
def test_physical_type_str(physical_type, expected_str):
    """Test using `str` on a `PhysicalType` instance."""
    assert str(physical_type) == expected_str


@pytest.mark.parametrize(
    "physical_type, expected_repr", physical_type_with_expected_repr
)
def physical_type_repr(physical_type, expected_repr):
    """Test using `repr` on a `PhysicalType` instance."""
    assert repr(physical_type) == expected_repr


def test_physical_type_hash():
    """Test that a `PhysicalType` instance can be used as a dict key."""
    dictionary = {length: 42}
    assert dictionary[length] == 42


@pytest.mark.parametrize("multiplicand", [[], 42, 0, -1])
def test_physical_type_multiplication(multiplicand):
    """
    Test that multiplication of a physical type returns `NotImplemented`
    when attempted for an invalid type.
    """
    with pytest.raises(TypeError):
        length * multiplicand


def test_unrecognized_unit_physical_type():
    """
    Test basic functionality for the physical type of an unrecognized
    unit.
    """
    unrecognized_unit = u.Unit("parrot", parse_strict="silent")
    physical_type = unrecognized_unit.physical_type
    assert isinstance(physical_type, physical.PhysicalType)
    assert physical_type == "unknown"


invalid_inputs = [(42,), ("valid input", 42)]


@pytest.mark.parametrize("invalid_input", invalid_inputs)
def test_invalid_physical_types(invalid_input):
    """
    Test that `PhysicalType` cannot be instantiated when one of the
    supplied names is not a string, while making sure that the physical
    type for the unit remains unknown.
    """
    obscure_unit = u.s**87
    with pytest.raises(ValueError):
        physical.PhysicalType(obscure_unit, invalid_input)
    assert obscure_unit.physical_type == "unknown"


class TestDefPhysType:
    weird_unit = u.m**99
    strange_unit = u.s**42

    def test_attempt_to_define_unknown_physical_type(self):
        """Test that a unit cannot be defined as unknown."""
        with pytest.raises(ValueError):
            physical.def_physical_type(self.weird_unit, "unknown")
        assert "unknown" not in physical._unit_physical_mapping

    def test_multiple_same_physical_type_names(self):
        """
        Test that `def_physical_type` raises an exception when it tries to
        set the physical type of a new unit as the name of an existing
        physical type.
        """
        with pytest.raises(ValueError):
            physical.def_physical_type(self.weird_unit, {"time", "something"})
        assert self.weird_unit.physical_type == "unknown"

    def test_expanding_names_for_physical_type(self):
        """
        Test that calling `def_physical_type` on an existing physical
        type adds a new physical type name.
        """
        weird_name = "weird name"
        strange_name = "strange name"

        try:
            physical.def_physical_type(self.weird_unit, weird_name)
            assert self.weird_unit.physical_type == weird_name, (
                f"unable to set physical type for {self.weird_unit}"
            )
        finally:  # cleanup added name
            physical._attrname_physical_mapping.pop(weird_name.replace(" ", "_"), None)
            physical._name_physical_mapping.pop(weird_name, None)

        # add both strange_name and weird_name
        try:
            physical.def_physical_type(self.weird_unit, strange_name)
            assert set((self.weird_unit).physical_type) == {
                weird_name,
                strange_name,
            }, "did not correctly append a new physical type name."
        finally:  # cleanup added names
            physical._attrname_physical_mapping.pop(
                strange_name.replace(" ", "_"), None
            )
            physical._name_physical_mapping.pop(strange_name, None)
            physical._attrname_physical_mapping.pop(weird_name.replace(" ", "_"), None)
            physical._name_physical_mapping.pop(weird_name, None)

    def test_redundant_physical_type(self):
        """
        Test that a physical type name already in use cannot be assigned
        for another unit (excluding `"unknown"`).
        """
        with pytest.raises(ValueError):
            physical.def_physical_type(self.weird_unit, "length")

    @staticmethod
    def _undef_physical_type(unit):
        """Reset the physical type of unit to "unknown"."""
        for name in list(unit.physical_type):
            del physical._unit_physical_mapping[name]
        del physical._physical_unit_mapping[unit._physical_type_id]
        assert unit.physical_type == "unknown"

    def teardown_method(self):
        """
        Remove the definitions of the physical types that were added
        using `def_physical_unit` for testing purposes.
        """
        for unit in [self.weird_unit, self.strange_unit]:
            physical_type = physical.get_physical_type(unit)
            if physical_type != "unknown":
                self._undef_physical_type(unit)
            assert unit.physical_type == "unknown", (
                f"the physical type for {unit}, which was added for"
                "testing, was not deleted."
            )


@pytest.mark.parametrize("ptype_name", ["length", "speed", "entropy"])
def test_pickling(ptype_name):
    # Regression test for #11685
    ptype = u.get_physical_type(ptype_name)
    pkl = pickle.dumps(ptype)
    other = pickle.loads(pkl)
    assert other == ptype


def test_physical_types_module_access():
    # all physical type names in dir
    assert set(dir(physical)).issuperset(physical._attrname_physical_mapping.keys())
    assert set(dir(physical)).issuperset(physical.__all__)

    # all physical type can be accessed by name
    for pname in physical._attrname_physical_mapping.keys():
        ptype = physical._attrname_physical_mapping[pname]
        assert hasattr(physical, pname)  # make sure works in lazy load
        assert getattr(physical, pname) is ptype

    # a failed access
    with pytest.raises(AttributeError, match="has no attribute"):
        physical.not_a_valid_physical_type_name
