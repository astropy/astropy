# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Defines the physical types that correspond to different units."""

import numbers

from . import core
from . import si
from . import astrophys
from . import cgs
from . import misc

__all__ = ["def_physical_type", "get_physical_type", "PhysicalType"]

_units_and_physical_types = [
    (core.dimensionless_unscaled, "dimensionless"),
    (si.m, "length"),
    (si.m ** 2, "area"),
    (si.m ** 3, "volume"),
    (si.s, "time"),
    (si.rad, "angle"),
    (si.sr, "solid angle"),
    (si.m / si.s, "speed"),
    (si.m / si.s ** 2, "acceleration"),
    (si.Hz, "frequency"),
    (si.g, "mass"),
    (si.mol, "amount of substance"),
    (si.K, "temperature"),
    (si.W * si.m ** -1 * si.K ** -1, "thermal conductivity"),
    (si.J * si.K ** -1, {"heat capacity", "entropy"}),
    (si.J * si.K ** -1 * si.kg ** -1, {"specific heat capacity", "specific entropy"}),
    (si.N, "force"),
    (si.J, {"energy", "work", "torque"}),
    (si.J * si.m ** -2 * si.s ** -1, {"energy flux", "irradiance"}),
    (si.Pa, {"pressure", "energy density", "stress"}),
    (si.W, {"power", "radiant flux"}),
    (si.kg * si.m ** -3, "mass density"),
    (si.m ** 3 / si.kg, "specific volume"),
    (si.mol / si.m ** 3, "molar concentration"),
    (si.m ** 3 / si.mol, "molar volume"),
    (si.kg * si.m / si.s, {"momentum", "impulse"}),
    (si.kg * si.m ** 2 / si.s, {"angular momentum", "action"}),
    (si.rad / si.s, "angular speed"),
    (si.rad / si.s ** 2, "angular acceleration"),
    (si.rad / si.m, "plate scale"),
    (si.g / (si.m * si.s), "dynamic viscosity"),
    (si.m ** 2 / si.s, {"diffusivity", "kinematic viscosity"}),
    (si.m ** -1, "wavenumber"),
    (si.m ** -2, "column density"),
    (si.A, "electrical current"),
    (si.C, "electrical charge"),
    (si.V, "electrical potential"),
    (si.Ohm, {"electrical resistance", "electrical impedance", "electrical reactance"}),
    (si.Ohm * si.m, "electrical resistivity"),
    (si.S, "electrical conductance"),
    (si.S / si.m, "electrical conductivity"),
    (si.F, "electrical capacitance"),
    (si.C * si.m, "electrical dipole moment"),
    (si.A / si.m ** 2, "electrical current density"),
    (si.V / si.m, "electrical field strength"),
    (si.C / si.m ** 2,
        {"electrical flux density", "surface charge density", "polarization density"},
    ),
    (si.C / si.m ** 3, "electrical charge density"),
    (si.F / si.m, "permittivity"),
    (si.Wb, "magnetic flux"),
    (si.T, "magnetic flux density"),
    (si.A / si.m, "magnetic field strength"),
    (si.m ** 2 * si.A, "magnetic moment"),
    (si.H / si.m, {"electromagnetic field strength", "permeability"}),
    (si.H, "inductance"),
    (si.cd, "luminous intensity"),
    (si.lm, "luminous flux"),
    (si.lx, {"luminous emittance", "illuminance"}),
    (si.W / si.sr, "radiant intensity"),
    (si.cd / si.m ** 2, "luminance"),
    (si.m ** -3 * si.s ** -1, "volumetric rate"),
    (astrophys.Jy, "spectral flux density"),
    (si.W * si.m ** 2 * si.Hz ** -1, "surface tension"),
    (si.J * si.m ** -3 * si.s ** -1, {"spectral flux density wav", "power density"}),
    (astrophys.photon / si.Hz / si.cm ** 2 / si.s, "photon flux density"),
    (astrophys.photon / si.AA / si.cm ** 2 / si.s, "photon flux density wav"),
    (astrophys.R, "photon flux"),
    (misc.bit, "data quantity"),
    (misc.bit / si.s, "bandwidth"),
    (cgs.Franklin, "electrical charge (ESU)"),
    (cgs.statampere, "electrical current (ESU)"),
    (cgs.Biot, "electrical current (EMU)"),
    (cgs.abcoulomb, "electrical charge (EMU)"),
    (si.m * si.s ** -3, {"jerk", "jolt"}),
    (si.m * si.s ** -4, {"snap", "jounce"}),
    (si.m * si.s ** -5, "crackle"),
    (si.m * si.s ** -6, {"pop", "pounce"}),
    (si.K / si.m, "temperature gradient"),
    (si.J / si.kg, "specific energy"),
    (si.mol * si.m ** -3 * si.s ** -1, "reaction rate"),
    (si.kg * si.m ** 2, "moment of inertia"),
    (si.mol / si.s, "catalytic activity"),
    (si.J * si.K ** -1 * si.mol ** -1, "molar heat capacity"),
    (si.mol / si.kg, "molality"),
    (si.m * si.s, {"absement", "sustained displacement"}),
    (si.m * si.s ** 2, "absity"),
    (si.m ** 3 / si.s, "volumetric flow rate"),
    (si.s ** -2, "frequency drift"),
    (si.Pa ** -1, "compressibility"),
    (astrophys.electron * si.m ** -3, "electron density"),
    (astrophys.electron * si.m ** -2 * si.s ** -1, "electron flux"),
    (si.kg / si.m ** 2, "surface mass density"),
    (si.W / si.m ** 2 / si.sr, "radiance"),
    (si.J / si.mol, "chemical potential"),
    (si.kg / si.m, "linear density"),
    (si.H ** -1, "magnetic reluctance"),
    (si.W / si.K, "thermal conductance"),
    (si.K / si.W, "thermal resistance"),
    (si.K * si.m / si.W, "thermal resistivity"),
    (si.N / si.s, "yank"),
    (si.S * si.m ** 2 / si.mol, "molar conductivity"),
    (si.m ** 2 / si.V / si.s, "electrical mobility"),
    (si.lumen / si.W, "luminous efficacy"),
    (si.m ** 2 / si.kg, {"opacity", "mass attenuation coefficient"}),
    (si.kg * si.m ** -2 * si.s ** -1, {"mass flux", "momentum density"}),
    (si.m ** -3, "number density"),
    (si.m ** -2 * si.s ** -1, "particle flux"),
]

_physical_unit_mapping = {}
_unit_physical_mapping = {}
_name_physical_mapping = {}


def _physical_type_from_str(name):
    """
    Return the `PhysicalType` instance associated with the name of a
    physical.
    """
    if name == "unknown":
        raise ValueError("cannot uniquely identify an unknown physical type.")
    if name in _name_physical_mapping:
        return _name_physical_mapping[name]
    else:
        raise ValueError(f"{name!r} is not a known physical type.")


def _replace_temperatures_with_kelvin(unit):
    """
    If a unit contains a temperature unit besides kelvin, then replace
    that unit with kelvin.

    Temperatures cannot be converted directly between K, °F, °C, and
    °Ra, in particular since there would be different conversions for
    T and ΔT.  However, each of these temperatures each represents the
    physical type.  Replacing the different temperature units with
    kelvin allows the physical type to be treated consistently.
    """
    physical_type_id = unit._get_physical_type_id()

    physical_type_id_components = []
    substitution_was_made = False

    for base, power in physical_type_id:
        if base in ["deg_F", "deg_C", "deg_R"]:
            base = "K"
            substitution_was_made = True
        physical_type_id_components.append((base, power))

    if substitution_was_made:
        return core.Unit._from_physical_type_id(tuple(physical_type_id_components))
    else:
        return unit


def _standardize_physical_type_names(physical_type_input):
    """
    Convert a string or `set` of strings into a `set` containing
    string representations of physical types.

    The strings provided in ``physical_type_input`` can each contain
    multiple physical types that are separated by a regular slash.
    Underscores are treated as spaces so that variable names could
    be identical to physical type names.
    """
    if isinstance(physical_type_input, str):
        physical_type_input = {physical_type_input}

    standardized_physical_types = set()

    for ptype_input in physical_type_input:
        if not isinstance(ptype_input, str):
            raise ValueError(f"expecting a string, but got {ptype_input}")
        input_set = set(ptype_input.split("/"))
        processed_set = {s.strip().replace("_", " ") for s in input_set}
        standardized_physical_types |= processed_set

    return standardized_physical_types


class PhysicalType:
    """
    Represents the physical type(s) that are dimensionally compatible
    with a set of units.

    Instances of this class should be accessed through either
    `get_physical_type` or by using the
    `~astropy.units.core.UnitBase.physical_type` attribute of units.
    This class is not intended to be instantiated directly in user code.

    Parameters
    ----------
    unit : `~astropy.units.Unit`
        The unit to be represented by the physical type.

    physical_types : `str` or `set` of `str`
        A `str` representing the name of the physical type of the unit,
        or a `set` containing strings that represent one or more names
        of physical types.

    Notes
    -----
    A physical type will be considered equal to an equivalent
    `PhysicalType` instance (recommended) or a string that contains a
    name of the physical type.  The latter method is not recommended
    in packages, as the names of some physical types may change in the
    future.

    To maintain backwards compatibility, two physical type names may be
    included in one string if they are separated with a slash (e.g.,
    ``"momentum/impulse"``).  String representations of physical types
    may include underscores instead of spaces.

    Examples
    --------
    `PhysicalType` instances may be accessed via the
    `~astropy.units.core.UnitBase.physical_type` attribute of units.

    >>> import astropy.units as u
    >>> u.meter.physical_type
    PhysicalType('length')

    `PhysicalType` instances may also be accessed by calling
    `get_physical_type`. This function will accept a unit, a string
    containing the name of a physical type, or the number one.

    >>> u.get_physical_type(u.m ** -3)
    PhysicalType('number density')
    >>> u.get_physical_type("volume")
    PhysicalType('volume')
    >>> u.get_physical_type(1)
    PhysicalType('dimensionless')

    Some units are dimensionally compatible with multiple physical types.
    A pascal is intended to represent pressure and stress, but the unit
    decomposition is equivalent to that of energy density.

    >>> pressure = get_physical_type("pressure")
    >>> pressure
    PhysicalType({'energy density', 'pressure', 'stress'})
    >>> 'energy density' in pressure
    True

    Physical types can be tested for equality against other physical
    type objects or against strings that may contain the name of a
    physical type.

    >>> area = (u.m ** 2).physical_type
    >>> area == u.barn.physical_type
    True
    >>> area == "area"
    True

    Multiplication, division, and exponentiation are enabled so that
    physical types may be used for dimensional analysis.

    >>> length = u.pc.physical_type
    >>> area = (u.cm ** 2).physical_type
    >>> length * area
    PhysicalType('volume')
    >>> area / length
    PhysicalType('length')
    >>> length ** 3
    PhysicalType('volume')

    may also be performed using a string that contains the name of a
    physical type.

    >>> "length" * area
    PhysicalType('volume')
    >>> "area" / length
    PhysicalType('length')

    Unknown physical types are labelled as ``"unknown"``.

    >>> (u.s ** 13).physical_type
    PhysicalType('unknown')

    Dimensional analysis may be performed for unknown physical types too.

    >>> length_to_19th_power = (u.m ** 19).physical_type
    >>> length_to_20th_power = (u.m ** 20).physical_type
    >>> length_to_20th_power / length_to_19th_power
    PhysicalType('length')
    """

    def __init__(self, unit, physical_types):
        self._unit = _replace_temperatures_with_kelvin(unit)
        self._physical_type_id = self._unit._get_physical_type_id()
        self._physical_type = _standardize_physical_type_names(physical_types)
        self._physical_type_list = sorted(self._physical_type)

    def __iter__(self):
        yield from self._physical_type_list

    def __eq__(self, other):
        """
        Return `True` if ``other`` represents a physical type that is
        consistent with the physical type of the `PhysicalType` instance.
        """
        if isinstance(other, PhysicalType):
            return self._physical_type_id == other._physical_type_id
        elif isinstance(other, str):
            other = _standardize_physical_type_names(other)
            return other.issubset(self._physical_type)
        else:
            return NotImplemented

    def __ne__(self, other):
        equality = self.__eq__(other)
        return not equality if isinstance(equality, bool) else NotImplemented

    def _name_string_as_ordered_set(self):
        return "{" + str(self._physical_type_list)[1:-1] + "}"

    def __repr__(self):
        if len(self._physical_type) == 1:
            names = "'" + self._physical_type_list[0] + "'"
        else:
            names = self._name_string_as_ordered_set()
        return f"PhysicalType({names})"

    def __str__(self):
        if len(self._physical_type) == 1:
            return self._physical_type_list[0]
        else:
            return self._name_string_as_ordered_set()

    @staticmethod
    def _dimensionally_compatible_unit(obj):
        """
        Return a unit that corresponds to the provided argument.

        If a unit is passed in, return that unit.  If a physical type
        (or a `str` with the name of a physical type) is passed in,
        return a unit that corresponds to that physical type.  If the
        number equal to ``1`` is passed in, return a dimensionless unit.
        Otherwise, return `NotImplemented`.
        """
        if isinstance(obj, core.UnitBase):
            return _replace_temperatures_with_kelvin(obj)
        elif isinstance(obj, PhysicalType):
            return obj._unit
        elif isinstance(obj, numbers.Real) and obj == 1:
            return core.dimensionless_unscaled
        elif isinstance(obj, str):
            return _physical_type_from_str(obj)._unit
        else:
            return NotImplemented

    def _dimensional_analysis(self, other, operation):
        other_unit = self._dimensionally_compatible_unit(other)
        if other_unit is NotImplemented:
            return NotImplemented
        other_unit = _replace_temperatures_with_kelvin(other_unit)
        new_unit = getattr(self._unit, operation)(other_unit)
        return new_unit.physical_type

    def __mul__(self, other):
        return self._dimensional_analysis(other, "__mul__")

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self._dimensional_analysis(other, "__truediv__")

    def __rtruediv__(self, other):
        other = self._dimensionally_compatible_unit(other)
        if other is NotImplemented:
            return NotImplemented
        return other.physical_type._dimensional_analysis(self, "__truediv__")

    def __pow__(self, power):
        return (self._unit ** power).physical_type

    def __hash__(self):
        return hash(self._physical_type_id)

    def __len__(self):
        return len(self._physical_type)

    # We need to prevent operations like where a Unit instance left
    # multiplies a PhysicalType instance from returning a `Quantity`
    # instance with a PhysicalType as the value.  We can do this by
    # preventing np.array from casting a PhysicalType instance as
    # an object array.
    __array__ = None


def _get_names_in_use_except_from(unit):
    """
    Get a `set` containing all of the physical type names in use, except
    for the names that correspond to ``unit``.
    """
    names_in_use = set(_unit_physical_mapping.keys())
    id = unit._get_physical_type_id()
    known_physical_type = id in _physical_unit_mapping.keys()
    names_for_unit = set(_physical_unit_mapping[id]) if known_physical_type else set()
    return names_in_use - names_for_unit


def def_physical_type(unit, name):
    """
    Add a mapping between a unit and the corresponding physical type(s).

    If a physical type already exists for a unit, add new physical type
    names so long as those names are not already in use for other
    physical types.

    Parameters
    ----------
    unit : `~astropy.units.Unit`
        The unit to be represented by the physical type.

    name : `str` or `set` of `str`
        A `str` representing the name of the physical type of the unit,
        or a `set` containing strings that represent one or more names
        of physical types.

    Raises
    ------
    ValueError
        If a physical type name is already in use for another unit, or
        if attempting to name a unit as ``"unknown"``.
    """
    physical_type_names = _standardize_physical_type_names(name)

    if "unknown" in physical_type_names:
        raise ValueError("cannot uniquely define an unknown physical type")

    names_for_other_units = _get_names_in_use_except_from(unit)

    names_already_in_use = physical_type_names & names_for_other_units
    if names_already_in_use:
        raise ValueError(
            f"the following physical type names are already in use: "
            f"{names_already_in_use}."
        )

    physical_type_id = unit._get_physical_type_id()
    unit_already_in_use = physical_type_id in _physical_unit_mapping
    if unit_already_in_use:
        physical_type = _physical_unit_mapping[physical_type_id]
        physical_type_names |= set(physical_type)
        physical_type.__init__(unit, physical_type_names)
    else:
        physical_type = PhysicalType(unit, physical_type_names)
        _physical_unit_mapping[physical_type_id] = physical_type

    for ptype in physical_type:
        _unit_physical_mapping[ptype] = physical_type_id

    for ptype_name in physical_type_names:
        _name_physical_mapping[ptype_name] = physical_type


def get_physical_type(obj):
    """
    Return the physical type that corresponds to a unit (or another
    physical type representation).

    Parameters
    ----------
    obj : `~astropy.units.UnitBase`, `str`, number, or `~astropy.units.Quantity`
        The object for which to retreive the physical type.  This object
        should be a unit, the name of a physical type, the number one, a
        `~astropy.units.Quantity`, or an `object` that has a ``unit``
        attribute that is a unit.

    Returns
    -------
    physical_type : `~astropy.units.PhysicalType`
        A representation of the physical type(s) of the unit.

    Examples
    --------
    The physical type may be retrieved from a unit or a
    `~astropy.units.Quantity`.

    >>> import astropy.units as u
    >>> u.get_physical_type(u.meter ** -2)
    PhysicalType('column density')
    >>> u.get_physical_type(0.62 * u.barn * u.Mpc)
    PhysicalType('volume')

    The physical type may also be retrieved by providing a `str` that
    contains the name of a physical type.

    >>> u.get_physical_type("energy")
    PhysicalType({'energy', 'torque', 'work'})

    The number one corresponds to a dimensionless physical type.

    >>> u.get_physical_type(1)
    PhysicalType('dimensionless')
    """
    if isinstance(obj, str):
        return _physical_type_from_str(obj)

    if hasattr(obj, "unit"):
        unit = obj.unit
    elif isinstance(obj, numbers.Real) and obj == 1:
        unit = core.dimensionless_unscaled
    else:
        unit = obj

    unit = _replace_temperatures_with_kelvin(unit)
    physical_type_id = unit._get_physical_type_id()
    unit_has_known_physical_type = physical_type_id in _physical_unit_mapping
    if unit_has_known_physical_type:
        return _physical_unit_mapping[physical_type_id]
    else:
        return PhysicalType(unit, "unknown")


for unit, physical_type in _units_and_physical_types:
    def_physical_type(unit, physical_type)

del unit, physical_type
