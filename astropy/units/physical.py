# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Defines the physical types that correspond to different units.

This module is not intended for use by user code directly.  Instead,
the physical unit name(s) of a `Unit` can be obtained using its
``physical_type`` property.
"""

import numbers
import re
import warnings

from . import core
from . import si
from . import astrophys
from . import cgs
from . import misc

__all__ = ["def_physical_type", "get_physical_type", "PhysicalType"]

_units_and_physical_types = [
    (core.Unit(1), "dimensionless"),
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
    (si.kg / si.m ** 3, "mass density"),
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
    (si.Ohm, {"electrical resistance", "impedance", "reactance"}),
    (si.Ohm * si.m, "electrical resistivity"),
    (si.S, "electrical conductance"),
    (si.S / si.m, "electrical conductivity"),
    (si.F, "electrical capacitance"),
    (si.C * si.m, "electrical dipole moment"),
    (si.A / si.m ** 2, "electrical current density"),
    (si.V / si.m, "electrical field strength"),
    (
        si.C / si.m ** 2,
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
    (astrophys.count, "number of counts"),
    (misc.pixel, "number of elements on 2D regular grid"),
    (misc.voxel, "number of elements on 3D regular grid"),
    (astrophys.electron, "number of electrons"),
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
    (si.kg * si.m ** -2 * si.s ** -1, "momentum density"),
]

_physical_unit_mapping = {}
_unit_physical_mapping = {}


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


class PhysicalType:
    """
    Represents and provides information on the physical type(s) that are
    associated with a set of units.

    Instances of this class should be accessed using the ``physical_type``
    property of `Unit` instances.  This class is not intended to be
    instantiated directly in user code.

    Parameters
    ----------
    unit : u.Unit
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
    The preferred method to access a physical type is by accessing the
    ``physical_type`` attribute of a unit.

    >>> import astropy.units as u
    >>> u.meter.physical_type
    'length'

    Some units correspond to multiple physical types.

    >>> pressure = u.Pa.physical_type
    >>> print(pressure)
    {'energy density', 'pressure', 'stress'}
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
    'volume'
    >>> area / length
    'length'
    >>> length ** 3
    'volume'

    Unknown physical types are labelled as ``"unknown"``.

    >>> (u.s ** 13).physical_type
    'unknown'

    Dimensional analysis may be performed for unknown physical types too.

    >>> length_to_19th_power = (u.m ** 19).physical_type
    >>> length_to_20th_power = (u.m ** 20).physical_type
    >>> length_to_19th_power, length_to_20th_power
    ('unknown', 'unknown')
    >>> length_to_20th_power / length_to_19th_power
    'length'
    """

    @staticmethod
    def _standardize_physical_type_names(physical_type_input):
        """
        Convert a string or `set` of strings into a `set` containing
        string representations of physical types.

        The strings provided in ``physical_type_input`` can each contain
        multiple physical types that are separated by a regular slash.
        Underscores are treated as spaces so that variable names could
        be identical to physical type names.
        """

        if not isinstance(physical_type_input, (set, str)):
            raise TypeError(
                f"expecting a string or a set containing strings that "
                f"represent a physical type, not {physical_type_input}"
            )

        if isinstance(physical_type_input, str):
            physical_type_input = {physical_type_input}

        standardized_physical_types = set()

        for ptype_input in physical_type_input:
            if not isinstance(ptype_input, str):
                raise TypeError(f"expecting a string, but got {ptype_input}")
            input_set = set(ptype_input.split("/"))
            processed_set = {s.strip().replace("_", " ") for s in input_set}
            standardized_physical_types |= processed_set

        return standardized_physical_types

    def __init__(self, unit, physical_types):

        if not isinstance(unit, core.UnitBase):
            raise TypeError(f"Expecting a unit, not {unit}.")

        self._unit = _replace_temperatures_with_kelvin(unit)
        self._physical_type_id = self._unit._get_physical_type_id()
        self._physical_type = self._standardize_physical_type_names(physical_types)
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
            other = self._standardize_physical_type_names(other)
            return other.issubset(self._physical_type)
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return repr(str(self)) if len(self._physical_type) == 1 else str(self)

    def __str__(self):
        if len(self._physical_type) == 1:
            return list(self._physical_type)[0]
        else:
            return "{" + str(sorted(self._physical_type))[1:-1] + "}"

    @staticmethod
    def _identify_unit_from_unit_or_physical_type(obj):
        """
        If a unit is passed in, return that unit.  If a physical type is
        passed in, return a unit that corresponds to that physical type.
        """
        if isinstance(obj, core.UnitBase):
            return obj
        elif isinstance(obj, PhysicalType):
            return obj._unit
        else:
            raise TypeError("Expecting a unit or a physical type")

    def __mul__(self, other):
        other_unit = self._identify_unit_from_unit_or_physical_type(other)
        other_unit = _replace_temperatures_with_kelvin(other_unit)
        new_unit = self._unit * other_unit
        return new_unit.physical_type

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other_unit = self._identify_unit_from_unit_or_physical_type(other)
        other_unit = _replace_temperatures_with_kelvin(other_unit)
        new_unit = self._unit / other_unit
        return new_unit.physical_type

    def __rtruediv__(self, other):
        other_unit = self._identify_unit_from_unit_or_physical_type(other)
        other_unit = _replace_temperatures_with_kelvin(other_unit)
        new_unit = other_unit / self._unit
        return new_unit.physical_type

    def __pow__(self, power):
        if not isinstance(power, numbers.Real):
            raise TypeError(f"{power} is not a real number")
        return (self._unit ** power).physical_type


def def_physical_type(unit, name):
    """
    Add a mapping between a unit and the corresponding physical type(s).

    This function is not intended to be called directly in user code.

    Parameters
    ----------
    unit : u.Unit
        The unit to be represented by the physical type.

    name : `str` or `set` of `str`
        A `str` representing the name of the physical type of the unit,
        or a `set` containing strings that represent one or more names
        of physical types.

    Raises
    ------
    ValueError
        If the unit had previously been defined.
    """

    if name in ("unknown", {"unknown"}) or "unknown" in name:
        raise ValueError("unable to uniquely define an unknown physical type")

    physical_type_id = unit._get_physical_type_id()
    if physical_type_id in _physical_unit_mapping:
        raise ValueError(
            f"{physical_type_id!r} ({name!r}) already defined as "
            f"{_physical_unit_mapping[physical_type_id]!r}"
        )

    physical_type = PhysicalType(unit, name)

    for ptype in physical_type:
        if ptype in _unit_physical_mapping.keys():
            raise ValueError(f"{ptype} has already been defined")
        _unit_physical_mapping[ptype] = physical_type_id

    _physical_unit_mapping[physical_type_id] = physical_type


def get_physical_type(unit):
    """
    Return a representation of the physical type of a unit.

    This function is not intended to be called directly in user code.
    Physical types should be accessed using the ``physical_type``
    property of `Unit` instances.

    Parameters
    ----------
    unit : `~astropy.units.UnitBase` instance
        The unit for which to get the physical type.

    Returns
    -------
    physical : PhysicalType
        A representation of the physical type(s) of the unit.
    """
    if not isinstance(unit, core.UnitBase):
        raise TypeError("the input to get_physical_type must be a unit.")

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
