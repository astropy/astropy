# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numbers
from collections.abc import Iterator
from typing import Final, Union

from astropy.utils.compat import COPY_IF_NEEDED

from . import core, quantity
from .typing import PhysicalTypeID, QuantityLike, UnitPowerLike

__all__: Final = ["PhysicalType", "def_physical_type", "get_physical_type"]


class PhysicalType:
    """
    Represents the physical type(s) that are dimensionally compatible
    with a set of units.

    Instances of this class should be accessed through either
    `get_physical_type` or by using the
    `~astropy.units.core.UnitBase.physical_type` attribute of units.
    This class is not intended to be instantiated directly in user code.

    For a list of physical types, see `astropy.units.physical`.

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

    >>> pressure = u.get_physical_type("pressure")
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

    def __init__(self, unit: core.UnitBase, physical_types: str | set[str]) -> None:
        self._unit = _replace_temperatures_with_kelvin(unit)
        self._physical_type = sorted(_standardize_physical_type_names(physical_types))

    def __iter__(self) -> Iterator[str]:
        yield from self._physical_type

    def __eq__(self, other: object) -> bool:
        """
        Return `True` if ``other`` represents a physical type that is
        consistent with the physical type of the `PhysicalType` instance.
        """
        if isinstance(other, PhysicalType):
            return self._unit._physical_type_id == other._unit._physical_type_id
        elif isinstance(other, str):
            other = _standardize_physical_type_names(other)
            return other.issubset(self._physical_type)
        else:
            return NotImplemented

    def __repr__(self) -> str:
        if len(self._physical_type) == 1:
            names = "'" + self._physical_type[0] + "'"
        else:
            names = "{" + str(self._physical_type)[1:-1] + "}"
        return f"PhysicalType({names})"

    def __str__(self) -> str:
        return "/".join(self._physical_type)

    @staticmethod
    def _dimensionally_compatible_unit(obj: object) -> core.UnitBase | None:
        """
        Return a unit that corresponds to the provided argument.
        """
        if isinstance(obj, core.UnitBase):
            return _replace_temperatures_with_kelvin(obj)
        elif isinstance(obj, PhysicalType):
            return obj._unit
        elif isinstance(obj, numbers.Real) and obj == 1:
            return core.dimensionless_unscaled
        elif isinstance(obj, str):
            return _physical_type_from_str(obj)._unit
        return None

    def __mul__(
        self, other: Union["PhysicalType", core.UnitBase, numbers.Real, str]
    ) -> "PhysicalType":
        if other_unit := self._dimensionally_compatible_unit(other):
            return (self._unit * other_unit).physical_type
        return NotImplemented

    def __rmul__(
        self, other: Union["PhysicalType", core.UnitBase, str]
    ) -> "PhysicalType":
        return self.__mul__(other)

    def __truediv__(
        self, other: Union["PhysicalType", core.UnitBase, numbers.Real, str]
    ) -> "PhysicalType":
        if other_unit := self._dimensionally_compatible_unit(other):
            return (self._unit / other_unit).physical_type
        return NotImplemented

    def __rtruediv__(
        self, other: Union["PhysicalType", core.UnitBase, numbers.Real, str]
    ) -> "PhysicalType":
        if other_unit := self._dimensionally_compatible_unit(other):
            return (other_unit / self._unit).physical_type
        return NotImplemented

    def __pow__(self, power: UnitPowerLike) -> "PhysicalType":
        return (self._unit**power).physical_type

    def __hash__(self) -> int:
        return hash(self._unit._physical_type_id)

    def __len__(self) -> int:
        return len(self._physical_type)

    # We need to prevent operations like where a Unit instance left
    # multiplies a PhysicalType instance from returning a `Quantity`
    # instance with a PhysicalType as the value.  We can do this by
    # preventing np.array from casting a PhysicalType instance as
    # an object array.
    __array__: Final = None


_physical_unit_mapping: Final[dict[PhysicalTypeID, PhysicalType]] = {}
_unit_physical_mapping: Final[dict[str, PhysicalTypeID]] = {}
_name_physical_mapping: Final[dict[str, PhysicalType]] = {}
# mapping from attribute-accessible name (no spaces, etc.) to the actual name.
_attrname_physical_mapping: Final[dict[str, PhysicalType]] = {}


def _physical_type_from_str(name: str) -> PhysicalType:
    """
    Return the `PhysicalType` instance associated with the name of a
    physical type.
    """
    if name == "unknown":
        raise ValueError("cannot uniquely identify an 'unknown' physical type.")

    elif name in _attrname_physical_mapping:
        return _attrname_physical_mapping[name]  # convert attribute-accessible
    elif name in _name_physical_mapping:
        return _name_physical_mapping[name]
    else:
        raise ValueError(f"{name!r} is not a known physical type.")


def _replace_temperatures_with_kelvin(unit: core.UnitBase) -> core.UnitBase:
    """Replace °F, and °C in the bases of `unit` with K.

    The Kelvin, Celsius and Fahrenheit scales have different zero points,
    which is a problem for the unit conversion machinery (without the
    `temperature` equivalency). Replacing °F, and °C with kelvin allows the
    physical type to be treated consistently. The Rankine scale has the
    same zero point as the Kelvin scale, so degrees Rankine do not have to
    be special-cased.
    """
    physical_type_id = unit._physical_type_id

    physical_type_id_components = []
    substitution_was_made = False

    for base, power in physical_type_id:
        if base in ["deg_F", "deg_C"]:
            base = "K"
            substitution_was_made = True
        physical_type_id_components.append((base, power))

    if substitution_was_made:
        return core.Unit._from_physical_type_id(tuple(physical_type_id_components))
    else:
        return unit


def _standardize_physical_type_names(physical_type_input: str | set[str]) -> set[str]:
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


def def_physical_type(unit: core.UnitBase, name: str | set[str]) -> None:
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

    Notes
    -----
    For a list of physical types, see `astropy.units.physical`.
    """
    physical_type_id = unit._physical_type_id
    physical_type_names = _standardize_physical_type_names(name)

    if "unknown" in physical_type_names:
        raise ValueError("cannot uniquely define an unknown physical type")

    names_for_other_units = set(_unit_physical_mapping.keys()).difference(
        _physical_unit_mapping.get(physical_type_id, {})
    )
    names_already_in_use = physical_type_names & names_for_other_units
    if names_already_in_use:
        raise ValueError(
            "the following physical type names are already in use: "
            f"{names_already_in_use}."
        )

    unit_already_in_use = physical_type_id in _physical_unit_mapping
    if unit_already_in_use:
        physical_type = _physical_unit_mapping[physical_type_id]
        physical_type._physical_type = sorted(physical_type_names | set(physical_type))
    else:
        physical_type = PhysicalType(unit, physical_type_names)
        _physical_unit_mapping[physical_type_id] = physical_type

    for ptype in physical_type:
        _unit_physical_mapping[ptype] = physical_type_id
        _name_physical_mapping[ptype] = physical_type
        # attribute-accessible name
        attr_name = ptype.replace(" ", "_").replace("(", "").replace(")", "")
        _attrname_physical_mapping[attr_name] = physical_type


def get_physical_type(
    obj: PhysicalType | str | core.UnitBase | QuantityLike,
) -> PhysicalType:
    """
    Return the physical type that corresponds to a unit (or another
    physical type representation).

    Parameters
    ----------
    obj : quantity-like or `~astropy.units.PhysicalType`-like
        An object that (implicitly or explicitly) has a corresponding
        physical type. This object may be a unit, a
        `~astropy.units.Quantity`, an object that can be converted to a
        `~astropy.units.Quantity` (such as a number or array), a string
        that contains a name of a physical type, or a
        `~astropy.units.PhysicalType` instance.

    Returns
    -------
    `~astropy.units.PhysicalType`
        A representation of the physical type(s) of the unit.

    Notes
    -----
    For a list of physical types, see `astropy.units.physical`.

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

    Numbers and arrays of numbers correspond to a dimensionless physical
    type.

    >>> u.get_physical_type(1)
    PhysicalType('dimensionless')
    """
    if isinstance(obj, PhysicalType):
        return obj

    if isinstance(obj, str):
        return _physical_type_from_str(obj)

    if isinstance(obj, core.UnitBase):
        unit = obj
    else:
        try:
            unit = quantity.Quantity(obj, copy=COPY_IF_NEEDED).unit
        except TypeError as exc:
            raise TypeError(f"{obj} does not correspond to a physical type.") from exc

    unit = _replace_temperatures_with_kelvin(unit)
    physical_type_id = unit._physical_type_id
    unit_has_known_physical_type = physical_type_id in _physical_unit_mapping

    if unit_has_known_physical_type:
        return _physical_unit_mapping[physical_type_id]
    else:
        return PhysicalType(unit, "unknown")
