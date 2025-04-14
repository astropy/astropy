"""Typing module for supporting type annotations related to :mod:`~astropy.units`."""

from __future__ import annotations

__all__ = [
    "QuantityLike",
    "UnitLike",
    "UnitPower",
    "UnitPowerLike",
    "UnitScale",
    "UnitScaleLike",
]


from fractions import Fraction
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from astropy.units import Quantity, UnitBase

if TYPE_CHECKING:
    from typing import TypeAlias


UnitLike: TypeAlias = UnitBase | str | Quantity
"""Type alias for input that can be converted to a Unit.

See :term:`unit-like`. Note that this includes only scalar quantities.
"""

# Note: Quantity is technically covered by npt.ArrayLike, but we want to
# explicitly include it here so that it is clear that we are also including
# Quantity objects in the definition of QuantityLike.
QuantityLike: TypeAlias = Quantity | npt.ArrayLike
"""Type alias for a quantity-like object.

This is an object that can be converted to a :class:`~astropy.units.Quantity` object
using the :func:`~astropy.units.Quantity` constructor.

Examples
--------
We assume the following imports:

    >>> from astropy import units as u

This is a non-exhaustive list of examples of quantity-like objects:

Integers and floats:

    >>> u.Quantity(1, u.meter)
    <Quantity 1.0 m>

    >>> u.Quantity(1.0, u.meter)
    <Quantity 1.0 m>

Lists and tuples:

    >>> u.Quantity([1.0, 2.0], u.meter)
    <Quantity [1., 2.] m>

    >>> u.Quantity((1.0, 2.0), u.meter)
    <Quantity [1., 2.] m>

Numpy arrays:

    >>> u.Quantity(np.array([1.0, 2.0]), u.meter)

:class:`~astropy.units.Quantity` objects:

    >>> u.Quantity(u.Quantity(1.0, u.meter))
    <Quantity 1.0 m>

Strings:

    >>> u.Quantity('1.0 m')
    <Quantity 1.0 m>

For more examples see the :mod:`numpy.typing` definition of
:obj:`numpy.typing.ArrayLike`.
"""


UnitPower: TypeAlias = int | float | Fraction
"""Alias for types that can be powers of the components of a
`~astropy.units.UnitBase` instance"""
UnitPowerLike: TypeAlias = UnitPower | np.integer | np.floating
"""Alias for types that can be used to create powers of the components of a
`~astropy.units.UnitBase` instance"""
UnitScale: TypeAlias = float | complex
"Alias for types that can be scale factors of a `~astropy.units.CompositeUnit`"
UnitScaleLike: TypeAlias = UnitScale | int | Fraction | np.number
"""Alias for types that can be used to create scale factors of a
`~astropy.units.CompositeUnit`"""
PhysicalTypeID: TypeAlias = tuple[tuple[str, UnitPower], ...]
