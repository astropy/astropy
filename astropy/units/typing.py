"""Typing module for supporting type annotations related to :mod:`~astropy.units`."""

from __future__ import annotations

__all__ = ["QuantityLike", "UnitPower"]


from fractions import Fraction
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from astropy.units import Quantity

if TYPE_CHECKING:
    from typing import TypeAlias

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


# The classes from the standard library `numbers` module are not suitable for
# type checking (https://github.com/python/mypy/issues/3186). For now we define
# our own number types, but if a good definition becomes available upstream
# then we should switch to that.
Real: TypeAlias = int | float | Fraction | np.integer | np.floating
Complex: TypeAlias = Real | complex | np.complexfloating

UnitPower: TypeAlias = int | float | Fraction
"""A type alias for the possible powers of a unit."""
