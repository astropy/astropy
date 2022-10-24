# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

# STDLIB
import warnings
from typing import Literal

# THIRD-PARTY
import numpy as np

# LOCAL
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .core import Attribute

if HAS_SCIPY:
    # THIRD-PARTY
    from scipy.spatial.transform import Rotation as ScipyRotation
else:  # make no-op version
    warnings.warn("`scipy` must be installed for the Rotation attribute")


__all__ = ["RotationAttribute"]


class RotationAttribute(Attribute):
    """
    A frame attribute that is a 3D :class:`~scipy.spatial.transform.Rotation`.

    Parameters
    ----------
    default : ndarray or None, optional
        Default value for the attribute if the user does not supply one.
        If None, taken to be the I[``shape``] identity matrix.
    secondary_attribute : str, optional
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    unit : unit-like or None, optional
        Name of a unit that the input will be converted into. If None, no
        unit-checking or conversion is performed
    """

    def __init__(self, default: np.ndarray | None = None, secondary_attribute: str="") -> None:
        if default is None and not secondary_attribute:
            default = None
        else:
            default = self.convert_input(default)[0]
        super().__init__(default, secondary_attribute)

    def convert_input(self, value: object) -> tuple[ScipyRotation | None, Literal[False]]:
        """
        Checks that the input is a `~scipy.spatial.transform.Rotation`.

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : ndarray, bool
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        TypeError
            If the input is not valid for this attribute.
        """
        if not HAS_SCIPY:
            raise ModuleNotFoundError("No module named 'scipy'")

        if value is None:
            converted = False
        elif isinstance(value, ScipyRotation):
            converted = False
        else:
            raise TypeError(
                f"value must be a `scipy.spatial.transform.Rotation`, not {value}"
            )

        return value, converted
