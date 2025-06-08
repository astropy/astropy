# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Global Curvature.

This is private API. See `~astropy.cosmology.traits` for public API.

"""

__all__ = ["CurvatureComponent"]

__doctest_requires__ = {"CurvatureComponent.Ok0": ["numpy>=2.0"]}

from typing import Any

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class CurvatureComponent:
    """The object has attributes and methods related to the global curvature.

    This is a trait class; it is not meant to be instantiated directly, but
    instead to be used as a mixin to other classes.

    """

    Ok0: float | np.floating[Any]
    """Omega curvature; the effective curvature density/critical density at z=0."""

    @property
    @abc.abstractmethod
    def is_flat(self) -> bool:
        """Return `bool`; `True` if the cosmology is globally flat."""
        raise NotImplementedError

    @deprecated_keywords("z", since="7.0")
    def Ok(self, z: Quantity | ArrayLike) -> npt.NDArray | float:
        """Return the equivalent density parameter for curvature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Ok : ndarray or float
            The equivalent density parameter for curvature at each redshift.
            Returns `float` if the input is scalar.
        """
        z = aszarr(z)
        if self.Ok0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros(z.shape) if hasattr(z, "shape") else 0.0
        return self.Ok0 * (z + 1.0) ** 2 * self.inv_efunc(z) ** 2
