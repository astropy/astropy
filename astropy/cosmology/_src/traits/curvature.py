# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Global Curvature.

This is private API. See `~astropy.cosmology.traits` for public API.

"""

__all__ = ["CurvatureComponent"]

import abc

import numpy as np
from numpy.typing import ArrayLike, NDArray

from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class CurvatureComponent:
    """The object has attributes and methods related to the global curvature.

    This is a trait class; it is not meant to be instantiated directly, but
    instead to be used as a mixin to other classes.

    """

    @property
    @abc.abstractmethod
    def Ok0(self) -> float | np.floating:
        """Omega curvature; the effective curvature density/critical density at z=0."""
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def is_flat(self) -> bool:
        """Return `bool`; `True` if the cosmology is globally flat."""
        raise NotImplementedError

    @deprecated_keywords("z", since="7.0")
    def Ok(self, z: Quantity | ArrayLike) -> NDArray[np.floating]:
        """Return the equivalent density parameter for curvature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Ok : ndarray
            The equivalent density parameter for curvature at each redshift.

            .. versionchanged:: 7.2
                Always returns a numpy object, never a `float`.

        Examples
        --------
        >>> import numpy as np
        >>> from astropy.cosmology import Planck18, units as cu

        >>> Planck18.Ok(2)
        array(0.)

        >>> Planck18.Ok([1, 2])
        array([0., 0.])

        >>> Planck18.Ok(np.array([2]))
        array([0.])

        >>> Planck18.Ok(2 * cu.redshift)
        array(0.)

        >>> cosmo = Planck18.clone(Ode0=0.71, to_nonflat=True)

        >>> cosmo.Ok0
        np.float64(-0.021153694455455927)

        >>> cosmo.Ok(100)
        np.float64(-0.0006557825253017665)

        """
        z = aszarr(z)
        if self.Ok0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros(getattr(z, "shape", ()))
        return self.Ok0 * (z + 1.0) ** 2 * self.inv_efunc(z) ** 2
