# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ("TotalComponent",)

from abc import abstractmethod

import numpy as np
from numpy.typing import ArrayLike

from astropy.cosmology._src.typing import FArray
from astropy.units import Quantity


class TotalComponent:
    """The cosmology has attributes and methods for the total density.

    This trait has the abstract ``Otot`` method that returns the total density
    parameter at redshift ``z``. It should be the sum of all other components.

    Examples
    --------
    For an example of a real cosmology that implements this trait, see
    :class:`~astropy.cosmology.LambdaCDM`. Here we will define an illustrative example
    class that meets the minimum API requirements, but is not cosmologically meaningful:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.cosmology.traits import TotalComponent, MatterComponent, DarkEnergyComponent
    >>> import dataclasses


    >>> @dataclasses.dataclass(frozen=True)
    ... class ExampleHasTotalDensity(TotalComponent, MatterComponent, DarkEnergyComponent):
    ...     Om0: float = 0.3
    ...     Ode0: float = 0.7
    ...     def Otot0(self): return self.Om0 + self.Ode0
    ...     def Otot(self, z): return self.Om(z) + self.Ode(z)

    """

    @property
    @abstractmethod
    def Otot0(self) -> float | np.floating:
        """Omega total; the total density/critical density at z=0."""
        raise NotImplementedError  # pragma: no cover

    @abstractmethod
    def Otot(self, z: Quantity | ArrayLike, /) -> FArray:
        """The total density parameter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshifts.

        Returns
        -------
        Otot : array
            The total density relative to the critical density at each
            redshift.
        """
        raise NotImplementedError  # pragma: no cover
