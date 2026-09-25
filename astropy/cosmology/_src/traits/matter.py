"""Matter component."""

import numpy as np
from numpy.typing import ArrayLike

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr
from astropy.units import Quantity

__all__ = ("MatterComponent",)


class MatterComponent:
    """The cosmology has attributes and methods for the matter density.

    This is a trait class; it is not meant to be instantiated directly, but instead to
    be used as a mixin to other classes.

    Examples
    --------
    For an example of a real cosmology that implements this trait, see
    :class:`~astropy.cosmology.LambdaCDM`. Here we will define an illustrative example
    class that meets the minimum API requirements, but is not cosmologically meaningful:

    >>> from astropy.cosmology.traits import MatterComponent
    >>> import dataclasses

    >>> @dataclasses.dataclass(frozen=True)
    ... class ExampleHasMatter(MatterComponent):
    ...     Om0: float
    ...     def inv_efunc(self, z): return 1.0  # necessary for Om(z)

    >>> cosmo = ExampleHasMatter(Om0=0.3)
    >>> cosmo.Om0
    0.3

    >>> cosmo.Om([0.0, 1.0, 2.0])
    array([0.3, 2.4, 8.1])

    """

    Om0: float | np.floating
    """Omega matter; matter density/critical density at z=0."""

    def Om(self, z: Quantity | ArrayLike, /) -> FArray:
        """Return the density parameter for non-relativistic matter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

            .. versionchanged:: 8.0
               z must be a positional argument.

        Returns
        -------
        Om : ndarray
            The density of non-relativistic matter relative to the critical
            density at each redshift.

        Notes
        -----
        This does not include neutrinos, even if non-relativistic at the
        redshift of interest.
        """
        z = aszarr(z)
        return self.Om0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2
