# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""CMB Temperature.

This is private API. See `~astropy.cosmology.traits` for public API.

"""

__all__ = ["TemperatureCMB"]


from numpy.typing import ArrayLike

from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class TemperatureCMB:
    """The trait for computing the cosmological background temperature."""

    Tcmb0: Quantity
    """Temperature of the CMB at z=0."""

    @deprecated_keywords("z", since="7.0")
    def Tcmb(self, z: Quantity | ArrayLike) -> Quantity:
        """Compute the CMB temperature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Tcmb : Quantity ['temperature']
            The temperature of the CMB.

        Examples
        --------
        >>> import astropy.units as u
        >>> from astropy.cosmology import Planck18, units as cu

        >>> Planck18.Tcmb(u.Quantity([0.5, 1.0], cu.redshift))
        <Quantity [4.08825, 5.451  ] K>

        >>> Planck18.Tcmb(u.Quantity(0.5, ''))
        <Quantity 4.08825 K>

        >>> Planck18.Tcmb(0.5)
        <Quantity 4.08825 K>

        >>> Planck18.Tcmb([0.5, 1.0])
        <Quantity [4.08825, 5.451  ] K>
        """
        return self.Tcmb0 * (aszarr(z) + 1.0)
