# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""CMB Temperature.

This is private API. See `~astropy.cosmology.parts` for public API.

"""

__all__ = ["_TemperatureCMB"]


from numpy.typing import ArrayLike

from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class _TemperatureCMB:
    """The object has attributes and methods for computing the cosmological background temperature.

    Attributes
    ----------
    Tcmb0 : |Quantity|
        Temperature of the CMB at redshift 0.

    Methods
    -------
    Tcmb
        Compute the CMB temperature at a given redshift.
    """

    Tcmb0: Quantity
    """Temperature of the CMB as |Quantity| at z=0."""

    @deprecated_keywords("z", since="7.0")
    def Tcmb(self, z: Quantity | ArrayLike) -> Quantity:
        """Return the CMB temperature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Tcmb : Quantity ['temperature']
            The temperature of the CMB in K.
        """
        return self.Tcmb0 * (aszarr(z) + 1.0)
