# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""CMB Temperature.

This is private API. See `~astropy.cosmology.traits` for public API.

"""

__all__ = ("TemperatureCMB",)


from numpy.typing import ArrayLike

from astropy.cosmology._src.utils import aszarr
from astropy.units import Quantity


class TemperatureCMB:
    """The trait for computing the cosmological background temperature.

    This is a trait class; it is not meant to be instantiated directly, but instead to
    be used as a mixin to other classes.

    Examples
    --------
    For an example of a real cosmology that implements this trait, see
    :class:`~astropy.cosmology.LambdaCDM`. Here we will define an illustrative example
    class that meets the minimum API requirements, but is not cosmologically meaningful:

    >>> import astropy.units as u
    >>> from astropy.cosmology.traits import TemperatureCMB
    >>> import dataclasses

    >>> @dataclasses.dataclass(frozen=True)
    ... class ExampleHasTcmb(TemperatureCMB):
    ...     Tcmb0: u.Quantity

    >>> cosmo = ExampleHasTcmb(Tcmb0=2.7255 * u.K)
    >>> cosmo.Tcmb0
    <Quantity 2.7255 K>

    >>> cosmo.Tcmb([0.0, 1.0, 2.0])
    <Quantity [2.7255, 5.451 , 8.1765] K>

    """

    Tcmb0: Quantity
    """Temperature of the CMB at z=0."""

    def Tcmb(self, z: Quantity | ArrayLike, /) -> Quantity:
        """Compute the CMB temperature at redshift ``z``.

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
