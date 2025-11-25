# Licensed under a 3-clause BSD style license - see LICENSE.rst
r"""Neutrino component trait.

This is private API. See `~astropy.cosmology.traits` for public API.
"""

__all__ = ["NeutrinoComponent"]

from abc import abstractmethod
from collections.abc import Callable
from functools import cached_property
from typing import Final

import numpy as np
from numpy.typing import ArrayLike

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity

# Physics constants for neutrino calculations
TEMP_NEUTRINO: Final = 0.7137658555036082  # (4/11)^1/3
NEUTRINO_FERMI_DIRAC_CORRECTION: Final = 0.22710731766  # 7/8 (4/11)^4/3


##############################################################################


class NeutrinoComponent:
    r"""The cosmology has attributes and methods for the neutrino density.

    This trait handles both massless neutrinos (relativistic, radiation-like)
    and massive neutrinos (with complex evolution).

    This is an abstract trait. Subclasses must implement:

    - `has_massive_nu` (property): Whether there are massive neutrinos
    - `Onu0` (property): Neutrino density parameter at z=0
    - `nu_relative_density` (method): Neutrino-to-photon density ratio at redshift z

    The parent class must provide ``Tcmb0`` (CMB temperature) and ``Ogamma``
    (method to compute photon density at redshift z).

    Notes
    -----
    The density in neutrinos is given by:

    .. math::

        \rho_{\nu} \left(a\right) = 0.2271 \, N_{eff} \,
        f\left(m_{\nu} a / T_{\nu 0} \right) \,
        \rho_{\\gamma} \left( a \right)

    where

    .. math::

        f \left(y\right) = \frac{120}{7 \pi^4}
        \int_0^{\\infty} \, dx \frac{x^2 \\sqrt{x^2 + y^2}}
        {e^x + 1}

    assuming that all neutrino species have the same mass.

    If they have different masses, a similar term is calculated for each
    one.

    Note that ``f`` has the asymptotic behavior :math:`f(0) = 1`. This
    method returns :math:`0.2271 f` using an analytical fitting formula
    (Komatsu et al., 2011), ApJS, 192, 18.

    The neutrino density evolution depends on whether neutrinos are massive or massless:

    - **Massless neutrinos**: Behave like radiation with density scaling as (1+z)^4.
      The density is simply proportional to the photon density with a constant ratio
      determined by Neff and Fermi-Dirac statistics.

    - **Massive neutrinos**: Have complex evolution that transitions from relativistic
      (radiation-like) at early times to non-relativistic (matter-like) at late times.
      The implementation typically uses the Komatsu fitting formula (Komatsu
      et al., 2011) for computational efficiency.

    References
    ----------
    Komatsu et al. (2011), "Seven-Year Wilkinson Microwave Anisotropy Probe
    (WMAP) Observations: Cosmological Interpretation", ApJS, 192, 18.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.cosmology.traits import NeutrinoComponent
    >>> NEUTRINO_FERMI_DIRAC_CORRECTION = 0.22710731766  # 7/8 (4/11)^4/3
    >>>
    >>> class ExampleNeutrinoCosmology(NeutrinoComponent):
    ...     def __init__(self):
    ...         self.Tcmb0 = 2.7255 * u.K
    ...         self.Neff = 3.046
    ...         self.Ogamma0 = 5e-5
    ...     @property
    ...     def has_massive_nu(self):
    ...         return False
    ...     @property
    ...     def Onu0(self):
    ...         return NEUTRINO_FERMI_DIRAC_CORRECTION * self.Neff * self.Ogamma0
    ...     def nu_relative_density(self, z):
    ...         return NEUTRINO_FERMI_DIRAC_CORRECTION * self.Neff * np.ones_like(np.asarray(z))
    ...     def Ogamma(self, z):
    ...         return self.Ogamma0 * (np.asarray(z) + 1.0) ** 4
    """

    # Type annotations for dependencies (provided by parent class)
    Tcmb0: Quantity
    Ogamma: Callable[[ArrayLike], FArray]

    @property
    @abstractmethod
    def has_massive_nu(self) -> bool:
        """Does this cosmology have at least one massive neutrino species?

        Returns
        -------
        has_massive_nu : bool
            True if at least one neutrino species has non-zero mass.

        Notes
        -----
        Subclasses must implement this property.
        """

    @property
    @abstractmethod
    def Onu0(self) -> float:
        """Omega nu; the density/critical density of neutrinos at z=0.

        Returns
        -------
        Onu0 : float
            The density parameter for neutrinos at z=0.

        Notes
        -----
        Subclasses must implement this property.
        """

    @abstractmethod
    def nu_relative_density(self, z: Quantity | ArrayLike) -> FArray:
        r"""Neutrino density function relative to the energy density in photons.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        f : array
            The neutrino density scaling factor relative to the density in
            photons at each redshift.

        Notes
        -----
        Subclasses must implement this method. For massless neutrinos, this
        should return a constant. For massive neutrinos, this should use an
        appropriate fitting formula (e.g., Komatsu et al. 2011).
        """

    @cached_property
    def Tnu0(self) -> Quantity:
        """Temperature of the neutrino background as |Quantity| at z=0.

        Returns
        -------
        Tnu0 : Quantity ['temperature']
            The neutrino temperature at z=0 in Kelvin.

        Notes
        -----
        The neutrino temperature is related to the CMB temperature by:

        .. math::

            T_{\\nu 0} = \\left(\\frac{4}{11}\\right)^{1/3} T_{CMB}

        This comes from the decoupling of neutrinos before electron-positron
        annihilation. See Weinberg 'Cosmology' p 154 eq (3.1.21).
        """
        # The constant in front is (4/11)^1/3 -- see any cosmology book for an
        # explanation -- for example, Weinberg 'Cosmology' p 154 eq (3.1.21).
        return TEMP_NEUTRINO * self.Tcmb0

    @deprecated_keywords("z", since="7.0")
    def Onu(self, z: Quantity | ArrayLike) -> FArray:
        r"""Return the density parameter for neutrinos at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Onu : ndarray
            The energy density of neutrinos relative to the critical density at
            each redshift. Note that this includes their kinetic energy (if
            they have mass), so it is not equal to the commonly used
            :math:`\sum \frac{m_{\nu}}{94 eV}`, which does not include
            kinetic energy.

        Notes
        -----
        The neutrino density parameter evolves with redshift according to:

        .. math::

            \\Omega_{\\nu}(z) = \\Omega_{\\gamma}(z) \\times f(z)

        where f(z) is the neutrino-to-photon density ratio computed by
        nu_relative_density(z). For massless neutrinos, f(z) is constant.
        For massive neutrinos, f(z) evolves as neutrinos transition from
        relativistic to non-relativistic.
        """
        z = aszarr(z)
        if self.Onu0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros_like(z)
        return self.Ogamma(z) * self.nu_relative_density(z)

    @deprecated_keywords("z", since="7.0")
    def Tnu(self, z: Quantity | ArrayLike) -> Quantity:
        """Return the neutrino temperature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Tnu : Quantity ['temperature']
            The temperature of the cosmic neutrino background in K.

        Notes
        -----
        The neutrino temperature scales with redshift as:

        .. math::

            T_{\\nu}(z) = T_{\\nu 0} (1 + z)

        This simple scaling applies to both massless and massive neutrinos,
        as the temperature depends only on the expansion of the universe.
        """
        return self.Tnu0 * (aszarr(z) + 1.0)
