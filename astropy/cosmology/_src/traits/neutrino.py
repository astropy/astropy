# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Neutrino component trait.

This is private API. See `~astropy.cosmology.traits` for public API.
"""

__all__ = ["NeutrinoComponent"]

from collections.abc import Callable
from functools import cached_property
from typing import Any, Final

import numpy as np
from numpy.typing import ArrayLike, NDArray

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity

# Physics constants for neutrino calculations
# See Komatsu et al. 2011, eq 26 and the surrounding discussion for an explanation of
# what this does. However, this is modified to handle multiple neutrino masses by
# computing the above for each mass, then summing
NEUTRINO_FERMI_DIRAC_CORRECTION: Final = 0.22710731766  # 7/8 (4/11)^4/3

# These are purely fitting constants -- see the Komatsu paper
KOMATSU_P: Final = 1.83
KOMATSU_INVP: Final = 0.54644808743  # 1.0 / p
KOMATSU_K: Final = 0.3173

TEMP_NEUTRINO: Final = 0.7137658555036082  # (4/11)^1/3


class NeutrinoComponent:
    """The cosmology has attributes and methods for the neutrino density.

    This trait handles both massless neutrinos (relativistic, radiation-like)
    and massive neutrinos (with complex evolution using Komatsu fitting formula).

    This trait provides neutrino-related properties and methods for cosmologies.
    The parent class must provide the required dependencies listed below.

    Required Dependencies from Parent
    ----------------------------------
    Neff : float
        Number of effective neutrino species.
    m_nu : Quantity or None
        Neutrino masses in eV.
    Tcmb0 : Quantity
        CMB temperature at z=0.
    Ogamma0 : float
        Photon density parameter (typically from PhotonComponent trait).
    inv_efunc : callable
        Inverse E-function from parent cosmology.
    _nu_info : NeutrinoInfo
        Internal neutrino information structure (set in parent's __post_init__).

    Notes
    -----
    The neutrino density evolution depends on whether neutrinos are massive or massless:

    - **Massless neutrinos**: Behave like radiation with density scaling as (1+z)^4.
      The density is simply proportional to the photon density with a constant ratio
      determined by Neff and Fermi-Dirac statistics.

    - **Massive neutrinos**: Have complex evolution that transitions from relativistic
      (radiation-like) at early times to non-relativistic (matter-like) at late times.
      The implementation uses the Komatsu fitting formula [1]_ for computational
      efficiency.

    References
    ----------
    .. [1] Komatsu et al. (2011). Seven-Year Wilkinson Microwave Anisotropy Probe
           (WMAP) Observations: Cosmological Interpretation. ApJS, 192, 18.
    """

    # Type annotations for dependencies (provided by parent class)
    Neff: float
    m_nu: Quantity | None
    Tcmb0: Quantity
    Ogamma0: float
    inv_efunc: Callable[[NDArray[Any]], NDArray[Any]]
    _nu_info: Any  # NeutrinoInfo NamedTuple

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

    @property
    def has_massive_nu(self) -> bool:
        """Does this cosmology have at least one massive neutrino species?

        Returns
        -------
        has_massive_nu : bool
            True if at least one neutrino species has non-zero mass.

        Notes
        -----
        This returns False if Tcmb0 is zero, as neutrinos are turned off in
        that case regardless of the m_nu parameter.
        """
        if self.Tnu0.value == 0:
            return False
        return self._nu_info.has_massive_nu

    @cached_property
    def Onu0(self) -> float:
        """Omega nu; the density/critical density of neutrinos at z=0.

        Returns
        -------
        Onu0 : float
            The density parameter for neutrinos at z=0.

        Notes
        -----
        For massless neutrinos, this is calculated directly using the Fermi-Dirac
        correction factor:

        .. math::

            \\Omega_{\\nu 0} = 0.2271 \\, N_{eff} \\, \\Omega_{\\gamma 0}

        where 0.2271 = 7/8 (4/11)^(4/3) accounts for the temperature difference
        and Fermi-Dirac vs. Bose-Einstein statistics.

        For massive neutrinos, the calculation uses the full nu_relative_density
        method evaluated at z=0.
        """
        if self._nu_info.has_massive_nu:
            return self.Ogamma0 * self.nu_relative_density(0)
        else:
            # This case is particularly simple, so do it directly The 0.2271...
            # is 7/8 (4/11)^(4/3) -- the temperature bit ^4 (blackbody energy
            # density) times 7/8 for FD vs. BE statistics.
            return NEUTRINO_FERMI_DIRAC_CORRECTION * self.Neff * self.Ogamma0

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

    @deprecated_keywords("z", since="7.0")
    def nu_relative_density(self, z: Quantity | ArrayLike) -> FArray:
        r"""Neutrino density function relative to the energy density in photons.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        f : array
            The neutrino density scaling factor relative to the density in
            photons at each redshift.

        Notes
        -----
        The density in neutrinos is given by

        .. math::

           \rho_{\nu} \left(a\right) = 0.2271 \, N_{eff} \,
           f\left(m_{\nu} a / T_{\nu 0} \right) \,
           \rho_{\gamma} \left( a \right)

        where

        .. math::

           f \left(y\right) = \frac{120}{7 \pi^4}
           \int_0^{\infty} \, dx \frac{x^2 \sqrt{x^2 + y^2}}
           {e^x + 1}

        assuming that all neutrino species have the same mass.
        If they have different masses, a similar term is calculated for each
        one. Note that ``f`` has the asymptotic behavior :math:`f(0) = 1`. This
        method returns :math:`0.2271 f` using an analytical fitting formula
        given in Komatsu et al. 2011, ApJS 192, 18.

        For massless neutrinos (or at high redshift when all neutrinos are
        relativistic), this returns a constant value proportional to Neff.
        For massive neutrinos, this uses the Komatsu fitting formula to
        efficiently approximate the integral above.

        References
        ----------
        .. [1] Komatsu et al. (2011). Seven-Year Wilkinson Microwave Anisotropy
               Probe (WMAP) Observations: Cosmological Interpretation.
               ApJS, 192, 18.
        """
        # Note that there is also a scalar-z-only cython implementation of
        # this in scalar_inv_efuncs.pyx, so if you find a problem in this
        # you need to update there too.

        # The massive and massless contribution must be handled separately
        # But check for common cases first
        z = aszarr(z)
        if not self._nu_info.has_massive_nu:
            return NEUTRINO_FERMI_DIRAC_CORRECTION * self.Neff * np.ones_like(z)

        curr_nu_y = self._nu_info.nu_y / (1.0 + np.expand_dims(z, axis=-1))
        rel_mass_per = (1.0 + (KOMATSU_K * curr_nu_y) ** KOMATSU_P) ** KOMATSU_INVP
        rel_mass = rel_mass_per.sum(-1) + self._nu_info.n_massless_nu

        return NEUTRINO_FERMI_DIRAC_CORRECTION * self._nu_info.neff_per_nu * rel_mass
