# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

from math import pi, sqrt
from typing import Any, Protocol, runtime_checkable

import numpy as np

import astropy.constants as const
import astropy.units as u
from astropy.cosmology.core import CosmologyAPI
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.utils import aszarr

__all__ = ["FLRWAPI"]

__doctest_requires__ = {"*": ["scipy"]}

#####################ar#########################################################
# Parameters

# Some conversion constants -- useful to compute them once here and reuse in
# the initialization rather than have every object do them.
_H0units = u.km / u.s / u.Mpc
_critdens_const = (3 / (8 * pi * const.G)).to(u.g * u.s**2 / u.cm**3)


##############################################################################


@runtime_checkable
class FLRWAPI(CosmologyAPI, Protocol):
    """API for :class:`~astropy.cosmology.FLRW`."""

    H0: u.Quantity | Parameter
    Om0: Any | Parameter
    Ode0: Any | Parameter
    Tcmb0: u.Quantity | Parameter
    Neff: Any | Parameter
    m_nu: u.Quantity | Parameter
    Ob0: Any | Parameter

    # --- Derived parameters ---

    @property
    def is_flat(self) -> bool:
        return bool((self.Ok0 == 0.0) and (self.Otot0 == 1.0))

    @property
    def Otot0(self) -> Any:
        """Omega total; the total density/critical density at z=0."""
        return self.Om0 + self.Ogamma0 + self.Onu0 + self.Ode0 + self.Ok0

    @property
    def Odm0(self) -> Any:
        """Omega dark matter; dark matter density/critical density at z=0."""
        ...

    @property
    def Ok0(self) -> Any:
        """Omega curvature; the effective curvature density/critical density at z=0."""
        ...

    @property
    def Tnu0(self) -> u.Quantity:
        """Temperature of the neutrino background at z=0."""
        ...

    @property
    def has_massive_nu(self) -> bool:
        """Whether there are (at least one) massive neutrino species."""
        ...

    @property
    def h(self) -> Any:
        """Dimensionless Hubble constant: h = H_0 / 100 [km/sec/Mpc]."""
        return self.H0.to_value(_H0units) / 100.0

    @property
    def hubble_time(self) -> u.Quantity:
        """Hubble time."""
        return 1 / self.H0.to(1 / u.Gyr)

    @property
    def hubble_distance(self) -> u.Quantity:
        """Hubble distance."""
        return (const.c / self.H0).to(u.Mpc)

    @property
    def critical_density0(self) -> u.Quantity:
        """Critical density at z = 0."""
        return (_critdens_const * self.H0) << u.g / u.cm**3

    @property
    def Ogamma0(self) -> Any:
        """Omega gamma; the density/critical density of photons at z=0."""
        ...

    @property
    def Onu0(self) -> Any:
        """Omega nu; the density/critical density of neutrinos at z=0."""
        ...

    # --- Methods ---

    def w(self, z, /):
        r"""The dark energy equation of state.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        w : ndarray or float
            The dark energy equation of state.
            `float` if scalar input.

        Notes
        -----
        The dark energy equation of state is defined as
        :math:`w(z) = P(z)/\rho(z)`, where :math:`P(z)` is the pressure at
        redshift z and :math:`\rho(z)` is the density at redshift z, both in
        units where c=1.

        This must be overridden by subclasses.
        """
        ...

    def Otot(self, z, /):
        """The total density parameter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts.

        Returns
        -------
        Otot : ndarray or float
            The total density relative to the critical density at each redshift.
            Returns float if input scalar.
        """
        return self.Om(z) + self.Ogamma(z) + self.Onu(z) + self.Ode(z) + self.Ok(z)

    def Om(self, z, /):
        """
        Return the density parameter for non-relativistic matter
        at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        Om : ndarray or float
            The density of non-relativistic matter relative to the critical
            density at each redshift.
            Returns `float` if the input is scalar.

        Notes
        -----
        This does not include neutrinos, even if non-relativistic at the
        redshift of interest; see `Onu`.
        """
        z = aszarr(z)
        return self.Om0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2

    def Ob(self, z, /):
        """Return the density parameter for baryonic matter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        Ob : ndarray or float
            The density of baryonic matter relative to the critical density at
            each redshift.
            Returns `float` if the input is scalar.

        Raises
        ------
        ValueError
            If ``Ob0`` is `None`.
        """
        if self.Ob0 is None:
            raise ValueError("Baryon density not set for this cosmology")
        z = aszarr(z)
        return self.Ob0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2

    def Odm(self, z, /):
        """Return the density parameter for dark matter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        Odm : ndarray or float
            The density of non-relativistic dark matter relative to the
            critical density at each redshift.
            Returns `float` if the input is scalar.

        Raises
        ------
        ValueError
            If ``Ob0`` is `None`.

        Notes
        -----
        This does not include neutrinos, even if non-relativistic at the
        redshift of interest.
        """
        if self.Odm0 is None:
            raise ValueError(
                "Baryonic density not set for this cosmology, "
                "unclear meaning of dark matter density"
            )
        z = aszarr(z)
        return self.Odm0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2

    def Ok(self, z, /):
        """
        Return the equivalent density parameter for curvature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

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

    def Ode(self, z, /):
        """Return the density parameter for dark energy at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        Ode : ndarray or float
            The density of non-relativistic matter relative to the critical
            density at each redshift.
            Returns `float` if the input is scalar.
        """

        z = aszarr(z)
        if self.Ode0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros(z.shape) if hasattr(z, "shape") else 0.0
        return self.Ode0 * self.de_density_scale(z) * self.inv_efunc(z) ** 2

    def Ogamma(self, z, /):
        """Return the density parameter for photons at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        Ogamma : ndarray or float
            The energy density of photons relative to the critical density at
            each redshift.
            Returns `float` if the input is scalar.
        """
        z = aszarr(z)
        return self.Ogamma0 * (z + 1.0) ** 4 * self.inv_efunc(z) ** 2

    def Onu(self, z, /):

        r"""Return the density parameter for neutrinos at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        Onu : ndarray or float
            The energy density of neutrinos relative to the critical density at
            each redshift. Note that this includes their kinetic energy (if
            they have mass), so it is not equal to the commonly used
            :math:`\sum \frac{m_{\nu}}{94 eV}`, which does not include
            kinetic energy.
            Returns `float` if the input is scalar.
        """
        z = aszarr(z)
        if self.Onu0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros(z.shape) if hasattr(z, "shape") else 0.0
        return self.Ogamma(z) * self.nu_relative_density(z)

    def Tcmb(self, z, /) -> u.Quantity:
        """Return the CMB temperature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        Tcmb : `~astropy.units.Quantity` ['temperature']
            The temperature of the CMB in K.
        """
        return self.Tcmb0 * (aszarr(z) + 1.0)

    def Tnu(self, z, /) -> u.Quantity:
        """Return the neutrino temperature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        Tnu : `~astropy.units.Quantity` ['temperature']
            The temperature of the cosmic neutrino background in K.
        """
        return self.Tnu0 * (aszarr(z) + 1.0)

    def nu_relative_density(self, z, /):
        r"""Neutrino density function relative to the energy density in photons.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        f : ndarray or float
            The neutrino density scaling factor relative to the density in
            photons at each redshift.
            Only returns `float` if z is scalar.

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
        """
        ...

    def de_density_scale(self, z, /):
        r"""Evaluates the redshift dependence of the dark energy density.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        I : ndarray or float
            The scaling of the energy density of dark energy with redshift.
            Returns `float` if the input is scalar.

        Notes
        -----
        The scaling factor, I, is defined by :math:`\rho(z) = \rho_0 I`,
        and is given by

        .. math::

           I = \exp \left( 3 \int_{a}^1 \frac{ da^{\prime} }{ a^{\prime} }
                          \left[ 1 + w\left( a^{\prime} \right) \right] \right)

        The actual integral used is rewritten from [1]_ to be in terms of z.

        It will generally helpful for subclasses to overload this method if
        the integral can be done analytically for the particular dark
        energy equation of state that they implement.

        References
        ----------
        .. [1] Linder, E. (2003). Exploring the Expansion History of the
               Universe. Phys. Rev. Lett., 90, 091301.
        """
        ...

    def efunc(self, z, /):
        """Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        E : ndarray or float
            The redshift scaling of the Hubble constant.
            Returns `float` if the input is scalar.
            Defined such that :math:`H(z) = H_0 E(z)`.
        """
        ...

    def inv_efunc(self, z, /):
        """Inverse of ``efunc``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        E : ndarray or float
            The redshift scaling of the inverse Hubble constant.
            Returns `float` if the input is scalar.
        """
        ...

    def H(self, z, /) -> u.Quantity:
        """Hubble parameter (km/s/Mpc) at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        H : `~astropy.units.Quantity` ['frequency']
            Hubble parameter at each input redshift.
        """
        return self.H0 * self.efunc(z)

    def scale_factor(self, z, /) -> u.Quantity:
        """Scale factor at redshift ``z``.

        The scale factor is defined as :math:`a = 1 / (1 + z)`.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        a : ndarray or float
            Scale factor at each input redshift.
            Returns `float` if the input is scalar.
        """
        return 1.0 / (aszarr(z) + 1.0)

    def lookback_time(self, z, /) -> u.Quantity:
        """Lookback time in Gyr to redshift ``z``.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            Lookback time in Gyr to each input redshift.

        See Also
        --------
        z_at_value : Find the redshift corresponding to a lookback time.
        """
        ...

    def lookback_distance(self, z, /) -> u.Quantity:
        """
        The lookback distance is the light travel time distance to a given
        redshift. It is simply c * lookback_time. It may be used to calculate
        the proper distance between two redshifts, e.g. for the mean free path
        to ionizing radiation.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Lookback distance in Mpc
        """
        return (self.lookback_time(z) * const.c).to(u.Mpc)

    def age(self, z, /) -> u.Quantity:

        """Age of the universe in Gyr at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            The age of the universe in Gyr at each input redshift.

        See Also
        --------
        z_at_value : Find the redshift corresponding to an age.
        """
        ...

    def critical_density(self, z, /) -> u.Quantity:
        """Critical density in grams per cubic cm at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        rho : `~astropy.units.Quantity`
            Critical density in g/cm^3 at each input redshift.
        """
        return self.critical_density0 * self.efunc(z) ** 2

    def comoving_distance(self, z, /) -> u.Quantity:
        """Comoving line-of-sight distance in Mpc at a given redshift.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving distance in Mpc to each input redshift.
        """
        ...

    def comoving_transverse_distance(self, z, /) -> u.Quantity:
        r"""Comoving transverse distance in Mpc at a given redshift.

        This value is the transverse comoving distance at redshift ``z``
        corresponding to an angular separation of 1 radian. This is the same as
        the comoving distance if :math:`\Omega_k` is zero (as in the current
        concordance Lambda-CDM model).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving transverse distance in Mpc at each input redshift.

        Notes
        -----
        This quantity is also called the 'proper motion distance' in some texts.
        """
        Ok0 = self.Ok0
        dc = self.comoving_distance(z)
        if Ok0 == 0:
            return dc
        sqrtOk0 = sqrt(abs(Ok0))
        dh = self.hubble_distance
        if Ok0 > 0:
            return dh / sqrtOk0 * np.sinh(sqrtOk0 * dc.value / dh.value)
        else:
            return dh / sqrtOk0 * np.sin(sqrtOk0 * dc.value / dh.value)

    def angular_diameter_distance(self, z, /) -> u.Quantity:
        """Angular diameter distance in Mpc at a given redshift.

        This gives the proper (sometimes called 'physical') transverse
        distance corresponding to an angle of 1 radian for an object
        at redshift ``z`` ([1]_, [2]_, [3]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Angular diameter distance in Mpc at each input redshift.

        References
        ----------
        .. [1] Weinberg, 1972, pp 420-424; Weedman, 1986, pp 421-424.
        .. [2] Weedman, D. (1986). Quasar astronomy, pp 65-67.
        .. [3] Peebles, P. (1993). Principles of Physical Cosmology, pp 325-327.
        """
        z = aszarr(z)
        return self.comoving_transverse_distance(z) / (z + 1.0)

    def luminosity_distance(self, z, /) -> u.Quantity:
        """Luminosity distance in Mpc at redshift ``z``.

        This is the distance to use when converting between the bolometric flux
        from an object at redshift ``z`` and its bolometric luminosity [1]_.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Luminosity distance in Mpc at each input redshift.

        See Also
        --------
        z_at_value : Find the redshift corresponding to a luminosity distance.

        References
        ----------
        .. [1] Weinberg, 1972, pp 420-424; Weedman, 1986, pp 60-62.
        """
        z = aszarr(z)
        return (z + 1.0) * self.comoving_transverse_distance(z)

    def angular_diameter_distance_z1z2(self, z1, z2, /) -> u.Quantity:
        """Angular diameter distance between objects at 2 redshifts.

        Useful for gravitational lensing, for example computing the angular
        diameter distance between a lensed galaxy and the foreground lens.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts. For most practical applications such as
            gravitational lensing, ``z2`` should be larger than ``z1``. The
            method will work for ``z2 < z1``; however, this will return
            negative distances.

        Returns
        -------
        d : `~astropy.units.Quantity`
            The angular diameter distance between each input redshift pair.
            Returns scalar if input is scalar, array else-wise.
        """
        ...

    def absorption_distance(self, z, /) -> u.Quantity:
        """Absorption distance at redshift ``z``.

        This is used to calculate the number of objects with some cross section
        of absorption and number density intersecting a sightline per unit
        redshift path ([1]_, [2]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        d : float or ndarray
            Absorption distance (dimensionless) at each input redshift.
            Returns `float` if input scalar, `~numpy.ndarray` otherwise.

        References
        ----------
        .. [1] Hogg, D. (1999). Distance measures in cosmology, section 11.
               arXiv e-prints, astro-ph/9905116.
        .. [2] Bahcall, John N. and Peebles, P.J.E. 1969, ApJ, 156L, 7B
        """
        ...

    def distmod(self, z, /) -> u.Quantity:
        """Distance modulus at redshift ``z``.

        The distance modulus is defined as the (apparent magnitude - absolute
        magnitude) for an object at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        distmod : `~astropy.units.Quantity` ['length']
            Distance modulus at each input redshift, in magnitudes.

        See Also
        --------
        z_at_value : Find the redshift corresponding to a distance modulus.
        """
        # Remember that the luminosity distance is in Mpc Abs is necessary because in
        # certain obscure closed cosmologies the distance modulus can be negative --
        # which is okay because it enters as the square.
        val = 5.0 * np.log10(abs(self.luminosity_distance(z).value)) + 25.0
        return u.Quantity(val, u.mag)

    def comoving_volume(self, z, /) -> u.Quantity:
        r"""Comoving volume in cubic Mpc at redshift ``z``.

        This is the volume of the universe encompassed by redshifts less than
        ``z``. For the case of :math:`\Omega_k = 0` it is a sphere of radius
        `comoving_distance` but it is less intuitive if :math:`\Omega_k` is not.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        V : `~astropy.units.Quantity`
            Comoving volume in :math:`Mpc^3` at each input redshift.
        """
        Ok0 = self.Ok0
        if Ok0 == 0:
            return 4.0 / 3.0 * pi * self.comoving_distance(z) ** 3

        dh = self.hubble_distance.value  # .value for speed
        dm = self.comoving_transverse_distance(z).value
        term1 = 4.0 * pi * dh**3 / (2.0 * Ok0) * u.Mpc**3
        term2 = dm / dh * np.sqrt(1 + Ok0 * (dm / dh) ** 2)
        term3 = sqrt(abs(Ok0)) * dm / dh

        if Ok0 > 0:
            return term1 * (term2 - 1.0 / sqrt(abs(Ok0)) * np.arcsinh(term3))
        else:
            return term1 * (term2 - 1.0 / sqrt(abs(Ok0)) * np.arcsin(term3))

    def differential_comoving_volume(self, z, /) -> u.Quantity:
        """Differential comoving volume at redshift z.

        Useful for calculating the effective comoving volume.
        For example, allows for integration over a comoving volume that has a
        sensitivity function that changes with redshift. The total comoving
        volume is given by integrating ``differential_comoving_volume`` to
        redshift ``z`` and multiplying by a solid angle.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        dV : `~astropy.units.Quantity`
            Differential comoving volume per redshift per steradian at each
            input redshift.
        """
        dm = self.comoving_transverse_distance(z)
        return self.hubble_distance * dm**2.0 / (self.efunc(z) << u.steradian)

    # TODO: deprecate? trivially derived from the above.
    # - kpc_comoving_per_arcmin
    # - kpc_proper_per_arcmin
    # - arcsec_per_kpc_comoving
    # - arcsec_per_kpc_proper
