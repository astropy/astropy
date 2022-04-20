# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

import warnings
from abc import abstractmethod
from math import exp, floor, log, pi, sqrt
from numbers import Number
from typing import Any, Mapping, Optional, TypeVar

import numpy as np
from numpy import inf, sin

import astropy.constants as const
import astropy.units as u
from astropy.cosmology.core import Cosmology, FlatCosmologyMixin
from astropy.cosmology.parameter import Parameter, _validate_non_negative, _validate_with_unit
from astropy.cosmology.utils import aszarr, vectorize_redshift_method
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.decorators import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning

# isort: split
if HAS_SCIPY:
    from scipy.integrate import quad
else:
    def quad(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'scipy.integrate'")


__all__ = ["FLRW", "FlatFLRWMixin"]

__doctest_requires__ = {'*': ['scipy']}


##############################################################################
# Parameters

# Some conversion constants -- useful to compute them once here and reuse in
# the initialization rather than have every object do them.
_H0units_to_invs = (u.km / (u.s * u.Mpc)).to(1.0 / u.s)
_sec_to_Gyr = u.s.to(u.Gyr)
# const in critical density in cgs units (g cm^-3)
_critdens_const = (3 / (8 * pi * const.G)).cgs.value
# angle conversions
_radian_in_arcsec = (1 * u.rad).to(u.arcsec)
_radian_in_arcmin = (1 * u.rad).to(u.arcmin)
# Radiation parameter over c^2 in cgs (g cm^-3 K^-4)
_a_B_c2 = (4 * const.sigma_sb / const.c ** 3).cgs.value
# Boltzmann constant in eV / K
_kB_evK = const.k_B.to(u.eV / u.K)


# typing
_FLRWT = TypeVar("_FLRWT", bound="FLRW")
_FlatFLRWMixinT = TypeVar("_FlatFLRWMixinT", bound="FlatFLRWMixin")

##############################################################################


class FLRW(Cosmology):
    """
    A class describing an isotropic and homogeneous
    (Friedmann-Lemaitre-Robertson-Walker) cosmology.

    This is an abstract base class -- you cannot instantiate examples of this
    class, but must work with one of its subclasses, such as
    :class:`~astropy.cosmology.LambdaCDM` or :class:`~astropy.cosmology.wCDM`.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0.  If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0. Note that this does not include massive
        neutrinos.

    Ode0 : float
        Omega dark energy: density of dark energy in units of the critical
        density at z=0.

    Tcmb0 : float or scalar quantity-like ['temperature'], optional
        Temperature of the CMB z=0. If a float, must be in [K]. Default: 0 [K].
        Setting this to zero will turn off both photons and neutrinos
        (even massive ones).

    Neff : float, optional
        Effective number of Neutrino species. Default 3.04.

    m_nu : quantity-like ['energy', 'mass'] or array-like, optional
        Mass of each neutrino species in [eV] (mass-energy equivalency enabled).
        If this is a scalar Quantity, then all neutrino species are assumed to
        have that mass. Otherwise, the mass of each species. The actual number
        of neutrino species (and hence the number of elements of m_nu if it is
        not scalar) must be the floor of Neff. Typically this means you should
        provide three neutrino masses unless you are considering something like
        a sterile neutrino.

    Ob0 : float or None, optional
        Omega baryons: density of baryonic matter in units of the critical
        density at z=0.  If this is set to None (the default), any computation
        that requires its value will raise an exception.

    name : str or None (optional, keyword-only)
        Name for this cosmological object.

    meta : mapping or None (optional, keyword-only)
        Metadata for the cosmology, e.g., a reference.

    Notes
    -----
    Class instances are immutable -- you cannot change the parameters' values.
    That is, all of the above attributes (except meta) are read only.

    For details on how to create performant custom subclasses, see the
    documentation on :ref:`astropy-cosmology-fast-integrals`.
    """

    H0 = Parameter(doc="Hubble constant as an `~astropy.units.Quantity` at z=0.",
                   unit="km/(s Mpc)", fvalidate="scalar")
    Om0 = Parameter(doc="Omega matter; matter density/critical density at z=0.",
                    fvalidate="non-negative")
    Ode0 = Parameter(doc="Omega dark energy; dark energy density/critical density at z=0.",
                     fvalidate="float")
    Tcmb0 = Parameter(doc="Temperature of the CMB as `~astropy.units.Quantity` at z=0.",
                      unit="Kelvin", fvalidate="scalar")
    Neff = Parameter(doc="Number of effective neutrino species.", fvalidate="non-negative")
    m_nu = Parameter(doc="Mass of neutrino species.",
                     unit="eV", equivalencies=u.mass_energy())
    Ob0 = Parameter(doc="Omega baryon; baryonic matter density/critical density at z=0.")

    def __init__(self, H0, Om0, Ode0, Tcmb0=0.0*u.K, Neff=3.04, m_nu=0.0*u.eV,
                 Ob0=None, *, name=None, meta=None):
        super().__init__(name=name, meta=meta)

        # Assign (and validate) Parameters
        self.H0 = H0
        self.Om0 = Om0
        self.Ode0 = Ode0
        self.Tcmb0 = Tcmb0
        self.Neff = Neff
        self.m_nu = m_nu  # (reset later, this is just for unit validation)
        self.Ob0 = Ob0  # (must be after Om0)

        # Derived quantities:
        # Dark matter density; matter - baryons, if latter is not None.
        self._Odm0 = None if Ob0 is None else (self._Om0 - self._Ob0)

        # 100 km/s/Mpc * h = H0 (so h is dimensionless)
        self._h = self._H0.value / 100.0
        # Hubble distance
        self._hubble_distance = (const.c / self._H0).to(u.Mpc)
        # H0 in s^-1
        H0_s = self._H0.value * _H0units_to_invs
        # Hubble time
        self._hubble_time = (_sec_to_Gyr / H0_s) << u.Gyr

        # Critical density at z=0 (grams per cubic cm)
        cd0value = _critdens_const * H0_s ** 2
        self._critical_density0 = cd0value << u.g / u.cm ** 3

        # Compute photon density from Tcmb
        self._Ogamma0 = _a_B_c2 * self._Tcmb0.value ** 4 / self._critical_density0.value

        # Compute Neutrino temperature:
        # The constant in front is (4/11)^1/3 -- see any cosmology book for an
        # explanation -- for example, Weinberg 'Cosmology' p 154 eq (3.1.21).
        self._Tnu0 = 0.7137658555036082 * self._Tcmb0

        # Compute neutrino parameters:
        if self._m_nu is None:
            self._nneutrinos = 0
            self._neff_per_nu = None
            self._massivenu = False
            self._massivenu_mass = None
            self._nmassivenu = self._nmasslessnu = None
        else:
            self._nneutrinos = floor(self._Neff)

            # We are going to share Neff between the neutrinos equally. In
            # detail this is not correct, but it is a standard assumption
            # because properly calculating it is a) complicated b) depends on
            # the details of the massive neutrinos (e.g., their weak
            # interactions, which could be unusual if one is considering
            # sterile neutrinos).
            self._neff_per_nu = self._Neff / self._nneutrinos

            # Now figure out if we have massive neutrinos to deal with, and if
            # so, get the right number of masses. It is worth keeping track of
            # massless ones separately (since they are easy to deal with, and a
            # common use case is to have only one massive neutrino).
            massive = np.nonzero(self._m_nu.value > 0)[0]
            self._massivenu = massive.size > 0
            self._nmassivenu = len(massive)
            self._massivenu_mass = self._m_nu[massive].value if self._massivenu else None
            self._nmasslessnu = self._nneutrinos - self._nmassivenu

        # Compute Neutrino Omega and total relativistic component for massive
        # neutrinos. We also store a list version, since that is more efficient
        # to do integrals with (perhaps surprisingly! But small python lists
        # are more efficient than small NumPy arrays).
        if self._massivenu:  # (`_massivenu` set in `m_nu`)
            nu_y = self._massivenu_mass / (_kB_evK * self._Tnu0)
            self._nu_y = nu_y.value
            self._nu_y_list = self._nu_y.tolist()
            self._Onu0 = self._Ogamma0 * self.nu_relative_density(0)
        else:
            # This case is particularly simple, so do it directly The 0.2271...
            # is 7/8 (4/11)^(4/3) -- the temperature bit ^4 (blackbody energy
            # density) times 7/8 for FD vs. BE statistics.
            self._Onu0 = 0.22710731766 * self._Neff * self._Ogamma0
            self._nu_y = self._nu_y_list = None

        # Compute curvature density
        self._Ok0 = 1.0 - self._Om0 - self._Ode0 - self._Ogamma0 - self._Onu0

        # Subclasses should override this reference if they provide
        #  more efficient scalar versions of inv_efunc.
        self._inv_efunc_scalar = self.inv_efunc
        self._inv_efunc_scalar_args = ()

    # ---------------------------------------------------------------
    # Parameter details

    @Ob0.validator
    def Ob0(self, param, value):
        """Validate baryon density to None or positive float > matter density."""
        if value is None:
            return value

        value = _validate_non_negative(self, param, value)
        if value > self.Om0:
            raise ValueError("baryonic density can not be larger than total matter density.")
        return value

    @m_nu.validator
    def m_nu(self, param, value):
        """Validate neutrino masses to right value, units, and shape.

        There are no neutrinos if floor(Neff) or Tcmb0 are 0.
        The number of neutrinos must match floor(Neff).
        Neutrino masses cannot be negative.
        """
        # Check if there are any neutrinos
        if (nneutrinos := floor(self._Neff)) == 0 or self._Tcmb0.value == 0:
            return None  # None, regardless of input

        # Validate / set units
        value = _validate_with_unit(self, param, value)

        # Check values and data shapes
        if value.shape not in ((), (nneutrinos,)):
            raise ValueError("unexpected number of neutrino masses — "
                             f"expected {nneutrinos}, got {len(value)}.")
        elif np.any(value.value < 0):
            raise ValueError("invalid (negative) neutrino mass encountered.")

        # scalar -> array
        if value.isscalar:
            value = np.full_like(value, value, shape=nneutrinos)

        return value

    # ---------------------------------------------------------------
    # properties

    @property
    def is_flat(self):
        """Return bool; `True` if the cosmology is flat."""
        return bool((self._Ok0 == 0.0) and (self.Otot0 == 1.0))

    @property
    def Otot0(self):
        """Omega total; the total density/critical density at z=0."""
        return self._Om0 + self._Ogamma0 + self._Onu0 + self._Ode0 + self._Ok0

    @property
    def Odm0(self):
        """Omega dark matter; dark matter density/critical density at z=0."""
        return self._Odm0

    @property
    def Ok0(self):
        """Omega curvature; the effective curvature density/critical density at z=0."""
        return self._Ok0

    @property
    def Tnu0(self):
        """Temperature of the neutrino background as `~astropy.units.Quantity` at z=0."""
        return self._Tnu0

    @property
    def has_massive_nu(self):
        """Does this cosmology have at least one massive neutrino species?"""
        if self._Tnu0.value == 0:
            return False
        return self._massivenu

    @property
    def h(self):
        """Dimensionless Hubble constant: h = H_0 / 100 [km/sec/Mpc]."""
        return self._h

    @property
    def hubble_time(self):
        """Hubble time as `~astropy.units.Quantity`."""
        return self._hubble_time

    @property
    def hubble_distance(self):
        """Hubble distance as `~astropy.units.Quantity`."""
        return self._hubble_distance

    @property
    def critical_density0(self):
        """Critical density as `~astropy.units.Quantity` at z=0."""
        return self._critical_density0

    @property
    def Ogamma0(self):
        """Omega gamma; the density/critical density of photons at z=0."""
        return self._Ogamma0

    @property
    def Onu0(self):
        """Omega nu; the density/critical density of neutrinos at z=0."""
        return self._Onu0

    # ---------------------------------------------------------------

    @abstractmethod
    def w(self, z):
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
        raise NotImplementedError("w(z) is not implemented")

    def Otot(self, z):
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

    def Om(self, z):
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
        return self._Om0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2

    def Ob(self, z):
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
        if self._Ob0 is None:
            raise ValueError("Baryon density not set for this cosmology")
        z = aszarr(z)
        return self._Ob0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2

    def Odm(self, z):
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
        if self._Odm0 is None:
            raise ValueError("Baryonic density not set for this cosmology, "
                             "unclear meaning of dark matter density")
        z = aszarr(z)
        return self._Odm0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2

    def Ok(self, z):
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
        if self._Ok0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros(z.shape) if hasattr(z, "shape") else 0.0
        return self._Ok0 * (z + 1.0) ** 2 * self.inv_efunc(z) ** 2

    def Ode(self, z):
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
        if self._Ode0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros(z.shape) if hasattr(z, "shape") else 0.0
        return self._Ode0 * self.de_density_scale(z) * self.inv_efunc(z) ** 2

    def Ogamma(self, z):
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
        return self._Ogamma0 * (z + 1.0) ** 4 * self.inv_efunc(z) ** 2

    def Onu(self, z):
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
        if self._Onu0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros(z.shape) if hasattr(z, "shape") else 0.0
        return self.Ogamma(z) * self.nu_relative_density(z)

    def Tcmb(self, z):
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
        return self._Tcmb0 * (aszarr(z) + 1.0)

    def Tnu(self, z):
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
        return self._Tnu0 * (aszarr(z) + 1.0)

    def nu_relative_density(self, z):
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
        # Note that there is also a scalar-z-only cython implementation of
        # this in scalar_inv_efuncs.pyx, so if you find a problem in this
        # you need to update there too.

        # See Komatsu et al. 2011, eq 26 and the surrounding discussion
        # for an explanation of what we are doing here.
        # However, this is modified to handle multiple neutrino masses
        # by computing the above for each mass, then summing
        prefac = 0.22710731766  # 7/8 (4/11)^4/3 -- see any cosmo book

        # The massive and massless contribution must be handled separately
        # But check for common cases first
        z = aszarr(z)
        if not self._massivenu:
            return prefac * self._Neff * (np.ones(z.shape) if hasattr(z, "shape") else 1.0)

        # These are purely fitting constants -- see the Komatsu paper
        p = 1.83
        invp = 0.54644808743  # 1.0 / p
        k = 0.3173

        curr_nu_y = self._nu_y / (1. + np.expand_dims(z, axis=-1))
        rel_mass_per = (1.0 + (k * curr_nu_y) ** p) ** invp
        rel_mass = rel_mass_per.sum(-1) + self._nmasslessnu

        return prefac * self._neff_per_nu * rel_mass

    def _w_integrand(self, ln1pz):
        """Internal convenience function for w(z) integral (eq. 5 of [1]_).

        Parameters
        ----------
        ln1pz : `~numbers.Number` or scalar ndarray
            Assumes scalar input, since this should only be called inside an
            integral.

        References
        ----------
        .. [1] Linder, E. (2003). Exploring the Expansion History of the
               Universe. Phys. Rev. Lett., 90, 091301.
        """
        return 1.0 + self.w(exp(ln1pz) - 1.0)

    def de_density_scale(self, z):
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
        # This allows for an arbitrary w(z) following eq (5) of
        # Linder 2003, PRL 90, 91301.  The code here evaluates
        # the integral numerically.  However, most popular
        # forms of w(z) are designed to make this integral analytic,
        # so it is probably a good idea for subclasses to overload this
        # method if an analytic form is available.
        z = aszarr(z)
        if not isinstance(z, (Number, np.generic)):  # array/Quantity
            ival = np.array([quad(self._w_integrand, 0, log(1 + redshift))[0]
                             for redshift in z])
            return np.exp(3 * ival)
        else:  # scalar
            ival = quad(self._w_integrand, 0, log(z + 1.0))[0]
            return exp(3 * ival)

    def efunc(self, z):
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

        Notes
        -----
        It is not necessary to override this method, but if de_density_scale
        takes a particularly simple form, it may be advantageous to.
        """
        Or = self._Ogamma0 + (self._Onu0 if not self._massivenu
                              else self._Ogamma0 * self.nu_relative_density(z))
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return np.sqrt(zp1 ** 2 * ((Or * zp1 + self._Om0) * zp1 + self._Ok0) +
                       self._Ode0 * self.de_density_scale(z))

    def inv_efunc(self, z):
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
        # Avoid the function overhead by repeating code
        Or = self._Ogamma0 + (self._Onu0 if not self._massivenu
                              else self._Ogamma0 * self.nu_relative_density(z))
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return (zp1 ** 2 * ((Or * zp1 + self._Om0) * zp1 + self._Ok0) +
                self._Ode0 * self.de_density_scale(z))**(-0.5)

    def _lookback_time_integrand_scalar(self, z):
        """Integrand of the lookback time (equation 30 of [1]_).

        Parameters
        ----------
        z : float
            Input redshift.

        Returns
        -------
        I : float
            The integrand for the lookback time.

        References
        ----------
        .. [1] Hogg, D. (1999). Distance measures in cosmology, section 11.
               arXiv e-prints, astro-ph/9905116.
        """
        return self._inv_efunc_scalar(z, *self._inv_efunc_scalar_args) / (z + 1.0)

    def lookback_time_integrand(self, z):
        """Integrand of the lookback time (equation 30 of [1]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        I : float or array
            The integrand for the lookback time.

        References
        ----------
        .. [1] Hogg, D. (1999). Distance measures in cosmology, section 11.
               arXiv e-prints, astro-ph/9905116.
        """
        z = aszarr(z)
        return self.inv_efunc(z) / (z + 1.0)

    def _abs_distance_integrand_scalar(self, z):
        """Integrand of the absorption distance [1]_.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        X : float
            The integrand for the absorption distance.

        References
        ----------
        .. [1] Hogg, D. (1999). Distance measures in cosmology, section 11.
               arXiv e-prints, astro-ph/9905116.
        """
        args = self._inv_efunc_scalar_args
        return (z + 1.0) ** 2 * self._inv_efunc_scalar(z, *args)

    def abs_distance_integrand(self, z):
        """Integrand of the absorption distance [1]_.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        X : float or array
            The integrand for the absorption distance.

        References
        ----------
        .. [1] Hogg, D. (1999). Distance measures in cosmology, section 11.
               arXiv e-prints, astro-ph/9905116.
        """
        z = aszarr(z)
        return (z + 1.0) ** 2 * self.inv_efunc(z)

    def H(self, z):
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
        return self._H0 * self.efunc(z)

    def scale_factor(self, z):
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

    def lookback_time(self, z):
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
        return self._lookback_time(z)

    def _lookback_time(self, z):
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
        """
        return self._hubble_time * self._integral_lookback_time(z)

    @vectorize_redshift_method
    def _integral_lookback_time(self, z, /):
        """Lookback time to redshift ``z``. Value in units of Hubble time.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : float or ndarray
            Lookback time to each input redshift in Hubble time units.
            Returns `float` if input scalar, `~numpy.ndarray` otherwise.
        """
        return quad(self._lookback_time_integrand_scalar, 0, z)[0]

    def lookback_distance(self, z):
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

    def age(self, z):
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
        return self._age(z)

    def _age(self, z):
        """Age of the universe in Gyr at redshift ``z``.

        This internal function exists to be re-defined for optimizations.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            The age of the universe in Gyr at each input redshift.
        """
        return self._hubble_time * self._integral_age(z)

    @vectorize_redshift_method
    def _integral_age(self, z, /):
        """Age of the universe at redshift ``z``. Value in units of Hubble time.

        Calculated using explicit integration.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : float or ndarray
            The age of the universe at each input redshift in Hubble time units.
            Returns `float` if input scalar, `~numpy.ndarray` otherwise.

        See Also
        --------
        z_at_value : Find the redshift corresponding to an age.
        """
        return quad(self._lookback_time_integrand_scalar, z, inf)[0]

    def critical_density(self, z):
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

        return self._critical_density0 * (self.efunc(z)) ** 2

    def comoving_distance(self, z):
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
        return self._comoving_distance_z1z2(0, z)

    def _comoving_distance_z1z2(self, z1, z2):
        """
        Comoving line-of-sight distance in Mpc between objects at redshifts
        ``z1`` and ``z2``.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving distance in Mpc between each input redshift.
        """
        return self._integral_comoving_distance_z1z2(z1, z2)

    @vectorize_redshift_method(nin=2)
    def _integral_comoving_distance_z1z2_scalar(self, z1, z2, /):
        """
        Comoving line-of-sight distance between objects at redshifts ``z1`` and
        ``z2``. Value in Mpc.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts.

        Returns
        -------
        d : float or ndarray
            Comoving distance in Mpc between each input redshift.
            Returns `float` if input scalar, `~numpy.ndarray` otherwise.
        """
        return quad(self._inv_efunc_scalar, z1, z2, args=self._inv_efunc_scalar_args)[0]

    def _integral_comoving_distance_z1z2(self, z1, z2):
        """
        Comoving line-of-sight distance in Mpc between objects at redshifts
        ``z1`` and ``z2``. The comoving distance along the line-of-sight
        between two objects remains constant with time for objects in the
        Hubble flow.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'] or array-like
            Input redshifts.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving distance in Mpc between each input redshift.
        """
        return self._hubble_distance * self._integral_comoving_distance_z1z2_scalar(z1, z2)

    def comoving_transverse_distance(self, z):
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
        return self._comoving_transverse_distance_z1z2(0, z)

    def _comoving_transverse_distance_z1z2(self, z1, z2):
        r"""Comoving transverse distance in Mpc between two redshifts.

        This value is the transverse comoving distance at redshift ``z2`` as
        seen from redshift ``z1`` corresponding to an angular separation of
        1 radian. This is the same as the comoving distance if :math:`\Omega_k`
        is zero (as in the current concordance Lambda-CDM model).

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving transverse distance in Mpc between input redshift.

        Notes
        -----
        This quantity is also called the 'proper motion distance' in some texts.
        """
        Ok0 = self._Ok0
        dc = self._comoving_distance_z1z2(z1, z2)
        if Ok0 == 0:
            return dc
        sqrtOk0 = sqrt(abs(Ok0))
        dh = self._hubble_distance
        if Ok0 > 0:
            return dh / sqrtOk0 * np.sinh(sqrtOk0 * dc.value / dh.value)
        else:
            return dh / sqrtOk0 * sin(sqrtOk0 * dc.value / dh.value)

    def angular_diameter_distance(self, z):
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

    def luminosity_distance(self, z):
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

    def angular_diameter_distance_z1z2(self, z1, z2):
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
        z1, z2 = aszarr(z1), aszarr(z2)
        if np.any(z2 < z1):
            warnings.warn(f"Second redshift(s) z2 ({z2}) is less than first "
                          f"redshift(s) z1 ({z1}).", AstropyUserWarning)
        return self._comoving_transverse_distance_z1z2(z1, z2) / (z2 + 1.0)

    @vectorize_redshift_method
    def absorption_distance(self, z, /):
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
        return quad(self._abs_distance_integrand_scalar, 0, z)[0]

    def distmod(self, z):
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
        # Remember that the luminosity distance is in Mpc
        # Abs is necessary because in certain obscure closed cosmologies
        #  the distance modulus can be negative -- which is okay because
        #  it enters as the square.
        val = 5. * np.log10(abs(self.luminosity_distance(z).value)) + 25.0
        return u.Quantity(val, u.mag)

    def comoving_volume(self, z):
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
        Ok0 = self._Ok0
        if Ok0 == 0:
            return 4.0 / 3.0 * pi * self.comoving_distance(z) ** 3

        dh = self._hubble_distance.value  # .value for speed
        dm = self.comoving_transverse_distance(z).value
        term1 = 4.0 * pi * dh ** 3 / (2.0 * Ok0) * u.Mpc ** 3
        term2 = dm / dh * np.sqrt(1 + Ok0 * (dm / dh) ** 2)
        term3 = sqrt(abs(Ok0)) * dm / dh

        if Ok0 > 0:
            return term1 * (term2 - 1. / sqrt(abs(Ok0)) * np.arcsinh(term3))
        else:
            return term1 * (term2 - 1. / sqrt(abs(Ok0)) * np.arcsin(term3))

    def differential_comoving_volume(self, z):
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
        return self._hubble_distance * (dm ** 2.0) / (self.efunc(z) << u.steradian)

    def kpc_comoving_per_arcmin(self, z):
        """
        Separation in transverse comoving kpc corresponding to an arcminute at
        redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            The distance in comoving kpc corresponding to an arcmin at each
            input redshift.
        """
        return self.comoving_transverse_distance(z).to(u.kpc) / _radian_in_arcmin

    def kpc_proper_per_arcmin(self, z):
        """
        Separation in transverse proper kpc corresponding to an arcminute at
        redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            The distance in proper kpc corresponding to an arcmin at each input
            redshift.
        """
        return self.angular_diameter_distance(z).to(u.kpc) / _radian_in_arcmin

    def arcsec_per_kpc_comoving(self, z):
        """
        Angular separation in arcsec corresponding to a comoving kpc at
        redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        theta : `~astropy.units.Quantity` ['angle']
            The angular separation in arcsec corresponding to a comoving kpc at
            each input redshift.
        """
        return _radian_in_arcsec / self.comoving_transverse_distance(z).to(u.kpc)

    def arcsec_per_kpc_proper(self, z):
        """
        Angular separation in arcsec corresponding to a proper kpc at redshift
        ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        theta : `~astropy.units.Quantity` ['angle']
            The angular separation in arcsec corresponding to a proper kpc at
            each input redshift.
        """
        return _radian_in_arcsec / self.angular_diameter_distance(z).to(u.kpc)


class FlatFLRWMixin(FlatCosmologyMixin):
    """
    Mixin class for flat FLRW cosmologies. Do NOT instantiate directly.
    Must precede the base class in the multiple-inheritance so that this
    mixin's ``__init__`` proceeds the base class'.
    Note that all instances of ``FlatFLRWMixin`` are flat, but not all
    flat cosmologies are instances of ``FlatFLRWMixin``. As example,
    ``LambdaCDM`` **may** be flat (for the a specific set of parameter values),
    but ``FlatLambdaCDM`` **will** be flat.
    """

    Ode0 = FLRW.Ode0.clone(derived=True)  # same as FLRW, but now a derived param.

    def __init_subclass__(cls):
        super().__init_subclass__()
        if "Ode0" in cls._init_signature.parameters:
            raise TypeError("subclasses of `FlatFLRWMixin` cannot have `Ode0` in `__init__`")

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)  # guaranteed not to have `Ode0`
        # Do some twiddling after the fact to get flatness
        self._Ok0 = 0.0
        self._Ode0 = 1.0 - (self._Om0 + self._Ogamma0 + self._Onu0 + self._Ok0)

    @lazyproperty
    def nonflat(self: _FlatFLRWMixinT) -> _FLRWT:
        # Create BoundArgument to handle args versus kwargs.
        # This also handles all errors from mismatched arguments
        ba = self._nonflat_cls_._init_signature.bind_partial(**self._init_arguments,
                                                             Ode0=self.Ode0)
        # Make new instance, respecting args vs kwargs
        inst = self._nonflat_cls_(*ba.args, **ba.kwargs)
        # Because of machine precision, make sure parameters exactly match
        for n in inst.__all_parameters__ + ("Ok0", ):
            setattr(inst, "_" + n, getattr(self, n))

        return inst

    def clone(self, *, meta: Optional[Mapping] = None, to_nonflat: bool = None, **kwargs: Any):
        """Returns a copy of this object with updated parameters, as specified.

        This cannot be used to change the type of the cosmology, except for
        changing to the non-flat version of this cosmology.

        Parameters
        ----------
        meta : mapping or None (optional, keyword-only)
            Metadata that will update the current metadata.
        to_nonflat : bool or None, optional keyword-only
            Whether to change to the non-flat version of this cosmology.
        **kwargs
            Cosmology parameter (and name) modifications. If any parameter is
            changed and a new name is not given, the name will be set to "[old
            name] (modified)".

        Returns
        -------
        newcosmo : `~astropy.cosmology.Cosmology` subclass instance
            A new instance of this class with updated parameters as specified.
            If no arguments are given, then a reference to this object is
            returned instead of copy.

        Examples
        --------
        To make a copy of the ``Planck13`` cosmology with a different matter
        density (``Om0``), and a new name:

            >>> from astropy.cosmology import Planck13
            >>> Planck13.clone(name="Modified Planck 2013", Om0=0.35)
            FlatLambdaCDM(name="Modified Planck 2013", H0=67.77 km / (Mpc s),
              Om0=0.35, ...

        If no name is specified, the new name will note the modification.

            >>> Planck13.clone(Om0=0.35).name
            'Planck13 (modified)'

        The keyword 'to_nonflat' can be used to clone on the non-flat equivalent
        cosmology.

            >>> Planck13.clone(to_nonflat=True)
            LambdaCDM(name="Planck13", ...

            >>> Planck13.clone(H0=70, to_nonflat=True)
            LambdaCDM(name="Planck13 (modified)", H0=70.0 km / (Mpc s), ...

        With 'to_nonflat' `True`, ``Ode0`` can be modified.

            >>> Planck13.clone(to_nonflat=True, Ode0=1)
            LambdaCDM(name="Planck13 (modified)", H0=67.77 km / (Mpc s),
                      Om0=0.30712, Ode0=1.0, ...
        """
        return super().clone(meta=meta, to_nonflat=to_nonflat, **kwargs)

    @property
    def Otot0(self):
        """Omega total; the total density/critical density at z=0."""
        return 1.0

    def Otot(self, z):
        """The total density parameter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts.

        Returns
        -------
        Otot : ndarray or float
            Returns float if input scalar. Value of 1.
        """
        return 1.0 if isinstance(z, (Number, np.generic)) else np.ones_like(z, subok=False)

    def __equiv__(self, other):
        """flat-FLRW equivalence. Use ``.is_equivalent()`` for actual check!

        Parameters
        ----------
        other : `~astropy.cosmology.FLRW` subclass instance
            The object in which to compare.

        Returns
        -------
        bool or `NotImplemented`
            `True` if 'other' is of the same class / non-flat class (e.g.
            ``FlatLambdaCDM`` and ``LambdaCDM``) has matching parameters
            and parameter values. `False` if 'other' is of the same class but
            has different parameters. `NotImplemented` otherwise.
        """
        # check if case (1): same class & parameters
        if isinstance(other, FlatFLRWMixin):
            return super().__equiv__(other)

        # check cases (3, 4), if other is the non-flat version of this class
        # this makes the assumption that any further subclass of a flat cosmo
        # keeps the same physics.
        comparable_classes = [c for c in self.__class__.mro()[1:]
                              if (issubclass(c, FLRW) and c is not FLRW)]
        if other.__class__ not in comparable_classes:
            return NotImplemented

        # check if have equivalent parameters
        # check all parameters in other match those in 'self' and 'other' has
        # no extra parameters (case (2)) except for 'Ode0' and that other
        params_eq = (
            set(self.__all_parameters__) == set(other.__all_parameters__)  # no extra
            and all(np.all(getattr(self, k) == getattr(other, k))  # equal
                    for k in self.__parameters__)
            and other.is_flat
        )

        return params_eq
