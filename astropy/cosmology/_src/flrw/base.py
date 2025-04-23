# Licensed under a 3-clause BSD style license - see LICENSE.rst
# ruff: noqa: RUF009

from __future__ import annotations

__all__ = ["FLRW", "FlatFLRWMixin"]

import inspect
import warnings
from abc import abstractmethod
from dataclasses import field
from functools import cached_property
from inspect import signature
from math import exp, floor, log, pi, sqrt
from numbers import Number
from typing import TYPE_CHECKING, TypeVar, overload

import numpy as np
from numpy import inf, sin

import astropy.constants as const
import astropy.units as u
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.decorators import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning

# isort: split
from astropy.cosmology._src.core import (
    Cosmology,
    FlatCosmologyMixin,
    dataclass_decorator,
)
from astropy.cosmology._src.parameter import (
    Parameter,
    validate_non_negative,
    validate_with_unit,
)
from astropy.cosmology._src.traits import (
    ScaleFactor,
    TemperatureCMB,
    _BaryonComponent,
    _CriticalDensity,
)
from astropy.cosmology._src.utils import (
    aszarr,
    deprecated_keywords,
    vectorize_redshift_method,
)

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Self

    import astropy.units


# isort: split
if HAS_SCIPY:
    from scipy.integrate import quad
else:

    def quad(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'scipy.integrate'")


__doctest_requires__ = {"*": ["scipy"]}
_InputT = TypeVar("_InputT", bound=u.Quantity | np.ndarray | np.generic | Number)


##############################################################################
# Parameters

# Some conversion constants -- useful to compute them once here and reuse in
# the initialization rather than have every object do them.
_H0units_to_invs = (u.km / (u.s * u.Mpc)).to(1.0 / u.s)
_sec_to_Gyr = u.s.to(u.Gyr)
# angle conversions
_radian_in_arcsec = (1 * u.rad).to(u.arcsec)
_radian_in_arcmin = (1 * u.rad).to(u.arcmin)
# Radiation parameter over c^2 in cgs (g cm^-3 K^-4)
_a_B_c2 = (4 * const.sigma_sb / const.c**3).cgs.value
# Boltzmann constant in eV / K
_kB_evK = const.k_B.to(u.eV / u.K)


# typing
_FLRWT = TypeVar("_FLRWT", bound="FLRW")
_FlatFLRWMixinT = TypeVar("_FlatFLRWMixinT", bound="FlatFLRWMixin")

##############################################################################


ParameterOde0 = Parameter(
    doc="Omega dark energy; dark energy density/critical density at z=0.",
    fvalidate="float",
)


@dataclass_decorator
class FLRW(Cosmology, ScaleFactor, TemperatureCMB, _CriticalDensity, _BaryonComponent):
    """An isotropic and homogeneous (Friedmann-Lemaitre-Robertson-Walker) cosmology.

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

    Ob0 : float, optional
        Omega baryons: density of baryonic matter in units of the critical
        density at z=0.

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

    H0: Parameter = Parameter(
        doc="Hubble constant at z=0.",
        unit="km/(s Mpc)",
        fvalidate="scalar",
    )
    Om0: Parameter = Parameter(
        doc="Omega matter; matter density/critical density at z=0.",
        fvalidate="non-negative",
    )
    Ode0: Parameter = ParameterOde0.clone()
    Tcmb0: Parameter = Parameter(
        default=0.0 * u.K,
        doc="Temperature of the CMB at z=0.",
        unit="Kelvin",
        fvalidate="scalar",
    )
    Neff: Parameter = Parameter(
        default=3.04,
        doc="Number of effective neutrino species.",
        fvalidate="non-negative",
    )
    m_nu: Parameter = Parameter(
        default=0.0 * u.eV,
        doc="Mass of neutrino species.",
        unit="eV",
        equivalencies=u.mass_energy(),
    )
    Ob0: Parameter = Parameter(
        default=0.0,
        doc="Omega baryon; baryonic matter density/critical density at z=0.",
    )

    def __post_init__(self):
        # Compute neutrino parameters:
        if self.m_nu is None:
            nneutrinos = 0
            neff_per_nu = None
            massivenu = False
            massivenu_mass = None
            nmassivenu = nmasslessnu = None
        else:
            nneutrinos = floor(self.Neff)

            # We are going to share Neff between the neutrinos equally. In
            # detail this is not correct, but it is a standard assumption
            # because properly calculating it is a) complicated b) depends on
            # the details of the massive neutrinos (e.g., their weak
            # interactions, which could be unusual if one is considering
            # sterile neutrinos).
            neff_per_nu = self.Neff / nneutrinos

            # Now figure out if we have massive neutrinos to deal with, and if
            # so, get the right number of masses. It is worth keeping track of
            # massless ones separately (since they are easy to deal with, and a
            # common use case is to have only one massive neutrino).
            massive = np.nonzero(self.m_nu.value > 0)[0]
            massivenu = massive.size > 0
            nmassivenu = len(massive)
            massivenu_mass = self.m_nu[massive].value if massivenu else None
            nmasslessnu = nneutrinos - nmassivenu

        object.__setattr__(self, "_nneutrinos", nneutrinos)
        object.__setattr__(self, "_neff_per_nu", neff_per_nu)
        object.__setattr__(self, "_massivenu", massivenu)
        object.__setattr__(self, "_massivenu_mass", massivenu_mass)
        object.__setattr__(self, "_nmassivenu", nmassivenu)
        object.__setattr__(self, "_nmasslessnu", nmasslessnu)

        # Compute Neutrino Omega and total relativistic component for massive
        # neutrinos. We also store a list version, since that is more efficient
        # to do integrals with (perhaps surprisingly! But small python lists
        # are more efficient than small NumPy arrays).
        if self._massivenu:  # (`_massivenu` set in `m_nu`)
            nu_y = (self._massivenu_mass / (_kB_evK * self.Tnu0)).value
            nu_y_list = nu_y.tolist()
        else:
            nu_y = nu_y_list = None
        object.__setattr__(self, "_nu_y", nu_y)
        object.__setattr__(self, "_nu_y_list", nu_y_list)

        # Subclasses should override this reference if they provide
        #  more efficient scalar versions of inv_efunc.
        object.__setattr__(self, "_inv_efunc_scalar", self.inv_efunc)
        object.__setattr__(self, "_inv_efunc_scalar_args", ())

    # ---------------------------------------------------------------
    # Parameter details

    @Ob0.validator
    def Ob0(self, param, value):
        """Validate baryon density to a non-negative float > matter density."""
        if value is None:
            warnings.warn(
                "Ob0=None is deprecated, use Ob0=0 instead, "
                "which never causes methods to raise exceptions.",
                category=DeprecationWarning,
                stacklevel=2,
            )
            return 0.0

        value = validate_non_negative(self, param, value)
        if value > self.Om0:
            raise ValueError(
                "baryonic density can not be larger than total matter density."
            )
        return value

    @m_nu.validator
    def m_nu(self, param, value):
        """Validate neutrino masses to right value, units, and shape.

        There are no neutrinos if floor(Neff) or Tcmb0 are 0. The number of
        neutrinos must match floor(Neff). Neutrino masses cannot be
        negative.
        """
        # Check if there are any neutrinos
        if (nneutrinos := floor(self.Neff)) == 0 or self.Tcmb0.value == 0:
            return None  # None, regardless of input

        # Validate / set units
        value = validate_with_unit(self, param, value)

        # Check values and data shapes
        if value.shape not in ((), (nneutrinos,)):
            raise ValueError(
                "unexpected number of neutrino masses — "
                f"expected {nneutrinos}, got {len(value)}."
            )
        elif np.any(value.value < 0):
            raise ValueError("invalid (negative) neutrino mass encountered.")

        # scalar -> array
        if value.isscalar:
            value = np.full_like(value, value, shape=nneutrinos)

        return value

    # ---------------------------------------------------------------
    # properties

    @property
    def is_flat(self) -> bool:
        """Return bool; `True` if the cosmology is flat."""
        return bool((self.Ok0 == 0.0) and (self.Otot0 == 1.0))

    @property
    def Otot0(self) -> float:
        """Omega total; the total density/critical density at z=0."""
        return self.Om0 + self.Ogamma0 + self.Onu0 + self.Ode0 + self.Ok0

    @cached_property
    def Odm0(self) -> float:
        """Omega dark matter; dark matter density/critical density at z=0."""
        return self.Om0 - self.Ob0

    @cached_property
    def Ok0(self) -> float:
        """Omega curvature; the effective curvature density/critical density at z=0."""
        return 1.0 - self.Om0 - self.Ode0 - self.Ogamma0 - self.Onu0

    @cached_property
    def Tnu0(self) -> u.Quantity:
        """Temperature of the neutrino background as |Quantity| at z=0."""
        # The constant in front is (4/11)^1/3 -- see any cosmology book for an
        # explanation -- for example, Weinberg 'Cosmology' p 154 eq (3.1.21).
        return 0.7137658555036082 * self.Tcmb0

    @property
    def has_massive_nu(self) -> bool:
        """Does this cosmology have at least one massive neutrino species?"""
        if self.Tnu0.value == 0:
            return False
        return self._massivenu

    @cached_property
    def h(self) -> float:
        """Dimensionless Hubble constant: h = H_0 / 100 [km/sec/Mpc]."""
        return self.H0.value / 100.0

    @cached_property
    def hubble_time(self) -> u.Quantity:
        """Hubble time."""
        return (_sec_to_Gyr / (self.H0.value * _H0units_to_invs)) << u.Gyr

    @cached_property
    def hubble_distance(self) -> u.Quantity:
        """Hubble distance."""
        return (const.c / self.H0).to(u.Mpc)

    @cached_property
    def critical_density0(self) -> u.Quantity:
        r"""Critical mass density at z=0.

        The critical density is the density of the Universe at which the Universe is
        flat. It is defined as :math:`\rho_{\text{crit}} = 3 H_0^2 / (8 \pi G)`.

        """
        return (3 * self.H0**2 / (8 * pi * const.G)).cgs

    @cached_property
    def Ogamma0(self) -> float:
        """Omega gamma; the density/critical density of photons at z=0."""
        # photon density from Tcmb
        return _a_B_c2 * self.Tcmb0.value**4 / self.critical_density0.value

    @cached_property
    def Onu0(self) -> float:
        """Omega nu; the density/critical density of neutrinos at z=0."""
        if self._massivenu:  # (`_massivenu` set in `m_nu`)
            return self.Ogamma0 * self.nu_relative_density(0)
        else:
            # This case is particularly simple, so do it directly The 0.2271...
            # is 7/8 (4/11)^(4/3) -- the temperature bit ^4 (blackbody energy
            # density) times 7/8 for FD vs. BE statistics.
            return 0.22710731766 * self.Neff * self.Ogamma0

    # ---------------------------------------------------------------

    @abstractmethod
    @deprecated_keywords("z", since="7.0")
    def w(self, z):
        r"""The dark energy equation of state.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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

    @deprecated_keywords("z", since="7.0")
    def Otot(self, z):
        """The total density parameter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshifts.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Otot : ndarray or float
            The total density relative to the critical density at each redshift.
            Returns float if input scalar.
        """
        return self.Om(z) + self.Ogamma(z) + self.Onu(z) + self.Ode(z) + self.Ok(z)

    @deprecated_keywords("z", since="7.0")
    def Om(self, z):
        """Return the density parameter for non-relativistic matter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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

    @deprecated_keywords("z", since="7.0")
    def Odm(self, z):
        """Return the density parameter for dark matter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Odm : ndarray or float
            The density of non-relativistic dark matter relative to the
            critical density at each redshift.
            Returns `float` if the input is scalar.

        Notes
        -----
        This does not include neutrinos, even if non-relativistic at the
        redshift of interest.
        """
        z = aszarr(z)
        return self.Odm0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2

    @deprecated_keywords("z", since="7.0")
    def Ok(self, z):
        """Return the equivalent density parameter for curvature at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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

    @deprecated_keywords("z", since="7.0")
    def Ode(self, z):
        """Return the density parameter for dark energy at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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

    @deprecated_keywords("z", since="7.0")
    def Ogamma(self, z):
        """Return the density parameter for photons at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Ogamma : ndarray or float
            The energy density of photons relative to the critical density at
            each redshift.
            Returns `float` if the input is scalar.
        """
        z = aszarr(z)
        return self.Ogamma0 * (z + 1.0) ** 4 * self.inv_efunc(z) ** 2

    @deprecated_keywords("z", since="7.0")
    def Onu(self, z):
        r"""Return the density parameter for neutrinos at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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

    @deprecated_keywords("z", since="7.0")
    def Tnu(self, z):
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
        """
        return self.Tnu0 * (aszarr(z) + 1.0)

    @deprecated_keywords("z", since="7.0")
    def nu_relative_density(self, z):
        r"""Neutrino density function relative to the energy density in photons.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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
            return (
                prefac * self.Neff * (np.ones(z.shape) if hasattr(z, "shape") else 1.0)
            )

        # These are purely fitting constants -- see the Komatsu paper
        p = 1.83
        invp = 0.54644808743  # 1.0 / p
        k = 0.3173

        curr_nu_y = self._nu_y / (1.0 + np.expand_dims(z, axis=-1))
        rel_mass_per = (1.0 + (k * curr_nu_y) ** p) ** invp
        rel_mass = rel_mass_per.sum(-1) + self._nmasslessnu

        return prefac * self._neff_per_nu * rel_mass

    def _w_integrand(self, ln1pz, /):
        """Internal convenience function for w(z) integral (eq. 5 of [1]_).

        Parameters
        ----------
        ln1pz : `~numbers.Number` or scalar ndarray, positional-only
            Assumes scalar input, since this should only be called inside an
            integral.

            .. versionchanged:: 7.0
                The argument is positional-only.

        References
        ----------
        .. [1] Linder, E. (2003). Exploring the Expansion History of the
               Universe. Phys. Rev. Lett., 90, 091301.
        """
        return 1.0 + self.w(exp(ln1pz) - 1.0)

    @deprecated_keywords("z", since="7.0")
    def de_density_scale(self, z):
        r"""Evaluates the redshift dependence of the dark energy density.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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
            ival = np.array(
                [quad(self._w_integrand, 0, log(1 + redshift))[0] for redshift in z]
            )
            return np.exp(3 * ival)
        else:  # scalar
            ival = quad(self._w_integrand, 0, log(z + 1.0))[0]
            return exp(3 * ival)

    @deprecated_keywords("z", since="7.0")
    def efunc(self, z):
        """Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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
        Or = self.Ogamma0 + (
            self.Onu0
            if not self._massivenu
            else self.Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return np.sqrt(
            zp1**2 * ((Or * zp1 + self.Om0) * zp1 + self.Ok0)
            + self.Ode0 * self.de_density_scale(z)
        )

    @deprecated_keywords("z", since="7.0")
    def inv_efunc(self, z):
        """Inverse of ``efunc``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        E : ndarray or float
            The redshift scaling of the inverse Hubble constant.
            Returns `float` if the input is scalar.
        """
        # Avoid the function overhead by repeating code
        Or = self.Ogamma0 + (
            self.Onu0
            if not self._massivenu
            else self.Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return (
            zp1**2 * ((Or * zp1 + self.Om0) * zp1 + self.Ok0)
            + self.Ode0 * self.de_density_scale(z)
        ) ** (-0.5)

    def _lookback_time_integrand_scalar(self, z, /):
        """Integrand of the lookback time (equation 30 of [1]_).

        Parameters
        ----------
        z : float, positional-only
            Input redshift.

            .. versionchanged:: 7.0
                The argument is positional-only.

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

    @deprecated_keywords("z", since="7.0")
    def lookback_time_integrand(self, z):
        """Integrand of the lookback time (equation 30 of [1]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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

    def _abs_distance_integrand_scalar(self, z, /):
        """Integrand of the absorption distance (eq. 4, [1]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, positional-only
            Input redshift.

            .. versionchanged:: 7.0
                The argument is positional-only.

        Returns
        -------
        dX : float
            The integrand for the absorption distance (dimensionless).

        References
        ----------
        .. [1] Bahcall, John N. and Peebles, P.J.E. 1969, ApJ, 156L, 7B
        """
        return (z + 1.0) ** 2 * self._inv_efunc_scalar(z, *self._inv_efunc_scalar_args)

    @deprecated_keywords("z", since="7.0")
    def abs_distance_integrand(self, z):
        """Integrand of the absorption distance (eq. 4, [1]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        dX : float or array
            The integrand for the absorption distance (dimensionless).

        References
        ----------
        .. [1] Bahcall, John N. and Peebles, P.J.E. 1969, ApJ, 156L, 7B
        """
        z = aszarr(z)
        return (z + 1.0) ** 2 * self.inv_efunc(z)

    @deprecated_keywords("z", since="7.0")
    def H(self, z):
        """Hubble parameter (km/s/Mpc) at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        H : Quantity ['frequency']
            Hubble parameter at each input redshift.
        """
        return self.H0 * self.efunc(z)

    @deprecated_keywords("z", since="7.0")
    def lookback_time(self, z):
        """Lookback time in Gyr to redshift ``z``.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        t : Quantity ['time']
            Lookback time in Gyr to each input redshift.

        See Also
        --------
        z_at_value : Find the redshift corresponding to a lookback time.
        """
        return self._lookback_time(z)

    def _lookback_time(self, z, /):
        """Lookback time in Gyr to redshift ``z``.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, positional-only
            Input redshift.

            .. versionchanged:: 7.0
                The argument is positional-only.

        Returns
        -------
        t : Quantity ['time']
            Lookback time in Gyr to each input redshift.
        """
        return self.hubble_time * self._integral_lookback_time(z)

    @vectorize_redshift_method
    def _integral_lookback_time(self, z, /):
        """Lookback time to redshift ``z``. Value in units of Hubble time.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, positional-only
            Input redshift.

            .. versionchanged:: 7.0
                The argument is positional-only.

        Returns
        -------
        t : float or ndarray
            Lookback time to each input redshift in Hubble time units.
            Returns `float` if input scalar, `~numpy.ndarray` otherwise.
        """
        return quad(self._lookback_time_integrand_scalar, 0, z)[0]

    @deprecated_keywords("z", since="7.0")
    def lookback_distance(self, z):
        """The lookback distance is the light travel time distance to a given redshift.

        It is simply c * lookback_time. It may be used to calculate
        the proper distance between two redshifts, e.g. for the mean free path
        to ionizing radiation.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : Quantity ['length']
            Lookback distance in Mpc
        """
        return (self.lookback_time(z) * const.c).to(u.Mpc)

    @deprecated_keywords("z", since="7.0")
    def age(self, z):
        """Age of the universe in Gyr at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        t : Quantity ['time']
            The age of the universe in Gyr at each input redshift.

        See Also
        --------
        z_at_value : Find the redshift corresponding to an age.
        """
        return self._age(z)

    def _age(self, z, /):
        """Age of the universe in Gyr at redshift ``z``.

        This internal function exists to be re-defined for optimizations.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, positional-only
            Input redshift.

            .. versionchanged:: 7.0
                The argument is positional-only.

        Returns
        -------
        t : Quantity ['time']
            The age of the universe in Gyr at each input redshift.
        """
        return self.hubble_time * self._integral_age(z)

    @vectorize_redshift_method
    def _integral_age(self, z, /):
        """Age of the universe at redshift ``z``. Value in units of Hubble time.

        Calculated using explicit integration.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, positional-only
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

    # ---------------------------------------------------------------
    # Comoving distance

    @overload
    def comoving_distance(self, z: _InputT) -> astropy.units.Quantity: ...

    @overload
    def comoving_distance(self, z: _InputT, z2: _InputT) -> astropy.units.Quantity: ...

    @deprecated_keywords("z2", since="7.1")
    @deprecated_keywords("z", since="7.0")
    def comoving_distance(
        self, z: _InputT, z2: _InputT | None = None
    ) -> astropy.units.Quantity:
        r"""Comoving line-of-sight distance :math:`d_c(z1, z2)` in Mpc.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        Parameters
        ----------
        z, z2 : Quantity ['redshift'], positional-only
            Input redshifts. If one argument ``z`` is given, the distance
            :math:`d_c(0, z)` is returned. If two arguments ``z1, z2`` are
            given, the distance :math:`d_c(z_1, z_2)` is returned.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Quantity ['length']
            Comoving distance in Mpc between each input redshift.
        """
        z1, z2 = (0.0, z) if z2 is None else (z, z2)
        return self._comoving_distance_z1z2(z1, z2)

    def _comoving_distance_z1z2(self, z1, z2, /):
        """Comoving line-of-sight distance in Mpc between redshifts ``z1`` and ``z2``.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, positional-only
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : Quantity ['length']
            Comoving distance in Mpc between each input redshift.
        """
        return self._integral_comoving_distance_z1z2(z1, z2)

    def _integral_comoving_distance_z1z2(self, z1, z2, /):
        """Comoving line-of-sight distance (Mpc) between objects at redshifts z1 and z2.

        The comoving distance along the line-of-sight between two objects remains
        constant with time for objects in the Hubble flow.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'] or array-like
            Input redshifts.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        |Quantity| ['length']
            Comoving distance in Mpc between each input redshift.
        """
        return self.hubble_distance * self._integral_comoving_distance_z1z2_scalar(z1, z2)  # fmt: skip

    @vectorize_redshift_method(nin=2)
    def _integral_comoving_distance_z1z2_scalar(self, z1, z2, /):
        """Comoving line-of-sight distance in Mpc between objects at redshifts ``z1`` and ``z2``.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, positional-only
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : float or ndarray
            Comoving distance in Mpc between each input redshift.
            Returns `float` if input scalar, `~numpy.ndarray` otherwise.
        """
        return quad(self._inv_efunc_scalar, z1, z2, args=self._inv_efunc_scalar_args)[0]

    # ---------------------------------------------------------------

    @deprecated_keywords("z", since="7.0")
    def comoving_transverse_distance(self, z):
        r"""Comoving transverse distance in Mpc at a given redshift.

        This value is the transverse comoving distance at redshift ``z``
        corresponding to an angular separation of 1 radian. This is the same as
        the comoving distance if :math:`\Omega_k` is zero (as in the current
        concordance Lambda-CDM model).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : Quantity ['length']
            Comoving transverse distance in Mpc at each input redshift.

        Notes
        -----
        This quantity is also called the 'proper motion distance' in some texts.
        """
        return self._comoving_transverse_distance_z1z2(0, z)

    def _comoving_transverse_distance_z1z2(self, z1, z2, /):
        r"""Comoving transverse distance in Mpc between two redshifts.

        This value is the transverse comoving distance at redshift ``z2`` as
        seen from redshift ``z1`` corresponding to an angular separation of
        1 radian. This is the same as the comoving distance if :math:`\Omega_k`
        is zero (as in the current concordance Lambda-CDM model).

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, positional-only
            Input redshifts.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : Quantity ['length']
            Comoving transverse distance in Mpc between input redshift.

        Notes
        -----
        This quantity is also called the 'proper motion distance' in some texts.
        """
        Ok0 = self.Ok0
        dc = self._comoving_distance_z1z2(z1, z2)
        if Ok0 == 0:
            return dc
        sqrtOk0 = sqrt(abs(Ok0))
        dh = self.hubble_distance
        if Ok0 > 0:
            return dh / sqrtOk0 * np.sinh(sqrtOk0 * dc.value / dh.value)
        else:
            return dh / sqrtOk0 * sin(sqrtOk0 * dc.value / dh.value)

    @deprecated_keywords("z", since="7.0")
    def angular_diameter_distance(self, z):
        """Angular diameter distance in Mpc at a given redshift.

        This gives the proper (sometimes called 'physical') transverse
        distance corresponding to an angle of 1 radian for an object
        at redshift ``z`` ([1]_, [2]_, [3]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : Quantity ['length']
            Angular diameter distance in Mpc at each input redshift.

        References
        ----------
        .. [1] Weinberg, 1972, pp 420-424; Weedman, 1986, pp 421-424.
        .. [2] Weedman, D. (1986). Quasar astronomy, pp 65-67.
        .. [3] Peebles, P. (1993). Principles of Physical Cosmology, pp 325-327.
        """
        z = aszarr(z)
        return self.comoving_transverse_distance(z) / (z + 1.0)

    @deprecated_keywords("z", since="7.0")
    def luminosity_distance(self, z):
        """Luminosity distance in Mpc at redshift ``z``.

        This is the distance to use when converting between the bolometric flux
        from an object at redshift ``z`` and its bolometric luminosity [1]_.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : Quantity ['length']
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
        z1, z2 : Quantity-like ['redshift'], array-like
            Input redshifts. For most practical applications such as
            gravitational lensing, ``z2`` should be larger than ``z1``. The
            method will work for ``z2 < z1``; however, this will return
            negative distances.

        Returns
        -------
        d : Quantity ['length']
            The angular diameter distance between each input redshift pair.
            Returns scalar if input is scalar, array else-wise.
        """
        z1, z2 = aszarr(z1), aszarr(z2)
        if np.any(z2 < z1):
            warnings.warn(
                f"Second redshift(s) z2 ({z2}) is less than first "
                f"redshift(s) z1 ({z1}).",
                AstropyUserWarning,
            )
        return self._comoving_transverse_distance_z1z2(z1, z2) / (z2 + 1.0)

    @vectorize_redshift_method
    def absorption_distance(self, z, /):
        """Absorption distance at redshift ``z`` (eq. 4, [1]_).

        This is used to calculate the number of objects with some cross section
        of absorption and number density intersecting a sightline per unit
        redshift path [1]_.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, positional-only
            Input redshift.

        Returns
        -------
        X : float or ndarray
            Absorption distance (dimensionless) at each input redshift.
            Returns `float` if input scalar, `~numpy.ndarray` otherwise.

        References
        ----------
        .. [1] Bahcall, John N. and Peebles, P.J.E. 1969, ApJ, 156L, 7B
        """
        return quad(self._abs_distance_integrand_scalar, 0, z)[0]

    @deprecated_keywords("z", since="7.0")
    def distmod(self, z):
        """Distance modulus at redshift ``z``.

        The distance modulus is defined as the (apparent magnitude - absolute
        magnitude) for an object at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        distmod : Quantity ['length']
            Distance modulus at each input redshift, in magnitudes.

        See Also
        --------
        z_at_value : Find the redshift corresponding to a distance modulus.
        """
        # Remember that the luminosity distance is in Mpc
        # Abs is necessary because in certain obscure closed cosmologies
        #  the distance modulus can be negative -- which is okay because
        #  it enters as the square.
        val = 5.0 * np.log10(abs(self.luminosity_distance(z).value)) + 25.0
        return u.Quantity(val, u.mag)

    @deprecated_keywords("z", since="7.0")
    def comoving_volume(self, z):
        r"""Comoving volume in cubic Mpc at redshift ``z``.

        This is the volume of the universe encompassed by redshifts less than
        ``z``. For the case of :math:`\Omega_k = 0` it is a sphere of radius
        `comoving_distance` but it is less intuitive if :math:`\Omega_k` is not.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        V : Quantity ['volume']
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

    @deprecated_keywords("z", since="7.0")
    def differential_comoving_volume(self, z):
        """Differential comoving volume at redshift z.

        Useful for calculating the effective comoving volume.
        For example, allows for integration over a comoving volume that has a
        sensitivity function that changes with redshift. The total comoving
        volume is given by integrating ``differential_comoving_volume`` to
        redshift ``z`` and multiplying by a solid angle.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        dV : Quantity
            Differential comoving volume per redshift per steradian at each
            input redshift.
        """
        dm = self.comoving_transverse_distance(z)
        return self.hubble_distance * (dm**2.0) / (self.efunc(z) << u.steradian)

    @deprecated_keywords("z", since="7.0")
    def kpc_comoving_per_arcmin(self, z):
        """Separation in transverse comoving kpc equal to an arcmin at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : Quantity ['length']
            The distance in comoving kpc corresponding to an arcmin at each
            input redshift.
        """
        return self.comoving_transverse_distance(z).to(u.kpc) / _radian_in_arcmin

    @deprecated_keywords("z", since="7.0")
    def kpc_proper_per_arcmin(self, z):
        """Separation in transverse proper kpc equal to an arcminute at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        d : Quantity ['length']
            The distance in proper kpc corresponding to an arcmin at each input
            redshift.
        """
        return self.angular_diameter_distance(z).to(u.kpc) / _radian_in_arcmin

    @deprecated_keywords("z", since="7.0")
    def arcsec_per_kpc_comoving(self, z):
        """Angular separation in arcsec equal to a comoving kpc at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        theta : Quantity ['angle']
            The angular separation in arcsec corresponding to a comoving kpc at
            each input redshift.
        """
        return _radian_in_arcsec / self.comoving_transverse_distance(z).to(u.kpc)

    @deprecated_keywords("z", since="7.0")
    def arcsec_per_kpc_proper(self, z):
        """Angular separation in arcsec corresponding to a proper kpc at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        theta : Quantity ['angle']
            The angular separation in arcsec corresponding to a proper kpc at
            each input redshift.
        """
        return _radian_in_arcsec / self.angular_diameter_distance(z).to(u.kpc)


@dataclass_decorator
class FlatFLRWMixin(FlatCosmologyMixin):
    """Mixin class for flat FLRW cosmologies.

    Do NOT instantiate directly. Must precede the base class in the
    multiple-inheritance so that this mixin's ``__init__`` proceeds the
    base class'. Note that all instances of ``FlatFLRWMixin`` are flat, but
    not all flat cosmologies are instances of ``FlatFLRWMixin``. As
    example, ``LambdaCDM`` **may** be flat (for the a specific set of
    parameter values), but ``FlatLambdaCDM`` **will** be flat.
    """

    Ode0: Parameter = field(  # now a derived param.
        default=ParameterOde0.clone(default=0, derived=True),
        init=False,
        repr=False,
    )

    def __init_subclass__(cls):
        super().__init_subclass__()

        # Check that Ode0 is not in __init__
        if (
            getattr(
                vars(cls).get("Ode0", cls.__dataclass_fields__.get("Ode0")),
                "init",
                True,
            )
            or "Ode0" in signature(cls.__init__).parameters
        ):
            msg = "subclasses of `FlatFLRWMixin` cannot have `Ode0` in `__init__`"
            raise TypeError(msg)

    def __post_init__(self):
        self.__dict__["Ode0"] = 0
        super().__post_init__()
        # Do some twiddling after the fact to get flatness
        self.__dict__["Ok0"] = 0.0
        self.__dict__["Ode0"] = 1.0 - (self.Om0 + self.Ogamma0 + self.Onu0 + self.Ok0)

    @lazyproperty
    def nonflat(self: _FlatFLRWMixinT) -> _FLRWT:  # noqa: PYI019
        # Create BoundArgument to handle args versus kwargs.
        # This also handles all errors from mismatched arguments
        ba = inspect.signature(self.__nonflatclass__).bind_partial(
            **dict(self.parameters), Ode0=self.Ode0, name=self.name
        )
        # Make new instance, respecting args vs kwargs
        inst = self.__nonflatclass__(*ba.args, **ba.kwargs)
        # Because of machine precision, make sure parameters exactly match
        inst.__dict__.update(inst.parameters)
        inst.__dict__.update(inst._derived_parameters)
        inst.__dict__["Ok0"] = self.Ok0

        return inst

    @property
    def Otot0(self):
        """Omega total; the total density/critical density at z=0."""
        return 1.0

    @deprecated_keywords("z", since="7.0")
    def Otot(self, z):
        """The total density parameter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Otot : ndarray or float
            Returns float if input scalar. Value of 1.
        """
        return (
            1.0 if isinstance(z, (Number, np.generic)) else np.ones_like(z, subok=False)
        )

    def clone(
        self, *, meta: Mapping | None = None, to_nonflat: bool = False, **kwargs
    ) -> Self:
        if not to_nonflat and kwargs.get("Ode0") is not None:
            msg = "Cannot set 'Ode0' in clone unless 'to_nonflat=True'. "
            raise ValueError(msg)
        return super().clone(meta=meta, to_nonflat=to_nonflat, **kwargs)
