# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
from abc import abstractmethod
from math import acos, cos, exp, floor, inf, log, pi, sin, sqrt
from numbers import Number

import numpy as np

import astropy.constants as const
import astropy.units as u
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.decorators import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning

from . import scalar_inv_efuncs
from . import units as cu
from .core import Cosmology, FlatCosmologyMixin, Parameter
from .parameter import _validate_non_negative, _validate_with_unit
from .utils import aszarr, vectorize_redshift_method

# isort: split
if HAS_SCIPY:
    from scipy.integrate import quad
    from scipy.special import ellipkinc, hyp2f1
else:
    def quad(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'scipy.integrate'")

    def ellipkinc(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'scipy.special'")

    def hyp2f1(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'scipy.special'")


__all__ = ["FLRW", "LambdaCDM", "FlatLambdaCDM", "wCDM", "FlatwCDM",
           "w0waCDM", "Flatw0waCDM", "wpwaCDM", "w0wzCDM", "FlatFLRWMixin"]

__doctest_requires__ = {'*': ['scipy']}


# Some conversion constants -- useful to compute them once here and reuse in
# the initialization rather than have every object do them.
H0units_to_invs = (u.km / (u.s * u.Mpc)).to(1.0 / u.s)
sec_to_Gyr = u.s.to(u.Gyr)
# const in critical density in cgs units (g cm^-3)
critdens_const = (3 / (8 * pi * const.G)).cgs.value
arcsec_in_radians = pi / (3600. * 180)
arcmin_in_radians = pi / (60. * 180)
# Radiation parameter over c^2 in cgs (g cm^-3 K^-4)
a_B_c2 = (4 * const.sigma_sb / const.c ** 3).cgs.value
# Boltzmann constant in eV / K
kB_evK = const.k_B.to(u.eV / u.K)


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
                      unit="Kelvin", fmt="0.4g", fvalidate="scalar")
    Neff = Parameter(doc="Number of effective neutrino species.", fvalidate="non-negative")
    m_nu = Parameter(doc="Mass of neutrino species.",
                     unit="eV", equivalencies=u.mass_energy(), fmt="")
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

        # Derived quantities
        self._Odm0 = None if Ob0 is None else (self._Om0 - self._Ob0)

        # 100 km/s/Mpc * h = H0 (so h is dimensionless)
        self._h = self._H0.value / 100.
        # Hubble distance
        self._hubble_distance = (const.c / self._H0).to(u.Mpc)
        # H0 in s^-1; don't use units for speed
        H0_s = self._H0.value * H0units_to_invs
        # Hubble time; again, avoiding units package for speed
        self._hubble_time = u.Quantity(sec_to_Gyr / H0_s, u.Gyr)

        # critical density at z=0 (grams per cubic cm)
        cd0value = critdens_const * H0_s ** 2
        self._critical_density0 = u.Quantity(cd0value, u.g / u.cm ** 3)

        # -------------------
        # neutrinos

        # Load up neutrino masses.
        self._nneutrinos = floor(self._Neff)

        # We are going to share Neff between the neutrinos equally.
        # In detail this is not correct, but it is a standard assumption
        # because properly calculating it is a) complicated b) depends
        # on the details of the massive neutrinos (e.g., their weak
        # interactions, which could be unusual if one is considering sterile
        # neutrinos)
        self._massivenu = False
        if self._nneutrinos > 0 and self._Tcmb0.value > 0:
            self._neff_per_nu = self._Neff / self._nneutrinos

            # Now, figure out if we have massive neutrinos to deal with,
            # and, if so, get the right number of masses
            # It is worth the effort to keep track of massless ones separately
            # (since they are quite easy to deal with, and a common use case
            # is to set only one neutrino to have mass)
            m_nu = self._m_nu
            if m_nu.isscalar:
                # Assume all neutrinos have the same mass
                if m_nu.value == 0:
                    self._nmasslessnu = self._nneutrinos
                    self._nmassivenu = 0
                else:
                    self._massivenu = True
                    self._nmasslessnu = 0
                    self._nmassivenu = self._nneutrinos
                    self._massivenu_mass = (m_nu.value *
                                            np.ones(self._nneutrinos))
            else:
                # Make sure we have the right number of masses
                # -unless- they are massless, in which case we cheat a little
                if m_nu.value.min() < 0:
                    raise ValueError("Invalid (negative) neutrino mass"
                                     " encountered")
                if m_nu.value.max() == 0:
                    self._nmasslessnu = self._nneutrinos
                    self._nmassivenu = 0
                else:
                    self._massivenu = True
                    if len(m_nu) != self._nneutrinos:
                        errstr = "Unexpected number of neutrino masses"
                        raise ValueError(errstr)
                    # Segregate out the massless ones
                    self._nmasslessnu = len(np.nonzero(m_nu.value == 0)[0])
                    self._nmassivenu = self._nneutrinos - self._nmasslessnu
                    w = np.nonzero(m_nu.value > 0)[0]
                    self._massivenu_mass = m_nu[w]

        # Compute photon density, Tcmb, neutrino parameters
        # Tcmb0=0 removes both photons and neutrinos, is handled
        # as a special case for efficiency
        if self._Tcmb0.value > 0:
            # Compute photon density from Tcmb
            self._Ogamma0 = a_B_c2 * self._Tcmb0.value ** 4 / self._critical_density0.value

            # Compute Neutrino temperature
            # The constant in front is (4/11)^1/3 -- see any
            #  cosmology book for an explanation -- for example,
            #  Weinberg 'Cosmology' p 154 eq (3.1.21)
            self._Tnu0 = 0.7137658555036082 * self._Tcmb0

            # Compute Neutrino Omega and total relativistic component
            # for massive neutrinos.  We also store a list version,
            # since that is more efficient to do integrals with (perhaps
            # surprisingly!  But small python lists are more efficient
            # than small numpy arrays).
            if self._massivenu:
                nu_y = self._massivenu_mass / (kB_evK * self._Tnu0)
                self._nu_y = nu_y.value
                self._nu_y_list = self._nu_y.tolist()
                self._Onu0 = self._Ogamma0 * self.nu_relative_density(0)
            else:
                # This case is particularly simple, so do it directly
                # The 0.2271... is 7/8 (4/11)^(4/3) -- the temperature
                # bit ^4 (blackbody energy density) times 7/8 for
                # FD vs. BE statistics.
                self._Onu0 = 0.22710731766 * self._Neff * self._Ogamma0

        else:
            self._Ogamma0 = 0.0
            self._Tnu0 = 0.0 * u.K
            self._Onu0 = 0.0

        # now set m_nu Parameter
        if self._nneutrinos == 0 or self._Tnu0.value == 0:
            self._m_nu = None
        else:
            if not self._massivenu:  # only massless
                m = np.zeros(self._nmasslessnu)
            elif self._nmasslessnu == 0:  # only massive
                m = self._massivenu_mass
            else:  # a mix -- the most complicated case
                m = np.append(np.zeros(self._nmasslessnu),
                              self._massivenu_mass.value)
            self._m_nu = m << self._m_nu.unit

        # -------------------

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

    # ---------------------------------------------------------------
    # properties

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
        if self._massivenu:
            Or = self._Ogamma0 * (1 + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
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
        if self._massivenu:
            Or = self._Ogamma0 * (1 + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
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
        return quad(self._lookback_time_integrand_scalar, z, np.inf)[0]

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
            return dh / sqrtOk0 * np.sin(sqrtOk0 * dc.value / dh.value)

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
            return 4. / 3. * pi * self.comoving_distance(z) ** 3

        dh = self._hubble_distance.value  # .value for speed
        dm = self.comoving_transverse_distance(z).value
        term1 = 4. * pi * dh ** 3 / (2. * Ok0) * u.Mpc ** 3
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
        return (self.comoving_transverse_distance(z).to(u.kpc) *
                arcmin_in_radians / u.arcmin)

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
        return (self.angular_diameter_distance(z).to(u.kpc) *
                arcmin_in_radians / u.arcmin)

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
        return u.arcsec / (self.comoving_transverse_distance(z).to(u.kpc) *
                           arcsec_in_radians)

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
        return u.arcsec / (self.angular_diameter_distance(z).to(u.kpc) *
                           arcsec_in_radians)


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

    Ode0 = Parameter(doc="Omega dark energy; dark energy density/critical density at z=0.",
                     derived=True)  # no longer a Parameter

    def __init_subclass__(cls):
        super().__init_subclass__()
        if "Ode0" in cls._init_signature.parameters:
            raise TypeError("subclasses of `FlatFLRWMixin` cannot have `Ode0` in `__init__`")

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)  # guaranteed not to have `Ode0`
        # Do some twiddling after the fact to get flatness
        self._Ode0 = 1.0 - self._Om0 - self._Ogamma0 - self._Onu0
        self._Ok0 = 0.0

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
            set(self.__all_parameters__) == set(other.__all_parameters__) # no extra params
            and all(np.all(getattr(self, k) == getattr(other, k))  # params equal
                    for k in self.__all_parameters__)
            # flatness conditions
            and other.Ok0 == 0.0  # `Ode0` is checked in __all_parameters__
        )

        return params_eq


class LambdaCDM(FLRW):
    """FLRW cosmology with a cosmological constant and curvature.

    This has no additional attributes beyond those of FLRW.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0.  If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

    Ode0 : float
        Omega dark energy: density of the cosmological constant in units of
        the critical density at z=0.

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

    Examples
    --------
    >>> from astropy.cosmology import LambdaCDM
    >>> cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, Ode0, Tcmb0=0.0*u.K, Neff=3.04, m_nu=0.0*u.eV,
                 Ob0=None, *, name=None, meta=None):
        super().__init__(H0=H0, Om0=Om0, Ode0=Ode0, Tcmb0=Tcmb0, Neff=Neff,
                         m_nu=m_nu, Ob0=Ob0, name=name, meta=meta)

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.lcdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0)
            if self._Ok0 == 0:
                self._optimize_flat_norad()
            else:
                self._comoving_distance_z1z2 = self._elliptic_comoving_distance_z1z2
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.lcdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0 + self._Onu0)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.lcdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list)

    def _optimize_flat_norad(self):
        """Set optimizations for flat LCDM cosmologies with no radiation."""
        # Call out the Om0=0 (de Sitter) and Om0=1 (Einstein-de Sitter)
        # The dS case is required because the hypergeometric case
        #    for Omega_M=0 would lead to an infinity in its argument.
        # The EdS case is three times faster than the hypergeometric.
        if self._Om0 == 0:
            self._comoving_distance_z1z2 = self._dS_comoving_distance_z1z2
            self._age = self._dS_age
            self._lookback_time = self._dS_lookback_time
        elif self._Om0 == 1:
            self._comoving_distance_z1z2 = self._EdS_comoving_distance_z1z2
            self._age = self._EdS_age
            self._lookback_time = self._EdS_lookback_time
        else:
            self._comoving_distance_z1z2 = self._hypergeometric_comoving_distance_z1z2
            self._age = self._flat_age
            self._lookback_time = self._flat_lookback_time

    def w(self, z):
        r"""Returns dark energy equation of state at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        w : ndarray or float
            The dark energy equation of state.
            Returns `float` if the input is scalar.

        Notes
        -----
        The dark energy equation of state is defined as
        :math:`w(z) = P(z)/\rho(z)`, where :math:`P(z)` is the pressure at
        redshift z and :math:`\rho(z)` is the density at redshift z, both in
        units where c=1. Here this is :math:`w(z) = -1`.
        """
        z = aszarr(z)
        return -1.0 * (np.ones(z.shape) if hasattr(z, "shape") else 1.0)

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
        and in this case is given by :math:`I = 1`.
        """
        z = aszarr(z)
        return np.ones(z.shape) if hasattr(z, "shape") else 1.0

    def _elliptic_comoving_distance_z1z2(self, z1, z2):
        r"""Comoving transverse distance in Mpc between two redshifts.

        This value is the transverse comoving distance at redshift ``z``
        corresponding to an angular separation of 1 radian. This is the same as
        the comoving distance if :math:`\Omega_k` is zero.

        For :math:`\Omega_{rad} = 0` the comoving distance can be directly
        calculated as an elliptic integral [1]_.

        Not valid or appropriate for flat cosmologies (Ok0=0).

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving distance in Mpc between each input redshift.

        References
        ----------
        .. [1] Kantowski, R., Kao, J., & Thomas, R. (2000). Distance-Redshift
               in Inhomogeneous FLRW. arXiv e-prints, astro-ph/0002334.
        """
        try:
            z1, z2 = np.broadcast_arrays(z1, z2)
        except ValueError as e:
            raise ValueError("z1 and z2 have different shapes") from e

        # The analytic solution is not valid for any of Om0, Ode0, Ok0 == 0.
        # Use the explicit integral solution for these cases.
        if self._Om0 == 0 or self._Ode0 == 0 or self._Ok0 == 0:
            return self._integral_comoving_distance_z1z2(z1, z2)

        b = -(27. / 2) * self._Om0**2 * self._Ode0 / self._Ok0**3
        kappa = b / abs(b)
        if (b < 0) or (2 < b):
            def phi_z(Om0, Ok0, kappa, y1, A, z):
                return np.arccos(((z + 1.0) * Om0 / abs(Ok0) + kappa * y1 - A) /
                                 ((z + 1.0) * Om0 / abs(Ok0) + kappa * y1 + A))

            v_k = pow(kappa * (b - 1) + sqrt(b * (b - 2)), 1. / 3)
            y1 = (-1 + kappa * (v_k + 1 / v_k)) / 3
            A = sqrt(y1 * (3 * y1 + 2))
            g = 1 / sqrt(A)
            k2 = (2 * A + kappa * (1 + 3 * y1)) / (4 * A)

            phi_z1 = phi_z(self._Om0, self._Ok0, kappa, y1, A, z1)
            phi_z2 = phi_z(self._Om0, self._Ok0, kappa, y1, A, z2)
        # Get lower-right 0<b<2 solution in Om0, Ode0 plane.
        # Fot the upper-left 0<b<2 solution the Big Bang didn't happen.
        elif (0 < b) and (b < 2) and self._Om0 > self._Ode0:
            def phi_z(Om0, Ok0, y1, y2, z):
                return np.arcsin(np.sqrt((y1 - y2) /
                                         ((z + 1.0) * Om0 / abs(Ok0) + y1)))

            yb = cos(acos(1 - b) / 3)
            yc = sqrt(3) * sin(acos(1 - b) / 3)
            y1 = (1. / 3) * (-1 + yb + yc)
            y2 = (1. / 3) * (-1 - 2 * yb)
            y3 = (1. / 3) * (-1 + yb - yc)
            g = 2 / sqrt(y1 - y2)
            k2 = (y1 - y3) / (y1 - y2)
            phi_z1 = phi_z(self._Om0, self._Ok0, y1, y2, z1)
            phi_z2 = phi_z(self._Om0, self._Ok0, y1, y2, z2)
        else:
            return self._integral_comoving_distance_z1z2(z1, z2)

        prefactor = self._hubble_distance / sqrt(abs(self._Ok0))
        return prefactor * g * (ellipkinc(phi_z1, k2) - ellipkinc(phi_z2, k2))

    def _dS_comoving_distance_z1z2(self, z1, z2):
        r"""
        Comoving line-of-sight distance in Mpc between objects at redshifts
        ``z1`` and ``z2`` in a flat, :math:`\Omega_{\Lambda}=1` cosmology
        (de Sitter).

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        The de Sitter case has an analytic solution.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts. Must be 1D or scalar.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving distance in Mpc between each input redshift.
        """
        try:
            z1, z2 = np.broadcast_arrays(z1, z2)
        except ValueError as e:
            raise ValueError("z1 and z2 have different shapes") from e

        return self._hubble_distance * (z2 - z1)

    def _EdS_comoving_distance_z1z2(self, z1, z2):
        r"""
        Comoving line-of-sight distance in Mpc between objects at redshifts
        ``z1`` and ``z2`` in a flat, :math:`\Omega_M=1` cosmology
        (Einstein - de Sitter).

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        For :math:`\Omega_M=1`, :math:`\Omega_{rad}=0` the comoving distance
        has an analytic solution.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts. Must be 1D or scalar.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving distance in Mpc between each input redshift.
        """
        try:
            z1, z2 = np.broadcast_arrays(z1, z2)
        except ValueError as e:
            raise ValueError("z1 and z2 have different shapes") from e

        prefactor = 2 * self._hubble_distance
        return prefactor * ((z1 + 1.0)**(-1./2) - (z2 + 1.0)**(-1./2))

    def _hypergeometric_comoving_distance_z1z2(self, z1, z2):
        r"""
        Comoving line-of-sight distance in Mpc between objects at redshifts
        ``z1`` and ``z2``.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        For :math:`\Omega_{rad} = 0` the comoving distance can be directly
        calculated as a hypergeometric function [1]_.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshifts.

        Returns
        -------
        d : `~astropy.units.Quantity` ['length']
            Comoving distance in Mpc between each input redshift.

        References
        ----------
        .. [1] Baes, M., Camps, P., & Van De Putte, D. (2017). Analytical
               expressions and numerical evaluation of the luminosity distance
               in a flat cosmology. MNRAS, 468(1), 927-930.
        """
        try:
            z1, z2 = np.broadcast_arrays(z1, z2)
        except ValueError as e:
            raise ValueError("z1 and z2 have different shapes") from e

        s = ((1 - self._Om0) / self._Om0) ** (1./3)
        # Use np.sqrt here to handle negative s (Om0>1).
        prefactor = self._hubble_distance / np.sqrt(s * self._Om0)
        return prefactor * (self._T_hypergeometric(s / (z1 + 1.0)) -
                            self._T_hypergeometric(s / (z2 + 1.0)))

    def _T_hypergeometric(self, x):
        r"""Compute value using Gauss Hypergeometric function 2F1.

        .. math::

           T(x) = 2 \sqrt(x) _{2}F_{1}\left(\frac{1}{6}, \frac{1}{2};
                                            \frac{7}{6}; -x^3 \right)

        Notes
        -----
        The :func:`scipy.special.hyp2f1` code already implements the
        hypergeometric transformation suggested by Baes et al. [1]_ for use in
        actual numerical evaulations.

        References
        ----------
        .. [1] Baes, M., Camps, P., & Van De Putte, D. (2017). Analytical
           expressions and numerical evaluation of the luminosity distance
           in a flat cosmology. MNRAS, 468(1), 927-930.
        """
        return 2 * np.sqrt(x) * hyp2f1(1./6, 1./2, 7./6, -x**3)

    def _dS_age(self, z):
        """Age of the universe in Gyr at redshift ``z``.

        The age of a de Sitter Universe is infinite.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            The age of the universe in Gyr at each input redshift.
        """
        t = (inf if isinstance(z, Number) else np.full_like(z, inf, dtype=float))
        return self._hubble_time * t

    def _EdS_age(self, z):
        r"""Age of the universe in Gyr at redshift ``z``.

        For :math:`\Omega_{rad} = 0` (:math:`T_{CMB} = 0`; massless neutrinos)
        the age can be directly calculated as an elliptic integral [1]_.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            The age of the universe in Gyr at each input redshift.

        References
        ----------
        .. [1] Thomas, R., & Kantowski, R. (2000). Age-redshift relation for
               standard cosmology. PRD, 62(10), 103507.
        """
        return (2./3) * self._hubble_time * (aszarr(z) + 1.0) ** (-1.5)

    def _flat_age(self, z):
        r"""Age of the universe in Gyr at redshift ``z``.

        For :math:`\Omega_{rad} = 0` (:math:`T_{CMB} = 0`; massless neutrinos)
        the age can be directly calculated as an elliptic integral [1]_.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            The age of the universe in Gyr at each input redshift.

        References
        ----------
        .. [1] Thomas, R., & Kantowski, R. (2000). Age-redshift relation for
               standard cosmology. PRD, 62(10), 103507.
        """
        # Use np.sqrt, np.arcsinh instead of math.sqrt, math.asinh
        # to handle properly the complex numbers for 1 - Om0 < 0
        prefactor = (2./3) * self._hubble_time / np.emath.sqrt(1 - self._Om0)
        arg = np.arcsinh(np.emath.sqrt((1 / self._Om0 - 1 + 0j) / (aszarr(z) + 1.0)**3))
        return (prefactor * arg).real

    def _EdS_lookback_time(self, z):
        r"""Lookback time in Gyr to redshift ``z``.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        For :math:`\Omega_{rad} = 0` (:math:`T_{CMB} = 0`; massless neutrinos)
        the age can be directly calculated as an elliptic integral.
        The lookback time is here calculated based on the ``age(0) - age(z)``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            Lookback time in Gyr to each input redshift.
        """
        return self._EdS_age(0) - self._EdS_age(z)

    def _dS_lookback_time(self, z):
        r"""Lookback time in Gyr to redshift ``z``.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        For :math:`\Omega_{rad} = 0` (:math:`T_{CMB} = 0`; massless neutrinos)
        the age can be directly calculated.

        .. math::

           a = exp(H * t) \  \text{where t=0 at z=0}

           t = (1/H) (ln 1 - ln a) = (1/H) (0 - ln (1/(1+z))) = (1/H) ln(1+z)

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            Lookback time in Gyr to each input redshift.
        """
        return self._hubble_time * np.log(aszarr(z) + 1.0)

    def _flat_lookback_time(self, z):
        r"""Lookback time in Gyr to redshift ``z``.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        For :math:`\Omega_{rad} = 0` (:math:`T_{CMB} = 0`; massless neutrinos)
        the age can be directly calculated.
        The lookback time is here calculated based on the ``age(0) - age(z)``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        t : `~astropy.units.Quantity` ['time']
            Lookback time in Gyr to each input redshift.
        """
        return self._flat_age(0) - self._flat_age(z)

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
        """
        # We override this because it takes a particularly simple
        # form for a cosmological constant
        if self._massivenu:
            Or = self._Ogamma0 * (1. + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return np.sqrt(zp1 ** 2 * ((Or * zp1 + self._Om0) * zp1 + self._Ok0) + self._Ode0)

    def inv_efunc(self, z):
        r"""Function used to calculate :math:`\frac{1}{H_z}`.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        E : ndarray or float
            The inverse redshift scaling of the Hubble constant.
            Returns `float` if the input is scalar.
            Defined such that :math:`H_z = H_0 / E`.
        """
        if self._massivenu:
            Or = self._Ogamma0 * (1 + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return (zp1 ** 2 * ((Or * zp1 + self._Om0) * zp1 + self._Ok0) + self._Ode0)**(-0.5)


class FlatLambdaCDM(FlatFLRWMixin, LambdaCDM):
    """FLRW cosmology with a cosmological constant and no curvature.

    This has no additional attributes beyond those of FLRW.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

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

    Examples
    --------
    >>> from astropy.cosmology import FlatLambdaCDM
    >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, Tcmb0=0.0*u.K, Neff=3.04, m_nu=0.0*u.eV,
                 Ob0=None, *, name=None, meta=None):
        super().__init__(H0=H0, Om0=Om0, Ode0=0.0, Tcmb0=Tcmb0, Neff=Neff,
                         m_nu=m_nu, Ob0=Ob0, name=name, meta=meta)

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.flcdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0)
            # Repeat the optimization reassignments here because the init
            # of the LambaCDM above didn't actually create a flat cosmology.
            # That was done through the explicit tweak setting self._Ok0.
            self._optimize_flat_norad()
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.flcdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0 + self._Onu0)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.flcdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list)

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
        """
        # We override this because it takes a particularly simple
        # form for a cosmological constant
        if self._massivenu:
            Or = self._Ogamma0 * (1 + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return np.sqrt(zp1 ** 3 * (Or * zp1 + self._Om0) + self._Ode0)

    def inv_efunc(self, z):
        r"""Function used to calculate :math:`\frac{1}{H_z}`.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        E : ndarray or float
            The inverse redshift scaling of the Hubble constant.
            Returns `float` if the input is scalar.
            Defined such that :math:`H_z = H_0 / E`.
        """
        if self._massivenu:
            Or = self._Ogamma0 * (1. + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])
        return (zp1 ** 3 * (Or * zp1 + self._Om0) + self._Ode0)**(-0.5)


class wCDM(FLRW):
    """
    FLRW cosmology with a constant dark energy equation of state and curvature.

    This has one additional attribute beyond those of FLRW.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

    Ode0 : float
        Omega dark energy: density of dark energy in units of the critical
        density at z=0.

    w0 : float, optional
        Dark energy equation of state at all redshifts. This is
        pressure/density for dark energy in units where c=1. A cosmological
        constant has w0=-1.0.

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

    Examples
    --------
    >>> from astropy.cosmology import wCDM
    >>> cosmo = wCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)
    """

    w0 = Parameter(doc="Dark energy equation of state.", fvalidate="float")

    def __init__(self, H0, Om0, Ode0, w0=-1.0, Tcmb0=0.0*u.K, Neff=3.04,
                 m_nu=0.0*u.eV, Ob0=None, *, name=None, meta=None):
        super().__init__(H0=H0, Om0=Om0, Ode0=Ode0, Tcmb0=Tcmb0, Neff=Neff,
                         m_nu=m_nu, Ob0=Ob0, name=name, meta=meta)
        self.w0 = w0

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.wcdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._w0)
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.wcdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0 + self._Onu0,
                                           self._w0)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.wcdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list, self._w0)

    def w(self, z):
        r"""Returns dark energy equation of state at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        w : ndarray or float
            The dark energy equation of state
            Returns `float` if the input is scalar.

        Notes
        -----
        The dark energy equation of state is defined as
        :math:`w(z) = P(z)/\rho(z)`, where :math:`P(z)` is the pressure at
        redshift z and :math:`\rho(z)` is the density at redshift z, both in
        units where c=1. Here this is :math:`w(z) = w_0`.
        """
        z = aszarr(z)
        return self._w0 * (np.ones(z.shape) if hasattr(z, "shape") else 1.0)

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
        and in this case is given by
        :math:`I = \left(1 + z\right)^{3\left(1 + w_0\right)}`
        """
        return (aszarr(z) + 1.0) ** (3.0 * (1. + self._w0))

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
        """
        if self._massivenu:
            Or = self._Ogamma0 * (1. + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return np.sqrt(zp1 ** 2 * ((Or * zp1 + self._Om0) * zp1 + self._Ok0) +
                       self._Ode0 * zp1 ** (3. * (1. + self._w0)))

    def inv_efunc(self, z):
        r"""Function used to calculate :math:`\frac{1}{H_z}`.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        E : ndarray or float
            The inverse redshift scaling of the Hubble constant.
            Returns `float` if the input is scalar.
            Defined such that :math:`H_z = H_0 / E`.
        """
        if self._massivenu:
            Or = self._Ogamma0 * (1. + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return (zp1 ** 2 * ((Or * zp1 + self._Om0) * zp1 + self._Ok0) +
                self._Ode0 * zp1 ** (3. * (1. + self._w0)))**(-0.5)


class FlatwCDM(FlatFLRWMixin, wCDM):
    """
    FLRW cosmology with a constant dark energy equation of state and no spatial
    curvature.

    This has one additional attribute beyond those of FLRW.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

    w0 : float, optional
        Dark energy equation of state at all redshifts. This is
        pressure/density for dark energy in units where c=1. A cosmological
        constant has w0=-1.0.

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

    Examples
    --------
    >>> from astropy.cosmology import FlatwCDM
    >>> cosmo = FlatwCDM(H0=70, Om0=0.3, w0=-0.9)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, w0=-1.0, Tcmb0=0.0*u.K, Neff=3.04, m_nu=0.0*u.eV,
                 Ob0=None, *, name=None, meta=None):
        super().__init__(H0=H0, Om0=Om0, Ode0=0.0, w0=w0, Tcmb0=Tcmb0,
                         Neff=Neff, m_nu=m_nu, Ob0=Ob0, name=name, meta=meta)

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.fwcdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._w0)
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.fwcdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0 + self._Onu0,
                                           self._w0)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.fwcdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list, self._w0)

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
        """
        if self._massivenu:
            Or = self._Ogamma0 * (1. + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return np.sqrt(zp1 ** 3 * (Or * zp1 + self._Om0) +
                       self._Ode0 * zp1 ** (3. * (1 + self._w0)))

    def inv_efunc(self, z):
        r"""Function used to calculate :math:`\frac{1}{H_z}`.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        E : ndarray or float
            The inverse redshift scaling of the Hubble constant.
            Returns `float` if the input is scalar.
            Defined such that :math:`H(z) = H_0 E(z)`.
        """
        if self._massivenu:
            Or = self._Ogamma0 * (1. + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return (zp1 ** 3 * (Or * zp1 + self._Om0) +
                self._Ode0 * zp1 ** (3. * (1. + self._w0)))**(-0.5)


class w0waCDM(FLRW):
    r"""FLRW cosmology with a CPL dark energy equation of state and curvature.

    The equation for the dark energy equation of state uses the
    CPL form as described in Chevallier & Polarski [1]_ and Linder [2]_:
    :math:`w(z) = w_0 + w_a (1-a) = w_0 + w_a z / (1+z)`.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

    Ode0 : float
        Omega dark energy: density of dark energy in units of the critical
        density at z=0.

    w0 : float, optional
        Dark energy equation of state at z=0 (a=1). This is pressure/density
        for dark energy in units where c=1.

    wa : float, optional
        Negative derivative of the dark energy equation of state with respect
        to the scale factor. A cosmological constant has w0=-1.0 and wa=0.0.

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

    Examples
    --------
    >>> from astropy.cosmology import w0waCDM
    >>> cosmo = w0waCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9, wa=0.2)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)

    References
    ----------
    .. [1] Chevallier, M., & Polarski, D. (2001). Accelerating Universes with
           Scaling Dark Matter. International Journal of Modern Physics D,
           10(2), 213-223.
    .. [2] Linder, E. (2003). Exploring the Expansion History of the
           Universe. Phys. Rev. Lett., 90, 091301.
    """

    w0 = Parameter(doc="Dark energy equation of state at z=0.", fvalidate="float")
    wa = Parameter(doc="Negative derivative of dark energy equation of state w.r.t. a.",
                   fvalidate="float")

    def __init__(self, H0, Om0, Ode0, w0=-1.0, wa=0.0, Tcmb0=0.0*u.K, Neff=3.04,
                 m_nu=0.0*u.eV, Ob0=None, *, name=None, meta=None):
        super().__init__(H0=H0, Om0=Om0, Ode0=Ode0, Tcmb0=Tcmb0, Neff=Neff,
                         m_nu=m_nu, Ob0=Ob0, name=name, meta=meta)
        self.w0 = w0
        self.wa = wa

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wacdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._w0, self._wa)
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wacdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0 + self._Onu0,
                                           self._w0, self._wa)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wacdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list, self._w0,
                                           self._wa)

    def w(self, z):
        r"""Returns dark energy equation of state at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        w : ndarray or float
            The dark energy equation of state
            Returns `float` if the input is scalar.

        Notes
        -----
        The dark energy equation of state is defined as
        :math:`w(z) = P(z)/\rho(z)`, where :math:`P(z)` is the pressure at
        redshift z and :math:`\rho(z)` is the density at redshift z, both in
        units where c=1. Here this is
        :math:`w(z) = w_0 + w_a (1 - a) = w_0 + w_a \frac{z}{1+z}`.
        """
        z = aszarr(z)
        return self._w0 + self._wa * z / (z + 1.0)

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
        and in this case is given by

        .. math::

           I = \left(1 + z\right)^{3 \left(1 + w_0 + w_a\right)}
                     \exp \left(-3 w_a \frac{z}{1+z}\right)
        """
        z = aszarr(z)
        zp1 = z + 1.0  # (converts z [unit] -> z [dimensionless])
        return zp1 ** (3 * (1 + self._w0 + self._wa)) * np.exp(-3 * self._wa * z / zp1)


class Flatw0waCDM(FlatFLRWMixin, w0waCDM):
    """FLRW cosmology with a CPL dark energy equation of state and no
    curvature.

    The equation for the dark energy equation of state uses the CPL form as
    described in Chevallier & Polarski [1]_ and Linder [2]_:
    :math:`w(z) = w_0 + w_a (1-a) = w_0 + w_a z / (1+z)`.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

    w0 : float, optional
        Dark energy equation of state at z=0 (a=1). This is pressure/density
        for dark energy in units where c=1.

    wa : float, optional
        Negative derivative of the dark energy equation of state with respect
        to the scale factor. A cosmological constant has w0=-1.0 and wa=0.0.

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

    Examples
    --------
    >>> from astropy.cosmology import Flatw0waCDM
    >>> cosmo = Flatw0waCDM(H0=70, Om0=0.3, w0=-0.9, wa=0.2)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)

    References
    ----------
    .. [1] Chevallier, M., & Polarski, D. (2001). Accelerating Universes with
           Scaling Dark Matter. International Journal of Modern Physics D,
           10(2), 213-223.
    .. [2] Linder, E. (2003). Exploring the Expansion History of the
           Universe. Phys. Rev. Lett., 90, 091301.
    """

    def __init__(self, H0, Om0, w0=-1.0, wa=0.0, Tcmb0=0.0*u.K, Neff=3.04,
                 m_nu=0.0*u.eV, Ob0=None, *, name=None, meta=None):
        super().__init__(H0=H0, Om0=Om0, Ode0=0.0, w0=w0, wa=wa, Tcmb0=Tcmb0,
                         Neff=Neff, m_nu=m_nu, Ob0=Ob0, name=name, meta=meta)

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.fw0wacdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._w0, self._wa)
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.fw0wacdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0 + self._Onu0,
                                           self._w0, self._wa)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.fw0wacdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list, self._w0,
                                           self._wa)


class wpwaCDM(FLRW):
    r"""
    FLRW cosmology with a CPL dark energy equation of state, a pivot redshift,
    and curvature.

    The equation for the dark energy equation of state uses the CPL form as
    described in Chevallier & Polarski [1]_ and Linder [2]_, but modified to
    have a pivot redshift as in the findings of the Dark Energy Task Force
    [3]_: :math:`w(a) = w_p + w_a (a_p - a) = w_p + w_a( 1/(1+zp) - 1/(1+z) )`.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

    Ode0 : float
        Omega dark energy: density of dark energy in units of the critical
        density at z=0.

    wp : float, optional
        Dark energy equation of state at the pivot redshift zp. This is
        pressure/density for dark energy in units where c=1.

    wa : float, optional
        Negative derivative of the dark energy equation of state with respect
        to the scale factor. A cosmological constant has wp=-1.0 and wa=0.0.

    zp : float or quantity-like ['redshift'], optional
        Pivot redshift -- the redshift where w(z) = wp

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

    Examples
    --------
    >>> from astropy.cosmology import wpwaCDM
    >>> cosmo = wpwaCDM(H0=70, Om0=0.3, Ode0=0.7, wp=-0.9, wa=0.2, zp=0.4)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)

    References
    ----------
    .. [1] Chevallier, M., & Polarski, D. (2001). Accelerating Universes with
           Scaling Dark Matter. International Journal of Modern Physics D,
           10(2), 213-223.
    .. [2] Linder, E. (2003). Exploring the Expansion History of the
           Universe. Phys. Rev. Lett., 90, 091301.
    .. [3] Albrecht, A., Amendola, L., Bernstein, G., Clowe, D., Eisenstein,
           D., Guzzo, L., Hirata, C., Huterer, D., Kirshner, R., Kolb, E., &
           Nichol, R. (2009). Findings of the Joint Dark Energy Mission Figure
           of Merit Science Working Group. arXiv e-prints, arXiv:0901.0721.
    """

    wp = Parameter(doc="Dark energy equation of state at the pivot redshift zp.", fvalidate="float")
    wa = Parameter(doc="Negative derivative of dark energy equation of state w.r.t. a.",
                   fvalidate="float")
    zp = Parameter(doc="The pivot redshift, where w(z) = wp.", unit=cu.redshift)

    def __init__(self, H0, Om0, Ode0, wp=-1.0, wa=0.0, zp=0.0 * cu.redshift,
                 Tcmb0=0.0*u.K, Neff=3.04, m_nu=0.0*u.eV, Ob0=None, *,
                 name=None, meta=None):
        super().__init__(H0=H0, Om0=Om0, Ode0=Ode0, Tcmb0=Tcmb0, Neff=Neff,
                         m_nu=m_nu, Ob0=Ob0, name=name, meta=meta)
        self.wp = wp
        self.wa = wa
        self.zp = zp

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        apiv = 1.0 / (1.0 + self._zp.value)
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.wpwacdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._wp, apiv, self._wa)
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.wpwacdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0 + self._Onu0,
                                           self._wp, apiv, self._wa)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.wpwacdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list, self._wp,
                                           apiv, self._wa)

    def w(self, z):
        r"""Returns dark energy equation of state at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        w : ndarray or float
            The dark energy equation of state
            Returns `float` if the input is scalar.

        Notes
        -----
        The dark energy equation of state is defined as
        :math:`w(z) = P(z)/\rho(z)`, where :math:`P(z)` is the pressure at
        redshift z and :math:`\rho(z)` is the density at redshift z, both in
        units where c=1. Here this is :math:`w(z) = w_p + w_a (a_p - a)` where
        :math:`a = 1/1+z` and :math:`a_p = 1 / 1 + z_p`.
        """
        apiv = 1.0 / (1.0 + self._zp.value)
        return self._wp + self._wa * (apiv - 1.0 / (aszarr(z) + 1.0))

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
        and in this case is given by

        .. math::

           a_p = \frac{1}{1 + z_p}

           I = \left(1 + z\right)^{3 \left(1 + w_p + a_p w_a\right)}
                     \exp \left(-3 w_a \frac{z}{1+z}\right)
        """
        z = aszarr(z)
        zp1 = z + 1.0  # (converts z [unit] -> z [dimensionless])
        apiv = 1. / (1. + self._zp.value)
        return zp1 ** (3. * (1. + self._wp + apiv * self._wa)) * \
            np.exp(-3. * self._wa * z / zp1)


class w0wzCDM(FLRW):
    """
    FLRW cosmology with a variable dark energy equation of state and curvature.

    The equation for the dark energy equation of state uses the simple form:
    :math:`w(z) = w_0 + w_z z`.

    This form is not recommended for z > 1.

    Parameters
    ----------
    H0 : float or scalar quantity-like ['frequency']
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc].

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

    Ode0 : float
        Omega dark energy: density of dark energy in units of the critical
        density at z=0.

    w0 : float, optional
        Dark energy equation of state at z=0. This is pressure/density for
        dark energy in units where c=1.

    wz : float, optional
        Derivative of the dark energy equation of state with respect to z.
        A cosmological constant has w0=-1.0 and wz=0.0.

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

    Examples
    --------
    >>> from astropy.cosmology import w0wzCDM
    >>> cosmo = w0wzCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9, wz=0.2)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)
    """

    w0 = Parameter(doc="Dark energy equation of state at z=0.", fvalidate="float")
    wz = Parameter(doc="Derivative of the dark energy equation of state w.r.t. z.", fvalidate="float")

    def __init__(self, H0, Om0, Ode0, w0=-1.0, wz=0.0, Tcmb0=0.0*u.K, Neff=3.04,
                 m_nu=0.0*u.eV, Ob0=None, *, name=None, meta=None):
        super().__init__(H0=H0, Om0=Om0, Ode0=Ode0, Tcmb0=Tcmb0, Neff=Neff,
                         m_nu=m_nu, Ob0=Ob0, name=name, meta=meta)
        self.w0 = w0
        self.wz = wz

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wzcdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._w0, self._wz)
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wzcdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0 + self._Onu0,
                                           self._w0, self._wz)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wzcdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list, self._w0,
                                           self._wz)

    def w(self, z):
        r"""Returns dark energy equation of state at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like, or `~numbers.Number`
            Input redshift.

        Returns
        -------
        w : ndarray or float
            The dark energy equation of state.
            Returns `float` if the input is scalar.

        Notes
        -----
        The dark energy equation of state is defined as
        :math:`w(z) = P(z)/\rho(z)`, where :math:`P(z)` is the pressure at
        redshift z and :math:`\rho(z)` is the density at redshift z, both in
        units where c=1. Here this is given by :math:`w(z) = w_0 + w_z z`.
        """
        return self._w0 + self._wz * aszarr(z)

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
        and in this case is given by

        .. math::

           I = \left(1 + z\right)^{3 \left(1 + w_0 - w_z\right)}
                     \exp \left(-3 w_z z\right)
        """
        z = aszarr(z)
        zp1 = z + 1.0  # (converts z [unit] -> z [dimensionless])
        return zp1 ** (3. * (1. + self._w0 - self._wz)) * np.exp(-3. * self._wz * z)
