# Licensed under a 3-clause BSD style license - see LICENSE.rst

from math import acos, cos, inf, sin, sqrt
from numbers import Number

import numpy as np
from numpy import log

import astropy.units as u
from astropy.cosmology.utils import aszarr
from astropy.utils.compat.optional_deps import HAS_SCIPY

from . import scalar_inv_efuncs
from .base import FLRW, FlatFLRWMixin

# isort: split
if HAS_SCIPY:
    from scipy.special import ellipkinc, hyp2f1
else:

    def ellipkinc(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'scipy.special'")

    def hyp2f1(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'scipy.special'")


__all__ = ["LambdaCDM", "FlatLambdaCDM"]

__doctest_requires__ = {"*": ["scipy"]}


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

    def __init__(
        self,
        H0,
        Om0,
        Ode0,
        Tcmb0=0.0 * u.K,
        Neff=3.04,
        m_nu=0.0 * u.eV,
        Ob0=None,
        *,
        name=None,
        meta=None
    ):
        super().__init__(
            H0=H0,
            Om0=Om0,
            Ode0=Ode0,
            Tcmb0=Tcmb0,
            Neff=Neff,
            m_nu=m_nu,
            Ob0=Ob0,
            name=name,
            meta=meta,
        )

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
            self._inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ok0,
                self._Ogamma0 + self._Onu0,
            )
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.lcdm_inv_efunc
            self._inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ok0,
                self._Ogamma0,
                self._neff_per_nu,
                self._nmasslessnu,
                self._nu_y_list,
            )

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

        b = -(27.0 / 2) * self._Om0**2 * self._Ode0 / self._Ok0**3
        kappa = b / abs(b)
        if (b < 0) or (2 < b):

            def phi_z(Om0, Ok0, kappa, y1, A, z):
                return np.arccos(
                    ((z + 1.0) * Om0 / abs(Ok0) + kappa * y1 - A)
                    / ((z + 1.0) * Om0 / abs(Ok0) + kappa * y1 + A)
                )

            v_k = pow(kappa * (b - 1) + sqrt(b * (b - 2)), 1.0 / 3)
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
                return np.arcsin(np.sqrt((y1 - y2) / ((z + 1.0) * Om0 / abs(Ok0) + y1)))

            yb = cos(acos(1 - b) / 3)
            yc = sqrt(3) * sin(acos(1 - b) / 3)
            y1 = (1.0 / 3) * (-1 + yb + yc)
            y2 = (1.0 / 3) * (-1 - 2 * yb)
            y3 = (1.0 / 3) * (-1 + yb - yc)
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
        return prefactor * ((z1 + 1.0) ** (-1.0 / 2) - (z2 + 1.0) ** (-1.0 / 2))

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

        s = ((1 - self._Om0) / self._Om0) ** (1.0 / 3)
        # Use np.sqrt here to handle negative s (Om0>1).
        prefactor = self._hubble_distance / np.sqrt(s * self._Om0)
        return prefactor * (
            self._T_hypergeometric(s / (z1 + 1.0))
            - self._T_hypergeometric(s / (z2 + 1.0))
        )

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
        return 2 * np.sqrt(x) * hyp2f1(1.0 / 6, 1.0 / 2, 7.0 / 6, -(x**3))

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
        t = inf if isinstance(z, Number) else np.full_like(z, inf, dtype=float)
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
        return (2.0 / 3) * self._hubble_time * (aszarr(z) + 1.0) ** (-1.5)

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
        prefactor = (2.0 / 3) * self._hubble_time / np.emath.sqrt(1 - self._Om0)
        arg = np.arcsinh(
            np.emath.sqrt((1 / self._Om0 - 1 + 0j) / (aszarr(z) + 1.0) ** 3)
        )
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
        return self._hubble_time * log(aszarr(z) + 1.0)

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
        Or = self._Ogamma0 + (
            self._Onu0
            if not self._massivenu
            else self._Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return np.sqrt(
            zp1**2 * ((Or * zp1 + self._Om0) * zp1 + self._Ok0) + self._Ode0
        )

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
        Or = self._Ogamma0 + (
            self._Onu0
            if not self._massivenu
            else self._Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return (zp1**2 * ((Or * zp1 + self._Om0) * zp1 + self._Ok0) + self._Ode0) ** (
            -0.5
        )


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

    To get an equivalent cosmology, but of type `astropy.cosmology.LambdaCDM`,
    use :attr:`astropy.cosmology.FlatFLRWMixin.nonflat`.

    >>> cosmo.nonflat
    LambdaCDM(H0=70.0 km / (Mpc s), Om0=0.3, ...
    """

    def __init__(
        self,
        H0,
        Om0,
        Tcmb0=0.0 * u.K,
        Neff=3.04,
        m_nu=0.0 * u.eV,
        Ob0=None,
        *,
        name=None,
        meta=None
    ):
        super().__init__(
            H0=H0,
            Om0=Om0,
            Ode0=0.0,
            Tcmb0=Tcmb0,
            Neff=Neff,
            m_nu=m_nu,
            Ob0=Ob0,
            name=name,
            meta=meta,
        )

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
            self._inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ogamma0 + self._Onu0,
            )
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.flcdm_inv_efunc
            self._inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ogamma0,
                self._neff_per_nu,
                self._nmasslessnu,
                self._nu_y_list,
            )

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
        Or = self._Ogamma0 + (
            self._Onu0
            if not self._massivenu
            else self._Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return np.sqrt(zp1**3 * (Or * zp1 + self._Om0) + self._Ode0)

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
        Or = self._Ogamma0 + (
            self._Onu0
            if not self._massivenu
            else self._Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])
        return (zp1**3 * (Or * zp1 + self._Om0) + self._Ode0) ** (-0.5)
