# Licensed under a 3-clause BSD style license - see LICENSE.rst
# ruff: noqa: RUF009


import numpy as np
from numpy import sqrt

from astropy.cosmology._src.core import dataclass_decorator
from astropy.cosmology._src.parameter import Parameter
from astropy.cosmology._src.utils import aszarr, deprecated_keywords

from . import scalar_inv_efuncs
from .base import FLRW, FlatFLRWMixin

__all__ = ["FlatwCDM", "wCDM"]

__doctest_requires__ = {"*": ["scipy"]}


@dataclass_decorator
class wCDM(FLRW):
    """FLRW cosmology with a constant dark energy EoS and curvature.

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

    w0: Parameter = Parameter(
        default=-1.0, doc="Dark energy equation of state.", fvalidate="float"
    )

    def __post_init__(self):
        super().__post_init__()

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self.Tcmb0.value == 0:
            inv_efunc_scalar = scalar_inv_efuncs.wcdm_inv_efunc_norel
            inv_efunc_scalar_args = (self.Om0, self.Ode0, self.Ok0, self.w0)
        elif not self._massivenu:
            inv_efunc_scalar = scalar_inv_efuncs.wcdm_inv_efunc_nomnu
            inv_efunc_scalar_args = (
                self.Om0,
                self.Ode0,
                self.Ok0,
                self.Ogamma0 + self.Onu0,
                self.w0,
            )
        else:
            inv_efunc_scalar = scalar_inv_efuncs.wcdm_inv_efunc
            inv_efunc_scalar_args = (
                self.Om0,
                self.Ode0,
                self.Ok0,
                self.Ogamma0,
                self._neff_per_nu,
                self._nmasslessnu,
                self._nu_y_list,
                self.w0,
            )

        object.__setattr__(self, "_inv_efunc_scalar", inv_efunc_scalar)
        object.__setattr__(self, "_inv_efunc_scalar_args", inv_efunc_scalar_args)

    @deprecated_keywords("z", since="7.0")
    def w(self, z):
        r"""Returns dark energy equation of state at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

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
        return self.w0 * (np.ones(z.shape) if hasattr(z, "shape") else 1.0)

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
        and in this case is given by
        :math:`I = \left(1 + z\right)^{3\left(1 + w_0\right)}`
        """
        return (aszarr(z) + 1.0) ** (3.0 * (1.0 + self.w0))

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
        """
        Or = self.Ogamma0 + (
            self.Onu0
            if not self._massivenu
            else self.Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return sqrt(
            zp1**2 * ((Or * zp1 + self.Om0) * zp1 + self.Ok0)
            + self.Ode0 * zp1 ** (3.0 * (1.0 + self.w0))
        )

    @deprecated_keywords("z", since="7.0")
    def inv_efunc(self, z):
        r"""Function used to calculate :math:`\frac{1}{H_z}`.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        E : ndarray or float
            The inverse redshift scaling of the Hubble constant.
            Returns `float` if the input is scalar.
            Defined such that :math:`H_z = H_0 / E`.
        """
        Or = self.Ogamma0 + (
            self.Onu0
            if not self._massivenu
            else self.Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return (
            zp1**2 * ((Or * zp1 + self.Om0) * zp1 + self.Ok0)
            + self.Ode0 * zp1 ** (3.0 * (1.0 + self.w0))
        ) ** (-0.5)


@dataclass_decorator
class FlatwCDM(FlatFLRWMixin, wCDM):
    """FLRW cosmology with a constant dark energy EoS and no spatial curvature.

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

    To get an equivalent cosmology, but of type `astropy.cosmology.wCDM`,
    use :attr:`astropy.cosmology.FlatFLRWMixin.nonflat`.

    >>> print(cosmo.nonflat)
    wCDM(H0=70.0 km / (Mpc s), Om0=0.3, Ode0=0.7, ...
    """

    def __post_init__(self):
        super().__post_init__()

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        if self.Tcmb0.value == 0:
            inv_efunc_scalar = scalar_inv_efuncs.fwcdm_inv_efunc_norel
            inv_efunc_scalar_args = (self.Om0, self.Ode0, self.w0)
        elif not self._massivenu:
            inv_efunc_scalar = scalar_inv_efuncs.fwcdm_inv_efunc_nomnu
            inv_efunc_scalar_args = (
                self.Om0,
                self.Ode0,
                self.Ogamma0 + self.Onu0,
                self.w0,
            )
        else:
            inv_efunc_scalar = scalar_inv_efuncs.fwcdm_inv_efunc
            inv_efunc_scalar_args = (
                self.Om0,
                self.Ode0,
                self.Ogamma0,
                self._neff_per_nu,
                self._nmasslessnu,
                self._nu_y_list,
                self.w0,
            )
        object.__setattr__(self, "_inv_efunc_scalar", inv_efunc_scalar)
        object.__setattr__(self, "_inv_efunc_scalar_args", inv_efunc_scalar_args)

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
        """
        Or = self.Ogamma0 + (
            self.Onu0
            if not self._massivenu
            else self.Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return sqrt(
            zp1**3 * (Or * zp1 + self.Om0) + self.Ode0 * zp1 ** (3.0 * (1 + self.w0))
        )

    @deprecated_keywords("z", since="7.0")
    def inv_efunc(self, z):
        r"""Function used to calculate :math:`\frac{1}{H_z}`.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        E : ndarray or float
            The inverse redshift scaling of the Hubble constant.
            Returns `float` if the input is scalar.
            Defined such that :math:`H(z) = H_0 E(z)`.
        """
        Or = self.Ogamma0 + (
            self.Onu0
            if not self._massivenu
            else self.Ogamma0 * self.nu_relative_density(z)
        )
        zp1 = aszarr(z) + 1.0  # (converts z [unit] -> z [dimensionless])

        return (
            zp1**3 * (Or * zp1 + self.Om0) + self.Ode0 * zp1 ** (3.0 * (1.0 + self.w0))
        ) ** (-0.5)
