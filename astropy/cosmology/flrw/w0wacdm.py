# Licensed under a 3-clause BSD style license - see LICENSE.rst

from numpy import exp

import astropy.units as u
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.utils import aszarr

from . import scalar_inv_efuncs
from .base import FLRW, FlatFLRWMixin

__all__ = ["w0waCDM", "Flatw0waCDM"]

__doctest_requires__ = {'*': ['scipy']}


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
        return zp1 ** (3 * (1 + self._w0 + self._wa)) * exp(-3 * self._wa * z / zp1)


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

    To get an equivalent cosmology, but of type `astropy.cosmology.w0waCDM`,
    use :attr:`astropy.cosmology.FlatFLRWMixin.nonflat`.

    >>> cosmo.nonflat
    w0waCDM(H0=70.0 km / (Mpc s), Om0=0.3, ...

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
