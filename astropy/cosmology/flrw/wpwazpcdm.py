# Licensed under a 3-clause BSD style license - see LICENSE.rst
# ruff: noqa: RUF009


from numpy import exp

from astropy.cosmology import units as cu
from astropy.cosmology._utils import aszarr
from astropy.cosmology.core import dataclass_decorator
from astropy.cosmology.parameter import Parameter

from . import scalar_inv_efuncs
from .base import FLRW, FlatFLRWMixin

__all__ = ["wpwaCDM", "FlatwpwaCDM"]

__doctest_requires__ = {"*": ["scipy"]}


@dataclass_decorator
class wpwaCDM(FLRW):
    r"""FLRW cosmology with a CPL dark energy EoS, a pivot redshift, and curvature.

    The equation for the dark energy equation of state (EoS) uses the CPL form as
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

    wp: Parameter = Parameter(
        default=-1.0,
        doc="Dark energy equation of state at the pivot redshift zp.",
        fvalidate="float",
    )
    wa: Parameter = Parameter(
        default=0.0,
        doc="Negative derivative of dark energy equation of state w.r.t. a.",
        fvalidate="float",
    )
    zp: Parameter = Parameter(
        default=0.0 * cu.redshift,
        doc="The pivot redshift, where w(z) = wp.",
        unit=cu.redshift,
    )

    def __post_init__(self):
        super().__post_init__()

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        apiv = 1.0 / (1.0 + self._zp.value)
        if self._Tcmb0.value == 0:
            inv_efunc_scalar = scalar_inv_efuncs.wpwacdm_inv_efunc_norel
            inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ok0,
                self._wp,
                apiv,
                self._wa,
            )
        elif not self._massivenu:
            inv_efunc_scalar = scalar_inv_efuncs.wpwacdm_inv_efunc_nomnu
            inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ok0,
                self._Ogamma0 + self._Onu0,
                self._wp,
                apiv,
                self._wa,
            )
        else:
            inv_efunc_scalar = scalar_inv_efuncs.wpwacdm_inv_efunc
            inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ok0,
                self._Ogamma0,
                self._neff_per_nu,
                self._nmasslessnu,
                self._nu_y_list,
                self._wp,
                apiv,
                self._wa,
            )
        object.__setattr__(self, "_inv_efunc_scalar", inv_efunc_scalar)
        object.__setattr__(self, "_inv_efunc_scalar_args", inv_efunc_scalar_args)

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
        apiv = 1.0 / (1.0 + self._zp.value)
        return zp1 ** (3.0 * (1.0 + self._wp + apiv * self._wa)) * exp(
            -3.0 * self._wa * z / zp1
        )


@dataclass_decorator
class FlatwpwaCDM(FlatFLRWMixin, wpwaCDM):
    r"""FLRW cosmology with a CPL dark energy EoS, a pivot redshift, and no curvature.

    The equation for the dark energy equation of state (EoS) uses the CPL form as
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
    >>> from astropy.cosmology import FlatwpwaCDM
    >>> cosmo = FlatwpwaCDM(H0=70, Om0=0.3, wp=-0.9, wa=0.2, zp=0.4)

    The comoving distance in Mpc at redshift z:

    >>> cosmo.comoving_distance(0.5)
    <Quantity 1868.68474438 Mpc>

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

    def __post_init__(self):
        super().__post_init__()

        # Please see :ref:`astropy-cosmology-fast-integrals` for discussion
        # about what is being done here.
        apiv = 1.0 / (1.0 + self._zp)
        if self._Tcmb0.value == 0:
            inv_efunc_scalar = scalar_inv_efuncs.fwpwacdm_inv_efunc_norel
            inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._wp,
                apiv,
                self._wa,
            )
        elif not self._massivenu:
            inv_efunc_scalar = scalar_inv_efuncs.fwpwacdm_inv_efunc_nomnu
            inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ogamma0 + self._Onu0,
                self._wp,
                apiv,
                self._wa,
            )
        else:
            inv_efunc_scalar = scalar_inv_efuncs.fwpwacdm_inv_efunc
            inv_efunc_scalar_args = (
                self._Om0,
                self._Ode0,
                self._Ogamma0,
                self._neff_per_nu,
                self._nmasslessnu,
                self._nu_y_list,
                self._wp,
                apiv,
                self._wa,
            )
        object.__setattr__(self, "_inv_efunc_scalar", inv_efunc_scalar)
        object.__setattr__(self, "_inv_efunc_scalar_args", inv_efunc_scalar_args)
