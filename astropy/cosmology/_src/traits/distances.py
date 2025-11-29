"""Distance measures."""

__all__ = ("DistanceMeasures",)

import warnings
from abc import abstractmethod
from numbers import Number
from typing import Final, TypeVar, overload

import numpy as np
from numpy.typing import ArrayLike

import astropy.constants as const
import astropy.units as u
from astropy.cosmology._src.scipy_compat import quad
from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr, vectorize_redshift_method
from astropy.utils.exceptions import AstropyUserWarning

from .hubble import HubbleParameter

_InputT = TypeVar("_InputT", bound=u.Quantity | np.ndarray | np.generic | Number)

# Some conversion constants -- useful to compute them once here and reuse in the
# initialization rather than have every object do them.

# angle conversions
RAD_IN_ARCSEC: Final = (1 * u.rad).to(u.arcsec)
RAD_IN_ARCMIN: Final = (1 * u.rad).to(u.arcmin)


class DistanceMeasures(HubbleParameter):
    """Distance measures for cosmological models."""

    # ===========================================
    # Lookback Time

    def lookback_time_integrand(self, z: u.Quantity | ArrayLike, /) -> FArray:
        """Integrand of the lookback time (equation 30 of [1]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        I : array
            The integrand for the lookback time.

        References
        ----------
        .. [1] Hogg, D. (1999). Distance measures in cosmology, section 11.
               arXiv e-prints, astro-ph/9905116.
        """
        z = aszarr(z)
        return self.inv_efunc(z) / (z + 1.0)

    @vectorize_redshift_method
    def lookback_time(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """Lookback time in Gyr to redshift ``z``.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        t : Quantity ['time']
            Lookback time in Gyr to each input redshift.

        See Also
        --------
        z_at_value : Find the redshift corresponding to a lookback time.
        """
        # NOTE: for performance, this method is overwritten in the FLRW subclass.
        return self.hubble_time * quad(self.lookback_time_integrand, 0, z)[0]

    # ===========================================
    # Lookback Distance

    def lookback_distance(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """The lookback distance is the light travel time distance to a given redshift.

        It is simply c * lookback_time. It may be used to calculate
        the proper distance between two redshifts, e.g. for the mean free path
        to ionizing radiation.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        d : Quantity ['length']
            Lookback distance in Mpc
        """
        return (self.lookback_time(z) * const.c).to(u.Mpc)

    # ===========================================
    # Age

    @vectorize_redshift_method
    def age(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """Age of the universe in Gyr at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        t : Quantity ['time']
            The age of the universe in Gyr at each input redshift.

        See Also
        --------
        z_at_value : Find the redshift corresponding to an age.
        """
        return self.hubble_time * quad(self.lookback_time_integrand, z, np.inf)[0]

    # ===========================================
    # Absorption Distance

    def abs_distance_integrand(self, z: u.Quantity | ArrayLike, /) -> FArray:
        """Integrand of the absorption distance (eq. 4, [1]_).

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        dX : array
            The integrand for the absorption distance (dimensionless).

        References
        ----------
        .. [1] Bahcall, John N. and Peebles, P.J.E. 1969, ApJ, 156L, 7B
        """
        z = aszarr(z)
        return (z + 1.0) ** 2 * self.inv_efunc(z)

    @vectorize_redshift_method
    def absorption_distance(self, z: u.Quantity | ArrayLike, /) -> FArray:
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
        X : array
            Absorption distance (dimensionless) at each input redshift.

        References
        ----------
        .. [1] Bahcall, John N. and Peebles, P.J.E. 1969, ApJ, 156L, 7B
        """
        return quad(self.abs_distance_integrand, 0, z)[0]

    # ===========================================
    # Comoving Distance

    @abstractmethod
    def comoving_distance(self, z: _InputT, z2: _InputT | None = None, /) -> u.Quantity:
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
        raise NotImplementedError  # pragma: no cover

    @abstractmethod
    def comoving_transverse_distance(
        self, z: _InputT, z2: _InputT | None = None, /
    ) -> u.Quantity:
        r"""Comoving transverse distance in Mpc at a given redshift.

        This value is the transverse comoving distance at redshift ``z``
        corresponding to an angular separation of 1 radian. This is the same as
        the comoving distance if :math:`\Omega_k` is zero (as in the current
        concordance Lambda-CDM model).

        Parameters
        ----------
        z, z2 : Quantity ['redshift'], positional-only
            Input redshifts. If one argument ``z`` is given, the distance
            :math:`d_M(0, z)` is returned. If two arguments ``z1, z2`` are
            given, the distance :math:`d_M(z_1, z_2)` is returned.

        Returns
        -------
        d : Quantity ['length']
            Comoving transverse distance in Mpc at each input redshift.

        Notes
        -----
        This quantity is also called the 'proper motion distance' in some texts.
        """
        raise NotImplementedError  # pragma: no cover

    # ===========================================
    # Angular Diameter Distance

    @overload
    def angular_diameter_distance(self, z: _InputT, /) -> u.Quantity: ...

    @overload
    def angular_diameter_distance(self, z: _InputT, z2: _InputT, /) -> u.Quantity: ...

    def angular_diameter_distance(
        self, z: _InputT, z2: _InputT | None = None, /
    ) -> u.Quantity:
        """Angular diameter distance between objects at 2 redshifts.

        When one redshift is given, this gives the proper (sometimes called 'physical')
        transverse distance corresponding to an angle of 1 radian for an object at
        redshift ``z`` ([1]_, [2]_, [3]_).

        The two redshift form is useful for e.g. gravitational lensing for
        computing the angular diameter distance between a lensed galaxy and the
        foreground lens.

        Parameters
        ----------
        z1, z2 : Quantity-like ['redshift'], array-like
            Input redshifts. For most practical applications such as gravitational
            lensing, ``z2`` should be larger than ``z1``. The method will work for ``z2
            < z1``; however, this will return negative distances.

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
        z1, z2 = (0.0, z) if z2 is None else (z, z2)
        z1, z2 = aszarr(z1), aszarr(z2)
        if np.any(z2 < z1):
            warnings.warn(
                f"Second redshift(s) z2 ({z2}) is less than first "
                f"redshift(s) z1 ({z1}).",
                AstropyUserWarning,
            )
        return self.comoving_transverse_distance(z1, z2) / (z2 + 1.0)

    # ===========================================
    # Luminosity Distance

    def luminosity_distance(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """Luminosity distance in Mpc at redshift ``z``.

        This is the distance to use when converting between the bolometric flux
        from an object at redshift ``z`` and its bolometric luminosity [1]_.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

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

    # ===========================================
    # Comoving Volume

    @abstractmethod
    def comoving_volume(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        r"""Comoving volume in cubic Mpc at redshift ``z``.

        This is the volume of the universe encompassed by redshifts less than
        ``z``. For the case of :math:`\Omega_k = 0` it is a sphere of radius
        `comoving_distance` but it is less intuitive if :math:`\Omega_k` is not.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        V : Quantity ['volume']
            Comoving volume in :math:`Mpc^3` at each input redshift.
        """
        raise NotImplementedError  # pragma: no cover

    def differential_comoving_volume(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
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

        Returns
        -------
        dV : Quantity
            Differential comoving volume per redshift per steradian at each
            input redshift.
        """
        dm = self.comoving_transverse_distance(z)
        return self.hubble_distance * (dm**2.0) / (self.efunc(z) << u.steradian)

    # ===========================================
    # Distance Modulus

    def distmod(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """Distance modulus at redshift ``z``.

        The distance modulus is defined as the (apparent magnitude - absolute
        magnitude) for an object at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

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

    # ===========================================
    # Angular Scale Conversions

    def kpc_comoving_per_arcmin(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """Separation in transverse comoving kpc equal to an arcmin at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        d : Quantity ['length']
            The distance in comoving kpc corresponding to an arcmin at each
            input redshift.
        """
        return self.comoving_transverse_distance(z).to(u.kpc) / RAD_IN_ARCMIN

    def kpc_proper_per_arcmin(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """Separation in transverse proper kpc equal to an arcminute at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        d : Quantity ['length']
            The distance in proper kpc corresponding to an arcmin at each input
            redshift.
        """
        return self.angular_diameter_distance(z).to(u.kpc) / RAD_IN_ARCMIN

    def arcsec_per_kpc_comoving(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """Angular separation in arcsec equal to a comoving kpc at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        theta : Quantity ['angle']
            The angular separation in arcsec corresponding to a comoving kpc at
            each input redshift.
        """
        return RAD_IN_ARCSEC / self.comoving_transverse_distance(z).to(u.kpc)

    def arcsec_per_kpc_proper(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        """Angular separation in arcsec corresponding to a proper kpc at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        theta : Quantity ['angle']
            The angular separation in arcsec corresponding to a proper kpc at
            each input redshift.
        """
        return RAD_IN_ARCSEC / self.angular_diameter_distance(z).to(u.kpc)
