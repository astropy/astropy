"""Matter-only cosmology."""

__all__ = ("MatterOnlyCosmology",)

from numbers import Number
from typing import TypeVar, overload

import numpy as np
from numpy.typing import ArrayLike

import astropy.units as u
from astropy.cosmology._src.core import dataclass_decorator
from astropy.cosmology._src.parameter import Parameter
from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr
from astropy.cosmology.traits import DistanceMeasures, HubbleParameter, MatterComponent

_InputT = TypeVar("_InputT", bound=u.Quantity | np.ndarray | np.generic | Number)


@dataclass_decorator
class MatterOnlyCosmology(DistanceMeasures, HubbleParameter, MatterComponent):
    r"""Matter-only cosmology.

    $$ E(z) = (1+z) \sqrt{1+\Omega_{m,0}z}} $$

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

    def Om(self, z: u.Quantity | ArrayLike, /) -> FArray:
        r"""Return the density parameter for non-relativistic matter at redshift ``z``.

        $$ \Omega_m(z) = \frac{\Omega_{m,0} (1 + z)}{1 + \Omega_{m,0} z} $$

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        Om : ndarray
            The density of non-relativistic matter relative to the critical
            density at each redshift.

        Examples
        --------
        >>> from astropy.cosmology import MatterOnlyCosmology

        >>> cosmo = MatterOnlyCosmology(H0=70, Om0=1.0)
        >>> cosmo.Om(0)
        np.float64(1.0)

        >>> cosmo.Om([0, 1, 2])
        array([1., 1., 1.])

        """
        z = aszarr(z)
        return self.Om0 * (1 + z) / (1 + self.Om0 * z)

    # ===============================================================
    # Hubble Parameter

    def efunc(self, z: u.Quantity | ArrayLike, /) -> FArray:
        r"""Return the dimensionless Hubble parameter at redshift ``z``.

        $$ E(z) = \sqrt{\Omega_{m,0} (1 + z)^3 + (1 - \Omega_{m,0}) (1 + z)^2} $$

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        E : ndarray
            Dimensionless Hubble parameter at each redshift.

        Examples
        --------
        >>> from astropy.cosmology import MatterOnlyCosmology
        >>> cosmo = MatterOnlyCosmology(H0=70, Om0=1.0)
        >>> cosmo.efunc(0)
        np.float64(1.0)

        >>> cosmo.efunc([0, 1, 2])
        array([1. , 2.82842712, 5.19615242])

        """
        z = aszarr(z)
        return (1 + z) * np.sqrt(1 + self.Om0 * z)

    def inv_efunc(self, z: u.Quantity | ArrayLike, /) -> FArray:
        r"""Return the inverse dimensionless Hubble parameter at redshift ``z``.

        $$ \frac{1}{E(z)} = \frac{1}{(1 + z) \sqrt{1 + \Omega_{m,0} z}} $$

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        1/E : ndarray
            Inverse dimensionless Hubble parameter at each redshift.

        Examples
        --------
        >>> from astropy.cosmology import MatterOnlyCosmology
        >>> cosmo = MatterOnlyCosmology(H0=70, Om0=1.0)
        >>> cosmo.inv_efunc(0)
        np.float64(1.0)

        >>> cosmo.inv_efunc([0, 1, 2])
        array([1. , 0.35355339, 0.19245009])

        """
        z = aszarr(z)
        return 1.0 / ((1 + z) * np.sqrt(1 + self.Om0 * z))

    # ===============================================================
    # Distance Measures

    def _comoving_distance(self, z1: np.ndarray, z2: np.ndarray, /) -> np.ndarray:
        """Comoving line-of-sight distance chi(z1, z2) in units of c/H0."""
        chi: FArray
        t1 = np.sqrt(1 + self.Om0 * z1)
        t2 = np.sqrt(1 + self.Om0 * z2)
        if self.Om0 == 1:  # Einstein-de Sitter
            chi = 2 * (1 / t1 - 1 / t2)
        elif self.Om0 < 1:  # Open
            a = np.sqrt(1 - self.Om0)
            chi = (1 / a) * np.log((t2 - a) * (t1 + a) / ((t2 + a) * (t1 - a)))
        else:  # Closed
            a = np.sqrt(self.Om0 - 1)
            chi = 2 / a * (np.arctan(t2 / a) - np.arctan(t1 / a))

        return chi

    @overload
    def comoving_distance(self, z: _InputT, /) -> u.Quantity: ...

    @overload
    def comoving_distance(self, z: _InputT, z2: _InputT, /) -> u.Quantity: ...

    def comoving_distance(self, z: _InputT, z2: _InputT | None = None, /) -> u.Quantity:
        r"""Comoving line-of-sight distance :math:`d_c(z1, z2)` in Mpc.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        $$ D_C(z_1, z_2)
            = \frac{c}{H_0}\int_{z_1}^{z_2} \frac{dz'}{E(z')},
            = \frac{c}{H_0}\,|\chi(z_2) - \chi(z_1)|.
        $$

        where, $E(z) = (1+z)\sqrt{1+\Omega_{m,0}z}$.

        There are three cases, based on $\Omega_k = 1 - \Omega_m$:

        1. `\Omega_k = 0`: The Einstein-de Sitter (flat) solution.

            $$ \chi(z) = 2 (1 - \frac{1}{\sqrt{1 + \Omega_m z}}) $$

        2. `\Omega_k > 0`: The open universe solution.

            For $tz = \sqrt{1 + \Omega_m z}$

            $$ \chi(z)
                = \frac{1}{\sqrt{\Omega_k}} \left[
                     \ln(tz - \sqrt{\Omega_k}) + \ln(1 + \sqrt{\Omega_k})
                    -\ln(tz + \sqrt{\Omega_k}) - \ln(1 - \sqrt{\Omega_k})
                \right]
            $$

        3. `\Omega_k < 0`: The closed universe solution.

            $$ \chi(z)
                = \frac{2}{\sqrt{-\Omega_k}} \left(
                      \arctan\left(\frac{tz}{\sqrt{-\Omega_k}}\right)
                    - \arctan\left(\frac{1}{\sqrt{-\Omega_k}}\right)
                \right)
            $$

        Parameters
        ----------
        z, z2 : Quantity ['redshift'], positional-only
            Input redshifts. If one argument ``z`` is given, the distance
            :math:`d_c(0, z)` is returned. If two arguments ``z1, z2`` are
            given, the distance :math:`d_c(z_1, z_2)` is returned.

        Returns
        -------
        Quantity ['length']
            Comoving distance in Mpc between each input redshift.
        """
        z1, z2 = (0.0, z) if z2 is None else (z, z2)
        return self.hubble_distance * self._comoving_distance(aszarr(z1), aszarr(z2))

    @overload
    def comoving_transverse_distance(self, z: _InputT, /) -> u.Quantity: ...

    @overload
    def comoving_transverse_distance(
        self, z: _InputT, z2: _InputT, /
    ) -> u.Quantity: ...

    def comoving_transverse_distance(
        self, z: _InputT, z2: _InputT | None = None, /
    ) -> u.Quantity:
        r"""Comoving transverse distance in Mpc at a given redshift.

        This value is the transverse comoving distance at redshift ``z``
        corresponding to an angular separation of 1 radian. This is the same as
        the comoving distance if :math:`\Omega_k` is zero (as in the current
        concordance Lambda-CDM model).

        Starting from the difference of comoving distances

        $$ \Delta\Xi = | \Xi(z_2) - \Xi(z_1) | $$

        We define the $S_k$ mapping

        $$ S_k(x) = \begin{cases}
            \sinh(x) / \sqrt{-\Omega_k},    & k = -1, \\
            x,                              & k=0, \\
            \sin(x) / \sqrt{\Omega_k},      & k=+1.
        \end{cases}
        $$

        And then the transverse comoving distance is given by:

        $$ D_M(z_1,z_2) = \frac{c}{H_0} S_k(\sqrt{|\Omega_k|} \Delta\chi) $$

        Parameters
        ----------
        z, z2 : Quantity ['redshift'], positional-only
            Input redshifts. If one argument ``z`` is given, the distance
            :math:`d_c(0, z)` is returned. If two arguments ``z1, z2`` are
            given, the distance :math:`d_c(z_1, z_2)` is returned.

        Returns
        -------
        d : Quantity ['length']
            Comoving transverse distance in Mpc at each input redshift.

        Notes
        -----
        This quantity is also called the 'proper motion distance' in some texts.
        """
        z1, z2 = (0.0, z) if z2 is None else (z, z2)
        delta_chi = self._comoving_distance(aszarr(z1), aszarr(z2))

        Ok0 = 1.0 - self.Om0
        if Ok0 == 0:  # Einstein-de Sitter
            DM = delta_chi

        elif Ok0 > 0:  # Open
            sqrtOk = np.sqrt(Ok0)
            DM = np.sinh(sqrtOk * delta_chi) / sqrtOk

        else:  # Closed
            sqrtOk = np.sqrt(-Ok0)
            DM = np.sin(sqrtOk * delta_chi) / sqrtOk

        return self.hubble_distance * DM

    def comoving_volume(self, z: u.Quantity | ArrayLike, /) -> u.Quantity:
        r"""Comoving volume in cubic Mpc at redshift ``z``.

        This is the volume of the universe encompassed by redshifts less than ``z``. For
        the case of $\Omega_k = 0$ it is a sphere of radius `comoving_distance`
        but it is less intuitive if $\Omega_k$ is not.

        $$ V_c(<z) = \begin{cases}
            \pi (c/(\sqrt{|\Omega_{k,0}|} H_0))^3 [2y-\sin(2y)], & k = -1, \\
            \frac{4\pi}{3}\,D_M^3(z), & k = 0, \\
            \pi (c/(\Omega_{k,0} H_0))^3 [\sinh(2y)-2y], & k = 1, \\
        \end{cases}
        $$

        where $y=\sqrt{|\Omega_{k,0}|}\,\chi(z)$.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        V : Quantity ['volume']
            Comoving volume in Mpc^3 at each input redshift.
        """
        z = aszarr(z)
        Ok0 = 1.0 - self.Om0

        if Ok0 == 0:  # Einstein-de Sitter
            cv = 4 * np.pi / 3 * self.comoving_transverse_distance(0, z) ** 3

        elif Ok0 > 0:  # Open
            sqrtOk = np.sqrt(Ok0)
            y = 2 * sqrtOk * self._comoving_distance(0, z)
            cv = np.pi / (sqrtOk**3) * self.hubble_distance**3 * (np.sinh(y) - y)

        else:  # Closed
            sqrtOk = np.sqrt(np.abs(Ok0))
            y = 2 * sqrtOk * self._comoving_distance(0, z)
            cv = np.pi / (sqrtOk**3) * self.hubble_distance**3 * (y - np.sin(y))

        return cv
