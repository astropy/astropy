from abc import abstractmethod
from math import exp, log

import numpy as np
from numpy.typing import ArrayLike

from astropy.cosmology._src.scipy_compat import quad
from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class DarkEnergyComponent:
    # Subclasses should use `Parameter` to make this a parameter of the cosmology.
    Ode0: float
    """Omega dark energy; dark energy density/critical density at z=0."""

    @abstractmethod
    @deprecated_keywords("z", since="7.0")
    def w(self, z: Quantity | ArrayLike) -> FArray:
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

    def _w_integrand(self, ln1pz: float | FArray, /) -> FArray:
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
    def de_density_scale(self, z: Quantity | ArrayLike) -> FArray:
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
        ival = (
            quad(self._w_integrand, 0, log(z + 1.0))[0]  # scalar
            if z.ndim == 0
            else np.asarray([quad(self._w_integrand, 0, log(1 + _z))[0] for _z in z])
        )
        return np.exp(3 * ival)

    @deprecated_keywords("z", since="7.0")
    def Ode(self, z: Quantity | ArrayLike) -> FArray:
        """Return the density parameter for dark energy at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Ode : ndarray
            The density of dark energy relative to the critical density at each
            redshift.
        """
        z = aszarr(z)
        if self.Ode0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros_like(z)
        # Ensure self.inv_efunc is implemented by the main class
        if not hasattr(self, "inv_efunc") or not callable(self.inv_efunc):
            msg = "The main class must implement an 'inv_efunc(z)' method."
            raise NotImplementedError(msg)
        return self.Ode0 * self.de_density_scale(z) * self.inv_efunc(z) ** 2
