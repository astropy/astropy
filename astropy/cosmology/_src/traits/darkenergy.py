from abc import abstractmethod
from math import exp, log
from numbers import Number

import numpy as np

from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.utils.compat.optional_deps import HAS_SCIPY

if HAS_SCIPY:
    from scipy.integrate import quad
else:

    def quad(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'scipy.integrate'")


class DarkEnergyComponent:
    # Subclasses should use `Parameter` to make this a parameter of the cosmology.
    Ode0: float
    """Omega dark energy; dark energy density/critical density at z=0."""

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
            The density of dark energy relative to the critical density at each
            redshift.
            Returns `float` if the input is scalar.
        """
        z = aszarr(z)
        if self.Ode0 == 0:  # Common enough to be worth checking explicitly
            return np.zeros(z.shape) if hasattr(z, "shape") else 0.0
        # Ensure self.inv_efunc is implemented by the main class
        if not hasattr(self, "inv_efunc") or not callable(self.inv_efunc):
            raise NotImplementedError(
                "The main class must implement an 'inv_efunc(z)' method."
            )
        return self.Ode0 * self.de_density_scale(z) * self.inv_efunc(z) ** 2
