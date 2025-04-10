"""Scale factor.

This is private API. See `~astropy.cosmology.parts` for public API.

"""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ["ScaleFactor"]


from numpy.typing import ArrayLike

import astropy.units as u
from astropy.cosmology._src.utils import aszarr, deprecated_keywords


class ScaleFactor:
    """The object has attributes and methods for computing the cosmological scale factor.

    The scale factor is defined as :math:`a = a_0 / (1 + z)`.

    """

    @property
    def scale_factor0(self) -> u.Quantity:
        r"""Scale factor at redshift 0.

        The scale factor is defined as :math:`a = \frac{a_0}{1 + z}`. The common
        convention is to set :math:`a_0 = 1`. However, in some cases, e.g. in
        some old CMB papers, :math:`a_0` is used to normalize `a` to be a
        convenient number at the redshift of interest for that paper. Explicitly
        using :math:`a_0` in both calculation and code avoids ambiguity.
        """
        return 1 << u.one

    @deprecated_keywords("z", since="7.0")
    def scale_factor(self, z: u.Quantity | ArrayLike) -> u.Quantity:
        """Scale factor at redshift ``z``.

        The scale factor is defined as :math:`a = a_0 / (1 + z)`.

        Parameters
        ----------
        z : Quantity-like ['redshift'] | array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        |Quantity| | ndarray | float
            Scale factor at each input redshift.
            Returns `float` if the input is scalar.
        """
        return self.scale_factor0 / (aszarr(z) + 1)
