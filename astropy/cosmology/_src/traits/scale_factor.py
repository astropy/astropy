# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Scale factor.

This is private API. See `~astropy.cosmology.traits` for public API.

"""

__all__ = ["ScaleFactor"]


from numpy.typing import ArrayLike

import astropy.units as u
from astropy.cosmology._src.utils import aszarr, deprecated_keywords


class ScaleFactor:
    """The trait for computing the cosmological scale factor.

    The scale factor is defined as :math:`a = a_0 / (1 + z)`.

    """

    @property
    def scale_factor0(self) -> u.Quantity:
        r"""Scale factor at redshift 0.

        The scale factor is defined as :math:`a = a_0 / (1 + z)`. The common convention
        is to set :math:`a_0 = 1`. However, in some cases, like in some old CMB papers,
        :math:`a_0` is used to normalize `a` to be a convenient number at the redshift
        of interest for that paper. Explicitly using :math:`a_0` in both calculation and
        code avoids ambiguity.
        """
        return 1 << u.one

    @deprecated_keywords("z", since="7.0")
    def scale_factor(self, z: u.Quantity | ArrayLike) -> u.Quantity:
        """Compute the scale factor at redshift ``z``.

        The scale factor is defined as :math:`a = a_0 / (1 + z)`.

        Parameters
        ----------
        z : Quantity-like ['redshift'] | array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        |Quantity|
            Scale factor at each input redshift.
        """
        return self.scale_factor0 / (aszarr(z) + 1)
