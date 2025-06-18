from numpy.typing import ArrayLike

from astropy.units import Quantity


class _MatterDensity:
    Om0: Quantity
    """Omega matter; matter density/critical density at z=0."""

    def Om(self, z: Quantity | ArrayLike) -> Quantity:
        """
        Calculate the matter density at redshift `z`.

        Parameters
        ----------
        z : Quantity or array-like
            Redshift.

        Returns
        -------
        Quantity
            Matter density at redshift `z`.
        """
        return self.Om0 * (1 + z) ** 3
