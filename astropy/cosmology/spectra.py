from abc import ABC, abstractmethod
from numbers import Integral, Real
import numpy as np


class PowerSpectrum(ABC):
    """
    Abstract base class for power spectrum computations.

    Attributes:
    scale (float): Scale factor for the spectrum.
    k_cutoff (float): Cutoff value for the wave number.
    spectrum_type (str): Type of the spectrum.
    """

    def __init__(
        self,
        ndims: Integral = 3,
        scale: Real = None,
        k_cutoff: Real = None,
    ) -> None:
        """
        Initializes the PowerSpectrum with the given scale and k_cutoff.

        Args:
            scale (float, optional): Scale factor for the spectrum. Defaults to 1.0.
            k_cutoff (float, optional): Cutoff value for the wave number.
        """
        self.ndims = ndims
        self.scale = scale if scale is not None else 1.
        self.k_cutoff = k_cutoff
        self.spectrum_type = None
        self.f = 1 / (2 * np.pi)**self.ndims

    @abstractmethod
    def _compute_spectrum(self) -> None:
        """
        Abstract method to compute the power spectrum.

        Args:
            k (float): Wave number.

        Raises:
            NotImplementedError: If the method is not overridden by a subclass.
        """
        raise NotImplementedError("Method must be overriden by child.")
    
    @abstractmethod
    def _compute_integral(self):
        raise NotImplementedError("Method must be overriden by child.")

    def _compute_cutoff(self, k: Real) -> Real:
        """
        Computes the cutoff function for a given wave number.

        Args:
            k (float): Wave number.

        Returns:
            float: Cutoff modifier.
        """
        return np.exp(- k * k / (self.k_cutoff * self.k_cutoff))

    def __call__(self, k) -> Real:
        """
        Evaluates the power spectrum with or without the cutoff.

        Args:
            k (float): Wave number.

        Returns:
            float: Evaluated power spectrum.
        """
        if self.k_cutoff is None:
            return self.f * self._compute_spectrum(k)
        else:
            return self.f * self._compute_spectrum(k) * self._compute_cutoff(k)
    
    def integrate(self, k) -> Real:
        return self.f * self._compute_integral(k)
    