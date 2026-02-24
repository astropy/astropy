# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains implementations of several power spectra for cosmology.
"""

__all__ = (
    "ScaleInvariantSpectrum",
    "PowerLawSpectrum",
    "RunningPowerLawSpectrum",
    "LogOscillationSpectrum",
    "BrokenPowerLawSpectrum",
)

def __dir__():
    """Directory, including lazily-imported objects."""
    return __all__
