# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains implementations of several power spectra for cosmology.
"""

__all__ = (
    "AnalyticalPowerSpectrum",
    "BrokenPowerLawSpectrum",
    "LogOscillationSpectrum",
    "PowerLawSpectrum",
    "PowerSpectrum",
    "RunningPowerLawSpectrum",
    "ScaleInvariantSpectrum",
)

from ._src.spectra import (
    AnalyticalPowerSpectrum,
    BrokenPowerLawSpectrum,
    LogOscillationSpectrum,
    PowerLawSpectrum,
    PowerSpectrum,
    RunningPowerLawSpectrum,
    ScaleInvariantSpectrum,
)
