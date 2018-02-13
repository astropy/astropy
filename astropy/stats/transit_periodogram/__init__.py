# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
transit_periodogram
===================

AstroPy-compatible reference implementation of the transit periorogram used
to discover transiting exoplanets.

"""

__all__ = ["TransitPeriodogram", "TransitPeriodogramResults"]

from .core import TransitPeriodogram, TransitPeriodogramResults
