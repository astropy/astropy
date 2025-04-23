# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Traits for building ``astropy`` :class:`~astropy.cosmology.Cosmology` classes."""

__all__ = [
    "ScaleFactor",
    "TemperatureCMB",
]

from ._src.traits import ScaleFactor, TemperatureCMB
