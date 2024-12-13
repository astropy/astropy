# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy Cosmology. **NOT public API**.

The public API is provided by `astropy.cosmology.parts`.
"""

__all__ = [
    "ScaleFactor",
    "_BaryonComponent",
    "_CriticalDensity",
    "_TemperatureCMB",
]

from .baryons import _BaryonComponent
from .rhocrit import _CriticalDensity
from .scale_factor import ScaleFactor
from .tcmb import _TemperatureCMB
