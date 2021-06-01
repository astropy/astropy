# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Read/Write methods for :mod:`astropy.cosmology`.

This module's namespace is flattened so that all public API functions are
here. No sub-module (e.g. ``ecsv`` or ``json``) is considered public API.

"""

from . import core
from .core import *

# Import the readers and writers to register them into the io registry.
# This is NOT public API
from . import ecsv, json  # noqa: F403

__all__ = core.__all__
