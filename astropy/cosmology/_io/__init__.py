# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Read/Write methods for :mod:`astropy.cosmology`. **NOT public API**.

"""

# Import the readers and writers to register them into the io registry.
from astropy.cosmology._io import mapping, table  # noqa: F403
from astropy.cosmology._io import ecsv, json  # noqa: F403
