# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines constants used in `astropy.samp`.
"""

from astropy.utils.data import get_pkg_data_filename

__all__ = [
    "SAFE_MTYPES",
    "SAMP_ICON",
    "SAMP_STATUS_ERROR",
    "SAMP_STATUS_OK",
    "SAMP_STATUS_WARNING",
]

__profile_version__ = "1.3"

#: General constant for samp.ok status string
SAMP_STATUS_OK = "samp.ok"
#: General constant for samp.warning status string
SAMP_STATUS_WARNING = "samp.warning"
#: General constant for samp.error status string
SAMP_STATUS_ERROR = "samp.error"

SAFE_MTYPES = [
    "samp.app.*",
    "samp.msg.progress",
    "table.*",
    "image.*",
    "coord.*",
    "spectrum.*",
    "bibcode.*",
    "voresource.*",
]

with open(get_pkg_data_filename("data/astropy_icon.png"), "rb") as f:
    SAMP_ICON = f.read()
