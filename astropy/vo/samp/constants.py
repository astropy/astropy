# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Defines constants used in `astropy.vo.samp`."""

import os

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

__all__ = ['SAMP_STATUS_OK', 'SAMP_STATUS_WARNING', 'SAMP_STATUS_ERROR',
           'SAMP_HUB_SINGLE_INSTANCE', 'SAMP_HUB_MULTIPLE_INSTANCE',
           'SAMP_RESTRICT_GROUP', 'SAMP_RESTRICT_OWNER',
           'SAFE_MTYPES', 'SAMP_ICON']

__profile_version__ = "1.3"

#: General constant for samp.ok status string
SAMP_STATUS_OK = "samp.ok"
#: General constant for samp.warning status string
SAMP_STATUS_WARNING = "samp.warning"
#: General constant for samp.error status string
SAMP_STATUS_ERROR = "samp.error"

#: General constant to specify single instance Hub running mode
SAMP_HUB_SINGLE_INSTANCE = "single"
#: General constant to specify multiple instance Hub running mode
SAMP_HUB_MULTIPLE_INSTANCE = "multiple"

#: General constant to specify the access restriction (through Basic Authentication) to the GROUP
SAMP_RESTRICT_GROUP = "GROUP"
#: General constant to specify the access restriction (through Basic Authentication) to the OWNER
SAMP_RESTRICT_OWNER = "OWNER"

SAFE_MTYPES = ["samp.app.*", "samp.msg.progress", "table.*", "image.*",
               "coord.*", "spectrum.*", "bibcode.*", "voresource.*"]

with open(os.path.join(DATA_DIR, 'astropy_icon.png'), 'rb') as f:
    SAMP_ICON = f.read()

# TODO: document this global variable. Is this the right place for it?
_THREAD_STARTED_COUNT = 0
