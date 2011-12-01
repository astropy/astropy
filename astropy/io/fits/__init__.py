# Licensed under a 3-clause BSD style license - see PYFITS.rst

"""
A package for reading and writing FITS files and manipulating their
contents.

A module for reading and writing Flexible Image Transport System
(FITS) files.  This file format was endorsed by the International
Astronomical Union in 1999 and mandated by NASA as the standard format
for storing high energy astrophysics data.  For details of the FITS
standard, see the NASA/Science Office of Standards and Technology
publication, NOST 100-2.0.
"""

import os

from astropy import setup_helpers
if setup_helpers.is_in_build_mode():
    pass
else:
    from . import py3compat

    # Public API compatibility imports
    from . import card
    from . import column
    from . import convenience
    from . import hdu
    from .card import *
    from .column import *
    from .convenience import *
    from .fitsrec import FITS_record, FITS_rec
    from .hdu import *

    from .hdu.groups import GroupData
    from .hdu.hdulist import fitsopen as open
    from .hdu.image import Section
    from .hdu.table import new_table
    from .header import Header


    # Additional imports used by the documentation (some of which should be
    # restructured at some point)
    from .verify import VerifyError


# Set module-global boolean variables--these variables can also get their
# values from environment variables
GLOBALS = [
     # Variable name                       # Default
    ('EXTENSION_NAME_CASE_SENSITIVE',      False),
    ('USE_MEMMAP',                         True),
    ('ENABLE_RECORD_VALUED_KEYWORD_CARDS', True)
]


for varname, default in GLOBALS:
    try:
        locals()[varname] = bool(int(os.environ.get('ASTROPY_FITS_' + varname,
                                                    default)))
    except ValueError:
        locals()[varname] = default


__all__ = (card.__all__ + column.__all__ + convenience.__all__ + hdu.__all__ +
          ['FITS_record', 'FITS_rec', 'GroupData', 'open', 'Section',
           'new_table', 'Header', 'VerifyError',
           'setExtensionNameCaseSensitive'] + [g[0] for g in GLOBALS])
