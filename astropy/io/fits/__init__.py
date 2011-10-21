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


    __all__ = (card.__all__ + column.__all__ + convenience.__all__ +
               hdu.__all__ +
               ['FITS_record', 'FITS_rec', 'GroupData', 'open', 'Section',
                'new_table', 'Header', 'VerifyError', 'USE_MEMMAP',
                'EXTENSION_NAME_CASE_SENSITIVE'])

# TODO: Hook these options into the config system
try:
    USE_MEMMAP = bool(int(os.environ.get('ASTROPY_FITS_USE_MEMMAP', 1)))
except ValueError:
    USE_MEMMAP = True

# Support case sensitive values for the value of a EXTNAME card in an extension
# header.  By default, Astropy converts the value of EXTNAME cards to upper
# case when reading from a file.  By setting EXTENSION_NAME_CASE_SENSITIVE to
# True the user may circumvent this process so that the EXTNAME value remains
# in the same case as it is in the file.
try:
    EXTENSION_NAME_CASE_SENSITIVE = \
        bool(int(os.environ.get('ASTROPY_FITS_EXTENSION_NAME_CASE_SENSITIVE',
                                0)))
except ValueError:
    EXTENSION_NAME_CASE_SENSITIVE = False
