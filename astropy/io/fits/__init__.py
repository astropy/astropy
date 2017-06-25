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


from . import py3compat

from ... import config as _config

# Set module-global boolean variables
# TODO: Make it possible to set these variables via environment variables
# again, once support for that is added to Astropy


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.io.fits`.
    """

    enable_record_valued_keyword_cards = _config.ConfigItem(
        True,
        'If True, enable support for record-valued keywords as described by '
        'FITS WCS distortion paper. Otherwise they are treated as normal '
        'keywords.',
        aliases=['astropy.io.fits.enabled_record_valued_keyword_cards'])
    extension_name_case_sensitive = _config.ConfigItem(
        False,
        'If True, extension names (i.e. the ``EXTNAME`` keyword) should be '
        'treated as case-sensitive.')
    strip_header_whitespace = _config.ConfigItem(
        True,
        'If True, automatically remove trailing whitespace for string values in '
        'headers.  Otherwise the values are returned verbatim, with all '
        'whitespace intact.')
    use_memmap = _config.ConfigItem(
        True,
        'If True, use memory-mapped file access to read/write the data in '
        'FITS files. This generally provides better performance, especially '
        'for large files, but may affect performance in I/O-heavy '
        'applications.')
    lazy_load_hdus = _config.ConfigItem(
        True,
        'If True, use lazy loading of HDUs when opening FITS files by '
        'default; that is fits.open() will only seek for and read HDUs on '
        'demand rather than reading all HDUs at once.  See the documentation '
        'for fits.open() for more datails.')
    enable_uint = _config.ConfigItem(
        True,
        'If True, default to recognizing the convention for representing '
        'unsigned integers in FITS--if an array has BITPIX > 0, BSCALE = 1, '
        'and BZERO = 2**BITPIX, represent the data as unsigned integers '
        'per this convention.')


conf = Conf()


# Public API compatibility imports
# These need to come after the global config variables, as some of the
# submodules use them
from . import card
from . import column
from . import convenience
from . import hdu
from .card import *
from .column import *
from .convenience import *
from .diff import *
from .fitsrec import FITS_record, FITS_rec
from .hdu import *

from .hdu.groups import GroupData
from .hdu.hdulist import fitsopen as open
from .hdu.image import Section
from .header import Header
from .verify import VerifyError


__all__ = (['Conf', 'conf'] + card.__all__ + column.__all__ +
           convenience.__all__ + hdu.__all__ +
           ['FITS_record', 'FITS_rec', 'GroupData', 'open', 'Section',
            'Header', 'VerifyError', 'conf'])
