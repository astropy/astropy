# -*- coding: utf-8 -*-

"""
Readers, Writers, and I/O Miscellany.
"""

__all__ = ["file_reader", "file_writer"]  # from `mypackage`

# e.g. a file reader and writer for ``myformat``
# this will be used in ``astropy_io.py``
from .core import file_reader, file_writer

# Register read and write methods into Astropy:
# determine if it is 1) installed and 2) the correct version (v5.0+)
try:
    import astropy
    from astropy.utils.introspection import minversion
except ImportError:
    ASTROPY_GE_5 = False
else:
    ASTROPY_GE_5 = minversion(astropy, "5.0")

if ASTROPY_GE_5:
    # Astropy is installed and v5.0+ so we import the following modules
    # to register "myformat" with Cosmology read/write and "mypackage"
    # with Cosmology to/from_format.
    from . import astropy_convert, astropy_io
