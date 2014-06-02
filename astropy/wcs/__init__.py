# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
.. _wcslib: http://www.atnf.csiro.au/~mcalabre/WCS/
.. _Paper IV: http://www.atnf.csiro.au/people/mcalabre/WCS/index.html
.. _SIP: http://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf
.. _ds9: http://hea-www.harvard.edu/RD/ds9/

Introduction
------------

`astropy.wcs` contains utilities for managing World Coordinate System
(WCS) transformations in FITS files.  These transformations map the
pixel locations in an image to their real-world units, such as their
position on the sky sphere.

It is at its base a wrapper around Mark Calabretta's `wcslib`_, but
also adds support for the Simple Imaging Polynomial (`SIP`_)
convention and table lookup distortions as defined in WCS `Paper IV`_.
Each of these transformations can be used independently or together in
a standard pipeline.
"""

from __future__ import division  # confidence high

try:
    # Not guaranteed available at setup time
    from .wcs import *
    from . import utils
except ImportError:
    if not _ASTROPY_SETUP_:
        raise


def get_include():
    """
    Get the path to astropy.wcs's C header files.
    """
    import os
    return os.path.join(os.path.dirname(__file__), "include")
