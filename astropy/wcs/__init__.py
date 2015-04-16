# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
.. _wcslib: http://www.atnf.csiro.au/~mcalabre/WCS/
.. _FITS WCS standard: http://www.atnf.csiro.au/people/mcalabre/WCS/index.html
.. _distortion paper: http://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf
.. _SIP: http://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf

`astropy.wcs` contains utilities for managing World Coordinate System
(WCS) transformations in FITS files.  These transformations map the
pixel locations in an image to their real-world units, such as their
position on the sky sphere.

It performs three separate classes of WCS transformations:

- Core WCS, as defined in the `FITS WCS standard`_, based on Mark
  Calabretta's `wcslib`_.  See `~astropy.wcs.Wcsprm`.
- Simple Imaging Polynomial (`SIP`_) convention.  See
  `~astropy.wcs.Sip`.
- table lookup distortions as defined in WCS `distortion paper`_.  See
  `~astropy.wcs.DistortionLookupTable`.

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
