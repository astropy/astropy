# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""World Coordinate System (WCS) transformations in FITS files.

.. _wcslib: https://www.atnf.csiro.au/people/mcalabre/WCS/wcslib/index.html
.. _distortion paper: https://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf
.. _SIP: https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf
.. _FITS WCS standard: https://fits.gsfc.nasa.gov/fits_wcs.html

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

from lazy_loader import attach_stub

__getattr__, __dir__, __all__ = attach_stub(__name__, __file__)
