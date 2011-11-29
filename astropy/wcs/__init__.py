"""
.. _wcslib: http://www.atnf.csiro.au/~mcalabre/WCS/
.. _pyfits: http://www.stsci.edu/resources/software_hardware/pyfits
.. _Paper IV: http://www.atnf.csiro.au/people/mcalabre/WCS/index.html
.. _SIP: http://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf
.. _ds9: http://hea-www.harvard.edu/RD/ds9/

astropy.wcs provides transformations following the `SIP`_ conventions,
`Paper IV`_ table lookup distortion, and the core WCS functionality
provided by `wcslib`_.  Each of these transformations can be used
independently or together in a standard pipeline.

The basic workflow is as follows:

    1. ``from astropy import wcs``

    2. Call the `wcs.WCS` constructor with a `pyfits`_ header
       and/or hdulist object.

    3. Optionally, if the FITS file uses any deprecated or
       non-standard features, you may need to call one of the
       `~astropy.wcs.WCS.fix` methods on the object.

    4. Use one of the following transformation methods:

       - `~WCS.all_pix2sky`: Perform all three transformations from
         pixel to sky coordinates.

       - `~WCS.wcs_pix2sky`: Perform just the core WCS transformation
         from pixel to sky coordinates.

       - `~WCS.wcs_sky2pix`: Perform just the core WCS transformation
         from sky to pixel coordinates.

       - `~WCS.sip_pix2foc`: Convert from pixel to focal plane
         coordinates using the `SIP`_ polynomial coefficients.

       - `~WCS.sip_foc2pix`: Convert from focal plane to pixel
         coordinates using the `SIP`_ polynomial coefficients.

       - `~WCS.p4_pix2foc`: Convert from pixel to focal plane
         coordinates using the table lookup distortion method
         described in `Paper IV`_.

       - `~WCS.det2im`: Convert from detector coordinates to image
         coordinates.  Commonly used for narrow column correction.
"""

from __future__ import division  # confidence high

import sys
from .wcs import *


class Wcsprm(Wcsprm):
    pass
