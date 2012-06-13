***************************************
World Coordinate System (`astropy.wcs`)
***************************************

.. _wcslib: http://www.atnf.csiro.au/~mcalabre/WCS/
.. _Paper IV: http://www.atnf.csiro.au/people/mcalabre/WCS/index.html
.. _SIP: http://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf
.. _ds9: http://hea-www.harvard.edu/RD/ds9/

Introduction
============

`astropy.wcs` contains utilities for managing World Coordinate System
(WCS) transformations in FITS files.  These transformations map the
pixel locations in an image to their real-world units, such as their
position on the sky sphere.

It is at its base a wrapper around Mark Calabretta's `wcslib`_, but
also adds support for the Simple Imaging Polynomial (`SIP`_)
convention and table lookup distortions as defined in WCS `Paper IV`_.
Each of these transformations can be used independently or together in
a standard pipeline.

Getting Started
===============

The basic workflow is as follows:

    1. ``from astropy import wcs``

    2. Call the `wcs.WCS` constructor with an `astropy.io.fits` header
       and/or hdulist object.

    3. Optionally, if the FITS file uses any deprecated or
       non-standard features, you may need to call one of the
       `~astropy.wcs.WCS.fix` methods on the object.

    4. Use one of the following transformation methods:

       - `~WCS.all_pix2world`: Perform all three transformations from
         pixel to world coordinates.

       - `~WCS.wcs_pix2world`: Perform just the core WCS
         transformation from pixel to world coordinates.

       - `~WCS.wcs_world2pix`: Perform just the core WCS
         transformation from world to pixel coordinates.

       - `~WCS.sip_pix2foc`: Convert from pixel to focal plane
         coordinates using the `SIP`_ polynomial coefficients.

       - `~WCS.sip_foc2pix`: Convert from focal plane to pixel
         coordinates using the `SIP`_ polynomial coefficients.

       - `~WCS.p4_pix2foc`: Convert from pixel to focal plane
         coordinates using the table lookup distortion method
         described in `Paper IV`_.

       - `~WCS.det2im`: Convert from detector coordinates to image
         coordinates.  Commonly used for narrow column correction.


Using `astropy.wcs`
===================

Loading WCS information from a FITS file
----------------------------------------

This example loads a FITS file (supplied on the commandline) and uses
the WCS cards in its primary header to transform.

.. literalinclude:: examples/from_file.py
   :language: python

Building a WCS structure programmatically
-----------------------------------------

This example, rather than starting from a FITS header, sets WCS values
programmatically, uses those settings to transform some points, and then
saves those settings to a new FITS header.

.. literalinclude:: examples/programmatic.py
   :language: python

Other information
=================

.. toctree::
   :maxdepth: 1

   relax
   units
   history



See Also
========

- `wcslib`_

Reference/API
=============

.. automodapi:: astropy.wcs


Acknowledgments and Licenses
============================

wcslib is licenced under the `GNU Lesser General Public License
<http://www.gnu.org/licenses/lgpl.html>`_.
