.. doctest-skip-all
.. _astropy-wcs:

***************************************
World Coordinate System (`astropy.wcs`)
***************************************

.. _wcslib: http://www.atnf.csiro.au/~mcalabre/WCS/
.. _FITS WCS standard: http://fits.gsfc.nasa.gov/fits_wcs.html
.. _distortion paper: http://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf
.. _SIP: http://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf

Introduction
============

`astropy.wcs` contains utilities for managing World Coordinate System
(WCS) transformations in FITS files.  These transformations map the
pixel locations in an image to their real-world units, such as their
position on the sky sphere.  These transformations can work both
forward (from pixel to sky) and backward (from sky to pixel).

It performs three separate classes of WCS transformations:

- Core WCS, as defined in the `FITS WCS standard`_, based on Mark
  Calabretta's `wcslib`_.  (Also includes ``TPV`` and ``TPD``
  distortion, but not ``SIP``).

- Simple Imaging Polynomial (`SIP`_) convention.

- table lookup distortions as defined in the FITS WCS `distortion
  paper`_.

Each of these transformations can be used independently or together in
a standard pipeline.

Getting Started
===============

The basic workflow is as follows:

    1. ``from astropy import wcs``

    2. Call the `~astropy.wcs.WCS` constructor with an
       `astropy.io.fits` `~astropy.io.fits.Header` and/or
       `~astropy.io.fits.HDUList` object.

    3. Optionally, if the FITS file uses any deprecated or
       non-standard features, you may need to call one of the
       `~astropy.wcs.wcs.WCS.fix` methods on the object.

    4. Use one of the following transformation methods:

       - From pixels to world coordinates:

         - `~astropy.wcs.wcs.WCS.all_pix2world`: Perform all three
           transformations in series (core WCS, SIP and table lookup
           distortions) from pixel to world coordinates.  Use this one
           if you're not sure which to use.

         - `~astropy.wcs.wcs.WCS.wcs_pix2world`: Perform just the core
           WCS transformation from pixel to world coordinates.

       - From world to pixel coordinates:

         - `~astropy.wcs.wcs.WCS.all_world2pix`: Perform all three
           transformations (core WCS, SIP and table lookup
           distortions) from world to pixel coordinates, using an
           iterative method if necessary.

         - `~astropy.wcs.wcs.WCS.wcs_world2pix`: Perform just the core
           WCS transformation from world to pixel coordinates.

       - Performing `SIP`_ transformations only:

         - `~astropy.wcs.wcs.WCS.sip_pix2foc`: Convert from pixel to
           focal plane coordinates using the `SIP`_ polynomial
           coefficients.

         - `~astropy.wcs.wcs.WCS.sip_foc2pix`: Convert from focal
           plane to pixel coordinates using the `SIP`_ polynomial
           coefficients.

       - Performing `distortion paper`_ transformations only:

         - `~astropy.wcs.wcs.WCS.p4_pix2foc`: Convert from pixel to
           focal plane coordinates using the table lookup distortion
           method described in the FITS WCS `distortion paper`_.

         - `~astropy.wcs.wcs.WCS.det2im`: Convert from detector
           coordinates to image coordinates.  Commonly used for narrow
           column correction.

For example, to convert pixel coordinates to world coordinates::

    >>> from astropy.wcs import WCS
    >>> w = WCS('image.fits')
    >>> lon, lat = w.all_pix2world(30, 40, 0)
    >>> print(lon, lat)


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

.. note::
    The members of the WCS object correspond roughly to the key/value
    pairs in the FITS header.  However, they are adjusted and
    normalized in a number of ways that make performing the WCS
    transformation easier.  Therefore, they can not be relied upon to
    get the original values in the header.  To build up a FITS header
    directly and specifically, use `astropy.io.fits.Header` directly.

.. _wcslint:

Validating the WCS keywords in a FITS file
------------------------------------------

Astropy includes a commandline tool, ``wcslint`` to check the WCS
keywords in a FITS file::

    > wcslint invalid.fits
    HDU 1:
      WCS key ' ':
        - RADECSYS= 'ICRS ' / Astrometric system
          RADECSYS is non-standard, use RADESYSa.
        - The WCS transformation has more axes (2) than the image it is
          associated with (0)
        - 'celfix' made the change 'PV1_5 : Unrecognized coordinate
          transformation parameter'.

    HDU 2:
      WCS key ' ':
        - The WCS transformation has more axes (3) than the image it is
          associated with (0)
        - 'celfix' made the change 'In CUNIT2 : Mismatched units type
          'length': have 'Hz', want 'm''.
        - 'unitfix' made the change 'Changed units: 'HZ      ' -> 'Hz''.

Bounds checking
---------------

Bounds checking is enabled by default, and any computed world
coordinates outside of [-180°, 180°] for longitude and [-90°, 90°] in
latitude are marked as invalid.  To disable this behavior, use
`astropy.wcs.Wcsprm.bounds_check`.

Supported projections
=====================

As `astropy.wcs` is based on `wcslib`_, it supports the standard
projections defined in the `FITS WCS standard`_.  These projection
codes are specified in the second part of the ``CTYPEn`` keywords
(accessible through `Wcsprm.ctype <astropy.wcs.Wcsprm.ctype>`), for
example, ``RA---TAN-SIP``.  The supported projection codes are:

- ``AZP``: zenithal/azimuthal perspective
- ``SZP``: slant zenithal perspective
- ``TAN``: gnomonic
- ``STG``: stereographic
- ``SIN``: orthographic/synthesis
- ``ARC``: zenithal/azimuthal equidistant
- ``ZPN``: zenithal/azimuthal polynomial
- ``ZEA``: zenithal/azimuthal equal area
- ``AIR``: Airy's projection
- ``CYP``: cylindrical perspective
- ``CEA``: cylindrical equal area
- ``CAR``: plate carrée
- ``MER``: Mercator's projection
- ``COP``: conic perspective
- ``COE``: conic equal area
- ``COD``: conic equidistant
- ``COO``: conic orthomorphic
- ``SFL``: Sanson-Flamsteed ("global sinusoid")
- ``PAR``: parabolic
- ``MOL``: Mollweide's projection
- ``AIT``: Hammer-Aitoff
- ``BON``: Bonne's projection
- ``PCO``: polyconic
- ``TSC``: tangential spherical cube
- ``CSC``: COBE quadrilateralized spherical cube
- ``QSC``: quadrilateralized spherical cube
- ``HPX``: HEALPix
- ``XPH``: HEALPix polar, aka "butterfly"

And, if built with wcslib 5.0 or later, the following polynomial
distortions are supported:

- ``TPV``: Polynomial distortion
- ``TUV``: Polynomial distortion

.. note::

    Though wcslib 5.4 and later handles ``SIP`` polynomial distortion,
    for backward compatibility, ``SIP`` is handled by astropy itself
    and methods exist to handle it specially.

Subsetting and Pixel Scales
===========================

WCS objects can be broken apart into their constituent axes using the
`~astropy.wcs.WCS.sub` function.  There is also a `~astropy.wcs.WCS.celestial`
convenience function that will return a WCS object with only the celestial axes
included.

The pixel scales of a celestial image or the pixel dimensions of a non-celestial
image can be extracted with the utility functions
`~astropy.wcs.utils.proj_plane_pixel_scales` and
`~astropy.wcs.utils.non_celestial_pixel_scales`. Likewise, celestial pixel
area can be extracted with the utility function
`~astropy.wcs.utils.proj_plane_pixel_area`.

Matplotlib plots with correct WCS projection
============================================

The `WCSAxes <http://wcsaxes.readthedocs.org>`_ affiliated package adds the
ability to use the :class:`~astropy.wcs.WCS` to define projections in
Matplotlib. More information on installing and using WCSAxes can be found `here
<http://wcsaxes.readthedocs.org>`__.

.. plot::
    :include-source:

    import numpy as np
    from matplotlib import rcParams, pyplot as plt
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.visualization import astropy_mpl_style
    rcParams.update(astropy_mpl_style)
    from astropy.utils.data import download_file

    fits_file = 'http://data.astropy.org/tutorials/FITS-images/HorseHead.fits'
    image_file = download_file(fits_file, cache=True )
    hdu = fits.open(image_file)[0]
    wcs = WCS(hdu.header)

    fig = plt.figure()
    fig.add_subplot(111, projection=wcs)
    plt.imshow(hdu.data, origin='lower', cmap='cubehelix')
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.show()

Other information
=================

.. toctree::
   :maxdepth: 1

   relax
   history

See Also
========

- `wcslib`_

Reference/API
=============

.. automodapi:: astropy.wcs

.. automodapi:: astropy.wcs.utils

Acknowledgments and Licenses
============================

`wcslib`_ is licenced under the `GNU Lesser General Public License
<http://www.gnu.org/licenses/lgpl.html>`_.
