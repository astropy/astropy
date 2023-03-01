.. include:: references.txt
.. _legacy_interface:

Legacy Interface
****************

astropy.wcs API
^^^^^^^^^^^^^^^

The ``Low Level API`` or ``Legacy Interface`` is the original `astropy.wcs` API.
It supports three types of transforms:

- Core WCS, as defined in the `FITS WCS standard`_, based on Mark
  Calabretta's `wcslib`_.  (Also includes ``TPV`` and ``TPD``
  distortion, but not ``SIP``).

- Simple Imaging Polynomial (`SIP`_) convention. (See :doc:`note about SIP in headers <note_sip>`.)

- Table lookup distortions as defined in the FITS WCS `distortion paper`_.

Each of these transformations can be used independently or together in a standard pipeline.
All methods support scalar and array inputs. Note, that all methods require an additional
positional argument which is the ``origin`` of the inputs. It has two possible values - ``0`` -
for zero-based coordinates like numpy arrays or ``1`` - for 1-based coordinates, like
the FITS standard, or those coming from ds9.

The basic workflow is to create a WCS object calling the WCS constructor with an
`~astropy.io.fits.Header` and/or `~astropy.io.fits.HDUList` object and calling
one of the methods below::

    >>> from astropy import wcs
    >>> from astropy.io import fits
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> fn = get_pkg_data_filename('data/j94f05bgq_flt.fits', package='astropy.wcs.tests')
    >>> f = fits.open(fn)
    >>> wcsobj = wcs.WCS(f[1].header)
    >>> f.close()

Optionally, if the FITS file uses any deprecated or non-standard features, you may need
to call one of the `~astropy.wcs.wcs.WCS.fix` methods on the object.

Use one of the following transformation methods.

1. Between pixels and world coordinates using all distortions:

  - `~astropy.wcs.wcs.WCS.all_pix2world`: Perform all three
    transformations in series (core WCS, SIP and table lookup
    distortions) from pixel to world coordinates.  Use this one
    if you're not sure which to use.

    >>> lon, lat = wcsobj.all_pix2world(30, 40, 0)
    >>> print(lon, lat)  # doctest: +FLOAT_CMP
    5.528442425094046 -72.05207808966726

  - `~astropy.wcs.wcs.WCS.all_world2pix`: Perform all three
     transformations (core WCS, SIP and table lookup
     distortions) from world to pixel coordinates, using an
     iterative method if necessary.

     >>> x, y = wcsobj.all_world2pix(lon, lat, 0)
     >>> print(x, y) # # doctest: +FLOAT_CMP
     30.00000214673885 39.999999958235094

 2. Performing `SIP`_ transformations only:

     - `~astropy.wcs.wcs.WCS.sip_pix2foc`: Convert from pixel to
        focal plane coordinates using the `SIP`_ polynomial
        coefficients.

        >>> xsip, ysip = wcsobj.sip_pix2foc(30, 40, 0)
        >>> print(xsip, ysip)  # doctest: +FLOAT_CMP
        -1985.8600487630586 -984.4223711273145

     - `~astropy.wcs.wcs.WCS.sip_foc2pix`: Convert from focal
        plane to pixel coordinates using the `SIP`_ polynomial
        coefficients. Note that this method only works if the
        inverse SIP distortion is specified in the header.

 3. Performing `distortion paper`_ transformations only:

     - `~astropy.wcs.wcs.WCS.p4_pix2foc`: Convert from pixel to
        focal plane coordinates using the table lookup distortion
        method described in the FITS WCS `distortion paper`_.

     - `~astropy.wcs.wcs.WCS.det2im`: Convert from detector
        coordinates to image coordinates.  Commonly used for narrow
        column correction.

Core wcslib API
^^^^^^^^^^^^^^^

The core wcslib API supports the FITS WCS standard defined in WCS
papers, I, II, III, IV. Note that distortions are not applied if
the functions in the core library are used.

1. From pixels to world coordinates:

    - `~astropy.wcs.wcs.WCS.wcs_pix2world`: Perform just the core WCS
       transformation from pixel to world coordinates.

        >>> lon, lat = wcsobj.wcs_pix2world(30, 40, 0)
        >>> print(lon, lat)  # doctest: +FLOAT_CMP
        5.527103615238458 -72.0522441352217

2. From world to pixel coordinates:

    - `~astropy.wcs.wcs.WCS.wcs_world2pix`: Perform the core WCS transformation
       from world to pixel coordinates.

        >>> x, y = wcsobj.wcs_world2pix(lon, lat, 0)
        >>> print(x, y)  # doctest: +FLOAT_CMP
        30.000000000223267 40.0000000003696
