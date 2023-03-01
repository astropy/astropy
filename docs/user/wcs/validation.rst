.. _validation:

Validation and Bounds checking
******************************

.. _wcslint:

Validating the WCS keywords in a FITS file
------------------------------------------

``astropy`` includes a command-line tool, ``wcslint``, to check the WCS
keywords in a FITS file. The example below shows it reporting back
results for a problematic file named ``invalid.fits``::

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

.. _wcs-bounds-check:

Bounds checking
---------------

Bounds checking is enabled by default, and any computed world
coordinates outside of [-180°, 180°] for longitude and [-90°, 90°] in
latitude are marked as invalid.  To disable this behavior, use
`astropy.wcs.Wcsprm.bounds_check`.
