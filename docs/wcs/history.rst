astropy.wcs History
===================

`astropy.wcs` began life as ``pywcs``.  Earlier version numbers refer to
that package.

pywcs Version 1.11
------------------

- Updated to wcslib version 4.8, which gives much more detailed error
  messages.

- Added functions get_pc() and get_cdelt().  These provide a way to
  always get the canonical representation of the linear transformation
  matrix, whether the header specified it in PC, CD or CROTA form.

- Long-running process will now release the Python GIL to better
  support Python multithreading.

- The dimensions of the `~astropy.wcs.Wcsprm.cd` and
  `~astropy.wcs.Wcsprm.pc` matrices were always returned as 2x2.  They
  now are sized according to naxis.

- Supports Python 3.x

- Builds on Microsoft Windows without severely patching wcslib.

- Lots of new unit tests

- ``pywcs`` will now run without ``pyfits``, though the SIP and distortion
  lookup table functionality is unavailable.

- Setting `~astropy.wcs.Wcsprm.cunit` will now verify that the values
  are valid unit strings.

pywcs Version 1.10
------------------

- Adds a ``UnitConversion`` class, which gives access to wcslib's unit
  conversion functionality.  Given two convertible unit strings, pywcs
  can convert arrays of values from one to the other.

- Now uses wcslib 4.7

- Changes to some wcs values would not always calculate secondary values.

pywcs Version 1.9
-----------------

- Support binary image arrays and pixel list format WCS by presenting
  a way to call wcslib's ``wcsbth()``

- Updated underlying wcslib to version 4.5, which fixes the following:

    - Fixed the interpretation of VELREF when translating
      AIPS-convention spectral types.  Such translation is now handled
      by a new special- purpose function, spcaips().  The wcsprm
      struct has been augmented with an entry for velref which is
      filled by wcspih() and wcsbth().  Previously, selection by
      VELREF of the radio or optical velocity convention for type VELO
      was not properly handled.

Bugs
````

- The `~astropy.wcs.Wcsprm.pc` member is now available with a default
  raw `~astropy.wcs.Wcsprm` object.

- Make properties that return arrays read-only, since modifying a
  (mutable) array could result in secondary values not being
  recomputed based on those changes.

- `float` properties can now be set using `int` values

pywcs Version 1.3a1
-------------------

Earlier versions of pywcs had two versions of every conversion method::

  X(...)      -- treats the origin of pixel coordinates at (0, 0)
  X_fits(...) -- treats the origin of pixel coordinates at (1, 1)

From version 1.3 onwards, there is only one method for each
conversion, with an 'origin' argument:

  - 0: places the origin at (0, 0), which is the C/Numpy convention.

  - 1: places the origin at (1, 1), which is the Fortran/FITS
    convention.

