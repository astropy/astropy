0.2 (unreleased)
----------------

New Features
^^^^^^^^^^^^

This is a brief overview of the new features included in Astropy 0.2--please
see the "What's New" section of the documentation for more details.

- ``astropy.cosmology``

  - Update to include cosmologies with variable dark energy equations of state.
    (This introduces some API incompatibilities with the older Cosmology
    objects).

  - Added parameters for relativistic species (photons, neutrinos) to the
    astropy.cosmology classes. The current treatment assumes that neutrinos are
    massless. [#365]

- ``astropy.table`` I/O infrastructure for custom readers/writers
  implemented. [#305]

  - Added support for reading/writing HDF5 files [#461]

- New ``astropy.time`` sub-package. [#332]

- New ``astropy.units`` sub-package.  This has the following effects on
  other subpackages:

  - In ``astropy.wcs``, the ``wcs.cunit`` list now takes and returns
    ``astropy.units.Unit`` objects.

- ``astropy.io.ascii``

  - Improved integration with the ``astropy.table`` Table class so that
    table and column metadata (e.g. keywords, units, description,
    formatting) are directly available in the output table object.  The
    CDS, DAOphot, and IPAC format readers now provide this type of
    integrated metadata.

  - Changed to using `astropy.table` masked tables instead of NumPy
    masked arrays for tables with missing values.

  - Added SExtractor table reader to ``astropy.io.ascii`` [#420]

  - Removed the Memory reader class which was used to convert data input
    passed to the ``write`` function into an internal table.  Instead
    ``write`` instantiates an astropy Table object using the data
    input to ``write``.

  - Removed the NumpyOutputter as the output of reading a table is now
    always a ``Table`` object.

  - Removed the option of supplying a function as a column output
    formatter.

  - Added a new ``strip_whitespace`` keyword argument to the ``write``
    function.  This controls whether whitespace is stripped from
    the left and right sides of table elements before writing.
    Default is True.

  - Fixed a bug in reading IPAC tables with null values.

- Added support for masked tables with missing or invalid data [#451]

- ``astropy.wcs``

  - From updating the the underlying wcslib 4.16:

    - When `astropy.wcs.WCS` constructs a default coordinate
      representation it will give it the special name "DEFAULTS", and
      will not report "Found one coordinate representation".

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Astropy doc themes moved into ``astropy.sphinx`` to allow affiliated packages
  to access them.

- Added expanded documentation for the ``astropy.cosmology`` sub-package.
  [#272]

- Added option to disable building of "legacy" packages (pyfits, vo, etc.).

- The value of the astronomical unit (au) has been updated to that adopted by
  IAU 2012 Resolution B2, and the values of the pc and kpc constants have been
  updated to reflect this. [#368]

- Added links to the documentation pages to directly edit the documentation on
  GitHub. [#347]

- Several updates merged from ``pywcs`` into ``astropy.wcs`` [#384]:

  - Improved the reading of distortion images.

  - Added a new option to choose whether or not to write SIP coefficients.

- Added HTML represention of tables in IPython notebook [#409]


Bug Fixes
^^^^^^^^^

- ``astropy.io.fits``

  - Fixed a bug that could cause a deadlock in the filesystem on OSX when
    reading the data from certain types of FITS files. This only occurred
    when used in conjunction with Numpy 1.7. [#369]

  - Fixed an issue where opening a FITS file containing a random group HDU in
    update mode could result in an unnecessary rewriting of the file even if
    no changes were made. This corresponds to PyFITS ticket 179.

  - Fixed a bug where opening a file containing compressed image HDUs in
    'update' mode and then immediately closing it without making any changes
    caused the file to be rewritten unncessarily.

  - Fixed two memory leaks that could occur when writing compressed image data,
    or in some cases when opening files containing compressed image HDUs in
    'update' mode.

  - Fixed a bug where ``ImageHDU.scale(option='old')`` wasn't working at
    all--it was not restoring the image to its original BSCALE and BZERO
    values.

  - Fixed a bug when writing out files containing zero-width table columns,
    where the TFIELDS keyword would be updated incorrectly, leaving the table
    largely unreadable.

  - Fixed a minor string formatting issue.

- ``astropy.io.votable``

  - The `Table` class now maintains a single array object which is a
    Numpy masked array.  For variable-length columns, the object that
    is stored there is also a Numpy masked array.

  - Changed the ``pedantic`` configuration option to be ``False`` by default
    due to the vast proliferation of non-compliant VO Tables. [#296]

  - Renamed `astropy.io.vo` to `astropy.io.votable`.

- ``astropy.table``

  - Added a workaround for an upstream bug in Numpy 1.6.2 that could cause
    a maxiumum recursion depth RuntimeError when printing table rows. [#341]

- ``astropy.wcs``

  - Updated to wcslib 4.15 [#418]

  - Fixed a problem with handling FITS headers on locales that do not use
    dot as a decimal separator. This required an upstream fix to wcslib which
    is included in wcslib 4.14. [#313]

- Fixed some tests that could fail due to missing/incorrect logging
  configuration--ensures that tests don't have any impact on the default log
  location or contents. [#291]

- Various minor documentation fixes [#293 and others]

- Fixed a bug where running the tests with the ``py.test`` command still tried
  to replace the system-installed pytest with the one bundled with Astropy.
  [#454]



0.1 (2012-06-19)
----------------

- Initial release.
