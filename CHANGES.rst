0.3 (unreleased)
----------------

- A configuration file with all options set to their defaults is now generated 
  when astropy is installed.  This file will be pulled in as the users' 
  astropy configuration file the first time they ``import astropy``.  [#498] 

- ``astropy.cosmology``

  - Add a WMAP9 object using the 9-year WMAP parameters from Hinshaw et al.
    Once this paper is accepted, this should be made the default, but for
    now WMAP7 remains the default. [#629]


0.2 (unreleased)
----------------


0.2b1 (2012-12-24)
------------------

New Features
^^^^^^^^^^^^

This is a brief overview of the new features included in Astropy 0.2--please
see the "What's New" section of the documentation for more details.

- ``astropy.coordinates``

  - This new subpackage contains a representation of celestial coordinates,
    and provides a wide range of related functionality.  While
    fully-functional, it is a work in progress and parts of the API may
    change in subsequent releases.

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

- New ``astropy.units`` sub-package that includes a class for units
  (``astropy.units.Unit``) and scalar quantities that have units
  (``astropy.units.Quantity``). [#370, #445]

  This has the following effects on other sub-packages:

  - In ``astropy.wcs``, the ``wcs.cunit`` list now takes and returns
    ``astropy.units.Unit`` objects. [#379]

  - In ``astropy.nddata``, units are now stored as ``astropy.units.Unit``
    objects. [#382]

  - In ``astropy.table``, units on columns are now stored as
    ``astropy.units.Unit`` objects. [#380]

  - In ``astropy.constants``, constants are now stored as
    ``astropy.units.Quantity`` objects. [#529]

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

- Added HTML representation of tables in IPython notebook [#409]

- Rewrote CFITSIO-based backend for handling tile compression of FITS files.
  It now uses a standard CFITSIO instead of heavily modified pieces of CFITSIO
  as before.  Astropy ships with its own copy of CFITSIO v3.30, but system
  packagers may choose instead to strip this out in favor of a
  system-installed version of CFITSIO.  This corresponds to PyFITS ticket 169.
  [#318]

- Moved `astropy.config.data` to `astropy.utils.data` and re-factored the I/O
  routines to separate out the generic I/O code that can be used to open any
  file or resource from the code used to access Astropy-related data. The
  'core' I/O routine is now `get_readable_fileobj`, which can be used to
  access any local as well as remote data, supports caching, and can
  decompress gzip and bzip2 files on-the-fly. [#425]

Bug Fixes
^^^^^^^^^

- ``astropy.io.fits``

  - Improved handling of scaled images and pseudo-unsigned integer images in
    compressed image HDUs.  They now work more transparently like normal image
    HDUs with support for the ``do_not_scale_image_data`` and ``uint`` options,
    as well as ``scale_back`` and ``save_backup``.  The ``.scale()`` method
    works better too. Corresponds to PyFITS ticket 88.

  - Permits non-string values for the EXTNAME keyword when reading in a file,
    rather than throwing an exception due to the malformatting.  Added
    verification for the format of the EXTNAME keyword when writing.
    Corresponds to PyFITS ticket 96.

  - Added support for EXTNAME and EXTVER in PRIMARY HDUs.  That is, if EXTNAME
    is specified in the header, it will also be reflected in the ``.name``
    attribute and in ``fits.info()``.  These keywords used to be verboten in
    PRIMARY HDUs, but the latest version of the FITS standard allows them.
    Corresponds to PyFITS ticket 151.

  - HCOMPRESS can again be used to compress data cubes (and higher-dimensional
    arrays) so long as the tile size is effectively 2-dimensional. In fact,
    compatible tile sizes will automatically be used even if they're not
    explicitly specified. Corresponds to PyFITS ticket 171.

  - Fixed a bug that could cause a deadlock in the filesystem on OSX when
    reading the data from certain types of FITS files. This only occurred
    when used in conjunction with Numpy 1.7. [#369]

  - Added support for the optional ``endcard`` parameter in the
    ``Header.fromtextfile()`` and ``Header.totextfile()`` methods.  Although
    ``endcard=False`` was a reasonable default assumption, there are still text
    dumps of FITS headers that include the END card, so this should have been
    more flexible. Corresponds to PyFITS ticket 176.

  - Fixed a crash when running fitsdiff on two empty (that is, zero row) tables.
    Corresponds to PyFITS ticket 178.

  - Fixed an issue where opening a FITS file containing a random group HDU in
    update mode could result in an unnecessary rewriting of the file even if
    no changes were made. This corresponds to PyFITS ticket 179.

  - Fixed a crash when generating diff reports from diffs using the
    ``ignore_comments`` options. Corresponds to PyFITS ticket 181.

  - Fixed some bugs with WCS Paper IV record-valued keyword cards:

    - Cards that looked kind of like RVKCs but were not intended to be were
      over-permissively treated as such--commentary keywords like COMMENT and
      HISTORY were particularly affected. Corresponds to PyFITS ticket 183.

    - Looking up a card in a header by its standard FITS keyword only should
      always return the raw value of that card.  That way cards containing
      values that happen to valid RVKCs but were not intended to be will still
      be treated like normal cards. Corresponds to PyFITS ticket 184.

    - Looking up a RVKC in a header with only part of the field-specifier (for
      example "DP1.AXIS" instead of "DP1.AXIS.1") was implicitly treated as a
      wildcard lookup. Corresponds to PyFITS ticket 184.

  - Fixed a crash when diffing two FITS files where at least one contains a
    compressed image HDU which was not recognized as an image instead of a
    table. Corresponds to PyFITS ticket 187.

  - Fixed a bug where opening a file containing compressed image HDUs in
    'update' mode and then immediately closing it without making any changes
    caused the file to be rewritten unnecessarily.

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

  - Fixed bugs in the backwards compatibility layer for the ``CardList.index``
    and ``CardList.count`` methods. Corresponds to PyFITS ticket 190.

  - Improved ``__repr__`` and text file representation of cards with long
    values that are split into CONTINUE cards. Corresponds to PyFITS ticket
    193.

  - Fixed a crash when trying to assign a long (> 72 character) value to blank
    ('') keywords. This also changed how blank keywords are represented--there
    are still exactly 8 spaces before any commentary content can begin; this
    *may* affect the exact display of header cards that assumed there could be
    fewer spaces in a blank keyword card before the content begins. However,
    the current approach is more in line with the requirements of the FITS
    standard. Corresponds to PyFITS ticket 194.

- ``astropy.io.votable``

  - The `Table` class now maintains a single array object which is a
    Numpy masked array.  For variable-length columns, the object that
    is stored there is also a Numpy masked array.

  - Changed the ``pedantic`` configuration option to be ``False`` by default
    due to the vast proliferation of non-compliant VO Tables. [#296]

  - Renamed `astropy.io.vo` to `astropy.io.votable`.

- ``astropy.table``

  - Added a workaround for an upstream bug in Numpy 1.6.2 that could cause
    a maximum recursion depth RuntimeError when printing table rows. [#341]

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
