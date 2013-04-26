0.3 (unreleased)
----------------

- `astropy.io.votable`

  - The format of the units of a VOTable file can be specified using
    the `unit_format` parameter.  Note that units are still always
    written out using the CDS format, to ensure compatibility with the
    standard.

- `astropy.time`

  - Add ``datetime`` format (:class:`~astropy.time.core.TimeDatetime`) which
    allows converting to and from standard library `~datetime.datetime` objects.

  - Add ``plot_date`` format (:class:`~astropy.time.core.TimePlotDate`) which
    allows converting to and from the date representation used when plotting
    dates with matplotlib via the `~matplotlib.pyplot.plot_date` function.

- `astropy.stats`

  - Added robust statistics functions `~astropy.stats.funcs.median_absolute_deviation`, 
    `~astropy.stats.funcs.biweight_location`, and 
    `~astropy.stats.funcs.biweight_midvariance`.

0.2.2 (unreleased)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.io.fits``

  - Improved round-tripping and preservation of manually assigned column
    attributes (``TNULLn``, ``TSCALn``, etc.) in table HDU headers. [#996]

- ``astropy.units``

  - Fixed an issue where the ``isiterable()`` utility returned ``True`` for
    quantities with scalar values.  Added an ``__iter__`` method for the
    ``Quantity`` class and fixed ``isiterable()`` to catch false positives.
    [#878]

- Misc

  - Fixed an unrelated error message that could occur when trying to import
    astropy from a source checkout without having build the extension modules
    first. This issue was claimed to be fixed in v0.2.1, but the fix itself had
    a bug. [#971]


0.2.1 (2013-04-03)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - Fixed encoding errors that could occur when formatting coordinate objects
    in code using ``from __future__ import unicode_literals``. [#817]

  - Fixed a bug where the minus sign was dropped when string formatting dms
    coordinates with -0 degrees. [#875]

- ``astropy.io.fits``

  - Properly supports the ZQUANTIZ keyword used to support quantization
    level--this includes working support for lossless GZIP compression of
    images.

  - Fixed support for opening gzipped FITS files in a writeable mode. [#256]

  - Added a more helpful exception message when trying to read invalid values
    from a table when the required ``TNULLn`` keyword is missing. [#309]

  - More refactoring of the tile compression handling to work around a
    potential memory access violation that was particularly prevalent on
    Windows. [#507]

  - Fixed an integer size mismatch in the compression module that could affect
    32-bit systems. [#786]

  - Fixed malformatting of the ``TFORMn`` keywords when writing compressed
    image tables (they ommitted the max array length parameter from the
    variable-length array format).

  - Fixed a crash that could occur when writing a table containing multi-
    dimensional array columns from an existing file into a new file.

  - Fixed a bug in fitsdiff that reported two header keywords contaning NaN
    as having different values.

- ``astropy.io.votable``

  - Fixed links to the ``astropy.io.votable`` documentation in the VOTable
    validator output. [#806]

  - When reading VOTables containing integers that are out of range for their
    column type, display a warning rather than raising an exception. [#825]

  - Changed the default string format for floating point values for better
    round-tripping. [#856]

  - Fixed opening VOTables through the ``Table.read()`` interface for tables
    that have no names. [#927]

  - Fixed creation of VOTables from an Astropy table that does not have a data
    mask. [#928]

  - Minor documentation fixes. [#932]

- ``astropy.nddata.convolution``

  - Added better handling of ``inf`` values to the ``convolve_fft`` family of
    functions. [#893]

- ``astropy.table``

  - Fixed silent failure to assign values to a row on multiple columns. [#764]

  - Fixed various buggy behavior when viewing a table after sorting by one of
    its columns. [#829]

  - Fixed using ``numpy.where()`` with table indexing. [#838]

  - Fixed a bug where opening a remote table with ``Table.read()`` could cause
    the entire table to be downloaded twice. [#845]

  - Fixed a bug where ``MaskedColumn`` no longer worked if the column being
    masked is renamed. [#916]

- ``astropy.units``

  - Added missing capability for array ``Quantity``\s to be initializable by
    a list of ``Quantity``\s. [#835]

  - Fixed the definition of year and lightyear to be in terms of Julian year
    per the IAU definition. [#861]

  - "degree" was removed from the list of SI base units. [#863]

- ``astropy.wcs``

  - Fixed ``TypeError`` when calling ``WCS.to_header_string()``. [#822]

- Misc

  - Fixed a minor issue when installing with ``./setup.py develop`` on a fresh
    git clone.  This is likely only of interest to developers on Astropy.
    [#725]

  - Fixes a crash with ``ImportError: No module named 'astropy.version'`` when
    running setup.py from a source checkout for the first time on OSX with
    Python 3.3. [#820]

  - Fixed an installation issue where running ``./setup.py install`` or when
    installing with pip the ``.astropy`` directory gets created in the home
    directory of the user running the command.  The user's ``.astropy``
    directory should only be created when they use astropy, not when they
    install it. [#867]

  - Fixed an exception when creating a ``ProgressBar`` with a "total" of 0.
    [#752]

  - Added better documentation of behavior that can occur when trying to import
    the astropy package from within a source checkout without first building
    the extension modules. [#795, #864]

  - Added link to the installation instructions in the README. [#797]

  - Catches segfaults in xmllint which can occur sometimes and is otherwise out
    of our control. [#803]

  - Minor changes to the documentation template. [#805]

  - Fixed a minor exception handling bug in ``download_file()``. [#808]

  - Added cleanup of any temporary files if an error occurs in
    ``download_file()``. [#857]

  - Filesystem free space is checked for before attempting to download a file
    with ``download_file()``. [#858]

  - Fixed package data locating to work across symlinks--required to work with
    some OS packaging layouts. [#827]

  - Fixed a bug when building Cython extensions where hidden files containing
    ``.pyx`` extensions could cause the build to crash. This can be an issue
    with software and filesystems that autogenerate hidden files. [#834]

  - Fixed bug that could cause a "script" called README.rst to be installed
    in a bin directory. [#852]

  - Fixed some miscellaneous and mostly rare reference leaks caught by
    cpychecker. [#914]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added logo and branding for Windows binary installers. [#741]

- Upgraded included version libexpat to 2.1.0. [#781]

- ~25% performance improvement in unit composition/decomposition. [#836]

- Added previously missing LaTeX formatting for ``L_sun`` and ``R_sun``. [#841]

- ``ConfigurationItem``\s now have a more useful and informative ``__repr__``
  and improved documentation for how to use them. [#855]

- Added a friendlier error message when trying to import astropy from a source
  checkout without first building the extension modules inplace. [#864]

- ``py.test`` now outputs more system information for help in debugging issues
  from users. [#869]

- Added unit definitions "mas" and "uas" for "milliarcsecond" and
  "microarcsecond" respectively. [#892]


0.2 (2013-02-19)
----------------

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

  - Add a WMAP9 object using the 9-year WMAP parameters from Hinshaw et al.
    Once this paper is accepted, this should be made the default, but for
    now WMAP7 remains the default. [#629]

- ``astropy.table`` I/O infrastructure for custom readers/writers
  implemented. [#305]

  - Added support for reading/writing HDF5 files [#461]

  - Added support for masked tables with missing or invalid data [#451]

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

- Generalized I/O infrastructure so that ``astropy.nddata`` can also have
  custom readers/writers [#659]

- ``astropy.wcs``

  - From updating the the underlying wcslib 4.16:

    - When `astropy.wcs.WCS` constructs a default coordinate
      representation it will give it the special name "DEFAULTS", and
      will not report "Found one coordinate representation".

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- A configuration file with all options set to their defaults is now generated 
  when astropy is installed.  This file will be pulled in as the users' 
  astropy configuration file the first time they ``import astropy``.  [#498] 

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

  - Uses the ``relax`` option by default so that non-standard keywords are
    allowed. [#585]


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

- Added a classmethod to 
  `astropy.coordinates.coordsystems.SphericalCoordinatesBase` that performs a 
  name resolve query using Sesame to retrieve coordinates for the requested
  object. This works for any subclass of `SphericalCoordinatesBase`, but 
  requires an internet connection. [#556]

- ``astropy.nddata.convolution`` removed requirement of PyFFTW3; uses Numpy's
  FFT by default instead with the added ability to specify an FFT
  implementation to use. [#660]


Bug Fixes
^^^^^^^^^

- ``astropy.io.ascii``

  - Fixed crash when pprinting a row with INDEF values. [#511]

  - Fixed failure when reading DAOphot files with empty keyword values [#666].

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

- Improved multiprocessing compatibility for file downloads. [#615]

- Fixed handling of Cython modules when building from a source checkout of a
  tagged release version. [#594]

- Added a workaround for a bug in Sphinx that could occur when using the
  ``:tocdepth:`` directive. [#595]

- Minor VOTable fixes [#596]

- Fixed how ``setup.py`` uses ``distribute_setup.py`` to prevent possible
  ``VersionConflict`` errors when an older version of distribute is already
  installed on the user's system. [#616][#640]

- Changed use of ``log.warn`` in the logging module to ``log.warning`` since
  the former is deprecated. [#624]


0.1 (2012-06-19)
----------------

- Initial release.
