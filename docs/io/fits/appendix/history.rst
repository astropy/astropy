.. doctest-skip-all

astropy.io.fits History
=======================

Prior to its inclusion in Astropy, the `astropy.io.fits` package was a stand-
alone package called `PyFITS`_.  Though for the time being active development
is continuing on PyFITS, that development is also being merged into Astropy.
This page documents the release history of PyFITS prior to its merge into
Astropy.

.. contents:: PyFITS Changelog
   :depth: 2
   :local:


3.3.0 (unreleased)
------------------

New Features
^^^^^^^^^^^^

- Added new verification options ``fix+ignore``, ``fix+warn``,
  ``fix+exception``, ``silentfix+ignore``, ``silentfix+warn``, and
  ``silentfix+exception`` which give more control over how to report fixable
  errors as opposed to unfixable errors.  See the "Verification" section in
  the PyFITS documentation for more details.

API Changes
^^^^^^^^^^^

- The ``pyfits.new_table`` function is now fully deprecated (though will not
  be removed for a long time, considering how widely it is used).

  Instead please use the more explicit ``pyfits.BinTableHDU.from_columns`` to
  create a new binary table HDU, and the similar
  ``pyfits.TableHDU.from_columns`` to create a new ASCII table.  These
  otherwise accept the same arguments as ``pyfits.new_table`` which is now
  just a wrapper for these.

- The ``.fromstring`` classmethod of each HDU type has been simplified such
  that, true to its namesake, it only initializes an HDU from a string
  containing its header *and* data. (spacetelescope/PyFITS#64)

- Fixed an issue where header wildcard matching (for example
  ``header['DATE*']``) can be used to match *any* characters that might appear
  in a keyword.  Previously this only matched keywords containing characters
  in the set ``[0-9A-Za-z_]``.  Now this can also match a hyphen ``-`` and any
  other characters, as some conventions like ``HIERARCH`` and record-valued
  keyword cards allow a wider range of valid characters than standard FITS
  keywords.

- This will be the *last* release to support the following APIs that have been
  marked deprecated since PyFITS v3.1:

  - The ``CardList`` class, which was part of the old header implementation.

  - The ``Card.key`` attribute.  Use ``Card.keyword`` instead.

  - The ``Card.cardimage`` and ``Card.ascardimage`` attributes.  Use simply
    ``Card.image`` or ``str(card)`` instead.

  - The ``create_card`` factory function.  Simply use the normal ``Card``
    constructor instead.

  - The ``create_card_from_string`` factory function.  Use ``Card.fromstring``
    instead.

  - The ``upper_key`` function.  Use ``Card.normalize_keyword`` method instead
    (this is not unlikely to be used outside of PyFITS itself, but it was
    technically public API).

  - The usage of ``Header.update`` with ``Header.update(keyword, value,
    comment)`` arguments.  ``Header.update`` should only be used analogously
    to ``dict.update``.  Use ``Header.set`` instead.

  - The ``Header.ascard`` attribute.  Use ``Header.cards`` instead for a list
    of all the ``Card`` objects in the header.

  - The ``Header.rename_key`` method.  Use ``Header.rename_keyword`` instead.

  - The ``Header.get_history`` method.  Use ``header['HISTORY']`` instead
    (normal keyword lookup).

  - The ``Header.get_comment`` method.  Use ``header['COMMENT']`` instead.

  - The ``Header.toTxtFile`` method.  Use ``header.totextfile`` instead.

  - The ``Header.fromTxtFile`` method.  Use ``Header.fromtextfile`` instead.

  - The ``pyfits.tdump`` and ``tcreate`` functions.  Use ``pyfits.tabledump``
    and ``pyfits.tableload`` respectively.

  - The ``BinTableHDU.tdump`` and ``tcreate`` methods.  Use
    ``BinTableHDU.dump`` and ``BinTableHDU.load`` respectively.

  - The ``txtfile`` argument to the ``Header`` constructor.  Use
    ``Header.fromfile`` instead.

  - The ``startColumn`` and ``endColumn`` arguments to the ``FITS_record``
    constructor.  These are unlikely to be used by any user code.

  These deprecated interfaces will be removed from the development version of
  PyFITS following the v3.3 release (they will still be available in any
  v3.3.x bugfix releases, however).

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- PyFITS has switched to a unified code base which supports Python 2.5 through
  3.4 simultaneously without translation.  This *shouldn't* have any
  significant performance impacts, but please report if anything seems
  noticeably slower.  As a reminder, support for Python 2.5 will be ended
  after PyFITS 3.3.x.

- Warnings for deprecated APIs in PyFITS are now always displayed by default.
  This is in line with a similar change made recently to Astropy:
  https://github.com/astropy/astropy/pull/1871
  To disable PyFITS deprecation warnings in scripts one may call
  ``pyfits.ignore_deprecation_warnings()`` after importing PyFITS.

- ``Card`` objects have a new ``is_blank`` attribute which returns ``True`` if
  the card represents a blank card (no keyword, value, or comment) and
  ``False`` otherwise.

Bug Fixes
^^^^^^^^^

- Fixed a regression where it was not possible to save an empty "compressed"
  image to a file (in this case there is nothing to compress, hence the
  quotes, but trying to do so caused a crash). (spacetelescope/PyFITS#69)

- Fixed a regression that may have been introduced in v3.2.1 with writing
  compressed image HDUs, particularly compressed images using a non-empty
  GZIP_COMPRESSED_DATA column. (spacetelescope/#71)


3.2.4 (unreleased)
------------------

- Fixed a regression where multiple consecutive calls of the ``writeto``
  method on the same HDU but to different files could lead to corrupt data or
  crashes on the subsequent calls after the first. (spacetelescope/PyFITS#40)


3.1.7 (unreleased)
------------------

- Nothing changed yet.


3.2.3 (2014-05-14)
------------------

- Nominal support for Python 3.4.

- Fixed a bug with using the ``tabledump`` and ``tableload`` functions with
  tables containing array columns (columns in which each element is an array
  instead of a single scalar value). (spacetelescope/PyFITS#22)

- Fixed an issue where PyFITS allowed newline characters in header values and
  comments. (spacetelescope/PyFITS#51)

- Fixed pickling of ``FITS_rec`` (table data) objects.
  (spacetelescope/PyFITS#53)

- Improved behavior when writing large compressed images on OSX by removing an
  unnecessary check for platform architecture. (spacetelescope/PyFITS#57)

- Allow reading FITS files from file-like objects that do not have a
  ``.closed`` attribute (and as such may not even have an "open" vs. "closed"
  concept). (spacetelescope/PyFITS#56)

- Fixed duplicate insertion of commentary keywords on compressed image
  headers. (spacetelescope/PyFITS#58)

- Fixed minor issue with comparison of header commentary card values.
  (spacetelescope/PyFITS#59)


3.1.6 (2014-05-14)
------------------

- Nominal support for Python 3.4.

- Fixed a bug with using the ``tabledump`` and ``tableload`` functions with
  tables containing array columns (columns in which each element is an array
  instead of a single scalar value). (Backported from 3.2.3)

- Fixed an issue where PyFITS allowed newline characters in header values and
  comments. (Backported from 3.2.3)

- Fixed pickling of ``FITS_rec`` (table data) objects.
  (Backported from 3.2.3)

- Improved behavior when writing large compressed images on OSX by removing an
  unnecessary check for platform architecture. (Backported from 3.2.3)

- Allow reading FITS files from file-like objects that do not have a
  ``.closed`` attribute (and as such may not even have an "open" vs. "closed"
  concept). (Backported from 3.2.3)

- Fixed minor issue with comparison of header commentary card values.
  (Backported from 3.2.3)


3.2.2 (2014-03-25)
------------------

- Fixed a regression on deletion of record-valued keyword cards using
  the Header wildcard syntax.  This was intended to be fixed before the
  v3.2.1 release.


3.1.5 (2014-03-25)
------------------

- Fixed a regression on deletion of record-valued keyword cards using
  the Header wildcard syntax.  This was intended to be fixed before the
  v3.1.4 release.


3.2.1 (2014-03-04)
------------------

- Nominal support for the upcoming Python 3.4.

- Added missing features from the ``Header.insert()`` method that were
  intended for inclusion in the original 3.1 release:  In addition to
  accepting an integer index as the first argument, it also supports supplying
  a keyword name as the first argument for insertion relative to a specific
  keyword.  It also now supports an optional ``after`` argument.  If
  ``after=True`` the insertion is made below the insertion point instead
  of above it. (spacetelescope/PyFITS#12)

- Fixed support for broadcasting of values assigned to table columns.
  (spacetelescope/PyFITS#48)

- A grab bag of minor performance improvements in headers.
  (spacetelescope/PyFITS#46)

- Fix an unrelated error that occurred when instantiating a ``ColDefs`` object
  with invalid input.

- Fixed an issue where opening an image containing pseudo-unsigned integers
  and immediately writing it to a new file using the ``writeto`` method would
  drop the scale factors that identified the data as unsigned.

- Fixed a bug where writing a file with ``checksum=True`` did not add the
  checksum on new files. (spacetelescope/PyFITS#8)

- Fixed an issue where validating an HDU's checksums removed the checksum from
  that HDU's header entirely (even if it was valid.)

- Fixed checksums on compressed images, so that the ``ZHECKSUM`` and
  ``ZDATASUM`` contain a checksum of the original image HDU, while
  ``CHECKSUM`` and ``DATASUM`` contain checksums of the compressed image HDU.
  This feature was supposed to be supported in 3.2, but the support was buggy.

- Fixed an issue where the size of the heap was sometimes not computed
  properly when writing an existing table containing variable-length array
  columns to a new FITS file.  This could result in corruption in the new FITS
  file. (spacetelescope/PyFITS#47)

- Fixed issue with updates to the header of ``CompImageHDU`` objects not being
  preserved on save. (spacetelescope/PyFITS#23)

- Fixed a bug where a boolean value of ``True`` in a header could not be
  replaced with the integer 1, and likewise for ``False`` and 0 and vice
  versa.

- Fixed an issue similar to the above one but for numeric values--now
  replacing a header value with an equivalent numeric value will up/downcast
  that value.  For example replacing '0' with '0.0' will write '0.0' to the
  header so that it is returned as a floating point value.  Likewise a float
  can be downcast to an integer. (spacetelescope/PyFITS#49)

- A handful of Python 3 compatibility fixes, especially for compatibility
  with the upcoming Python 3.4.

- Fixed unrelated crash when a header contains an invalid END card (for
  example "END = ").  This resulted in a cryptic traceback.  Now headers like
  this will detect "clearly intended" END cards and produce a warning about
  their invalidity and fix them. (#217)

- Allowed a sequence of ``Column`` objects to be passed in as the main
  argument to ``FITS_rec.from_columns`` as the documentation suggests should
  be possible.

- Fixed a display formatting issue with fitsdiff where sometimes it did not
  show the difference between two floating point numbers if they were the same
  up to some low number of digits. (spacetelescope/PyFITS#21)

- Fixed an issue where Python 2 sometimes allowed non-ASCII strings to be
  assigned as header values if they were assigned as old-style ``str`` objects
  and not ``unicode`` objects. (spacetelescope/PyFITS#37)


3.1.4 (2014-03-04)
------------------

- Added missing features from the ``Header.insert()`` method that were
  intended for inclusion in the original 3.1 release:  In addition to
  accepting an integer index as the first argument, it also supports supplying
  a keyword name as the first argument for insertion relative to a specific
  keyword.  It also now supports an optional ``after`` argument.  If
  ``after=True`` the insertion is made below the insertion point instead
  of above it. (Backported from 3.2.1)

- A grab bag of minor performance improvements in headers.
  (Backported from 3.2.1)

- Fixed an issue where opening an image containing pseudo-unsigned integers
  and immediately writing it to a new file using the ``writeto`` method would
  drop the scale factors that identified the data as unsigned.
  (Backported from 3.2.1)

- Fixed a bug where writing a file with ``checksum=True`` did not add the
  checksum on new files. (Backported from 3.2.1)

- Fixed an issue where validating an HDU's checksums removed the checksum from
  that HDU's header entirely (even if it was valid.)
  (Backported from 3.2.1)

- Fixed an issue where the size of the heap was sometimes not computed
  properly when writing an existing table containing variable-length array
  columns to a new FITS file.  This could result in corruption in the new FITS
  file. (Backported from 3.2.1)

- Fixed a bug where a boolean value of ``True`` in a header could not be
  replaced with the integer 1, and likewise for ``False`` and 0 and vice
  versa. (Backported from 3.2.1)

- Fixed an issue similar to the above one but for numeric values--now
  replacing a header value with an equivalent numeric value will up/downcast
  that value.  For example replacing '0' with '0.0' will write '0.0' to the
  header so that it is returned as a floating point value.  Likewise a float
  can be downcast to an integer. (Backported from 3.2.1)

- Fixed unrelated crash when a header contains an invalid END card (for
  example "END = ").  This resulted in a cryptic traceback.  Now headers like
  this will detect "clearly intended" END cards and produce a warning about
  their invalidity and fix them. (Backported from 3.2.1)

- Fixed a display formatting issue with fitsdiff where sometimes it did not
  show the difference between two floating point numbers if they were the same
  up to some low number of digits. (Backported from 3.2.1)

- Fixed an issue where Python 2 sometimes allowed non-ASCII strings to be
  assigned as header values if they were assigned as old-style ``str`` objects
  and not ``unicode`` objects. (Backported from 3.2.1)


3.0.13 (2014-03-04)
-------------------

- Fixed a bug where writing a file with ``checksum=True`` did not add the
  checksum on new files. (Backported from 3.2.1)

- Fixed an issue where validating an HDU's checksums removed the checksum from
  that HDU's header entirely (even if it was valid.)
  (Backported from 3.2.1)


3.2 (2013-11-26)
----------------

Highlights
^^^^^^^^^^

- Rewrote CFITSIO-based backend for handling tile compression of FITS files.
  It now uses a standard CFITSIO instead of heavily modified pieces of CFITSIO
  as before.  PyFITS ships with its own copy of CFITSIO v3.35 which supports
  the latest version of the Tiled Image Convention (v2.3), but system
  packagers may choose instead to strip this out in favor of a
  system-installed version of CFITSIO.  Earlier versions may work, but nothing
  earlier than 3.28 has been tested yet. (#169)

- Added support for reading and writing tables using the Q format for columns.
  The Q format is identical to the P format (variable-length arrays) except
  that it uses 64-bit integers for the data descriptors, allowing more than
  4 GB of variable-length array data in a single table. (#160)

- Added initial support for table columns containing pseudo-unsigned integers.
  This is currently enabled by using the ``uint=True`` option when opening
  files; any table columns with the correct BZERO value will be interpreted
  and returned as arrays of unsigned integers.

- Some refactoring of the table and ``FITS_rec`` modules in order to better
  separate the details of the FITS binary and ASCII table data structures from
  the HDU data structures that encapsulate them.  Most of these changes should
  not be apparent to users (but see API Changes below).


API Changes
^^^^^^^^^^^

- Assigning to values in ``ColDefs.names``, ``ColDefs.formats``,
  ``ColDefs.nulls`` and other attributes of ``ColDefs`` instances that return
  lists of column properties is no longer supported.  Assigning to those lists
  will no longer update the corresponding columns.  Instead, please just
  modify the ``Column`` instances directly (``Column.name``, ``Column.null``,
  etc.)

- The ``pyfits.new_table`` function is marked "pending deprecation".  This
  does not mean it will be removed outright or that its functionality has
  changed.  It will likely be replaced in the future for a function with
  similar, if not subtly different functionality.  A better, if not slightly
  more verbose approach is to use ``pyfits.FITS_rec.from_columns`` to create
  a new ``FITS_rec`` table--this has the same interface as
  ``pyfits.new_table``.  The difference is that it returns a plan ``FITS_rec``
  array, and not an HDU instance.  This ``FITS_rec`` object can then be used
  as the data argument in the constructors for ``BinTableHDU`` (for binary
  tables) or ``TableHDU`` (for ASCII tables).  This is analogous to creating
  an ``ImageHDU`` by passing in an image array.
  ``pyfits.FITS_rec.from_columns`` is just a simpler way of creating a
  FITS-compatible recarray from a FITS column specification.

- The ``updateHeader``, ``updateHeaderData``, and ``updateCompressedData``
  methods of the ``CompDataHDU`` class are pending deprecation and moved to
  internal methods.  The operation of these methods depended too much on
  internal state to be used safely by users; instead they are invoked
  automatically in the appropriate places when reading/writing compressed image
  HDUs.

- The ``CompDataHDU.compData`` attribute is pending deprecation in favor of
  the clearer and more PEP-8 compatible ``CompDataHDU.compressed_data``.

- The constructor for ``CompDataHDU`` has been changed to accept new keyword
  arguments.  The new keyword arguments are essentially the same, but are in
  underscore_separated format rather than camelCase format.  The old arguments
  are still pending deprecation.

- The internal attributes of HDU classes ``_hdrLoc``, ``_datLoc``, and
  ``_datSpan`` have been replaced with ``_header_offset``, ``_data_offset``,
  and ``_data_size`` respectively.  The old attribute names are still pending
  deprecation.  This should only be of interest to advanced users who have
  created their own HDU subclasses.

- The following previously deprecated functions and methods have been removed
  entirely: ``createCard``, ``createCardFromString``, ``upperKey``,
  ``ColDefs.data``, ``setExtensionNameCaseSensitive``, ``_File.getfile``,
  ``_TableBaseHDU.get_coldefs``, ``Header.has_key``, ``Header.ascardlist``.

  If you run your code with a previous version of PyFITS (>= 3.0, < 3.2) with
  the ``python -Wd`` argument, warnings for all deprecated interfaces still in
  use will be displayed.

- Interfaces that were pending deprecation are now fully deprecated.  These
  include: ``create_card``, ``create_card_from_string``, ``upper_key``,
  ``Header.get_history``, and ``Header.get_comment``.

- The ``.name`` attribute on HDUs is now directly tied to the HDU's header, so
  that if ``.header['EXTNAME']`` changes so does ``.name`` and vice-versa.

- The ``pyfits.file.PYTHON_MODES`` constant dict was renamed to
  ``pyfits.file.PYFITS_MODES`` which better reflects its purpose.  This is
  rarely used by client code, however.  Support for the old name will be
  removed by PyFITS 3.4.


Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- The new compression code also adds support for the ZQUANTIZ and ZDITHER0
  keywords added in more recent versions of this FITS Tile Compression spec.
  This includes support for lossless compression with GZIP. (#198) By default
  no dithering is used, but the ``SUBTRACTIVE_DITHER_1`` and
  ``SUBTRACTIVE_DITHER_2`` methods can be enabled by passing the correct
  constants to the ``quantize_method`` argument to the ``CompImageHDU``
  constructor.  A seed can be manually specified, or automatically generated
  using either the system clock or checksum-based methods via the
  ``dither_seed`` argument.  See the documentation for ``CompImageHDU`` for
  more details. (#198) (spacetelescope/PYFITS#32)

- Images compressed with the Tile Compression standard can now be larger than
  4 GB through support of the Q format. (#159)

- All HDUs now have a ``.ver`` ``.level`` attribute that returns the value of
  the EXTVAL and EXTLEVEL keywords from that HDU's header, if the exist.  This
  was added for consistency with the ``.name`` attribute which returns the
  EXTNAME value from the header.

- Then ``Column`` and ``ColDefs`` classes have new ``.dtype`` attributes
  which give the Numpy dtype for the column data in the first case, and the
  full Numpy compound dtype for each table row in the latter case.

- There was an issue where new tables created defaulted the values in all
  string columns to '0.0'.  Now string columns are filled with empty strings
  by default--this seems a less surprising default, but it may cause
  differences with tables created with older versions of PyFITS.

- Improved round-tripping and preservation of manually assigned column
  attributes (``TNULLn``, ``TSCALn``, etc.) in table HDU headers.
  (astropy/astropy#996)


Bug Fixes
^^^^^^^^^

- Binary tables containing compressed images may, optionally, contain other
  columns unrelated to the tile compression convention. Although this is an
  uncommon use case, it is permitted by the standard. (#159)

- Reworked some of the file I/O routines to allow simpler, more consistent
  mapping between OS-level file modes ('rb', 'wb', 'ab', etc.) and the more
  "PyFITS-specific" modes used by PyFITS like "readonly" and "update".
  That is, if reading a FITS file from an open file object, it doesn't matter
  as much what "mode" it was opened in so long as it has the right
  capabilities (read/write/etc.)  Also works around bugs in the Python io
  module in 2.6+ with regard to file modes. (spacetelescope/PyFITS#33)

- Fixed an obscure issue that can occur on systems that don't have flush to
  memory-mapped files implemented (namely GNU Hurd). (astropy/astropy#968)


3.1.3 (2013-11-26)
------------------

- Disallowed assigning NaN and Inf floating point values as header values,
  since the FITS standard does not define a way to represent them in. Because
  this is undefined, the previous behavior did not make sense and produced
  invalid FITS files. (spacetelescope/PyFITS#11)

- Added a workaround for a bug in 64-bit OSX that could cause truncation when
  writing files greater than 2^32 bytes in size. (spacetelescope/PyFITS#28)

- Fixed a long-standing issue where writing binary tables did not correctly
  write the TFORMn keywords for variable-length array columns (they ommitted
  the max array length parameter of the format).  This was thought fixed in
  v3.1.2, but it was only fixed there for compressed image HDUs and not for
  binary tables in general.

- Fixed an obscure issue that can occur on systems that don't have flush to
  memory-mapped files implemented (namely GNU Hurd). (Backported from 3.2)


3.0.12 (2013-11-26)
-------------------

- Disallowed assigning NaN and Inf floating point values as header values,
  since the FITS standard does not define a way to represent them in. Because
  this is undefined, the previous behavior did not make sense and produced
  invalid FITS files. (Backported from 3.1.3)

- Added a workaround for a bug in 64-bit OSX that could cause truncation when
  writing files greater than 2^32 bytes in size. (Backported from 3.1.3)

- Fixed a long-standing issue where writing binary tables did not correctly
  write the TFORMn keywords for variable-length array columns (they ommitted
  the max array length parameter of the format).  This was thought fixed in
  v3.1.2, but it was only fixed there for compressed image HDUs and not for
  binary tables in general. (Backported from 3.1.3)

- Fixed an obscure issue that can occur on systems that don't have flush to
  memory-mapped files implemented (namely GNU Hurd). (Backported from 3.2)


3.1.3 (unreleased)
------------------

- Disallowed assigning NaN and Inf floating point values as header values,
  since the FITS standard does not define a way to represent them in. Because
  this is undefined, the previous behavior did not make sense and produced
  invalid FITS files. (spacetelescope/PyFITS#11)


3.0.12 (unreleased)
-------------------

- Disallowed assigning NaN and Inf floating point values as header values,
  since the FITS standard does not define a way to represent them in. Because
  this is undefined, the previous behavior did not make sense and produced
  invalid FITS files. (Backported from 3.1.3)

- Added a workaround for a bug in 64-bit OSX that could cause truncation when
  writing files greater than 2^32 bytes in size. (Backported from 3.1.3)


3.1.2 (2013-04-22)
------------------

- When an error occurs opening a file in fitsdiff the exception message will
  now at least mention which file had the error. (#168)

- Fixed support for opening gzipped FITS files by filename in a writeable mode
  (PyFITS has supported writing to gzip files for some time now, but only
  enabled it when GzipFile objects were passed to ``pyfits.open()`` due to
  some legacy code preventing full gzip support. (#195)

- Added a more helpful error message in the case of malformatted FITS files
  that contain non-float NULL values in an ASCII table but are missing the
  required TNULLn keywords in the header. (#197)

- Fixed an (apparently long-standing) issue where writing compressed images
  did not correctly write the TFORMn keywords for variable-length array
  columns (they omitted the max array length parameter of the format). (#199)

- Slightly refactored how tables containing variable-length array columns are
  handled to add two improvements: Fixes an issue where accessing the data
  after a call to the ``pyfits.getdata`` convenience function caused an
  exception, and allows the VLA data to be read from an existing mmap of the
  FITS file. (#200)

- Fixed a bug that could occur when opening a table containing
  multi-dimensional columns (i.e. via the TDIMn keyword) and then writing it
  out to a new file. (#201)

- Added use of the console_scripts entry point to install the fitsdiff and
  fitscheck scripts, which if nothing else provides better Windows support.
  The generated scripts now override the ones explicitly defined in the
  scripts/ directory (which were just trivial stubs to begin with). (#202)

- Fixed a bug on Python 3 where attempting to open a non-existent file on
  Python 3 caused a seemingly unrelated traceback. (#203)

- Fixed a bug in fitsdiff that reported two header keywords containing NaN
  as value as different. (#204)

- Fixed an issue in the tests that caused some tests to fail if pyfits is
  installed with read-only permissions. (#208)

- Fixed a bug where instantiating a ``BinTableHDU`` from a numpy array
  containing boolean fields converted all the values to ``False``. (#215)

- Fixed an issue where passing an array of integers into the constructor of
  ``Column()`` when the column type is floats of the same byte width caused the
  column array to become garbled. (#218)

- Fixed inconsistent behavior in creating CONTINUE cards from byte strings
  versus Unicode strings in Python 2--CONTINUE cards can now be created
  properly from Unicode strings (so long as they are convertible to ASCII).
  (spacetelescope/PyFITS#1)

- Fixed a couple cases where creating a new table using TDIMn in some of the
  columns could caused a crash. (spacetelescope/PyFITS#3)

- Fixed a bug in parsing HIERARCH keywords that do not have a space after
  the first equals sign (before the value). (spacetelescope/PyFITS#5)

- Prevented extra leading whitespace on HIERARCH keywords from being treated
  as part of the keyword. (spacetelescope/PyFITS#6)

- Fixed a bug where HIERARCH keywords containing lower-case letters was
  mistakenly marked as invalid during header validation.
  (spacetelescope/PyFITS#7)

- Fixed an issue that was ancillary to (spacetelescope/PyFITS#7) where the
  ``Header.index()`` method did not work correctly with HIERARCH keywords
  containing lower-case letters.


3.0.11 (2013-04-17)
-------------------

- Fixed support for opening gzipped FITS files by filename in a writeable mode
  (PyFITS has supported writing to gzip files for some time now, but only
  enabled it when GzipFile objects were passed to ``pyfits.open()`` due to
  some legacy code preventing full gzip support. Backported from 3.1.2. (#195)

- Added a more helpful error message in the case of malformatted FITS files
  that contain non-float NULL values in an ASCII table but are missing the
  required TNULLn keywords in the header. Backported from 3.1.2. (#197)

- Fixed an (apparently long-standing) issue where writing compressed images did
  not correctly write the TFORMn keywords for variable-length array columns
  (they ommitted the max array length parameter of the format). Backported from
  3.1.2. (#199)

- Slightly refactored how tables containing variable-length array columns are
  handled to add two improvements: Fixes an issue where accessing the data
  after a call to the ``pyfits.getdata`` convenience function caused an
  exception, and allows the VLA data to be read from an existing mmap of the
  FITS file. Backported from 3.1.2. (#200)

- Fixed a bug that could occur when opening a table containing
  multi-dimensional columns (i.e. via the TDIMn keyword) and then writing it
  out to a new file. Backported from 3.1.2. (#201)

- Fixed a bug on Python 3 where attempting to open a non-existent file on
  Python 3 caused a seemingly unrelated traceback. Backported from 3.1.2.
  (#203)

- Fixed a bug in fitsdiff that reported two header keywords containing NaN
  as value as different. Backported from 3.1.2. (#204)

- Fixed an issue in the tests that caused some tests to fail if pyfits is
  installed with read-only permissions. Backported from 3.1.2. (#208)

- Fixed a bug where instantiating a ``BinTableHDU`` from a numpy array
  containing boolean fields converted all the values to ``False``. Backported
  from 3.1.2. (#215)

- Fixed an issue where passing an array of integers into the constructor of
  ``Column()`` when the column type is floats of the same byte width caused the
  column array to become garbled. Backported from 3.1.2. (#218)

- Fixed a couple cases where creating a new table using TDIMn in some of the
  columns could caused a crash. Backported from 3.1.2.
  (spacetelescope/PyFITS#3)


3.1.1 (2013-01-02)
------------------

This is a bug fix release for the 3.1.x series.

Bug Fixes
^^^^^^^^^

- Improved handling of scaled images and pseudo-unsigned integer images in
  compressed image HDUs.  They now work more transparently like normal image
  HDUs with support for the ``do_not_scale_image_data`` and ``uint`` options,
  as well as ``scale_back`` and ``save_backup``.  The ``.scale()`` method
  works better too. (#88)

- Permits non-string values for the EXTNAME keyword when reading in a file,
  rather than throwing an exception due to the malformatting.  Added
  verification for the format of the EXTNAME keyword when writing. (#96)

- Added support for EXTNAME and EXTVER in PRIMARY HDUs.  That is, if EXTNAME
  is specified in the header, it will also be reflected in the ``.name``
  attribute and in ``pyfits.info()``.  These keywords used to be verboten in
  PRIMARY HDUs, but the latest version of the FITS standard allows them.
  (#151)

- HCOMPRESS can again be used to compress data cubes (and higher-dimensional
  arrays) so long as the tile size is effectively 2-dimensional. In fact,
  PyFITS will automatically use compatible tile sizes even if they're not
  explicitly specified. (#171)

- Added support for the optional ``endcard`` parameter in the
  ``Header.fromtextfile()`` and ``Header.totextfile()`` methods.  Although
  ``endcard=False`` was a reasonable default assumption, there are still text
  dumps of FITS headers that include the END card, so this should have been
  more flexible. (#176)

- Fixed a crash when running fitsdiff on two empty (that is, zero row) tables.
  (#178)

- Fixed an issue where opening files containing random groups HDUs in update
  mode could cause an unnecessary rewrite of the file even if none of the
  data is modified. (#179)

- Fixed a bug that could caused a deadlock in the filesystem on OSX if PyFITS
  is used with Numpy 1.7 in some cases. (#180)

- Fixed a crash when generating diff reports from diffs using the
  ``ignore_comments`` options. (#181)

- Fixed some bugs with WCS Paper IV record-valued keyword cards:

  - Cards that looked kind of like RVKCs but were not intended to be were
    over-permissively treated as such--commentary keywords like COMMENT and
    HISTORY were particularly affected. (#183)

  - Looking up a card in a header by its standard FITS keyword only should
    always return the raw value of that card.  That way cards containing
    values that happen to valid RVKCs but were not intended to be will still
    be treated like normal cards. (#184)

  - Looking up a RVKC in a header with only part of the field-specifier (for
    example "DP1.AXIS" instead of "DP1.AXIS.1") was implicitly treated as a
    wildcard lookup. (#184)

- Fixed a crash when diffing two FITS files where at least one contains a
  compressed image HDU which was not recognized as an image instead of a
  table. (#187)

- Fixed bugs in the backwards compatibility layer for the ``CardList.index``
  and ``CardList.count`` methods. (#190)

- Improved ``__repr__`` and text file representation of cards with long values
  that are split into CONTINUE cards. (#193)

- Fixed a crash when trying to assign a long (> 72 character) value to blank
  ('') keywords. This also changed how blank keywords are represented--there
  are still exactly 8 spaces before any commentary content can begin; this
  *may* affect the exact display of header cards that assumed there could be
  fewer spaces in a blank keyword card before the content begins. However, the
  current approach is more in line with the requirements of the FITS standard.
  (#194)


3.0.10 (2013-01-02)
-------------------

- Improved handling of scaled images and pseudo-unsigned integer images in
  compressed image HDUs.  They now work more transparently like normal image
  HDUs with support for the ``do_not_scale_image_data`` and ``uint`` options,
  as well as ``scale_back`` and ``save_backup``.  The ``.scale()`` method
  works better too.  Backported from 3.1.1. (#88)

- Permits non-string values for the EXTNAME keyword when reading in a file,
  rather than throwing an exception due to the malformatting.  Added
  verification for the format of the EXTNAME keyword when writing.  Backported
  from 3.1.1. (#96)

- Added support for EXTNAME and EXTVER in PRIMARY HDUs.  That is, if EXTNAME
  is specified in the header, it will also be reflected in the ``.name``
  attribute and in ``pyfits.info()``.  These keywords used to be verbotten in
  PRIMARY HDUs, but the latest version of the FITS standard allows them.
  Backported from 3.1.1. (#151)

- HCOMPRESS can again be used to compress data cubes (and higher-dimensional
  arrays) so long as the tile size is effectively 2-dimensional. In fact,
  PyFITS will not automatically use compatible tile sizes even if they're not
  explicitly specified.  Backported from 3.1.1. (#171)

- Fixed a bug when writing out files containing zero-width table columns,
  where the TFIELDS keyword would be updated incorrectly, leaving the table
  largely unreadable.  Backported from 3.1.0. (#174)

- Fixed an issue where opening files containing random groups HDUs in update
  mode could cause an unnecessary rewrite of the file even if none of the
  data is modified.  Backported from 3.1.1. (#179)

- Fixed a bug that could caused a deadlock in the filesystem on OSX if PyFITS
  is used with Numpy 1.7 in some cases. Backported from 3.1.1. (#180)


3.1 (2012-08-08)
----------------

Highlights
^^^^^^^^^^

- The ``Header`` object has been significantly reworked, and ``CardList``
  objects are now deprecated (their functionality folded into the ``Header``
  class).  See API Changes below for more details.

- Memory maps are now used by default to access HDU data.  See API Changes
  below for more details.

- Now includes a new version of the ``fitsdiff`` program for comparing two
  FITS files, and a new FITS comparison API used by ``fitsdiff``.  See New
  Features below.

API Changes
^^^^^^^^^^^

- The ``Header`` class has been rewritten, and the ``CardList`` class is
  deprecated.  Most of the basic details of working with FITS headers are
  unchanged, and will not be noticed by most users.  But there are differences
  in some areas that will be of interest to advanced users, and to application
  developers.  For full details of the changes, see the "Header Interface
  Transition Guide" section in the PyFITS documentation.  See ticket #64 on
  the PyFITS Trac for further details and background. Some highlights are
  listed below:

  * The Header class now fully implements the Python dict interface, and can
    be used interchangeably with a dict, where the keys are header keywords.

  * New keywords can be added to the header using normal keyword assignment
    (previously it was necessary to use ``Header.update`` to add new
    keywords).  For example::

        >>> header['NAXIS'] = 2

    will update the existing 'FOO' keyword if it already exists, or add a new
    one if it doesn't exist, just like a dict.

  * It is possible to assign both a value and a comment at the same time using
    a tuple::

        >>> header['NAXIS'] = (2, 'Number of axes')

  * To add/update a new card and ensure it's added in a specific location, use
    ``Header.set()``::

        >>> header.set('NAXIS', 2, 'Number of axes', after='BITPIX')

    This works the same as the old ``Header.update()``.  ``Header.update()``
    still works in the old way too, but is deprecated.

  * Although ``Card`` objects still exist, it generally is not necessary to
    work with them directly.  ``Header.ascardlist()``/``Header.ascard`` are
    deprecated and should not be used.  To directly access the ``Card``
    objects in a header, use ``Header.cards``.

  * To access card comments, it is still possible to either go through the
    card itself, or through ``Header.comments``.  For example::

       >>> header.cards['NAXIS'].comment
       Number of axes
       >>> header.comments['NAXIS']
       Number of axes

  * ``Card`` objects can now be used interchangeably with
    ``(keyword, value, comment)`` 3-tuples.  They still have ``.value`` and
    ``.comment`` attributes as well.  The ``.key`` attribute has been renamed
    to ``.keyword`` for consistency, though ``.key`` is still supported (but
    deprecated).

- Memory mapping is now used by default to access HDU data.  That is,
  ``pyfits.open()`` uses ``memmap=True`` as the default.  This provides better
  performance in the majority of use cases--there are only some I/O intensive
  applications where it might not be desirable.  Enabling mmap by default also
  enabled finding and fixing a large number of bugs in PyFITS' handling of
  memory-mapped data (most of these bug fixes were backported to PyFITS
  3.0.5). (#85)

  * A new ``pyfits.USE_MEMMAP`` global variable was added.  Set
    ``pyfits.USE_MEMMAP = False`` to change the default memmap setting for
    opening files.  This is especially useful for controlling the behavior in
    applications where pyfits is deeply embedded.

  * Likewise, a new ``PYFITS_USE_MEMMAP`` environment variable is supported.
    Set ``PYFITS_USE_MEMMAP = 0`` in your environment to change the default
    behavior.

- The ``size()`` method on HDU objects is now a ``.size`` property--this
  returns the size in bytes of the data portion of the HDU, and in most cases
  is equivalent to ``hdu.data.nbytes`` (#83)

- ``BinTableHDU.tdump`` and ``BinTableHDU.tcreate`` are deprecated--use
  ``BinTableHDU.dump`` and ``BinTableHDU.load`` instead.  The new methods
  output the table data in a slightly different format from previous versions,
  which places quotes around each value.  This format is compatible with data
  dumps from previous versions of PyFITS, but not vice-versa due to a parsing
  bug in older versions.

- Likewise the ``pyfits.tdump`` and ``pyfits.tcreate`` convenience function
  versions of these methods have been renamed ``pyfits.tabledump`` and
  ``pyfits.tableload``.  The old deprecated, but currently retained for
  backwards compatibility. (r1125)

- A new global variable ``pyfits.EXTENSION_NAME_CASE_SENSITIVE`` was added.
  This serves as a replacement for ``pyfits.setExtensionNameCaseSensitive``
  which is not deprecated and may be removed in a future version.  To enable
  case-sensitivity of extension names (i.e. treat 'sci' as distinct from 'SCI')
  set ``pyfits.EXTENSION_NAME_CASE_SENSITIVE = True``.  The default is
  ``False``. (r1139)

- A new global configuration variable ``pyfits.STRIP_HEADER_WHITESPACE`` was
  added.  By default, if a string value in a header contains trailing
  whitespace, that whitespace is automatically removed when the value is read.
  Now if you set ``pyfits.STRIP_HEADER_WHITESPACE = False`` all whitespace is
  preserved. (#146)

- The old ``classExtensions`` extension mechanism (which was deprecated in
  PyFITS 3.0) is removed outright.  To our knowledge it was no longer used
  anywhere. (r1309)

- Warning messages from PyFITS issued through the Python warnings API are now
  output to stderr instead of stdout, as is the default.  PyFITS no longer
  modifies the default behavior of the warnings module with respect to which
  stream it outputs to. (r1319)

- The ``checksum`` argument to ``pyfits.open()`` now accepts a value of
  'remove', which causes any existing CHECKSUM/DATASUM keywords to be ignored,
  and removed when the file is saved.

New Features
^^^^^^^^^^^^

- Added support for the proposed "FITS" extension HDU type.  See
  http://listmgr.cv.nrao.edu/pipermail/fitsbits/2002-April/001094.html.  FITS
  HDUs contain an entire FITS file embedded in their data section.  ``FitsHDU``
  objects work like other HDU types in PyFITS.  Their ``.data`` attribute
  returns the raw data array.  However, they have a special ``.hdulist``
  attribute which processes the data as a FITS file and returns it as an
  in-memory HDUList object.  FitsHDU objects also support a
  ``FitsHDU.fromhdulist()`` classmethod which returns a new ``FitsHDU`` object
  that embeds the supplied HDUList. (#80)

- Added a new ``.is_image`` attribute on HDU objects, which is True if the HDU
  data is an 'image' as opposed to a table or something else.  Here the
  meaning of 'image' is fairly loose, and mostly just means a Primary or Image
  extension HDU, or possibly a compressed image HDU (#71)

- Added an ``HDUList.fromstring`` classmethod which can parse a FITS file
  already in memory and instantiate and ``HDUList`` object from it.  This
  could be useful for integrating PyFITS with other libraries that work on
  FITS file, such as CFITSIO.  It may also be useful in streaming
  applications.  The name is a slight misnomer, in that it actually accepts
  any Python object that implements the buffer interface, which includes
  ``bytes``, ``bytearray``, ``memoryview``, ``numpy.ndarray``, etc. (#90)

- Added a new ``pyfits.diff`` module which contains facilities for comparing
  FITS files.  One can use the ``pyfits.diff.FITSDiff`` class to compare two
  FITS files in their entirety.  There is also a ``pyfits.diff.HeaderDiff``
  class for just comparing two FITS headers, and other similar interfaces.
  See the PyFITS Documentation for more details on this interface.  The
  ``pyfits.diff`` module powers the new ``fitsdiff`` program installed with
  PyFITS.  After installing PyFITS, run ``fitsdiff --help`` for usage details.

- ``pyfits.open()`` now accepts a ``scale_back`` argument.  If set to
  ``True``, this automatically scales the data using the original BZERO and
  BSCALE parameters the file had when it was first opened, if any, as well as
  the original BITPIX.  For example, if the original BITPIX were 16, this
  would be equivalent to calling ``hdu.scale('int16', 'old')`` just before
  calling ``flush()`` or ``close()`` on the file.  This option applies to all
  HDUs in the file. (#120)

- ``pyfits.open()`` now accepts a ``save_backup`` argument.  If set to
  ``True``, this automatically saves a backup of the original file before
  flushing any changes to it (this of course only applies to update and append
  mode).  This may be especially useful when working with scaled image data.
  (#121)

Changes in Behavior
^^^^^^^^^^^^^^^^^^^

- Warnings from PyFITS are not output to stderr by default, instead of stdout
  as it has been for some time.  This is contrary to most users' expectations
  and makes it more difficult for them to separate output from PyFITS from the
  desired output for their scripts. (r1319)

Bug Fixes
^^^^^^^^^

- Fixed ``pyfits.tcreate()`` (now ``pyfits.tableload()``) to be more robust
  when encountering blank lines in a column definition file (#14)

- Fixed a fairly rare crash that could occur in the handling of CONTINUE cards
  when using Numpy 1.4 or lower (though 1.4 is the oldest version supported by
  PyFITS). (r1330)

- Fixed ``_BaseHDU.fromstring`` to actually correctly instantiate an HDU
  object from a string/buffer containing the header and data of that HDU.
  This allowed for the implementation of ``HDUList.fromstring`` described
  above. (#90)

- Fixed a rare corner case where, in some use cases, (mildly, recoverably)
  malformatted float values in headers were not properly returned as floats.
  (#137)

- Fixed a corollary to the previous bug where float values with a leading zero
  before the decimal point had the leading zero unnecessarily removed when
  saving changes to the file (eg. "0.001" would be written back as ".001" even
  if no changes were otherwise made to the file). (#137)

- When opening a file containing CHECKSUM and/or DATASUM keywords in update
  mode, the CHECKSUM/DATASUM are updated and preserved even if the file was
  opened with checksum=False.  This change in behavior prevents checksums from
  being unintentionally removed. (#148)

- Fixed a bug where ``ImageHDU.scale(option='old')`` wasn't working at all--it
  was not restoring the image to its original BSCALE and BZERO values. (#162)

- Fixed a bug when writing out files containing zero-width table columns,
  where the TFIELDS keyword would be updated incorrectly, leaving the table
  largely unreadable.  This fix will be backported to the 3.0.x series in
  version 3.0.10.  (#174)


3.0.9 (2012-08-06)
------------------

This is a bug fix release for the 3.0.x series.

Bug Fixes
^^^^^^^^^

- Fixed ``Header.values()``/``Header.itervalues()`` and ``Header.items()``/
  ``Header.iteritems()`` to correctly return the different values for
  duplicate keywords (particularly commentary keywords like HISTORY and
  COMMENT).  This makes the old Header implementation slightly more compatible
  with the new implementation in PyFITS 3.1. (#127)

  .. note::
      This fix did not change the existing behavior from earlier PyFITS
      versions where ``Header.keys()`` returns all keywords in the header with
      duplicates removed.  PyFITS 3.1 changes that behavior, so that
      ``Header.keys()`` includes duplicates.

- Fixed a bug where ``ImageHDU.scale(option='old')`` wasn't working at all--it
  was not restoring the image to its original BSCALE and BZERO values. (#162)

- Fixed a bug where opening a file containing compressed image HDUs in
  'update' mode and then immediately closing it without making any changes
  caused the file to be rewritten unnecessarily. (#167)

- Fixed two memory leaks that could occur when writing compressed image data,
  or in some cases when opening files containing compressed image HDUs in
  'update' mode. (#168)


3.0.8 (2012-06-04)
------------------

Changes in Behavior
^^^^^^^^^^^^^^^^^^^

- Prior to this release, image data sections did not work with scaled
  data--that is, images with non-trivial BSCALE and/or BZERO values.
  Previously, in order to read such images in sections, it was necessary to
  manually apply the BSCALE+BZERO to each section.  It's worth noting that
  sections *did* support pseudo-unsigned ints (flakily).  This change just
  extends that support for general BSCALE+BZERO values.

Bug Fixes
^^^^^^^^^

- Fixed a bug that prevented updates to values in boolean table columns from
  being saved.  This turned out to be a symptom of a deeper problem that could
  prevent other table updates from being saved as well. (#139)

- Fixed a corner case in which a keyword comment ending with the string "END"
  could, in some circumstances, cause headers (and the rest of the file after
  that point) to be misread. (#142)

- Fixed support for scaled image data and psuedo-unsigned ints in image data
  sections (``hdu.section``).  Previously this was not supported at all.  At
  some point support was supposedly added, but it was buggy and incomplete.
  Now the feature seems to work much better. (#143)

- Fixed the documentation to point out that image data sections *do* support
  non-contiguous slices (and have for a long time).  The documentation was
  never updated to reflect this, and misinformed users that only contiguous
  slices were supported, leading to some confusion. (#144)

- Fixed a bug where creating an ``HDUList`` object containing multiple PRIMARY
  HDUs caused an infinite recursion when validating the object prior to
  writing to a file. (#145)

- Fixed a rare but serious case where saving an update to a file that
  previously had a CHECKSUM and/or DATASUM keyword, but removed the checksum
  in saving, could cause the file to be slightly corrupted and unreadable.
  (#147)

- Fixed problems with reading "non-standard" FITS files with primary headers
  containing SIMPLE = F.  PyFITS has never made many guarantees as to how such
  files are handled.  But it should at least be possible to read their
  headers, and the data if possible.  Saving changes to such a file should not
  try to prepend an unwanted valid PRIMARY HDU. (#157)

- Fixed a bug where opening an image with ``disable_image_compression = True``
  caused compression to be disabled for all subsequent ``pyfits.open()`` calls.
  (r1651)


3.0.7 (2012-04-10)
------------------

Changes in Behavior
^^^^^^^^^^^^^^^^^^^

- Slices of GroupData objects now return new GroupData objects instead of
  extended multi-row _Group objects. This is analogous to how PyFITS 3.0 fixed
  FITS_rec slicing, and should have been fixed for GroupData at the same time.
  The old behavior caused bugs where functions internal to Numpy expected that
  slicing an ndarray would return a new ndarray.  As this is a rare usecase
  with a rare feature most users are unlikely to be affected by this change.

- The previously internal _Group object for representing individual group
  records in a GroupData object are renamed Group and are now a public
  interface.  However, there's almost no good reason to create Group objects
  directly, so it shouldn't be considered a "new feature".

- An annoyance from PyFITS 3.0.6 was fixed, where the value of the EXTEND
  keyword was always being set to F if there are not actually any extension
  HDUs.  It was unnecessary to modify this value.

Bug Fixes
^^^^^^^^^

- Fixed GroupData objects to return new GroupData objects when sliced instead
  of _Group record objects.  See "Changes in behavior" above for more details.

- Fixed slicing of Group objects--previously it was not possible to slice
  slice them at all.

- Made it possible to assign ``np.bool_`` objects as header values. (#123)

- Fixed overly strict handling of the EXTEND keyword; see "Changes in
  behavior" above. (#124)

- Fixed many cases where an HDU's header would be marked as "modified" by
  PyFITS and rewritten, even when no changes to the header are necessary.
  (#125)

- Fixed a bug where the values of the PTYPEn keywords in a random groups HDU
  were forced to be all lower-case when saving the file. (#130)

- Removed an unnecessary inline import in ``ExtensionHDU.__setattr__`` that was
  causing some slowdown when opening files containing a large number of
  extensions, plus a few other small (but not insignificant) performance
  improvements thanks to Julian Taylor. (#133)

- Fixed a regression where header blocks containing invalid end-of-header
  padding (i.e. null bytes instead of spaces) couldn't be parsed by PyFITS.
  Such headers can be parsed again, but a warning is raised, as such headers
  are not valid FITS. (#136)

- Fixed a memory leak where table data in random groups HDUs weren't being
  garbage collected. (#138)


3.0.6 (2012-02-29)
------------------

Highlights
^^^^^^^^^^

The main reason for this release is to fix an issue that was introduced in
PyFITS 3.0.5 where merely opening a file containing scaled data (that is, with
non-trivial BSCALE and BZERO keywords) in 'update' mode would cause the data
to be automatically rescaled--possibly converting the data from ints to
floats--as soon as the file is closed, even if the application did not touch
the data.  Now PyFITS will only rescale the data in an extension when the data
is actually accessed by the application.  So opening a file in 'update' mode
in order to modify the header or append new extensions will not cause any
change to the data in existing extensions.

This release also fixes a few Windows-specific bugs found through more
extensive Windows testing, and other miscellaneous bugs.

Bug Fixes
^^^^^^^^^

- More accurate error messages when opening files containing invalid header
  cards. (#109)

- Fixed a possible reference cycle/memory leak that was caught through more
  extensive testing on Windows. (#112)

- Fixed 'ostream' mode to open the underlying file in 'wb' mode instead of 'w'
  mode. (#112)

- Fixed a Windows-only issue where trying to save updates to a resized FITS
  file could result in a crash due to there being open mmaps on that file.
  (#112)

- Fixed a crash when trying to create a FITS table (i.e. with new_table())
  from a Numpy array containing bool fields. (#113)

- Fixed a bug where manually initializing an ``HDUList`` with a list of of
  HDUs wouldn't set the correct EXTEND keyword value on the primary HDU.
  (#114)

- Fixed a crash that could occur when trying to deepcopy a Header in Python <
  2.7. (#115)

- Fixed an issue where merely opening a scaled image in 'update' mode would
  cause the data to be converted to floats when the file is closed. (#119)


3.0.5 (2012-01-30)
------------------

- Fixed a crash that could occur when accessing image sections of files
  opened with memmap=True. (r1211)

- Fixed the inconsistency in the behavior of files opened in 'readonly' mode
  when memmap=True vs. when memmap=False.  In the latter case, although
  changes to array data were not saved to disk, it was possible to update the
  array data in memory.  On the other hand with memmap=True, 'readonly' mode
  prevented even in-memory modification to the data.  This is what
  'copyonwrite' mode was for, but difference in behavior was confusing.  Now
  'readonly' is equivalent to 'copyonwrite' when using memmap.  If the old
  behavior of denying changes to the array data is necessary, a new
  'denywrite' mode may be used, though it is only applicable to files opened
  with memmap. (r1275)

- Fixed an issue where files opened with memmap=True would return image data
  as a raw numpy.memmap object, which can cause some unexpected
  behaviors--instead memmap object is viewed as a numpy.ndarray. (r1285)

- Fixed an issue in Python 3 where a workaround for a bug in Numpy on Python 3
  interacted badly with some other software, namely to vo.table package (and
  possibly others). (r1320, r1337, and #110)

- Fixed buggy behavior in the handling of SIGINTs (i.e. Ctrl-C keyboard
  interrupts) while flushing changes to a FITS file.  PyFITS already prevented
  SIGINTs from causing an incomplete flush, but did not clean up the signal
  handlers properly afterwards, or reraise the keyboard interrupt once the
  flush was complete. (r1321)

- Fixed a crash that could occur in Python 3 when opening files with checksum
  checking enabled. (r1336)

- Fixed a small bug that could cause a crash in the ``StreamingHDU`` interface
  when using Numpy below version 1.5.

- Fixed a crash that could occur when creating a new ``CompImageHDU`` from an
  array of big-endian data. (#104)

- Fixed a crash when opening a file with extra zero padding at the end.
  Though FITS files should not have such padding, it's not explicitly forbidden
  by the format either, and PyFITS shouldn't stumble over it. (#106)

- Fixed a major slowdown in opening tables containing large columns of string
  values.  (#111)


3.0.4 (2011-11-22)
------------------

- Fixed a crash when writing HCOMPRESS compressed images that could happen on
  Python 2.5 and 2.6. (r1217)

- Fixed a crash when slicing an table in a file opened in 'readonly' mode with
  memmap=True. (r1230)

- Writing changes to a file or writing to a new file verifies the output in
  'fix' mode by default instead of 'exception'--that is, PyFITS will
  automatically fix common FITS format errors rather than raising an
  exception. (r1243)

- Fixed a bug where convenience functions such as getval() and getheader()
  crashed when specifying just 'PRIMARY' as the extension to use (r1263).

- Fixed a bug that prevented passing keyword arguments (beyond the standard
  data and header arguments) as positional arguments to the constructors of
  extension HDU classes.

- Fixed some tests that were failing on Windows--in this case the tests
  themselves failed to close some temp files and Windows refused to delete them
  while there were still open handles on them. (r1295)

- Fixed an issue with floating point formatting in header values on Python 2.5
  for Windows (and possibly other platforms).  The exponent was zero-padded to
  3 digits; although the FITS standard makes no specification on this, the
  formatting is now normalized to always pad the exponent to two digits.
  (r1295)

- Fixed a bug where long commentary cards (such as HISTORY and COMMENT) were
  broken into multiple CONTINUE cards.  However, commentary cards are not
  expected to be found in CONTINUE cards.  Instead these long cards are broken
  into multiple commentary cards. (#97)

- GZIP/ZIP-compressed FITS files can be detected and opened regardless of
  their filename extension. (#99)

- Fixed a serious bug where opening scaled images in 'update' mode and then
  closing the file without touching the data would cause the file to be
  corrupted. (#101)


3.0.3 (2011-10-05)
------------------

- Fixed several small bugs involving corner cases in record-valued keyword
  cards (#70)

- In some cases HDU creation failed if the first keyword value in the header
  was not a string value (#89)

- Fixed a crash when trying to compute the HDU checksum when the data array
  contains an odd number of bytes (#91)

- Disabled an unnecessary warning that was displayed on opening compressed
  HDUs with disable_image_compression = True (#92)

- Fixed a typo in code for handling HCOMPRESS compressed images.


3.0.2 (2011-09-23)
------------------

- The ``BinTableHDU.tcreate`` method and by extension the ``pyfits.tcreate``
  function don't get tripped up by blank lines anymore (#14)

- The presence, value, and position of the EXTEND keyword in Primary HDUs is
  verified when reading/writing a FITS file (#32)

- Improved documentation (in warning messages as well as in the handbook) that
  PyFITS uses zero-based indexing (as one would expect for C/Python code, but
  contrary to the PyFITS standard which was written with FORTRAN in mind)
  (#68)

- Fixed a bug where updating a header card comment could cause the value to be
  lost if it had not already been read from the card image string.

- Fixed a related bug where changes made directly to Card object in a header
  (i.e. assigning directly to card.value or card.comment) would not propagate
  when flushing changes to the file (#69) [Note: This and the bug above it
  were originally reported as being fixed in version 3.0.1, but the fix was
  never included in the release.]

- Improved file handling, particularly in Python 3 which had a few small file
  I/O-related bugs (#76)

- Fixed a bug where updating a FITS file would sometimes cause it to lose its
  original file permissions (#79)

- Fixed the handling of TDIMn keywords; 3.0 added support for them, but got
  the axis order backards (they were treated as though they were row-major)
  (#82)

- Fixed a crash when a FITS file containing scaled data is opened and
  immediately written to a new file without explicitly viewing the data first
  (#84)

- Fixed a bug where creating a table with columns named either 'names' or
  'formats' resulted in an infinite recursion (#86)


3.0.1 (2011-09-12)
------------------

- Fixed a bug where updating a header card comment could cause the value to be
  lost if it had not already been read from the card image string.

- Changed ``_TableBaseHDU.data`` so that if the data contain an empty table a
  ``FITS_rec`` object with zero rows is returned rather than ``None``.

- The ``.key`` attribute of ``RecordValuedKeywordCards`` now returns the full
  keyword+field-specifier value, instead of just the plain keyword (#46)

- Fixed a related bug where changes made directly to Card object in a header
  (i.e. assigning directly to card.value or card.comment) would not propagate
  when flushing changes to the file (#69)

- Fixed a bug where writing a table with zero rows could fail in some cases
  (#72)

- Miscellaneous small bug fixes that were causing some tests to fail,
  particularly on Python 3 (#74, #75)

- Fixed a bug where creating a table column from an array in non-native byte
  order would not preserve the byte order, thus interpreting the column array
  using the wrong byte order (#77)


3.0.0 (2011-08-23)
--------------------

- Contains major changes, bumping the version to 3.0

- Large amounts of refactoring and reorganization of the code; tried to
  preserve public API backwards-compatibility with older versions (private API
  has many changes and is not guaranteed to be backwards-compatible).  There
  are a few small public API changes to be aware of:

  * The pyfits.rec module has been removed completely.  If your version of
    numpy does not have the numpy.core.records module it is too old to be used
    with PyFITS.

  * The ``Header.ascardlist()`` method is deprecated--use the ``.ascard``
    attribute instead.

  * ``Card`` instances have a new ``.cardimage`` attribute that should be used
    rather than ``.ascardimage()``, which may become deprecated.

  * The ``Card.fromstring()`` method is now a classmethod.  It returns a new
    ``Card`` instance rather than modifying an existing instance.

  * The ``req_cards()`` method on HDU instances has changed:  The ``pos``
    argument is not longer a string.  It is either an integer value (meaning
    the card's position must match that value) or it can be a function that
    takes the card's position as it's argument, and returns True if the
    position is valid.  Likewise, the ``test`` argument no longer takes a
    string, but instead a function that validates the card's value and returns
    True or False.

  * The ``get_coldefs()`` method of table HDUs is deprecated.  Use the
    ``.columns`` attribute instead.

  * The ``ColDefs.data`` attribute is deprecated--use ``ColDefs.columns``
    instead (though in general you shouldn't mess with it directly--it might
    become internal at some point).

  * ``FITS_record`` objects take ``start`` and ``end`` as arguments instead of
    ``startColumn`` and ``endColumn`` (these are rarely created manually, so
    it's unlikely that this change will affect anyone).

  * ``BinTableHDU.tcreate()`` is now a classmethod, and returns a new
    ``BinTableHDU`` instance.

  * Use ``ExtensionHDU`` and ``NonstandardExtHDU`` for making new extension HDU
    classes.  They are now public interfaces, wheres previously they were
    private and prefixed with underscores.

  * Possibly others--please report if you find any changes that cause
    difficulties.

- Calls to deprecated functions will display a Deprecation warning.  However,
  in Python 2.7 and up Deprecation warnings are ignored by default, so run
  Python with the ``-Wd`` option to see if you're using any deprecated
  functions.  If we get close to actually removing any functions, we might
  make the Deprecation warnings display by default.

- Added basic Python 3 support

- Added support for multi-dimensional columns in tables as specified by the
  TDIMn keywords (#47)

- Fixed a major memory leak that occurred when creating new tables with the
  ``new_table()`` function (#49)
  be padded with zero-bytes) vs ASCII tables (where strings are padded with
  spaces) (#15)

- Fixed a bug in which the case of Random Access Group parameters names was not
  preserved when writing (#41)

- Added support for binary table fields with zero width (#42)

- Added support for wider integer types in ASCII tables; although this is non-
  standard, some GEIS images require it (#45)

- Fixed a bug that caused the index_of() method of HDULists to crash when the
  HDUList object is created from scratch (#48)

- Fixed the behavior of string padding in binary tables (where strings should
  be padded with nulls instead of spaces)

- Fixed a rare issue that caused excessive memory usage when computing
  checksums using a non-standard block size (see r818)

- Add support for forced uint data in image sections (#53)

- Fixed an issue where variable-length array columns were not extended when
  creating a new table with more rows than the original (#54)

- Fixed tuple and list-based indexing of FITS_rec objects (#55)

- Fixed an issue where BZERO and BSCALE keywords were appended to headers in
  the wrong location (#56)

- ``FITS_record`` objects (table rows) have full slicing support, including
  stepping, etc. (#59)

- Fixed a bug where updating multiple files simultaneously (such as when
  running parallel processes) could lead to a race condition with mktemp()
  (#61)

- Fixed a bug where compressed image headers were not in the order expected by
  the funpack utility (#62)


2.4.0 (2011-01-10)
--------------------
The following enhancements were added:

- Checksum support now correctly conforms to the FITS standard.  pyfits
  supports reading and writing both the old checksums and new
  standard-compliant checksums.  The ``fitscheck`` command-line utility is
  provided to verify and update checksums.

- Added a new optional keyword argument ``do_not_scale_image_data``
  to the ``pyfits.open`` convenience function.  When this argument
  is provided as True, and an ImageHDU is read that contains scaled
  data, the data is not automatically scaled when it is read.  This
  option may be used when opening a fits file for update, when you only
  want to update some header data.  Without the use of this argument, if
  the header updates required the size of the fits file to change, then
  when writing the updated information, the data would be read, scaled,
  and written back out in its scaled format (usually with a different
  data type) instead of in its non-scaled format.

- Added a new optional keyword argument ``disable_image_compression`` to the
  ``pyfits.open`` function.  When ``True``, any compressed image HDU's will
  be read in like they are binary table HDU's.

- Added a ``verify`` keyword argument to the ``pyfits.append`` function.  When
  ``False``, ``append`` will assume the existing FITS file is already valid
  and simply append new content to the end of the file, resulting in a large
  speed up appending to large files.

- Added HDU methods ``update_ext_name`` and ``update_ext_version`` for
  updating the name and version of an HDU.

- Added HDU method ``filebytes`` to calculate the number of bytes that will be
  written to the file associated with the HDU.

- Enhanced the section class to allow reading non-contiguous image data.
  Previously, the section class could only be used to read contiguous data.
  (CNSHD781626)

- Added method ``HDUList.fileinfo()`` that returns a dictionary with
  information about the location of header and data in the file associated
  with the HDU.

The following bugs were fixed:

- Reading in some malformed FITS headers would cause a ``NameError``
  exception, rather than information about the cause of the error.

- pyfits can now handle non-compliant ``CONTINUE`` cards produced by Java
  FITS.

- ``BinTable`` columns with ``TSCALn`` are now byte-swapped correctly.

- Ensure that floating-point card values are no longer than 20 characters.

- Updated ``flush`` so that when the data has changed in an HDU for a file
  opened in update mode, the header will be updated to match the changed data
  before writing out the HDU.

- Allow ``HIERARCH`` cards to contain a keyword and value whose total
  character length is 69 characters.  Previous length was limited at 68
  characters.

- Calls to ``FITS_rec['columnName']`` now return an ``ndarray``. exactly the
  same as a call to ``FITS_rec.field('columnName')`` or
  ``FITS_rec.columnName``.  Previously, ``FITS_rec['columnName']`` returned a
  much less useful ``fits_record`` object. (CNSHD789053)

- Corrected the ``append`` convenience function to eliminate the reading of
  the HDU data from the file that is being appended to.  (CNSHD794738)

- Eliminated common symbols between the pyfitsComp module and the cfitsio and
  zlib libraries.  These can cause problems on systems that use both PyFITS
  and cfitsio or zlib. (CNSHD795046)


2.3.1 (2010-06-03)
--------------------

The following bugs were fixed:

- Replaced code in the Compressed Image HDU extension which was covered under
  a GNU General Public License with code that is covered under a BSD License.
  This change allows the distribution of pyfits under a BSD License.


2.3 (2010-05-11)
------------------

The following enhancements were made:

- Completely eliminate support for numarray.

- Rework pyfits documentation to use Sphinx.

- Support python 2.6 and future division.

- Added a new method to get the file name associated with an HDUList object.
  The method HDUList.filename() returns the name of an associated file.  It
  returns None if no file is associated with the HDUList.

- Support the python 2.5 'with' statement when opening fits files.
  (CNSHD766308)  It is now possible to use the following construct:

    >>> from __future__ import with_statement import pyfits
    >>> with pyfits.open("input.fits") as hdul:
    ...    #process hdul
    >>>

- Extended the support for reading unsigned integer 16 values from an ImageHDU
  to include unsigned integer 32 and unsigned integer 64 values.  ImageHDU
  data is considered to be unsigned integer 16 when the data type is signed
  integer 16 and BZERO is equal to 2**15 (32784) and BSCALE is equal to 1.
  ImageHDU data is considered to be unsigned integer 32 when the data type is
  signed integer 32 and BZERO is equal to 2**31 and BSCALE is equal to 1.
  ImageHDU data is considered to be unsigned integer 64 when the data type is
  signed integer 64 and BZERO is equal to 2**63 and BSCALE is equal to 1.  An
  optional keyword argument (uint) was added to the open convenience function
  for this purpose.  Supplying a value of True for this argument will cause
  data of any of these types to be read in and scaled into the appropriate
  unsigned integer array (uint16, uint32, or uint64) instead of into the
  normal float 32 or float 64 array.  If an HDU associated with a file that
  was opened with the 'int' option and containing unsigned integer 16, 32, or
  64 data is written to a file, the data will be reverse scaled into a signed
  integer 16, 32, or 64 array and written out to the file along with the
  appropriate BSCALE/BZERO header cards.  Note that for backward
  compatibility, the 'uint16' keyword argument will still be accepted in the
  open function when handling unsigned integer 16 conversion.

- Provided the capability to access the data for a column of a fits table by
  indexing the table using the column name.  This is consistent with Record
  Arrays in numpy (array with fields).  (CNSHD763378)  The following example
  will illustrate this:

    >>> import pyfits
    >>> hdul = pyfits.open('input.fits')
    >>> table = hdul[1].data
    >>> table.names
    ['c1','c2','c3','c4']
    >>> print table.field('c2') # this is the data for column 2
    ['abc' 'xy']
    >>> print table['c2'] # this is also the data for column 2
    array(['abc', 'xy '], dtype='|S3')
    >>> print table[1] # this is the data for row 1
    (2, 'xy', 6.6999997138977054, True)

- Provided capabilities to create a BinaryTableHDU directly from a numpy
  Record Array (array with fields). The new capabilities include table
  creation, writing a numpy Record Array directly to a fits file using the
  pyfits.writeto and pyfits.append convenience functions.  Reading the data
  for a BinaryTableHDU from a fits file directly into a numpy Record Array
  using the pyfits.getdata convenience function.  (CNSHD749034)  Thanks to
  Erin Sheldon at Brookhaven National Laboratory for help with this.

  The following should illustrate these new capabilities:

    >>> import pyfits
    >>> import numpy
    >>> t=numpy.zeros(5,dtype=[('x','f4'),('y','2i4')]) \
    ... # Create a numpy Record Array with fields
    >>> hdu = pyfits.BinTableHDU(t) \
    ... # Create a Binary Table HDU directly from the Record Array
    >>> print hdu.data
    [(0.0, array([0, 0], dtype=int32))
     (0.0, array([0, 0], dtype=int32))
     (0.0, array([0, 0], dtype=int32))
     (0.0, array([0, 0], dtype=int32))
     (0.0, array([0, 0], dtype=int32))]
    >>> hdu.writeto('test1.fits',clobber=True) \
    ... # Write the HDU to a file
    >>> pyfits.info('test1.fits')
    Filename: test1.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU       4  ()            uint8
    1                BinTableHDU     12  5R x 2C       [E, 2J]
    >>> pyfits.writeto('test.fits', t, clobber=True) \
    ... # Write the Record Array directly to a file
    >>> pyfits.append('test.fits', t) \
    ... # Append another Record Array to the file
    >>> pyfits.info('test.fits')
    Filename: test.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU       4  ()            uint8
    1                BinTableHDU     12  5R x 2C       [E, 2J]
    2                BinTableHDU     12  5R x 2C       [E, 2J]
    >>> d=pyfits.getdata('test.fits',ext=1) \
    ... # Get the first extension from the file as a FITS_rec
    >>> print type(d)
    <class 'pyfits.core.FITS_rec'>
    >>> print d
    [(0.0, array([0, 0], dtype=int32))
     (0.0, array([0, 0], dtype=int32))
     (0.0, array([0, 0], dtype=int32))
     (0.0, array([0, 0], dtype=int32))
     (0.0, array([0, 0], dtype=int32))]
    >>> d=pyfits.getdata('test.fits',ext=1,view=numpy.ndarray) \
    ... # Get the first extension from the file as a numpy Record
          Array
    >>> print type(d)
    <type 'numpy.ndarray'>
    >>> print d
    [(0.0, [0, 0]) (0.0, [0, 0]) (0.0, [0, 0]) (0.0, [0, 0])
     (0.0, [0, 0])]
    >>> print d.dtype
    [('x', '>f4'), ('y', '>i4', 2)]
    >>> d=pyfits.getdata('test.fits',ext=1,upper=True,
    ...                  view=pyfits.FITS_rec) \
    ... # Force the Record Array field names to be in upper case
          regardless of how they are stored in the file
    >>> print d.dtype
    [('X', '>f4'), ('Y', '>i4', 2)]

- Provided support for writing fits data to file-like objects that do not
  support the random access methods seek() and tell().  Most pyfits functions
  or methods will treat these file-like objects as an empty file that cannot
  be read, only written.  It is also expected that the file-like object is in
  a writable condition (ie. opened) when passed into a pyfits function or
  method.  The following methods and functions will allow writing to a
  non-random access file-like object: HDUList.writeto(), HDUList.flush(),
  pyfits.writeto(), and pyfits.append().  The pyfits.open() convenience
  function may be used to create an HDUList object that is associated with the
  provided file-like object.  (CNSHD770036)

  An illustration of the new capabilities follows.  In this example fits data
  is written to standard output which is associated with a file opened in
  write-only mode:

    >>> import pyfits
    >>> import numpy as np
    >>> import sys
    >>>
    >>> hdu = pyfits.PrimaryHDU(np.arange(100,dtype=np.int32))
    >>> hdul = pyfits.HDUList()
    >>> hdul.append(hdu)
    >>> tmpfile = open('tmpfile.py','w')
    >>> sys.stdout = tmpfile
    >>> hdul.writeto(sys.stdout, clobber=True)
    >>> sys.stdout = sys.__stdout__
    >>> tmpfile.close()
    >>> pyfits.info('tmpfile.py')
    Filename: tmpfile.py
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU       5  (100,)        int32
    >>>

- Provided support for slicing a FITS_record object.  The FITS_record object
  represents the data from a row of a table.  Pyfits now supports the slice
  syntax to retrieve values from the row.  The following illustrates this new
  syntax:

    >>> hdul = pyfits.open('table.fits')
    >>> row = hdul[1].data[0]
    >>> row
    ('clear', 'nicmos', 1, 30, 'clear', 'idno= 100')
    >>> a, b, c, d, e = row[0:5]
    >>> a
    'clear'
    >>> b
    'nicmos'
    >>> c
    1
    >>> d
    30
    >>> e
    'clear'
    >>>

- Allow the assignment of a row value for a pyfits table using a tuple or a
  list as input.  The following example illustrates this new feature:

    >>> c1=pyfits.Column(name='target',format='10A')
    >>> c2=pyfits.Column(name='counts',format='J',unit='DN')
    >>> c3=pyfits.Column(name='notes',format='A10')
    >>> c4=pyfits.Column(name='spectrum',format='5E')
    >>> c5=pyfits.Column(name='flag',format='L')
    >>> coldefs=pyfits.ColDefs([c1,c2,c3,c4,c5])
    >>>
    >>> tbhdu=pyfits.new_table(coldefs, nrows = 5)
    >>>
    >>> # Assigning data to a table's row using a tuple
    >>> tbhdu.data[2] = ('NGC1',312,'A Note',
    ... num.array([1.1,2.2,3.3,4.4,5.5],dtype=num.float32),
    ... True)
    >>>
    >>> # Assigning data to a tables row using a list
    >>> tbhdu.data[3] = ['JIM1','33','A Note',
    ... num.array([1.,2.,3.,4.,5.],dtype=num.float32),True]

- Allow the creation of a Variable Length Format (P format) column from a list
  of data.  The following example illustrates this new feature:

    >>> a = [num.array([7.2e-20,7.3e-20]),num.array([0.0]),
    ... num.array([0.0])]
    >>> acol = pyfits.Column(name='testa',format='PD()',array=a)
    >>> acol.array
    _VLF([[  7.20000000e-20   7.30000000e-20], [ 0.], [ 0.]],
    dtype=object)
    >>>

- Allow the assignment of multiple rows in a table using the slice syntax. The
  following example illustrates this new feature:

    >>> counts = num.array([312,334,308,317])
    >>> names = num.array(['NGC1','NGC2','NGC3','NCG4'])
    >>> c1=pyfits.Column(name='target',format='10A',array=names)
    >>> c2=pyfits.Column(name='counts',format='J',unit='DN',
    ... array=counts)
    >>> c3=pyfits.Column(name='notes',format='A10')
    >>> c4=pyfits.Column(name='spectrum',format='5E')
    >>> c5=pyfits.Column(name='flag',format='L',array=[1,0,1,1])
    >>> coldefs=pyfits.ColDefs([c1,c2,c3,c4,c5])
    >>>
    >>> tbhdu1=pyfits.new_table(coldefs)
    >>>
    >>> counts = num.array([112,134,108,117])
    >>> names = num.array(['NGC5','NGC6','NGC7','NCG8'])
    >>> c1=pyfits.Column(name='target',format='10A',array=names)
    >>> c2=pyfits.Column(name='counts',format='J',unit='DN',
    ... array=counts)
    >>> c3=pyfits.Column(name='notes',format='A10')
    >>> c4=pyfits.Column(name='spectrum',format='5E')
    >>> c5=pyfits.Column(name='flag',format='L',array=[0,1,0,0])
    >>> coldefs=pyfits.ColDefs([c1,c2,c3,c4,c5])
    >>>
    >>> tbhdu=pyfits.new_table(coldefs)
    >>> tbhdu.data[0][3] = num.array([1.,2.,3.,4.,5.],
    ... dtype=num.float32)
    >>>
    >>> tbhdu2=pyfits.new_table(tbhdu1.data, nrows=9)
    >>>
    >>> # Assign the 4 rows from the second table to rows 5 thru
    ...   8 of the new table.  Note that the last row of the new
    ...   table will still be initialized to the default values.
    >>> tbhdu2.data[4:] = tbhdu.data
    >>>
    >>> print tbhdu2.data
    [ ('NGC1', 312, '0.0', array([ 0.,  0.,  0.,  0.,  0.],
    dtype=float32), True)
      ('NGC2', 334, '0.0', array([ 0.,  0.,  0.,  0.,  0.],
    dtype=float32), False)
      ('NGC3', 308, '0.0', array([ 0.,  0.,  0.,  0.,  0.],
    dtype=float32), True)
      ('NCG4', 317, '0.0', array([ 0.,  0.,  0.,  0.,  0.],
    dtype=float32), True)
      ('NGC5', 112, '0.0', array([ 1.,  2.,  3.,  4.,  5.],
    dtype=float32), False)
      ('NGC6', 134, '0.0', array([ 0.,  0.,  0.,  0.,  0.],
    dtype=float32), True)
      ('NGC7', 108, '0.0', array([ 0.,  0.,  0.,  0.,  0.],
    dtype=float32), False)
      ('NCG8', 117, '0.0', array([ 0.,  0.,  0.,  0.,  0.],
    dtype=float32), False)
      ('0.0', 0, '0.0', array([ 0.,  0.,  0.,  0.,  0.],
    dtype=float32), False)]
    >>>

The following bugs were fixed:

- Corrected bugs in HDUList.append and HDUList.insert to correctly handle the
  situation where you want to insert or append a Primary HDU as something
  other than the first HDU in an HDUList and the situation where you want to
  insert or append an Extension HDU as the first HDU in an HDUList.

- Corrected a bug involving scaled images (both compressed and not compressed)
  that include a BLANK, or ZBLANK card in the header.  When the image values
  match the BLANK or ZBLANK value, the value should be replaced with NaN after
  scaling.  Instead, pyfits was scaling the BLANK or ZBLANK value and
  returning it. (CNSHD766129)

- Corrected a byteswapping bug that occurs when writing certain column data.
  (CNSHD763307)

- Corrected a bug that occurs when creating a column from a chararray when one
  or more elements are shorter than the specified format length.  The bug
  wrote nulls instead of spaces to the file. (CNSHD695419)

- Corrected a bug in the HDU verification software to ensure that the header
  contains no NAXISn cards where n > NAXIS.

- Corrected a bug involving reading and writing compressed image data.  When
  written, the header keyword card ZTENSION will always have the value 'IMAGE'
  and when read, if the ZTENSION value is not 'IMAGE' the user will receive a
  warning, but the data will still be treated as image data.

- Corrected a bug that restricted the ability to create a custom HDU class and
  use it with pyfits.  The bug fix will allow something like this:

    >>> import pyfits
    >>> class MyPrimaryHDU(pyfits.PrimaryHDU):
    ...     def __init__(self, data=None, header=None):
    ...         pyfits.PrimaryHDU.__init__(self, data, header)
    ...     def _summary(self):
    ...         """
    ...         Reimplement a method of the class.
    ...         """
    ...         s = pyfits.PrimaryHDU._summary(self)
    ...         # change the behavior to suit me.
    ...         s1 = 'MyPRIMARY ' + s[11:]
    ...         return s1
    ...
    >>> hdul=pyfits.open("pix.fits",
    ... classExtensions={pyfits.PrimaryHDU: MyPrimaryHDU})
    >>> hdul.info()
    Filename: pix.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    MyPRIMARY  MyPrimaryHDU     59  (512, 512)    int16
    >>>

- Modified ColDefs.add_col so that instead of returning a new ColDefs object
  with the column added to the end, it simply appends the new column to the
  current ColDefs object in place.  (CNSHD768778)

- Corrected a bug in ColDefs.del_col which raised a KeyError exception when
  deleting a column from a ColDefs object.

- Modified the open convenience function so that when a file is opened in
  readonly mode and the file contains no HDU's an IOError is raised.

- Modified _TableBaseHDU to ensure that all locations where data is referenced
  in the object actually reference the same ndarray, instead of copies of the
  array.

- Corrected a bug in the Column class that failed to initialize data when the
  data is a boolean array.  (CNSHD779136)

- Corrected a bug that caused an exception to be raised when creating a
  variable length format column from character data (PA format).

- Modified installation code so that when installing on Windows, when a C++
  compiler compatible with the Python binary is not found, the installation
  completes with a warning that all optional extension modules failed to
  build.  Previously, an Error was issued and the installation stopped.


2.2.2 (2009-10-12)
--------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following bugs were fixed:

- Corrected a bug that caused an exception to be raised when creating a
  CompImageHDU using an initial header that does not match the image data in
  terms of the number of axis.


2.2.1 (2009-10-06)
--------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following bugs were fixed:

- Corrected a bug that prevented the opening of a fits file where a header
  contained a CHECKSUM card but no DATASUM card.

- Corrected a bug that caused NULLs to be written instead of blanks when an
  ASCII table was created using a numpy chararray in which the original data
  contained trailing blanks.  (CNSHD695419)


2.2 (2009-09-23)
------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following enhancements were made:

- Provide support for the FITS Checksum Keyword Convention.  (CNSHD754301)

- Adding the checksum=True keyword argument to the open convenience function
  will cause checksums to be verified on file open:

    >>> hdul=pyfits.open('in.fits', checksum=True)

- On output, CHECKSUM and DATASUM cards may be output to all HDU's in a fits
  file by using the keyword argument checksum=True in calls to the writeto
  convenience function, the HDUList.writeto method, the writeto methods of all
  of the HDU classes, and the append convenience function:

    >>> hdul.writeto('out.fits', checksum=True)

- Implemented a new insert method to the HDUList class that allows for the
  insertion of a HDU into a HDUList at a given index:

    >>> hdul.insert(2,hdu)

- Provided the capability to handle Unicode input for file names.

- Provided support for integer division required by Python 3.0.

The following bugs were fixed:

- Corrected a bug that caused an index out of bounds exception to be raised
  when iterating over the rows of a binary table HDU using the syntax  "for
  row in tbhdu.data:   ".  (CNSHD748609)

- Corrected a bug that prevented the use of the writeto convenience function
  for writing table data to a file.  (CNSHD749024)

- Modified the code to raise an IOError exception with the comment "Header
  missing END card." when pyfits can't find a valid END card for a header when
  opening a file.

  - This change addressed a problem with a non-standard fits file that
    contained several new-line characters at the end of each header and at the
    end of the file.  However, since some people want to be able to open these
    non-standard files anyway, an option was added to the open convenience
    function to allow these files to be opened without exception:

      >>> pyfits.open('infile.fits',ignore_missing_end=True)

- Corrected a bug that prevented the use of StringIO objects as fits files
  when reading and writing table data.  Previously, only image data was
  supported.  (CNSHD753698)

- Corrected a bug that caused a bus error to be generated when compressing
  image data using GZIP_1 under the Solaris operating system.

- Corrected bugs that prevented pyfits from properly reading Random Groups
  HDU's using numpy.  (CNSHD756570)

- Corrected a bug that can occur when writing a fits file.  (CNSHD757508)

  - If no default SIGINT signal handler has not been assigned, before the
    write, a TypeError exception is raised in the _File.flush() method when
    attempting to return the signal handler to its previous state.  Notably
    this occurred when using mod_python.  The code was changed to use SIG_DFL
    when no old handler was defined.

- Corrected a bug in CompImageHDU that prevented rescaling the image data
  using hdu.scale(option='old').


2.1.1 (2009-04-22)
-------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following bugs were fixed:

- Corrected a bug that caused an exception to be raised when closing a file
  opened for append, where an HDU was appended to the file, after data was
  accessed from the file.  This exception was only raised when running on a
  Windows platform.

- Updated the installation scripts, compression source code, and benchmark
  test scripts to properly install, build, and execute on a Windows platform.


2.1 (2009-04-14)
------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following enhancements were made:

- Added new tdump and tcreate capabilities to pyfits.

  - The new tdump convenience function allows the contents of a binary table
    HDU to be dumped to a set of three files in ASCII format.  One file will
    contain column definitions, the second will contain header parameters, and
    the third will contain header data.

  - The new tcreate convenience function allows the creation of a binary table
    HDU from the three files dumped by the tdump convenience function.

  - The primary use for the tdump/tcreate methods are to allow editing in a
    standard text editor of the binary table data and parameters.

- Added support for case sensitive values of the EXTNAME card in an extension
  header.  (CNSHD745784)

  - By default, pyfits converts the value of EXTNAME cards to upper case when
    reading from a file.  A new convenience function
    (setExtensionNameCaseSensitive) was implemented to allow a user to
    circumvent this behavior so that the EXTNAME value remains in the same
    case as it is in the file.

  - With the following function call, pyfits will maintain the case of all
    characters in the EXTNAME card values of all extension HDU's during the
    entire python session, or until another call to the function is made:

      >>> import pyfits
      >>> pyfits.setExtensionNameCaseSensitive()

  - The following function call will return pyfits to its default (all upper
    case) behavior:

      >>> pyfits.setExtensionNameCaseSensitive(False)


- Added support for reading and writing FITS files in which the value of the
  first card in the header is 'SIMPLE=F'.  In this case, the pyfits open
  function returns an HDUList object that contains a single HDU of the new
  type _NonstandardHDU.  The header for this HDU is like a normal header (with
  the exception that the first card contains SIMPLE=F instead of SIMPLE=T).
  Like normal HDU's the reading of the data is delayed until actually
  requested.  The data is read from the file into a string starting from the
  first byte after the header END card and continuing till the end of the
  file.  When written, the header is written, followed by the data string.  No
  attempt is made to pad the data string so that it fills into a standard 2880
  byte FITS block.  (CNSHD744730)

- Added support for FITS files containing  extensions with unknown XTENSION
  card values.  (CNSHD744730)  Standard FITS files support extension HDU's of
  types TABLE, IMAGE, BINTABLE, and A3DTABLE.  Accessing a nonstandard
  extension from a FITS file will now create a _NonstandardExtHDU object.
  Accessing the data of this object will cause the data to be read from the
  file into a string.  If the HDU is written back to a file the string data is
  written after the Header and padded to fill a standard 2880 byte FITS block.

The following bugs were fixed:

- Extensive changes were made to the tiled image compression code to support
  the latest enhancements made in CFITSIO version 3.13 to support this
  convention.

- Eliminated a memory leak in the tiled image compression code.

- Corrected a bug in the FITS_record.__setitem__ method which raised a
  NameError exception when attempting to set a value in a FITS_record object.
  (CNSHD745844)

- Corrected a bug that caused a TypeError exception to be raised when reading
  fits files containing large table HDU's (>2Gig).  (CNSHD745522)

- Corrected a bug that caused a TypeError exception to be raised for all calls
  to the warnings module when running under Python 2.6.  The formatwarning
  method in the warnings module was changed in Python 2.6 to include a new
  argument.  (CNSHD746592)

- Corrected the behavior of the membership (in) operator in the Header class
  to check against header card keywords instead of card values.  (CNSHD744730)

- Corrected the behavior of iteration on a Header object.  The new behavior
  iterates over the unique card keywords instead of the card values.


2.0.1 (2009-02-03)
--------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following bugs were fixed:

- Eliminated a memory leak when reading Table HDU's from a fits file.
  (CNSHD741877)


2.0 (2009-01-30)
------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following enhancements were made:

- Provide initial support for an image compression convention known as the
  "Tiled Image Compression Convention" `[1]`_.

  - The principle used in this convention is to first divide the n-dimensional
    image into a rectangular grid of subimages or "tiles".  Each tile is then
    compressed as a continuous block of data, and the resulting compressed
    byte stream is stored in a row of a variable length column in a FITS
    binary table.  Several commonly used algorithms for compressing image
    tiles are supported.  These include, GZIP, RICE, H-Compress and IRAF pixel
    list (PLIO).

  - Support for compressed image data is provided using the optional
    "pyfitsComp" module contained in a C shared library (pyfitsCompmodule.so).

  - The header of a compressed image HDU appears to the user like any image
    header.  The actual header stored in the FITS file is that of a binary
    table HDU with a set of special keywords, defined by the convention, to
    describe the structure of the compressed image.  The conversion between
    binary table HDU header and image HDU header is all performed behind the
    scenes.  Since the HDU is actually a binary table, it may not appear as a
    primary HDU in a FITS file.

  - The data of a compressed image HDU appears to the user as standard
    uncompressed image data.  The actual data is stored in the fits file as
    Binary Table data containing at least one column (COMPRESSED_DATA).  Each
    row of this variable-length column contains the byte stream that was
    generated as a result of compressing the corresponding image tile.
    Several optional columns may also appear.  These include,
    UNCOMPRESSED_DATA to hold the uncompressed pixel values for tiles that
    cannot be compressed, ZSCALE and ZZERO to hold the linear scale factor and
    zero point offset which may be needed to transform the raw uncompressed
    values back to the original image pixel values, and ZBLANK to hold the
    integer value used to represent undefined pixels (if any) in the image.

  - To create a compressed image HDU from scratch, simply construct a
    CompImageHDU object from an uncompressed image data array and its
    associated image header.  From there, the HDU can be treated just like any
    image HDU:

      >>> hdu=pyfits.CompImageHDU(imageData,imageHeader)
      >>> hdu.writeto('compressed_image.fits')

  - The signature for the CompImageHDU initializer method describes the
    possible options for constructing a CompImageHDU object::

      def __init__(self, data=None, header=None, name=None,
                   compressionType='RICE_1',
                   tileSize=None,
                   hcompScale=0.,
                   hcompSmooth=0,
                   quantizeLevel=16.):
          """
              data:            data of the image
              header:          header to be associated with the
                               image
              name:            the EXTNAME value; if this value
                               is None, then the name from the
                               input image header will be used;
                               if there is no name in the input
                               image header then the default name
                               'COMPRESSED_IMAGE' is used
              compressionType: compression algorithm 'RICE_1',
                               'PLIO_1', 'GZIP_1', 'HCOMPRESS_1'
              tileSize:        compression tile sizes default
                               treats each row of image as a tile
              hcompScale:      HCOMPRESS scale parameter
              hcompSmooth:     HCOMPRESS smooth parameter
              quantizeLevel:   floating point quantization level;
          """

- Added two new convenience functions.  The setval function allows the setting
  of the value of a single header card in a fits file.  The delval function
  allows the deletion of a single header card in a fits file.

- A modification was made to allow the reading of data from a fits file
  containing a Table HDU that has duplicate field names.  It is normally a
  requirement that the field names in a Table HDU be unique.  Prior to this
  change a ValueError was raised, when the data was accessed, to indicate that
  the HDU contained duplicate field names.  Now, a warning is issued and the
  field names are made unique in the internal record array.  This will not
  change the TTYPEn header card values.  You will be able to get the data from
  all fields using the field name, including the first field containing the
  name that is duplicated.  To access the data of the other fields with the
  duplicated names you will need to use the field number instead of the field
  name.  (CNSHD737193)

- An enhancement was made to allow the reading of unsigned integer 16 values
  from an ImageHDU when the data is signed integer 16 and BZERO is equal to
  32784 and BSCALE is equal to 1 (the standard way for scaling unsigned
  integer 16 data).  A new optional keyword argument (uint16) was added to the
  open convenience function.  Supplying a value of True for this argument will
  cause data of this type to be read in and scaled into an unsigned integer 16
  array, instead of a float 32 array.  If a HDU associated with a file that
  was opened with the uint16 option and containing unsigned integer 16 data is
  written to a file, the data will be reverse scaled into an integer 16 array
  and written out to the file and the BSCALE/BZERO header cards will be
  written with the values 1 and 32768 respectively.  (CHSHD736064) Reference
  the following example:

    >>> import pyfits
    >>> hdul=pyfits.open('o4sp040b0_raw.fits',uint16=1)
    >>> hdul[1].data
    array([[1507, 1509, 1505, ..., 1498, 1500, 1487],
           [1508, 1507, 1509, ..., 1498, 1505, 1490],
           [1505, 1507, 1505, ..., 1499, 1504, 1491],
           ...,
           [1505, 1506, 1507, ..., 1497, 1502, 1487],
           [1507, 1507, 1504, ..., 1495, 1499, 1486],
           [1515, 1507, 1504, ..., 1492, 1498, 1487]], dtype=uint16)
    >>> hdul.writeto('tmp.fits')
    >>> hdul1=pyfits.open('tmp.fits',uint16=1)
    >>> hdul1[1].data
    array([[1507, 1509, 1505, ..., 1498, 1500, 1487],
           [1508, 1507, 1509, ..., 1498, 1505, 1490],
           [1505, 1507, 1505, ..., 1499, 1504, 1491],
           ...,
           [1505, 1506, 1507, ..., 1497, 1502, 1487],
           [1507, 1507, 1504, ..., 1495, 1499, 1486],
           [1515, 1507, 1504, ..., 1492, 1498, 1487]], dtype=uint16)
    >>> hdul1=pyfits.open('tmp.fits')
    >>> hdul1[1].data
    array([[ 1507.,  1509.,  1505., ...,  1498.,  1500.,  1487.],
           [ 1508.,  1507.,  1509., ...,  1498.,  1505.,  1490.],
           [ 1505.,  1507.,  1505., ...,  1499.,  1504.,  1491.],
           ...,
           [ 1505.,  1506.,  1507., ...,  1497.,  1502.,  1487.],
           [ 1507.,  1507.,  1504., ...,  1495.,  1499.,  1486.],
           [ 1515.,  1507.,  1504., ...,  1492.,  1498.,  1487.]], dtype=float32)

- Enhanced the message generated when a ValueError exception is raised when
  attempting to access a header card with an unparsable value.  The message
  now includes the Card name.

The following bugs were fixed:

- Corrected a bug that occurs when appending a binary table HDU to a fits
  file.  Data was not being byteswapped on little endian machines.
  (CNSHD737243)

- Corrected a bug that occurs when trying to write an ImageHDU that is missing
  the required PCOUNT card in the header.  An UnboundLocalError exception
  complaining that the local variable 'insert_pos' was referenced before
  assignment was being raised in the method _ValidHDU.req_cards.  The code was
  modified so that it would properly issue a more meaningful ValueError
  exception with a description of what required card is missing in the header.

- Eliminated a redundant warning message about the PCOUNT card when validating
  an ImageHDU header with a PCOUNT card that is missing or has a value other
  than 0.

.. _[1]: http://fits.gsfc.nasa.gov/registry/tilecompression.html


1.4.1 (2008-11-04)
--------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following enhancements were made:

- Enhanced the way import errors are reported to provide more information.

The following bugs were fixed:

- Corrected a bug that occurs when a card value is a string and contains a
  colon but is not a record-valued keyword card.

- Corrected a bug where pyfits fails to properly handle a record-valued
  keyword card with values using exponential notation and trailing blanks.


1.4 (2008-07-07)
------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following enhancements were made:

- Added support for file objects and file like objects.

  - All convenience functions and class methods that take a file name will now
    also accept a file object or file like object.  File like objects
    supported are StringIO and GzipFile objects.  Other file like objects will
    work only if they implement all of the standard file object methods.

  - For the most part, file or file like objects may be either opened or
    closed at function call.  An opened object must be opened with the proper
    mode depending on the function or method called.  Whenever possible, if
    the object is opened before the method is called, it will remain open
    after the call.  This will not be possible when writing a HDUList that has
    been resized or when writing to a GzipFile object regardless of whether it
    is resized.  If the object is closed at the time of the function call,
    only the name from the object is used, not the object itself.  The pyfits
    code will extract the file name used by the object and use that to create
    an underlying file object on which the function will be performed.

- Added support for record-valued keyword cards as introduced in the "FITS WCS
  Paper IV proposal for representing a more general distortion model".

  - Record-valued keyword cards are string-valued cards where the string is
    interpreted as a definition giving a record field name, and its floating
    point value.  In a FITS header they have the following syntax::

      keyword= 'field-specifier: float'

    where keyword is a standard eight-character FITS keyword name, float is
    the standard FITS ASCII representation of a floating point number, and
    these are separated by a colon followed by a single blank.

    The grammar for field-specifier is::

      field-specifier:
          field
          field-specifier.field

      field:
          identifier
          identifier.index

    where identifier is a sequence of letters (upper or lower case),
    underscores, and digits of which the first character must not be a digit,
    and index is a sequence of digits.  No blank characters may occur in the
    field-specifier.  The index is provided primarily for defining array
    elements though it need not be used for that purpose.

    Multiple record-valued keywords of the same name but differing values may
    be present in a FITS header.  The field-specifier may be viewed as part of
    the keyword name.

    Some examples follow::

      DP1     = 'NAXIS: 2'
      DP1     = 'AXIS.1: 1'
      DP1     = 'AXIS.2: 2'
      DP1     = 'NAUX: 2'
      DP1     = 'AUX.1.COEFF.0: 0'
      DP1     = 'AUX.1.POWER.0: 1'
      DP1     = 'AUX.1.COEFF.1: 0.00048828125'
      DP1     = 'AUX.1.POWER.1: 1'

  - As with standard header cards, the value of a record-valued keyword card
    can be accessed using either the index of the card in a HDU's header or
    via the keyword name.  When accessing using the keyword name, the user may
    specify just the card keyword or the card keyword followed by a period
    followed by the field-specifier.  Note that while the card keyword is case
    insensitive, the field-specifier is not.  Thus, hdu['abc.def'],
    hdu['ABC.def'], or hdu['aBc.def'] are all equivalent but hdu['ABC.DEF'] is
    not.

  - When accessed using the card index of the HDU's header the value returned
    will be the entire string value of the card.  For example:

      >>> print hdr[10]
      NAXIS: 2
      >>> print hdr[11]
      AXIS.1: 1

  - When accessed using the keyword name exclusive of the field-specifier, the
    entire string value of the header card with the lowest index having that
    keyword name will be returned.  For example:

      >>> print hdr['DP1']
      NAXIS: 2

  - When accessing using the keyword name and the field-specifier, the value
    returned will be the floating point value associated with the
    record-valued keyword card.  For example:

      >>> print hdr['DP1.NAXIS']
      2.0

  - Any attempt to access a non-existent record-valued keyword card value will
    cause an exception to be raised (IndexError exception for index access or
    KeyError for keyword name access).

  - Updating the value of a record-valued keyword card can also be
    accomplished using either index or keyword name.  For example:

      >>> print hdr['DP1.NAXIS']
      2.0
      >>> hdr['DP1.NAXIS'] = 3.0
      >>> print hdr['DP1.NAXIS']
      3.0

  - Adding a new record-valued keyword card to an existing header is
    accomplished using the Header.update() method just like any other card.
    For example:

      >>> hdr.update('DP1', 'AXIS.3: 1', 'a comment', after='DP1.AXIS.2')

  - Deleting a record-valued keyword card from an existing header is
    accomplished using the standard list deletion syntax just like any other
    card.  For example:

      >>> del hdr['DP1.AXIS.1']

  - In addition to accessing record-valued keyword cards individually using a
    card index or keyword name, cards can be accessed in groups using a set of
    special pattern matching keys.  This access is made available via the
    standard list indexing operator providing a keyword name string that
    contains one or more of the special pattern matching keys.  Instead of
    returning a value, a CardList object will be returned containing shared
    instances of the Cards in the header that match the given keyword
    specification.

  - There are three special pattern matching keys.  The first key '*' will
    match any string of zero or more characters within the current level of
    the field-specifier.  The second key '?' will match a single character.
    The third key '...' must appear at the end of the keyword name string and
    will match all keywords that match the preceding pattern down all levels
    of the field-specifier.  All combinations of ?, \*, and ... are permitted
    (though ... is only permitted at the end).  Some examples follow:

      >>> cl=hdr['DP1.AXIS.*']
      >>> print cl
      DP1     = 'AXIS.1: 1'
      DP1     = 'AXIS.2: 2'
      >>> cl=hdr['DP1.*']
      >>> print cl
      DP1     = 'NAXIS: 2'
      DP1     = 'NAUX: 2'
      >>> cl=hdr['DP1.AUX...']
      >>> print cl
      DP1     = 'AUX.1.COEFF.0: 0'
      DP1     = 'AUX.1.POWER.0: 1'
      DP1     = 'AUX.1.COEFF.1: 0.00048828125'
      DP1     = 'AUX.1.POWER.1: 1'
      >>> cl=hdr['DP?.NAXIS']
      >>> print cl
      DP1     = 'NAXIS: 2'
      DP2     = 'NAXIS: 2'
      DP3     = 'NAXIS: 2'
      >>> cl=hdr['DP1.A*S.*']
      >>> print cl
      DP1     = 'AXIS.1: 1'
      DP1     = 'AXIS.2: 2'

  - The use of the special pattern matching keys for adding or updating header
    cards in an existing header is not allowed.  However, the deletion of
    cards from the header using the special keys is allowed.  For example:

      >>> del hdr['DP3.A*...']

- As noted above, accessing pyfits Header object using the special pattern
  matching keys will return a CardList object.  This CardList object can
  itself be searched in order to further refine the list of Cards.  For
  example:

      >>> cl=hdr['DP1...']
      >>> print cl
      DP1     = 'NAXIS: 2'
      DP1     = 'AXIS.1: 1'
      DP1     = 'AXIS.2: 2'
      DP1     = 'NAUX: 2'
      DP1     = 'AUX.1.COEFF.1: 0.000488'
      DP1     = 'AUX.2.COEFF.2: 0.00097656'
      >>> cl1=cl['*.*AUX...']
      >>> print cl1
      DP1     = 'NAUX: 2'
      DP1     = 'AUX.1.COEFF.1: 0.000488'
      DP1     = 'AUX.2.COEFF.2: 0.00097656'

  - The CardList keys() method will allow the retrieval of all of the key
    values in the CardList.  For example:

      >>> cl=hdr['DP1.AXIS.*']
      >>> print cl
      DP1     = 'AXIS.1: 1'
      DP1     = 'AXIS.2: 2'
      >>> cl.keys()
      ['DP1.AXIS.1', 'DP1.AXIS.2']

  - The CardList values() method will allow the retrieval of all of the values
    in the CardList.  For example:

      >>> cl=hdr['DP1.AXIS.*']
      >>> print cl
      DP1     = 'AXIS.1: 1'
      DP1     = 'AXIS.2: 2'
      >>> cl.values()
      [1.0, 2.0]

  - Individual cards can be retrieved from the list using standard list
    indexing.  For example:

      >>> cl=hdr['DP1.AXIS.*']
      >>> c=cl[0]
      >>> print c
      DP1     = 'AXIS.1: 1'
      >>> c=cl['DP1.AXIS.2']
      >>> print c
      DP1     = 'AXIS.2: 2'

  - Individual card values can be retrieved from the list using the value
    attribute of the card.  For example:

      >>> cl=hdr['DP1.AXIS.*']
      >>> cl[0].value
      1.0

  - The cards in the CardList are shared instances of the cards in the source
    header.  Therefore, modifying a card in the CardList also modifies it in
    the source header.  However, making an addition or a deletion to the
    CardList will not affect the source header.  For example:

      >>> hdr['DP1.AXIS.1']
      1.0
      >>> cl=hdr['DP1.AXIS.*']
      >>> cl[0].value = 4.0
      >>> hdr['DP1.AXIS.1']
      4.0
      >>> del cl[0]
      >>> print cl['DP1.AXIS.1']
      Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "NP_pyfits.py", line 977, in __getitem__
        return self.ascard[key].value
      File "NP_pyfits.py", line 1258, in __getitem__
        _key = self.index_of(key)
      File "NP_pyfits.py", line 1403, in index_of
        raise KeyError, 'Keyword %s not found.' % `key`
      KeyError: "Keyword 'DP1.AXIS.1' not found."
      >>> hdr['DP1.AXIS.1']
      4.0

  - A FITS header consists of card images.  In pyfits each card image is
    manifested by a Card object.  A pyfits Header object contains a list of
    Card objects in the form of a CardList object.  A record-valued keyword
    card image is represented in pyfits by a RecordValuedKeywordCard object.
    This object inherits from a Card object and has all of the methods and
    attributes of a Card object.

  - A new RecordValuedKeywordCard object is created with the
    RecordValuedKeywordCard constructor: RecordValuedKeywordCard(key, value,
    comment).  The key and value arguments may be specified in two ways.  The
    key value may be given as the 8 character keyword only, in which case the
    value must be a character string containing the field-specifier, a colon
    followed by a space, followed by the actual value.  The second option is
    to provide the key as a string containing the keyword and field-specifier,
    in which case the value must be the actual floating point value.  For
    example:

      >>> c1 = pyfits.RecordValuedKeywordCard('DP1', 'NAXIS: 2', 'Number of variables')
      >>> c2 = pyfits.RecordValuedKeywordCard('DP1.AXIS.1', 1.0, 'Axis number')

  - RecordValuedKeywordCards have attributes .key, .field_specifier, .value,
    and .comment.  Both .value and .comment can be changed but not .key or
    .field_specifier.  The constructor will extract the field-specifier from
    the input key or value, whichever is appropriate.  The .key attribute is
    the 8 character keyword.

  - Just like standard Cards, a RecordValuedKeywordCard may be constructed
    from a string using the fromstring() method or verified using the verify()
    method.  For example:

      >>> c1 = pyfits.RecordValuedKeywordCard().fromstring(
               "DP1     = 'NAXIS: 2' / Number of independent variables")
      >>> c2 = pyfits.RecordValuedKeywordCard().fromstring(
               "DP1     = 'AXIS.1: X' / Axis number")
      >>> print c1; print c2
      DP1     = 'NAXIS: 2' / Number of independent variables
      DP1     = 'AXIS.1: X' / Axis number
      >>> c2.verify()
      Output verification result:
      Card image is not FITS standard (unparsable value string).

  - A standard card that meets the criteria of a RecordValuedKeywordCard may
    be turned into a RecordValuedKeywordCard using the class method coerce.
    If the card object does not meet the required criteria then the original
    card object is just returned.

      >>> c1 = pyfits.Card('DP1','AUX: 1','comment')
      >>> c2 = pyfits.RecordValuedKeywordCard.coerce(c1)
      >>> print type(c2)
      <'pyfits.NP_pyfits.RecordValuedKeywordCard'>

  - Two other card creation methods are also available as
    RecordVauedKeywordCard class methods.  These are createCard() which will
    create the appropriate card object (Card or RecordValuedKeywordCard) given
    input key, value, and comment, and createCardFromString which will create
    the appropriate card object given an input string.  These two methods are
    also available as convenience functions:

      >>> c1 = pyfits.RecordValuedKeywordCard.createCard('DP1','AUX: 1','comment)

    or

      >>> c1 = pyfits.createCard('DP1','AUX: 1','comment)
      >>> print type(c1)
      <'pyfits.NP_pyfits.RecordValuedKeywordCard'>

      >>> c1 = pyfits.RecordValuedKeywordCard.createCard('DP1','AUX 1','comment)

    or

      >>> c1 = pyfits.createCard('DP1','AUX 1','comment)
      >>> print type(c1)
      <'pyfits.NP_pyfits.Card'>

      >>> c1 = pyfits.RecordValuedKeywordCard.createCardFromString \
               ("DP1 = 'AUX: 1.0' / comment")

    or

      >>> c1 = pyfits.createCardFromString("DP1     = 'AUX: 1.0' / comment")
      >>> print type(c1)
      <'pyfits.NP_pyfits.RecordValuedKeywordCard'>

The following bugs were fixed:

- Corrected a bug that occurs when writing a HDU out to a file.  During the
  write, any Keyboard Interrupts are trapped so that the write completes
  before the interrupt is handled.  Unfortunately, the Keyboard Interrupt was
  not properly reinstated after the write completed.  This was fixed.
  (CNSHD711138)

- Corrected a bug when using ipython, where temporary files created with the
  tempFile.NamedTemporaryFile method are not automatically removed.  This can
  happen for instance when opening a Gzipped fits file or when open a fits
  file over the internet.  The files will now be removed.  (CNSHD718307)

- Corrected a bug in the append convenience function's call to the writeto
  convenience function.  The classExtensions argument must be passed as a
  keyword argument.

- Corrected a bug that occurs when retrieving variable length character arrays
  from binary table HDUs (PA() format) and using slicing to obtain rows of
  data containing variable length arrays.  The code issued a TypeError
  exception.  The data can now be accessed with no exceptions. (CNSHD718749)

- Corrected a bug that occurs when retrieving data from a fits file opened in
  memory map mode when the file contains multiple image extensions or ASCII
  table or binary table HDUs.  The code issued a TypeError exception.  The
  data can now be accessed with no exceptions.  (CNSHD707426)

- Corrected a bug that occurs when attempting to get a subset of data from a
  Binary Table HDU and then use the data to create a new Binary Table HDU
  object.  A TypeError exception was raised.  The data can now be subsetted
  and used to create a new HDU.  (CNSHD723761)

- Corrected a bug that occurs when attempting to scale an Image HDU back to
  its original data type using the _ImageBaseHDU.scale method.  The code was
  not resetting the BITPIX header card back to the original data type.  This
  has been corrected.

- Changed the code to issue a KeyError exception instead of a NameError
  exception when accessing a non-existent field in a table.


1.3 (2008-02-22)
------------------

Updates described in this release are only supported in the NUMPY version of
pyfits.

The following enhancements were made:

- Provided support for a new extension to pyfits called *stpyfits*.

  - The *stpyfits* module is a wrapper around pyfits.  It provides all of the
    features and functions of pyfits along with some STScI specific features.
    Currently, the only new feature supported by stpyfits is the ability to
    read and write fits files that contain image data quality extensions with
    constant data value arrays.  See stpyfits `[2]`_ for more details on
    stpyfits.

- Added a new feature to allow trailing HDUs to be deleted from a fits file
  without actually reading the data from the file.

  - This supports a JWST requirement to delete a trailing HDU from a file
    whose primary Image HDU is too large to be read on a 32 bit machine.

- Updated pyfits to use the warnings module to issue warnings.  All warnings
  will still be issued to stdout, exactly as they were before, however, you
  may now suppress warnings with the -Wignore command line option.  For
  example, to run a script that will ignore warnings use the following command
  line syntax:

    python -Wignore yourscript.py

- Updated the open convenience function to allow the input of an already
  opened file object in place of a file name when opening a fits file.

- Updated the writeto convenience function to allow it to accept the
  output_verify option.

  - In this way, the user can use the argument output_verify='fix' to allow
    pyfits to correct any errors it encounters in the provided header before
    writing the data to the file.

- Updated the verification code to provide additional detail with a
  VerifyError exception.

- Added the capability to create a binary table HDU directly from a
  numpy.ndarray.  This may be done using either the new_table convenience
  function or the BinTableHDU constructor.


The following performance improvements were made:

- Modified the import logic to dramatically decrease the time it takes to
  import pyfits.

- Modified the code to provide performance improvements when copying and
  examining header cards.

The following bugs were fixed:

- Corrected a bug that occurs when reading the data from a fits file that
  includes BZERO/BSCALE scaling.  When the data is read in from the file,
  pyfits automatically scales the data using the BZERO/BSCALE values in the
  header.  In the previous release, pyfits created a 32 bit floating point
  array to hold the scaled data.  This could cause a problem when the value of
  BZERO is so large that the scaled value will not fit into the float 32.  For
  this release, when the input data is 32 bit integer, a 64 bit floating point
  array is used for the scaled data.

- Corrected a bug that caused an exception to be raised when attempting to
  scale image data using the ImageHDU.scale method.

- Corrected a bug in the new_table convenience function that occurred when a
  binary table was created using a ColDefs object as input and supplying an
  nrows argument for a number of rows that is greater than the number of rows
  present in the input ColDefs object.  The previous version of pyfits failed
  to allocate the necessary memory for the additional rows.

- Corrected a bug in the new_table convenience function that caused an
  exception to be thrown when creating an ASCII table.

- Corrected a bug in the new_table convenience function that will allow the
  input of a ColDefs object that was read from a file as a binary table with a
  data value equal to None.

- Corrected a bug in the construction of ASCII tables from Column objects that
  are created with noncontinuous start columns.

- Corrected bugs in a number of areas that would sometimes cause a failure to
  improperly raise an exception when an error occurred.

- Corrected a bug where attempting to open a non-existent fits file on a
  windows platform using a drive letter in the file specification caused a
  misleading IOError exception to be raised.

.. _[2]: http://stsdas.stsci.edu/stsci_python_sphinxdocs_2.13/tools/stpyfits.html


1.1 (2007-06-15)
------------------

- Modified to use either NUMPY or NUMARRAY.

- New file writing modes have been provided to allow streaming data to
  extensions without requiring the whole output extension image in memory. See
  documentation on StreamingHDU.

- Improvements to minimize byteswapping and memory usage by byteswapping in
  place.

- Now supports ':' characters in filenames.

- Handles keyboard interrupts during long operations.

- Preserves the byte order of the input image arrays.


1.0.1 (2006-03-24)
--------------------

The changes to PyFITS were primarily to improve the docstrings and to
reclassify some public functions and variables as private. Readgeis and
fitsdiff which were distributed with PyFITS in previous releases were moved to
pytools. This release of PyFITS is v1.0.1. The next release of PyFITS will
support both numarray and numpy (and will be available separately from
stsci_python, as are all the python packages contained within stsci_python).
An alpha release for PyFITS numpy support will be made around the time of this
stsci_python release.

- Updated docstrings for public functions.

- Made some previously public functions private.


1.0 (2005-11-01)
------------------

Major Changes since v0.9.6:

- Added support for the HIERARCH convention

- Added support for iteration and slicing for HDU lists

- PyFITS now uses the standard setup.py installation script

- Add utility functions at the module level, they include:

  - getheader
  - getdata
  - getval
  - writeto
  - append
  - update
  - info

Minor changes since v0.9.6:

- Fix a bug to make single-column ASCII table work.

- Fix a bug so a new table constructed from an existing table with X-formatted
  columns will work.

- Fix a problem in verifying HDUList right after the open statement.

- Verify that elements in an HDUList, besides the first one, are ExtensionHDU.

- Add output verification in methods flush() and close().

- Modify the design of the open() function to remove the output_verify
  argument.

- Remove the groups argument in GroupsHDU's constructor.

- Redesign the column definition class to make its column components more
  accessible.  Also to make it conducive for higher level functionalities,
  e.g. combining two column definitions.

- Replace the Boolean class with the Python Boolean type.  The old TRUE/FALSE
  will still work.

- Convert classes to the new style.

- Better format when printing card or card list.

- Add the optional argument clobber to all writeto() functions and methods.

- If adding a blank card, will not use existing blank card's space.

PyFITS Version 1.0 REQUIRES Python 2.3 or later.


0.9.6 (2004-11-11)
--------------------

Major changes since v0.9.3:

- Support for variable length array tables.

- Support for writing ASCII table extensions.

- Support for random groups, both reading and writing.

Some minor changes:

- Support for numbers with leading zeros in an ASCII table extension.

- Changed scaled columns' data type from Float32 to Float64 to preserve
  precision.

- Made Column constructor more flexible in accepting format specification.


0.9.3 (2004-07-02)
--------------------

Changes since v0.9.0:

- Lazy instanciation of full Headers/Cards for all HDU's when the file is
  opened.  At the open, only extracts vital info (e.g. NAXIS's) from the
  header parts.  This change will speed up the performance if the user only
  needs to access one extension in a multi-extension FITS file.

- Support the X format (bit flags) columns, both reading and writing, in a
  binary table.  At the user interface, they are converted to Boolean arrays
  for easy manipulation.  For example, if the column's TFORM is "11X",
  internally the data is stored in 2 bytes, but the user will see, at each row
  of this column, a Boolean array of 11 elements.

- Fix a bug such that when a table extension has no data, it will not try to
  scale the data when updating/writing the HDU list.


0.9 (2004-04-27)
------------------

Changes since v0.8.0:

- Rewriting of the Card class to separate the parsing and verification of
  header cards

- Restructure the keyword indexing scheme which speed up certain applications
  (update large number of new keywords and reading a header with larger
  numbers of cards) by a factor of 30 or more

- Change the default to be lenient FITS standard checking on input and strict
  FITS standard checking on output

- Support CONTINUE cards, both reading and writing

- Verification can now be performed at any of the HDUList, HDU, and Card
  levels

- Support (contiguous) subsection (attribute .section) of images to reduce
  memory usage for large images


0.8.0 (2003-08-19)
--------------------

**NOTE:** This version will only work with numarray Version 0.6.  In addition,
earlier versions of PyFITS will not work with numarray 0.6.  Therefore, both
must be updated simultaneously.

Changes since 0.7.6:

- Compatible with numarray 0.6/records 2.0

- For binary tables, now it is possible to update the original array if a
  scaled field is updated.

- Support of complex columns

- Modify the __getitem__ method in FITS_rec.  In order to make sure the scaled
  quantities are also viewing the same data as the original FITS_rec, all
  fields need to be "touched" when __getitem__ is called.

- Add a new attribute mmobject for HDUList, and close the memmap object when
  close HDUList object.  Earlier version does not close memmap object and can
  cause memory lockup.

- Enable 'update' as a legitimate memmap mode.

- Do not print message when closing an HDUList object which is not created
  from reading a FITS file.  Such message is confusing.

- remove the internal attribute "closed" and related method (__getattr__ in
  HDUList).  It is redundant.


0.7.6 (2002-11-22)

**NOTE:** This version will only work with numarray Version 0.4.

Changes since 0.7.5:

- Change x*=n to numarray.multiply(x, n, x) where n is a floating number, in
  order to make pyfits to work under Python 2.2. (2 occurrences)

- Modify the "update" method in the Header class to use the "fixed-format"
  card even if the card already exists.  This is to avoid the mis-alignment as
  shown below:

  After running drizzle on ACS images it creates a CD matrix whose elements
  have very many digits, *e.g.*:

    CD1_1   =  1.1187596304411E-05 / partial of first axis coordinate w.r.t. x
    CD1_2   = -8.502767249350019E-06 / partial of first axis coordinate w.r.t. y

  with pyfits, an "update" on these header items and write in new values which
  has fewer digits, *e.g.*:

    CD1_1   =        1.0963011E-05 / partial of first axis coordinate w.r.t. x
    CD1_2   =          -8.527229E-06 / partial of first axis coordinate w.r.t. y

- Change some internal variables to make their appearance more consistent:

    old name                new name

    __octalRegex            _octalRegex
    __readblock()           _readblock()
    __formatter()           _formatter().
    __value_RE              _value_RE
    __numr                  _numr
    __comment_RE            _comment_RE
    __keywd_RE              _keywd_RE
    __number_RE             _number_RE.
    tmpName()               _tmpName()
    dimShape                _dimShape
    ErrList                 _ErrList

- Move up the module description.  Move the copyright statement to the bottom
  and assign to the variable __credits__.

- change the following line:

    self.__dict__ = input.__dict__

  to

    self.__setstate__(input.__getstate__())

  in order for pyfits to run under numarray 0.4.

- edit _readblock to add the (optional) firstblock argument and raise IOError
  if the first 8 characters in the first block is not 'SIMPLE  ' or
  'XTENSION'.  Edit the function open to check for IOError to skip the last
  null filled block(s).  Edit readHDU to add the firstblock argument.


0.7.5 (2002-08-16)
--------------------

Changes since v0.7.3:

- Memory mapping now works for readonly mode, both for images and binary
  tables.

  Usage:  pyfits.open('filename', memmap=1)

- Edit the field method in FITS_rec class to make the column scaling for
  numbers use less temporary memory.  (does not work under 2.2, due to Python
  "bug" of array \*=)

- Delete bscale/bzero in the ImageBaseHDU constructor.

- Update bitpix in BaseImageHDU.__getattr__  after deleting bscale/bzero. (bug
  fix)

- In BaseImageHDU.__getattr__  point self.data to raw_data if float and if not
  memmap.  (bug fix).

- Change the function get_tbdata() to private: _get_tbdata().


0.7.3 (2002-07-12)
--------------------

Changes since v0.7.2:

- It will scale all integer image data to Float32, if BSCALE/BZERO != 1/0.  It
  will also expunge the BSCALE/BZERO keywords.

- Add the scale() method for ImageBaseHDU, so data can be scaled just before
  being written to the file.  It has the following arguments:

  type: destination data type (string), e.g. Int32, Float32, UInt8, etc.

  option: scaling scheme. if 'old', use the old BSCALE/BZERO values.  if
  'minmax', use the data range to fit into the full range of specified integer
  type.  Float destination data type will not be scaled for this option.

  bscale/bzero: user specifiable BSCALE/BZERO values.  They overwrite the
  "option".

- Deal with data area resizing in 'update' mode.

- Make the data scaling (both input and output) faster and use less memory.

- Bug fix to make column name change takes effect for field.

- Bug fix to avoid exception if the key is not present in the header already.
  This affects (fixes) add_history(), add_comment(), and add_blank().

- Bug fix in __getattr__() in Card class.  The change made in 0.7.2 to rstrip
  the comment must be string type to avoid exception.

0.7.2.1 (2002-06-25)
----------------------

A couple of bugs were addressed in this version.

- Fix a bug in _add_commentary(). Due to a change in index_of() during version
  0.6.5.5, _add_commentary needs to be modified to avoid exception if the key
  is not present in the header already. This affects (fixes) add_history(),
  add_comment(), and add_blank().

- Fix a bug in __getattr__() in Card class. The change made in 0.7.2 to rstrip
  the comment must be string type to avoid exception.


0.7.2 (2002-06-19)
--------------------

The two major improvements from Version 0.6.2 are:

- support reading tables  with "scaled" columns (e.g.  tscal/tzero, Boolean,
  and ASCII tables)

- a prototype output verification.

This version of PyFITS requires numarray version 0.3.4.

Other changes include:

- Implement the new HDU hierarchy proposed earlier this year.  This in turn
  reduces some of the redundant methods common to several HDU classes.

- Add 3 new methods to the Header class: add_history, add_comment, and
  add_blank.

- The table attributes _columns are now .columns and the attributes in ColDefs
  are now all without the underscores.  So, a user can get a list of column
  names by: hdu.columns.names.

- The "fill" argument in the new_table method now has a new meaning:<br> If
  set to true (=1), it will fill the entire new table with zeros/blanks.
  Otherwise (=0), just the extra rows/cells are filled with zeros/blanks.
  Fill values other than zero/blank are now not possible.

- Add the argument output_verify to the open method and writeto method.  Not
  in the flush or close methods yet, due to possible complication.

- A new copy method for tables, the copy is totally independent from the table
  it copies from.

- The tostring() call in writeHDUdata takes up extra space to store the string
  object.  Use tofile() instead, to save space.

- Make changes from _byteswap to _byteorder, following corresponding changes
  in numarray and recarray.

- Insert(update) EXTEND in PrimaryHDU only when header is None.

- Strip the trailing blanks for the comment value of a card.

- Add seek(0) right after the __buildin__.open(0), because for the 'ab+' mode,
  the pointer is at the end after open in Linux, but it is at the beginning in
  Solaris.

- Add checking of data against header, update header keywords (NAXIS's,
  BITPIX) when they don't agree with the data.

- change version to __version__.

There are also many other minor internal bug fixes and
technical changes.


0.6.2 (2002-02-12)
--------------------

This version requires numarray version 0.2.

Things not yet supported but are part of future development:

- Verification and/or correction of FITS objects being written to disk so that
  they are legal FITS. This is being added now and should be available in
  about a month.  Currently, one may construct FITS headers that are
  inconsistent with the data and write such FITS objects to disk. Future
  versions will provide options to either a) correct discrepancies and warn,
  b) correct discrepancies silently, c) throw a Python exception, or d) write
  illegal FITS (for test purposes!).

- Support for ascii tables or random groups format. Support for ASCII tables
  will be done soon (~1 month). When random group support is added is
  uncertain.

- Support for memory mapping FITS data (to reduce memory demands). We expect
  to provide this capability in about 3 months.

- Support for columns in binary tables having scaled values (e.g. BSCALE or
  BZERO) or boolean values. Currently booleans are stored as Int8 arrays and
  users must explicitly convert them into a boolean array. Likewise, scaled
  columns must be copied with scaling and offset by testing for those
  attributes explicitly. Future versions will produce such copies
  automatically.

- Support for tables with TNULL values. This awaits an enhancement to numarray
  to support mask arrays (planned).  (At least a couple of months off).

.. _PyFITS: http://www.stsci.edu/resources/software_hardware/pyfits
