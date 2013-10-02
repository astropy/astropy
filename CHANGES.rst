0.3 (unreleased)
----------------

New Features
^^^^^^^^^^^^

- ``astropy.coordinates``

  - Two classes `astropy.coordinates.Longitude` and `astropy.coordinates.Latitude`
    have been added.  These are derived from the new `Angle` class and used for
    all longitude-like (RA, azimuth, galactic L) and latitude-like coordinates
    (Dec, elevation, galactic B) respectively.  The `Longitude` class provides
    auto-wrapping capability and `Latitude` performs bounds checking.

  - `astropy.coordinates.Distance` supports conversion to and from distance
    modulii.

- ``astropy.io.votable``

  - The format of the units of a VOTable file can be specified using
    the `unit_format` parameter.  Note that units are still always
    written out using the CDS format, to ensure compatibility with the
    standard.

- ``astropy.time``

  - Update internal time manipulations so that arithmetic with Time and
    TimeDelta objects maintains sub-nanosecond precision over a time span
    longer than the age of the universe [#1189].

  - Use ``astropy.utils.iers`` to provide ``delta_ut1_utc``, so that
    automatic calculation of UT1 becomes possible [#1145].
  
  - Add ``datetime`` format which allows converting to and from standard
    library ``datetime.datetime`` objects.

  - Add ``plot_date`` format which allows converting to and from the date
    representation used when plotting dates with matplotlib via the
    ``matplotlib.pyplot.plot_date`` function.

  - Add ``gps`` format (seconds since 1980-01-01 00:00:00 UTC,
    including leap seconds)

  - Add array indexing to Time objects [#1132]

  - Allow for arithmetic of multi-element and single-element Time and TimeDelta
    objects [#1081].

  - Allow multiplication and division of TimeDelta objects by
    constants and arrays, as well as changing sign (negation) and
    taking the absolute value of TimeDelta objects [#1082].

  - Allow comparisons of Time and TimeDelta objects [#1171].

- ``astropy.stats``

  - Added robust statistics functions `~astropy.stats.funcs.median_absolute_deviation`,
    `~astropy.stats.funcs.biweight_location`, and
    `~astropy.stats.funcs.biweight_midvariance`.

  - Add `axis=int` option to `astropy.stats.funcs.sigma_clip` to allow clipping
    along a given axis for multidimensional data.  [#1083]

- ``astropy.table``

  - Table.read and Table.write now support reading and writing of FITS tables
    via the unified reading/writing interface [#591].

  - The 'units' and 'dtypes' attributes and keyword arguments in Column,
    MaskedColumn, Row, and Table are now deprecated in favor of the
    single-tense 'unit' and 'dtype' [#1174].

- ``astropy.convolution``

  - New class-based system for generating kernels, replacing `make_kernel`.
    [#1255]. The `astropy.nddata.convolution` sub-package has now been moved to
    `astropy.convolution`. [#1451]

- ``astropy.vo``

  - New package added to support Virtual Observatory Simple Cone Search query
    and service validation [#552].

- ``astropy.units``

  - Added new spectroscopic equivalencies for velocity conversions
    (relativistic, optical, and radio conventions are supported)

  - Added Brightness Temperature (antenna gain) equivalency for conversion
    between :math:`T_B` and flux density.

  - Added percent unit, and allowed any string containing just a number 
    to be interpreted as a scaled dimensionless unit.

- ``astropy.utils``

  - Added ``astropy.utils.iers`` which allows reading in of IERS A or
    IERS B bulletins and interpolation in UT1-UTC.

- ``astropy.extern.six`` 

  - Added `six <https://pypi.python.org/pypi/six/>`_ for python2/python3 
    compatibility

- Astropy now uses the ERFA library instead of the IAU SOFA library for
  fundamental time transformation routines.  
  The ERFA library is derived, with permission, from the IAU SOFA library but 
  is distributed under a BSD license. See ``license/ERFA.rst`` for details.

- ``astropy.logger``

  - The Astropy logger now no longer catches exceptions by default, and also
    only captures warnings emitted by Astropy itself (prior to this change,
    following an import of Astropy, any warning got re-directed through the
    Astropy logger). Logging to the Astropy log file has also been disabled by
    default. However, users of Astropy 0.2 will likely still see the previous
    behavior with Astropy 0.3 for exceptions and logging to file since the
    default configuration file installed by 0.2 set the exception logging to be
    on by default. To get the new behavior, set the ``log_exceptions`` and
    ``log_to_file`` configuration items to ``False`` in the ``astropy.cfg``
    file.

API Changes
^^^^^^^^^^^

- ``astropy.coordinates``

  - The `astropy.coordinates.Angle` class is now a subclass of
    `astropy.units.Quantity`. This means it has all of the methods of a
    `numpy.ndarray`.


  - The `astropy.coordinates.Distance` class is now a subclass of
    `astropy.units.Quantity`. This means it has all of the methods of a
    `numpy.ndarray`.

    - All angular units are now supported, not just `radian`, `degree`
      and `hour`, but now `arcsecond` and `arcminute` as well.  The
      object will retain its native unit, so when printing out a value
      initially provided in hours, its `to_string()` will, by default,
      also be expressed in hours.

    - The `Angle` class now supports arrays of angles.

    - To be consistent with `units.Unit`, `Angle.format` has been
      deprecated and renamed to `Angle.to_string`.

    - To be consistent with `astropy.units`, all plural forms of unit
      names have been removed.  Therefore, the following properties of
      `astropy.coordinates.Angle` should be renamed:

      - ``radians`` -> ``radian``
      - ``degrees`` -> ``degree``
      - ``hours`` -> ``hour``

    - Multiplication and division of two `Angle` objects used to raise
      `NotImplementedError`.  Now they raise `TypeError`.

  - The `astropy.coordinates.Angle` class no longer has a ``bounds`` attribute
    so there is no bounds-checking or auto-wrapping at this level.  This allows
    ``Angle`` objects to be used in arbitrary arithmetic expressions
    (e.g. coordinate distance computation).

  - The `astropy.coordinates.RA`and `astropy.coordinates.Dec` classes have
    been removed and replaced with `astropy.coordinates.Longitude` and
    `astropy.coordinates.Latitude` respectively.  These are now used for
    the components of Galactic and Horizontal (Alt-Az) coordinates as well
    instead of plain `Angle` objects.

  - `astropy.coordinates.angles.rotation_matrix` and
    `astropy.coordinates.angles.angle_axis` now take a `unit` kwarg
    instead of `degrees` kwarg to specify the units of the angles.
    `rotation_matrix` will also take the unit from the given `Angle`
    object if no unit is provided.

  - The `astropy.coordinates.AngularSeparation` class has been removed.  The output
    of the coordinates `separation()` method is now an `astropy.coordinates.Angle`.

- ``astropy.io.ascii``

  - In the ``read`` method of ``astropy.io.ascii``, empty column values in an ASCII table
    are now treated as missing values instead of the previous treatment as a zero-length
    string "".  This now corresponds to the behavior of other table readers like
    ``numpy.genfromtxt``.  To restore the previous behavior set ``fill_values=None`` in the
    call to ``ascii.read()``.

  - The ``read`` and ``write`` methods of ``astropy.io.ascii`` now have a ``format``
    argument for specifying the file format.  This is the preferred way to choose
    the format instead of the ``Reader`` and ``Writer`` arguments [#961].

- ``astropy.io.fits``

  - The ``updateHeader``, ``updateHeaderData``, and ``updateCompressedData``
    methods of the ``CompDataHDU`` class are pending deprecation and moved to
    internal methods.  The operation of these methods depended too much on
    internal state to be used safely by users; instead they are invoked
    automatically in the appropriate places when reading/writing compressed
    image HDUs.

  - The ``CompDataHDU.compData`` attribute is pending deprecation in favor of
    the clearer and more PEP-8 compatible ``CompDataHDU.compressed_data``.

  - The constructor for ``CompDataHDU`` has been changed to accept new keyword
    arguments.  The new keyword arguments are essentially the same, but are in
    underscore_separated format rather than camelCase format.  The old
    arguments are still pending deprecation.

  - The internal attributes of HDU classes ``_hdrLoc``, ``_datLoc``, and
    ``_datSpan`` have been replaced with ``_header_offset``, ``_data_offset``,
    and ``_data_size`` respectively.  The old attribute names are still pending
    deprecation.  This should only be of interest to advanced users who have
    created their own HDU subclasses.

  - The following previously deprecated functions and methods have been removed
    entirely: ``createCard``, ``createCardFromString``, ``upperKey``,
    ``ColDefs.data``, ``setExtensionNameCaseSensitive``, ``_File.getfile``,
    ``_TableBaseHDU.get_coldefs``, ``Header.has_key``, ``Header.ascardlist``.

  - Interfaces that were pending deprecation are now fully deprecated.  These
    include: ``create_card``, ``create_card_from_string``, ``upper_key``,
    ``Header.get_history``, and ``Header.get_comment``.

- ``astropy.nddata``

  - The `astropy.nddata.convolution` sub-package has now been moved to
    `astropy.convolution`. [#1451]

- ``astropy.stats.funcs``

  - For ``sigma_clip``, the ``maout`` optional parameter has been removed, and
    the function now always returns a masked array.  A new boolean parameter
    ``copy`` can be used to indicated whether the input data should be copied
    (``copy=True``, default) or used by reference (``copy=False``) in the
    output masked array.  [#1083]

- ``astropy.table``

  - The plural 'units' and 'dtypes' have been switched to 'unit' and 'dtype'
    where appropriate. The original attributes are still present in this version
    as deprecated attributes, but will be removed in the next version [#1174].

  - The ``copy`` methods of ``Column`` and ``MaskedColumn`` were changed so that
    the first argument is now ``order='C'``.  This is required for compatibility
    with Numpy 1.8 which is currently in development [#1250].

- ``astropy.time``

  - For consistency with ``Quantity``, the attributes ``val`` and
    ``is_scalar`` have been renamed to ``value`` and ``isscalar``,
    respectively, and the attribute ``vals`` has been dropped.

- ``astropy.units``

  - The ``Quantity`` class now inherits from the Numpy array class, and
    includes the following API changes [#929]:

    - Using ``float(...)``, ``int(...)``, and ``long(...)`` on a quantity will
      now only work if the quantity is dimensionless and unscaled.

    - All Numpy ufuncs should now treat units correctly (or raise an exception
      if not supported), rather than extract the value of quantities and
      operate on this, emitting a warning about the implicit loss of units.

    - When using relevant Numpy ufuncs on dimensionless quantities (e.g.
      ``np.exp(h * nu / (k_B * T))``), or combining dimensionless quantities
      with Python scalars or plain Numpy arrays ``1 + v / c``, the
      dimensionless Quantity will automatically be converted to an unscaled
      dimensionless Quantity.

    - When initializing a quantity from a value with no unit, it is now set to
      be dimensionless and unscaled by default. When initializing a Quantity
      from another Quantity and with no unit specified in the initializer, the
      unit is now taken from the unit of the Quantity being initialized from.

  - The exception ``astropy.units.UnitsException`` has been renamed to
    ``astropy.units.UnitsError`` to be more consistent with the naming
    of built-in Python exceptions.

Bug Fixes
^^^^^^^^^^

- ``astropy.io.ascii``

  - The ``write()`` function was ignoring the ``fill_values`` argument. [#910]

- ``astropy.units``

  - Fixed a bug that caused the order of multiplication/division of plain
    Numpy arrays with Quantities to matter (i.e. if the plain array comes
    first the units were not preserved in the output). [#899]


0.2.5 (unreleased)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - Fixed incorrect string formatting of Angles using ``precision=0``. [#1319]

  - Fixed string formatting of Angles using ``decimal=True`` which ignored the
    ``precision`` argument. [#1323]

- ``astropy.io.ascii``

  - Fixed a crash in the IPAC table reader when the ``include/exclude_names``
    option is set. [#1348]

  - Fixed writing AASTex tables to honor the ``tabletype`` option. [#1372]

- ``astropy.io.fits``

  - Fixed a bug that could cause a segfault when trying to decompress an
    compressed HDU whose contents are truncated (due to a corrupt file, for
    example). This still causes a Python traceback but better that than a
    segfault. [#1332]

  - Newly created ``CompImageHDU`` HDUs use the correct value of the
    ``DEFAULT_COMPRESSION_TYPE`` module-level constant instead of hard-coding
    "RICE_1" in the header.

  - Fixed a corner case where when extra memory is allocated to compress an
    image, it could lead to unnecessary in-memory copying of the compressed
    image data and a possible memory leak through Numpy.

  - Fixed a bug where assigning from an mmap'd array in one FITS file over
    the old (also mmap'd) array in another FITS file failed to update the
    destination file. Corresponds to PyFITS issue 25.

  - Some miscellaneous documentation fixes.

- ``astropy.io.votable``

  - Added a warning for when a VOTable 1.2 file contains no ``RESOURCES``
    elements (at least one should be present). [#1337]

  - Fixed a test failure specific to MIPS architecture caused by an errant
    floating point warning. [#1179]

- ``astropy.nddata.convolution``

  - Prevented in-place modification of the input arrays to ``convolve()``.
    [#1153]

- ``astropy.table``

  - Added HTML escaping for string values in tables when outputting the table
    as HTML. [#1347]

  - Added a workaround in a bug in Numpy that could cause a crash when
    accessing a table row in a masked table containing ``dtype=object``
    columns. [#1229]

  - Fixed an issue similar to the one in #1229, but specific to unmasked
    tables. [#1403]

- ``astropy.units``

  - Improved error handling for unparseable units and fixed parsing CDS units
    without mantissas in the exponent. [#1288]

  - Added a physical type for spectral flux density. [#1410]

  - Normalized conversions that should result in a scale of exactly 1.0 to
    round off slight floating point imprecisions. [#1407]

  - Added support in the CDS unit parser/formatter for unusual unit prefixes
    that are nonetheless required to be supported by that convention. [#1426]

- ``astropy.wcs``

  - When passing a single array to the wcs transformation functions,
    (`astropy.wcs.Wcs.all_pix2world`, etc.), its second dimension must
    now exactly match the number of dimensions in the
    transformation. [#1395]

  - Improved error message when incorrect arguments are passed to
    ``WCS.wcs_world2pix``. [#1394]

  - Fixed a crash when trying to read WCS from FITS headers on Python 3.3
    in Windows. [#1363]

- Misc

  - Fixed crash when the ``COLUMNS`` environment variable is set to a
    non-integer value. [#1291]

  - Fixed a bug in ``ProgressBar.map`` where ``multiprocess=True`` could cause
    it to hang on waiting for the process pool to be destroyed. [#1381]

  - Importing ``astropy`` from within a source checkout will no longer exit
    the Python interpreter if the extension modules have not been built
    in-place. [#1269]

  - Fixed a crash on Python 3.2 when affiliated packages try to use the
    ``astropy.utils.data.get_pkg_data_*`` functions. [#1256]

  - Fixed a minor path normalization issue that could occur on Windows in
    ``astropy.utils.data.get_pkg_data_filename``. [#1444]

  - Miscellaneous documentation fixes and improvements [#1308, #1317, #1377,
    #1393, #1362]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Astropy installation now requests setuptools >= 0.7 during build/installation
  if neither distribute or setuptools >= 0.7 is already installed.  In other
  words, if ``import setuptools`` fails, ``ez_setup.py`` is used to bootstrap
  the latest setuptools (rather than using ``distribute_setup.py`` to bootstrap
  the now obsolete distribute package). [#1197]

- When importing Astropy from a source checkout without having built the
  extension modules first an ``ImportError`` is raised rather than a
  ``SystemExit`` exception. [#1269]


0.2.4 (2013-07-24)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - Fixed the angle parser to support parsing the string "1 degree". [#1168]

- ``astropy.cosmology``

  - Fixed a crash in the ``comoving_volume`` method on non-flat cosmologies
    when passing it an array of redshifts.

- ``astropy.io.ascii``

  - Fixed a bug that prevented saving changes to the comment symbol when
    writing changes to a table. [#1167]

- ``astropy.io.fits``

  - Added a workaround for a bug in 64-bit OSX that could cause truncation when
    writing files greater than 2^32 bytes in size. [#839]

- ``astropy.io.votable``

  - Fixed incorrect reading of tables containing multiple ``<RESOURCE>``
    elements. [#1223]

- ``astropy.table``

  - Fixed a bug where ``Table.remove_column`` and ``Table.rename_column``
    could cause a maksed table to lose its masking. [#1120]

  - Fixed bugs where subclasses of ``Table`` did not preserver their class in
    certain operations. [#1142]

  - Fixed a bug where slicing a masked table did not preserve the mask. [#1187]

- ``astropy.units``

  - Fixed a bug where the ``.si`` and ``.cgs`` properties of dimensionless
    ``Quantity`` objects raised a ``ZeroDivisionError``. [#1150]

  - Fixed a bug where multiple subsequent calls to the ``.decompose()`` method
    on array quantities applied a scale factor each time. [#1163]

- Misc

  - Fixed an installation crash that could occur sometimes on Debian/Ubuntu
    and other \*NIX systems where ``pkg_resources`` can be installed without
    installing ``setuptools``. [#1150]

  - Updated the ``distribute_setup.py`` bootstrapper to use setuptools >= 0.7
    when installing on systems that don't already have an up to date version
    of distribute/setuptools. [#1180]

  - Changed the ``version.py`` template so that Astropy affiliated packages can
    (and they should) use their own ``cython_version.py`` and
    ``utils._compiler`` modules where appropriate. This issue only pertains to
    affiliated package maintainers. [#1198]

  - Fixed a corner case where the default config file generation could crash
    if building with matplotlib but *not* Sphinx installed in a virtualenv.
    [#1225]

  - Fixed a crash that could occur in the logging module on systems that
    don't have a default preferred encoding (in particular this happened
    in some versions of PyCharm). [#1244]

  - The Astropy log now supports passing non-string objects (and calling
    ``str()`` on them by default) to the logging methods, in line with Python's
    standard logging API. [#1267]

  - Minor documentation fixes [#582, #696, #1154, #1194, #1212, #1213, #1246,
    #1252]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``astropy.cosmology``

  - Added a new ``Plank13`` object representing the Plank 2013 results. [#895]

- ``astropy.units``

  - Performance improvements in initialization of ``Quantity`` objects with
    a large number of elements. [#1231]


0.2.3 (2013-05-30)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.time``

  - Fixed inaccurate handling of leap seconds when converting from UTC to UNIX
    timestamps. [#1118]

  - Tightened required accuracy in many of the time conversion tests. [#1121]

- Misc

  - Fixed a regression that was introduced in v0.2.2 by the fix to issue #992
    that was preventing installation of Astropy affiliated packages that use
    Astropy's setup framework. [#1124]


0.2.2 (2013-05-21)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.io``

  - Fixed issues in both the ``fits`` and ``votable`` sub-packages where array
    byte order was not being handled consistently, leading to possible crashes
    espcially on big-endian systems. [#1003]

- ``astropy.io.fits``

  - When an error occurs opening a file in fitsdiff the exception message will
    now at least mention which file had the error.

  - Fixed a couple cases where creating a new table using TDIMn in some of the
    columns could cause a crash.

  - Slightly refactored how tables containing variable-length array columns are
    handled to add two improvements: Fixes an issue where accessing the data
    after a call to the `astropy.io.fits.getdata` convenience function caused
    an exception, and allows the VLA data to be read from an existing mmap of
    the FITS file.

  - Fixed a bug on Python 3 where attempting to open a non-existent file on
    Python 3 caused a seemingly unrelated traceback.

  - Fixed an issue in the tests that caused some tests to fail if Astropy is
    installed with read-only permissions.

  - Fixed a bug where instantiating a ``BinTableHDU`` from a numpy array
    containing boolean fields converted all the values to ``False``.

  - Fixed an issue where passing an array of integers into the constructor of
    ``Column()`` when the column type is floats of the same byte width caused
    the column array to become garbled.

  - Fixed inconsistent behavior in creating CONTINUE cards from byte strings
    versus unicode strings in Python 2--CONTINUE cards can now be created
    properly from unicode strings (so long as they are convertable to ASCII).

  - Fixed a bug in parsing HIERARCH keywords that do not have a space after the
    first equals sign (before the value).

  - Prevented extra leading whitespace on HIERARCH keywords from being treated
    as part of the keyword.

  - Fixed a bug where HIERARCH keywords containing lower-case letters was
    mistakenly marked as invalid during header validation along with an
    ancillary issue where the ``Header.index()`` method id not work correctly
    with HIERARCH keywords containing lower-case letters.

  - Disallowed assigning NaN and Inf floating point values as header values,
    since the FITS standard does not define a way to represent them in. Because
    this is undefined, the previous behavior did not make sense and produced
    invalid FITS files. [#954]

  - Fixed an obscure issue that can occur on systems that don't have flush to
    memory-mapped files implemented (namely GNU Hurd). [#968]

  - Improved round-tripping and preservation of manually assigned column
    attributes (``TNULLn``, ``TSCALn``, etc.) in table HDU headers. [#996]

- ``astropy.io.votable``

  - Stopped deprecation warnings from the ``astropy.io.votable`` package that
    could occur during setup. [#970]

  - Fixed an issue where INFO elements were being incorrectly dropped when
    occurring inside a TABLE element. [#1000]

  - Fixed obscure test failures on MIPS platforms. [#1010]

- ``astropy.nddata.convolution``

  - Fixed an issue in ``make_kernel()`` when using an Airy function kernel.
    Also removed the superfluous 'brickwall' option. [#939]

- ``astropy.table``

  - Fixed a crash that could occur when adding a row to an empty (rowless)
    table with masked columns. [#973]

  - Made it possible to assign to one table row from the value of another row,
    effictively making it easier to copy rows, for example. [#1019]

- ``astropy.time``

  - Added appropriate ``__copy__`` and ``__deepcopy__`` behavior; this
    omission caused a seemingly unrelated error in FK5 coordinate separation.
    [#891]

- ``astropy.units``

  - Fixed an issue where the ``isiterable()`` utility returned ``True`` for
    quantities with scalar values.  Added an ``__iter__`` method for the
    ``Quantity`` class and fixed ``isiterable()`` to catch false positives.
    [#878]

  - Fixed previously undefined behavior when multiplying a unit by a string.
    [#949]

  - Added 'time' as a physical type--this was a simple omission. [#959]

  - Fixed issues with pickling unit objects so as to play nicer with the
    multiprocessing module. [#974]

  - Made it more difficult to accidentally override existing units with a new
    unit of the same name. [#1070]

  - Added several more physical types and units that were previously omitted,
    including 'mass density', 'specific volume', 'molar volume', 'momentum',
    'angular momentum', 'angular speed', 'angular acceleration', 'electric
    current', 'electric current density', 'electric field strength', 'electric
    flux density', 'electric charge density', 'permittivity', 'electromagnetic
    field strength', 'radiant intensity', 'data quantity', 'bandwidth'; and
    'knots', 'nautical miles', 'becquerels', and 'curies' respectively. [#1072]

- Misc

  - Fixed a permission error that could occur when running ``astropy.test()``
    on Python 3 when Astropy is installed as root. [#811]

  - Made it easier to filter warnings from the ``convolve()`` function and
    from ``Quantity`` objects. [#853]

  - Fixed a crash that could occur in Python 3 when generation of the default
    config file fails during setup. [#952]

  - Fixed an unrelated error message that could occur when trying to import
    astropy from a source checkout without having build the extension modules
    first. This issue was claimed to be fixed in v0.2.1, but the fix itself had
    a bug. [#971]

  - Fixed a crash that could occur when running the ``build_sphinx`` setup
    command in Python 3. [#977]

  - Added a more helpful error message when trying to run the
    ``setup.py build_sphinx`` command when Sphinx is not installed. [#1027]

  - Minor documentation fixes and restructuring.
    [#935, #967, #978, #1004, #1028, #1047]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Some performance improvements to the ``astropy.units`` package, in particular
  improving the time it takes to import the sub-package. [#1015]


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

  - Added new method `WCS.all_world2pix` for converting from world coordinates
    to pixel space, including inversion of the astrometric distortion
    correction. [#1066, #1281]


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

  - Add a WMAP9 object using the final (9-year) WMAP parameters from 
    Hinshaw et al. 2013. It has also been made the default cosmology
    [#629, #724].

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
