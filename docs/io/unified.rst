.. _table_io:

Unified File Read/Write Interface
*********************************

``astropy`` provides a unified interface for reading and writing data in
different formats. For many common cases this will streamline the process of
file I/O and reduce the need to master the separate details of all of the I/O
packages within ``astropy``. For details on the implementation see
:ref:`io_registry`.

Getting Started with Image I/O
==============================

Reading and writing image data in the unified I/O interface is supported
though the `~astropy.nddata.CCDData` class using FITS file format:

.. doctest-skip::

    >>> # Read CCD image
    >>> ccd = CCDData.read('image.fits')

.. doctest-skip::

    >>> # Write back CCD image
    >>> ccd.write('new_image.fits')

Note that the unit is stored in the ``BUNIT`` keyword in the header on saving,
and is read from the header if it is present.

Detailed help on the available keyword arguments for reading and writing
can be obtained via the ``help()`` method as follows:

.. doctest-skip::

    >>> CCDData.read.help('fits')  # Get help on the CCDData FITS reader
    >>> CCDData.writer.help('fits')  # Get help on the CCDData FITS writer

Getting Started with Table I/O
==============================

The :class:`~astropy.table.Table` class includes two methods,
:meth:`~astropy.table.Table.read` and
:meth:`~astropy.table.Table.write`, that make it possible to read from
and write to files. A number of formats are automatically supported (see
`Built-in table readers/writers`_) and new file formats and extensions can be
registered with the :class:`~astropy.table.Table` class (see
:ref:`io_registry`).

Examples
--------

..
  EXAMPLE START
  Reading a DAOPhot Table

To use this interface, first import the :class:`~astropy.table.Table` class,
then call the :class:`~astropy.table.Table`
:meth:`~astropy.table.Table.read` method with the name of the file and
the file format, for instance ``'ascii.daophot'``:

.. doctest-skip::

    >>> from astropy.table import Table
    >>> t = Table.read('photometry.dat', format='ascii.daophot')

..
  EXAMPLE END

..
  EXAMPLE START
  Reading a Table Directly from the Internet

It is possible to load tables directly from the Internet using URLs. For
example, download tables from Vizier catalogues in CDS format
(``'ascii.cds'``)::

    >>> t = Table.read("ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/snrs.dat",
    ...         readme="ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/ReadMe",
    ...         format="ascii.cds")  # doctest: +SKIP

For certain file formats the format can be automatically detected, for
example, from the filename extension::

    >>> t = Table.read('table.tex')  # doctest: +SKIP

..
  EXAMPLE END

..
  EXAMPLE START
  Writing a LaTeX Table

For writing a table, the format can be explicitly specified::

    >>> t.write(filename, format='latex')  # doctest: +SKIP

As for the :meth:`~astropy.table.Table.read` method, the format may
be automatically identified in some cases.

The underlying file handler will also automatically detect various
compressed data formats and transparently uncompress them as far as
supported by the Python installation (see
:meth:`~astropy.utils.data.get_readable_fileobj`).

For writing, you can also specify details about the `Table serialization
methods`_ via the ``serialize_method`` keyword argument. This allows
fine control of the way to write out certain columns, for instance
writing an ISO format Time column as a pair of JD1/JD2 floating
point values (for full resolution) or as a formatted ISO date string.

..
  EXAMPLE END

Getting Help on Readers and Writers
-----------------------------------

Each file format is handled by a specific reader or writer, and each of those
functions will have its own set of arguments. For examples of
this see the section `Built-in table readers/writers`_. This section also
provides the full list of choices for the ``format`` argument.

To get help on the available arguments for each format, use the ``help()``
method of the `~astropy.table.Table.read` or `~astropy.table.Table.write`
methods. Each of these calls prints a long help document which is divided
into two sections, the generic read/write documentation (common to any
call) and the format-specific documentation. For ASCII tables, the
format-specific documentation includes the generic `astropy.io.ascii` package
interface and then a description of the particular ASCII sub-format.

In the examples below we do not show the long output:

.. doctest-skip::

    >>> Table.read.help('fits')
    >>> Table.read.help('ascii')
    >>> Table.read.help('ascii.latex')
    >>> Table.write.help('hdf5')
    >>> Table.write.help('csv')

Command-Line Utility
--------------------

For convenience, the command-line tool ``showtable`` can be used to print the
content of tables for the formats supported by the unified I/O interface.

Example
^^^^^^^

..
  EXAMPLE START
  Viewing the Contents of a Table on the Command Line

To view the contents of a table on the command line::

    $ showtable astropy/io/fits/tests/data/table.fits

     target V_mag
    ------- -----
    NGC1001  11.1
    NGC1002  12.3
    NGC1003  15.2

To get full documentation on the usage and available options, do ``showtable
--help``.

..
  EXAMPLE END

.. _built_in_readers_writers:

Built-In Table Readers/Writers
==============================

The :class:`~astropy.table.Table` class has built-in support for various input
and output formats including :ref:`table_io_ascii`,
-:ref:`table_io_fits`, :ref:`table_io_hdf5`, :ref:`table_io_pandas`,
and :ref:`table_io_votable`.

A full list of the supported formats and corresponding classes is shown in the
table below. The ``Write`` column indicates those formats that support write
functionality, and the ``Suffix`` column indicates the filename suffix
indicating a particular format. If the value of ``Suffix`` is ``auto``, the
format is auto-detected from the file itself. Not all formats support auto-
detection.

===========================  =====  ======  ============================================================================================
           Format            Write  Suffix                                          Description
===========================  =====  ======  ============================================================================================
                      ascii    Yes          ASCII table in any supported format (uses guessing)
               ascii.aastex    Yes          :class:`~astropy.io.ascii.AASTex`: AASTeX deluxetable used for AAS journals
                ascii.basic    Yes          :class:`~astropy.io.ascii.Basic`: Basic table with custom delimiters
                  ascii.cds     No          :class:`~astropy.io.ascii.Cds`: CDS format table
     ascii.commented_header    Yes          :class:`~astropy.io.ascii.CommentedHeader`: Column names in a commented line
                  ascii.csv    Yes    .csv  :class:`~astropy.io.ascii.Csv`: Basic table with comma-separated values
              ascii.daophot     No          :class:`~astropy.io.ascii.Daophot`: IRAF DAOphot format table
                 ascii.ecsv    Yes   .ecsv  :class:`~astropy.io.ascii.Ecsv`: Basic table with Enhanced CSV (supporting metadata)
          ascii.fixed_width    Yes          :class:`~astropy.io.ascii.FixedWidth`: Fixed width
ascii.fixed_width_no_header    Yes          :class:`~astropy.io.ascii.FixedWidthNoHeader`: Fixed width with no header
 ascii.fixed_width_two_line    Yes          :class:`~astropy.io.ascii.FixedWidthTwoLine`: Fixed width with second header line
                 ascii.html    Yes   .html  :class:`~astropy.io.ascii.HTML`: HTML table
                 ascii.ipac    Yes          :class:`~astropy.io.ascii.Ipac`: IPAC format table
                ascii.latex    Yes    .tex  :class:`~astropy.io.ascii.Latex`: LaTeX table
            ascii.no_header    Yes          :class:`~astropy.io.ascii.NoHeader`: Basic table with no headers
                  ascii.rdb    Yes    .rdb  :class:`~astropy.io.ascii.Rdb`: Tab-separated with a type definition header line
                  ascii.rst    Yes    .rst  :class:`~astropy.io.ascii.RST`: reStructuredText simple format table
           ascii.sextractor     No          :class:`~astropy.io.ascii.SExtractor`: SExtractor format table
                  ascii.tab    Yes          :class:`~astropy.io.ascii.Tab`: Basic table with tab-separated values
                       fits    Yes    auto  :mod:`~astropy.io.fits`: Flexible Image Transport System file
                       hdf5    Yes    auto  HDF5_: Hierarchical Data Format binary file
                 pandas.csv    Yes          Wrapper around ``pandas.read_csv()`` and ``pandas.to_csv()``
                 pandas.fwf     No          Wrapper around ``pandas.read_fwf()`` (fixed width format)
                pandas.html    Yes          Wrapper around ``pandas.read_html()`` and ``pandas.to_html()``
                pandas.json    Yes          Wrapper around ``pandas.read_json()`` and ``pandas.to_json()``
                    votable    Yes    auto  :mod:`~astropy.io.votable`: Table format used by Virtual Observatory (VO) initiative
===========================  =====  ======  ============================================================================================

.. _table_io_ascii:

ASCII Formats
-------------

The :meth:`~astropy.table.Table.read` and
:meth:`~astropy.table.Table.write` methods can be used to read and write formats
supported by `astropy.io.ascii`.

Use ``format='ascii'`` in order to interface to the generic
:func:`~astropy.io.ascii.read` and :func:`~astropy.io.ascii.write`
functions from `astropy.io.ascii`. When reading a table, this means
that all supported ASCII table formats will be tried in order to successfully
parse the input.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading and Writing ASCII Formats

To read and write formats supported by `astropy.io.ascii`:

.. doctest-skip::

  >>> t = Table.read('astropy/io/ascii/tests/t/latex1.tex', format='ascii')
  >>> print(t)
  cola colb colc
  ---- ---- ----
     a    1    2
     b    3    4

When writing a table with ``format='ascii'`` the output is a basic
character-delimited file with a single header line containing the
column names.

All additional arguments are passed to the `astropy.io.ascii`
:func:`~astropy.io.ascii.read` and :func:`~astropy.io.ascii.write`
functions. Further details are available in the sections on
:ref:`io_ascii_read_parameters` and :ref:`io_ascii_write_parameters`. For
example, to change the column delimiter and the output format for the ``colc``
column use:

.. doctest-skip::

  >>> t.write(sys.stdout, format='ascii', delimiter='|', formats={'colc': '%0.2f'})
  cola|colb|colc
  a|1|2.00
  b|3|4.00


.. note::

   When specifying an ASCII table format using the unified interface, the
   format name is prefixed with ``ascii`` in order to identify the format as
   ASCII-based. Compare the table above to the `astropy.io.ascii` list of
   :ref:`supported formats <supported_formats>` where the prefix is not
   needed. Therefore the following are equivalent:

.. doctest-skip::

     >>> dat = ascii.read('file.dat', format='daophot')
     >>> dat = Table.read('file.dat', format='ascii.daophot')

   For compatibility with ``astropy`` version 0.2 and earlier, the following
   format values are also allowed in ``Table.read()``: ``daophot``, ``ipac``,
   ``html``, ``latex``, and ``rdb``.

.. attention:: **ECSV is recommended**

   For writing and reading tables to ASCII in a way that fully reproduces the
   table data, types, and metadata (i.e., the table will "round-trip"), we
   highly recommend using the :ref:`ecsv_format`. This writes the actual data
   in a space-delimited format (the ``basic`` format) that any ASCII table
   reader can parse, but also includes metadata encoded in a comment block that
   allows full reconstruction of the original columns. This includes support
   for :ref:`ecsv_format_mixin_columns` (such as
   `~astropy.coordinates.SkyCoord` or `~astropy.time.Time`) and
   :ref:`ecsv_format_masked_columns`.

..
  EXAMPLE END

.. _table_io_fits:

FITS
----

Reading and writing tables in `FITS <https://fits.gsfc.nasa.gov/>`_ format is
supported with ``format='fits'``. In most cases, existing FITS files should be
automatically identified as such based on the header of the file, but if not,
or if writing to disk, then the format should be explicitly specified.

Reading
^^^^^^^

If a FITS table file contains only a single table, then it can be read in
with:

.. doctest-skip::

    >>> from astropy.table import Table
    >>> t = Table.read('data.fits')

If more than one table is present in the file, you can select the HDU
as follows::

    >>> t = Table.read('data.fits', hdu=3)  # doctest: +SKIP

In this case if the ``hdu`` argument is omitted, then the first table found
will be read in and a warning will be emitted::

    >>> t = Table.read('data.fits')  # doctest: +SKIP
    WARNING: hdu= was not specified but multiple tables are present, reading in first available table (hdu=1) [astropy.io.fits.connect]

You can also read a table from the HDUs of an in-memory FITS file. This will
round-trip any :ref:`mixin_columns` that were written to that HDU, using the
header information to reconstruct them::

    >>> hdulist = astropy.io.fits.open('data.fits') # doctest: +SKIP
    >>> t = Table.read(hdulist[1])  # doctest: +SKIP

Writing
^^^^^^^

To write a table ``t`` to a new file::

    >>> t.write('new_table.fits')  # doctest: +SKIP

If the file already exists and you want to overwrite it, then set the
``overwrite`` keyword::

    >>> t.write('existing_table.fits', overwrite=True)  # doctest: +SKIP

At this time there is no support for appending an HDU to an existing
file or writing multi-HDU files using the Table interface. Instead, you
can use the convenience function
:func:`~astropy.io.fits.table_to_hdu` to create a single
binary table HDU and insert or append that to an existing
:class:`~astropy.io.fits.HDUList`.

As of ``astropy`` version 3.0 there is support for writing a table which
contains :ref:`mixin_columns` such as `~astropy.time.Time` or
`~astropy.coordinates.SkyCoord`. This uses FITS ``COMMENT`` cards to capture
additional information needed order to fully reconstruct the mixin columns when
reading back from FITS. The information is a Python `dict` structure which is
serialized using YAML.

Keywords
^^^^^^^^

The FITS keywords associated with an HDU table are represented in the ``meta``
ordered dictionary attribute of a :ref:`Table <astropy-table>`. After reading
a table you can view the available keywords in a readable format using:

.. doctest-skip::

  >>> for key, value in t.meta.items():
  ...     print('{0} = {1}'.format(key, value))

This does not include the "internal" FITS keywords that are required to specify
the FITS table properties (e.g., ``NAXIS``, ``TTYPE1``). ``HISTORY`` and
``COMMENT`` keywords are treated specially and are returned as a list of
values.

Conversely, the following shows examples of setting user keyword values for a
table ``t``:

.. doctest-skip::

  >>> t.meta['MY_KEYWD'] = 'my value'
  >>> t.meta['COMMENT'] = ['First comment', 'Second comment', 'etc']
  >>> t.write('my_table.fits', overwrite=True)

The keyword names (e.g., ``MY_KEYWD``) will be automatically capitalized prior
to writing.

At this time, the ``meta`` attribute of the :class:`~astropy.table.Table` class
is an ordered dictionary and does not fully represent the structure of a
FITS header (for example, keyword comments are dropped).

.. _fits_astropy_native:


TDISPn Keyword
^^^^^^^^^^^^^^

TDISPn FITS keywords will map to and from the `~astropy.table.Column` ``format``
attribute if the display format is convertible to and from a Python display
format. Below are the rules used for both conversion directions.

TDISPn to Python format string
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TDISPn format characters are defined in the table below.

============   ================================================================
   Format                              Description
============   ================================================================
Aw             Character
Lw             Logical
Iw.m           Integer
Bw.m           Binary, integers only
Ow.m           Octal, integers only
Zw.m           Hexadecimal, integers only
Fw.d           Floating-point, fixed decimal notation
Ew.dEe         Floating-point, exponential notation
ENw.d          Engineering; E format with exponent multiple of three
ESw.d          Scientific; same as EN but non-zero leading digit if not zero
Gw.dEe         General; appears as F if significance not lost, also E
Dw.dEe         Floating-point, exponential notation, double precision
============   ================================================================

Where w is the width in characters of displayed values, m is the minimum number
of digits displayed, d is the number of digits to the right of decimal, and e
is the number of digits in the exponent. The .m and Ee fields are optional.

The A (character), L (logical), F (floating point), and G (general) display
formats can be directly translated to Python format strings. The other formats
need to be modified to match Python display formats.

For the integer formats (I, B, O, and Z), the width (w) value is used to add
space padding to the left of the column value. The minimum number (m) value is
not used. For the E, G, D, EN, and ES formats (floating point exponential) the
width (w) and precision (d) are both used, but the exponential (e) is not used.

Python format string to TDISPn
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conversion from Python format strings back to TDISPn is slightly more
complicated.

Python strings map to the TDISP format A if the Python formatting string does
not contain right space padding. It will accept left space padding. The same
applies to the logical format L.

The integer formats (decimal integer, binary, octal, hexidecimal) map to the
I, B, O, and Z TDISP formats respectively. Integer formats do not accept a
zero padded format string or a format string with no left padding defined (a
width is required in the TDISP format standard for the Integer formats).

For all float and exponential values, zero padding is not accepted. There
must be at least a width or precision defined. If only a width is defined,
there is no precision set for the TDISPn format. If only a precision is
defined, the width is set to the precision plus an extra padding value
depending on format type, and both are set in the TDISPn format. Otherwise,
if both a width and precision are present they are both set in the TDISPn
format. A Python ``f`` or ``F`` map to TDISP F format. The Python ``g`` or
``G`` map to TDISP G format. The Python ``e`` and ``E`` map to TDISP E format.

Masked Columns
^^^^^^^^^^^^^^

Tables that contain `~astropy.table.MaskedColumn` columns can be written to
FITS. By default this will replace the masked data elements with certain
sentinel values according to the FITS standard:

- ``NaN`` for float columns.
- Value of ``TNULLn`` for integer columns, as defined by the column
  ``fill_value`` attribute.
- Null string for string columns (not currently implemented).

When the file is read back those elements are marked as masked in the returned
table, but see `issue #4708 <https://github.com/astropy/astropy/issues/4708>`_
for problems in all three cases.

The FITS standard has a few limitations:

- Not all data types are supported (e.g., logical / boolean).
- Integer columns require picking one value as the NULL indicator. If
  all possible values are represented in valid data (e.g., an unsigned
  int columns with all 256 possible values in valid data), then there
  is no way to represent missing data.
- The masked data values are permanently lost, precluding the possibility
  of later unmasking the values.

``astropy`` provides a work-around for this limitation that users can choose to
use. The key part is to use the ``serialize_method='data_mask'`` keyword
argument when writing the table. This tells the FITS writer to split each masked
column into two separate columns, one for the data and one for the mask.
When it gets read back that process is reversed and the two columns are
merged back into one masked column.

.. doctest-skip::

  >>> from astropy.table.table_helpers import simple_table
  >>> t = simple_table(masked=True)
  >>> t['d'] = [False, False, True]
  >>> t['d'].mask = [True, False, False]
  >>> t
  <Table masked=True length=3>
    a      b     c     d
  int64 float64 str1  bool
  ----- ------- ---- -----
     --     1.0    c    --
      2     2.0   -- False
      3      --    e  True

.. doctest-skip::

  >>> t.write('data.fits', serialize_method='data_mask', overwrite=True)
  >>> Table.read('data.fits')
  <Table masked=True length=3>
    a      b      c      d
  int64 float64 bytes1  bool
  ----- ------- ------ -----
     --     1.0      c    --
      2     2.0     -- False
      3      --      e  True

.. warning:: This option goes outside of the established FITS standard for
   representing missing data, so users should be careful about choosing this
   option, especially if other (non-``astropy``) users will be reading the
   file(s). Behind the scenes, ``astropy`` is converting the masked columns
   into two distinct data and mask columns, then writing metadata into
   ``COMMENT`` cards to allow reconstruction of the original data.

``astropy`` Native Objects (Mixin Columns)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to store not only standard `~astropy.table.Column` objects to a
FITS table HDU, but also any ``astropy`` native objects
(:ref:`mixin_columns`) within a `~astropy.table.Table` or
`~astropy.table.QTable`. This includes `~astropy.time.Time`,
`~astropy.units.Quantity`, `~astropy.coordinates.SkyCoord`, and many others.

In general, a mixin column may contain multiple data components as well as
object attributes beyond the standard Column attributes like ``format`` or
``description``. Abiding by the rules set by the FITS standard requires the
mapping of these data components and object attributes to the appropriate FITS
table columns and keywords. Thus, a well defined protocol has been developed
to allow the storage of these mixin columns in FITS while allowing the object to
"round-trip" through the file with no loss of data or attributes.

Quantity
~~~~~~~~

A `~astropy.units.Quantity` mixin column in a `~astropy.table.QTable` is
represented in a FITS table using the ``TUNITn`` FITS column keyword to
incorporate the unit attribute of Quantity. For example:

.. doctest-skip::

    >>> from astropy.table import QTable
    >>> import astropy.units as u
    >>> t = QTable([[1, 2] * u.angstrom)])
    >>> t.write('my_table.fits', overwrite=True)
    >>> qt = QTable.read('my_table.fits')
    >>> qt
    <QTable length=2>
      col0
    Angstrom
    float64
    --------
         1.0
         2.0

Time
~~~~

``astropy`` provides the following features for reading and writing ``Time``:

- Writing and reading `~astropy.time.Time` Table columns to and from FITS
  tables.
- Reading time coordinate columns in FITS tables (compliant with the time
  standard) as `~astropy.time.Time` Table columns.

Writing and reading ``astropy`` Time columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, a `~astropy.time.Time` mixin column within a `~astropy.table.Table`
or `~astropy.table.QTable` will be written to FITS in full precision. This will
be done using the FITS time standard by setting the necessary FITS header
keywords.

The default behavior for reading a FITS table into a `~astropy.table.Table`
has historically been to convert all FITS columns to `~astropy.table.Column`
objects, which have closely matching properties. For some columns, however,
closer native ``astropy`` representations are possible, and you can indicate
these should be used by passing ``astropy_native=True`` (for backwards
compatibility, this is not done by default). This will convert columns
conforming to the FITS time standard to `~astropy.time.Time` instances,
avoiding any loss of precision.

Example
~~~~~~~

..
  EXAMPLE START
  Writing and Reading Time Columns to/from FITS Tables

To read a FITS table into `~astropy.table.Table`:

.. doctest-skip::

    >>> from astropy.time import Time
    >>> from astropy.table import Table
    >>> from astropy.coordinates import EarthLocation
    >>> t = Table()
    >>> t['a'] = Time([100.0, 200.0], scale='tt', format='mjd',
    ...               location=EarthLocation(-2446354, 4237210, 4077985, unit='m'))
    >>> t.write('my_table.fits', overwrite=True)
    >>> tm = Table.read('my_table.fits', astropy_native=True)
    >>> tm['a']
    <Time object: scale='tt' format='jd' value=[ 2400100.5  2400200.5]>
    >>> tm['a'].location
    <EarthLocation (-2446354.,  4237210.,  4077985.) m>
    >>> all(tm['a'] == t['a'])
    True

The same will work with ``QTable``.

..
  EXAMPLE END

In addition to binary table columns, various global time informational FITS
keywords are treated specially with ``astropy_native=True``. In particular,
the keywords ``DATE``, ``DATE-*`` (ISO 8601 datetime strings), and the ``MJD-*``
(MJD date values) will be returned as ``Time`` objects in the Table ``meta``.
For more details regarding the FITS time paper and the implementation,
refer to :ref:`fits_time_column`.

Since not all FITS readers are able to use the FITS time standard, it is also
possible to store `~astropy.time.Time` instances using the `_time_format`.
For this case, none of the special header keywords associated with the
FITS time standard will be set. When reading this back into ``astropy``, the
column will be an ordinary Column instead of a `~astropy.time.Time` object.
See the `Details`_ section below for an example.

Reading FITS standard compliant time coordinate columns in binary tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reading FITS files which are compliant with the FITS time standard is supported
by ``astropy`` by following the multifarious rules and conventions set by the
standard. The standard was devised in order to describe time coordinates in
an unambiguous and comprehensive manner and also to provide flexibility for its
multiple use cases. Thus, while reading time coordinate columns in FITS-
compliant files, multiple aspects of the standard are taken into consideration.

Time coordinate columns strictly compliant with the two-vector JD subset of the
standard (described in the `Details`_ section below) can be read as native
`~astropy.time.Time` objects. The other subsets of the standard are also
supported by ``astropy``; a thorough examination of the FITS standard time-
related keywords is done and the time data is interpreted accordingly.

The standard describes the various components in the specification of time:

- Time coordinate frame
- Time unit
- Corrections, errors, etc.
- Durations

The keywords used to specify times define these components. Using these
keywords, time coordinate columns are identified and read as
`~astropy.time.Time` objects. Refer to :ref:`fits_time_column` for the
specification of these keywords and their description.

There are two aspects of the standard that require special attention due to the
subtleties involved while handling them. These are:

* Column named TIME with time unit

A common convention found in existing FITS files is that a FITS binary
table column with ``TTYPEn = ‘TIME’`` represents a time coordinate column.
Many astronomical data files, including official data products from major
observatories, follow this convention that predates the FITS standard.
The FITS time standard states that such a column will be controlled by
the global time reference frame keywords, and this will still be compliant
with the present standard.

Using this convention which has been incorporated into the standard, ``astropy``
can read time coordinate columns from all such FITS tables as native
`~astropy.time.Time` objects. Common examples of FITS files following
this convention are Chandra, XMM, and HST files.

Examples
~~~~~~~~

..
  EXAMPLE START
  Reading FITS Standard Compliant Time Coordinate Columns in Binary Tables

The following is an example of a Header extract of a Chandra event list:

.. parsed-literal::

    COMMENT      ---------- Globally valid key words ----------------
    DATE    = '2016-01-27T12:34:24' / Date and time of file creation
    TIMESYS = 'TT      '           / Time system
    MJDREF  =  5.0814000000000E+04 / [d] MJD zero point for times
    TIMEUNIT= 's       '           / Time unit
    TIMEREF = 'LOCAL   '           / Time reference (barycenter/local)

    COMMENT      ---------- Time Column -----------------------
    TTYPE1  = 'time    '           / S/C TT corresponding to mid-exposure
    TFORM1  = '1D      '           / format of field
    TUNIT1  = 's       '

When reading such a FITS table with ``astropy_native=True``, ``astropy`` checks
whether the name of a column is "TIME"/ "time" (``TTYPEn = ‘TIME’``) and
whether its unit is a FITS recognized time unit (``TUNITn`` is a time unit).

For example, reading a Chandra event list which has the above mentioned header
and the time coordinate column ``time`` as ``[1, 2]`` will give::

    >>> from astropy.table import Table
    >>> from astropy.time import Time, TimeDelta
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> chandra_events = get_pkg_data_filename('data/chandra_time.fits',
    ...                                        package='astropy.io.fits.tests')
    >>> native = Table.read(chandra_events, astropy_native=True)  # doctest: +IGNORE_WARNINGS
    >>> native['time']  # doctest: +FLOAT_CMP
    <Time object: scale='tt' format='mjd' value=[57413.76033393 57413.76033393]>
    >>> non_native = Table.read(chandra_events)
    >>> # MJDREF  =  5.0814000000000E+04, TIMESYS = 'TT'
    >>> ref_time = Time(non_native.meta['MJDREF'], format='mjd',
    ...                 scale=non_native.meta['TIMESYS'].lower())
    >>> # TTYPE1  = 'time', TUNIT1 = 's'
    >>> delta_time = TimeDelta(non_native['time'])
    >>> all(ref_time + delta_time == native['time'])
    True

By default, FITS table columns will be read as standard `~astropy.table.Column`
objects without taking the FITS time standard into consideration.

..
  EXAMPLE END

* String time column in ISO 8601 Datetime format

FITS uses a subset of ISO 8601 (which in itself does not imply a particular
timescale) for several time-related keywords, such as DATE-xxx. Following the
FITS standard, its values must be written as a character string in the
following ``datetime`` format:

.. parsed-literal::

    [+/-C]CCYY-MM-DD[Thh:mm:ss[.s...]]

A time coordinate column can be constructed using this representation of time.
The following is an example of an ISO 8601 ``datetime`` format time column:

.. parsed-literal::

    TIME
    ----
    1999-01-01T00:00:00
    1999-01-01T00:00:40
    1999-01-01T00:01:06
    .
    .
    .
    1999-01-20T01:10:00

The criteria for identifying a time coordinate column in ISO 8601 format is as
follows:

A time column is identified using the time coordinate frame keywords as
described in :ref:`fits_time_column`. Once it has been identified, its datatype
is checked in order to determine its representation format. Since ISO 8601
``datetime`` format is the only string representation of time, a time
coordinate column having string datatype will be automatically read as a
`~astropy.time.Time` object with ``format='fits'`` ('fits' represents the FITS
ISO 8601 format).

As this format does not imply a particular timescale, it is determined using
the timescale keywords in the header (``TCTYP`` or ``TIMESYS``) or their
defaults. The other time coordinate information is also determined in the same
way, using the time coordinate frame keywords. All ISO 8601 times are relative
to a globally accepted zero point (year 0 corresponds to 1 BCE) and are thus
not relative to the reference time keywords (MJDREF, JDREF, or DATEREF).
Hence, these keywords will be ignored while dealing with ISO 8601 time columns.

.. note::

   Reading FITS files with time coordinate columns *may* fail. ``astropy``
   supports a large subset of these files, but there are still some FITS files
   which are not compliant with any aspect of the standard. If you have such a
   file, please do not hesitate to let us know (by opening an issue in the
   `issue tracker <https://github.com/astropy/astropy/issues>`_).

   Also, reading a column having ``TTYPEn = ‘TIME’`` as `~astropy.time.Time`
   will fail if ``TUNITn`` for the column is not a FITS-recognized time unit.

Details
~~~~~~~

Time as a dimension in astronomical data presents challenges in its
representation in FITS files. The standard has therefore been extended to
describe rigorously the time coordinate in the ``World Coordinate System``
framework. Refer to `FITS WCS paper IV
<https://ui.adsabs.harvard.edu/abs/2015A%26A...574A..36R/>`_ for details.

Allowing ``Time`` columns to be written as time coordinate
columns in FITS tables thus involves storing time values in a way that
ensures retention of precision and mapping the associated metadata to the
relevant FITS keywords.

In accordance with the standard, which states that in binary tables one may use
pairs of doubles, the ``astropy`` Time column is written in such a table as a
vector of two doubles ``(TFORMn = ‘2D’) (jd1, jd2)`` where ``JD = jd1 + jd2``.
This reproduces the time values to double-double precision and is the
"lossless" version, exploiting the higher precision provided in binary tables.
Note that ``jd1`` is always a half-integer or integer, while ``abs(jd2) < 1``.
"Round-tripping" of ``astropy``-written FITS binary tables containing time
coordinate columns has been partially achieved by mapping selected metadata,
``scale`` and singular ``location`` of `~astropy.time.Time`, to corresponding
keywords. Note that the arbitrary metadata allowed in `~astropy.table.Table`
objects within the ``meta`` dict is not written and will be lost.

Examples
~~~~~~~~

..
  EXAMPLE START
  Time Columns in FITS Files

Consider the following Time column:

    >>> t['a'] = Time([100.0, 200.0], scale='tt', format='mjd')  # doctest: +SKIP

The FITS standard requires an additional translation layer back into
the desired format. The Time column ``t['a']`` will undergo the translation
``Astropy Time --> FITS --> Astropy Time`` which corresponds to the format
conversion ``mjd --> (jd1, jd2) --> jd``. Thus, the final conversion from
``(jd1, jd2)`` will require a software implementation which is fully compliant
with the FITS time standard.

Taking this into consideration, the functionality to read/write Time
from/to FITS can be explicitly turned off, by opting to store the time
representation values in the format specified by the ``format`` attribute
of the `~astropy.time.Time` column, instead of the ``(jd1, jd2)`` format, with
no extra metadata in the header. This is the "lossy" version, but can help
with portability. For the above example, the FITS column corresponding
to ``t['a']`` will then store ``[100.0 200.0]`` instead of
``[[ 2400100.5, 0. ], [ 2400200.5, 0. ]]``. This is done by setting the
`Table serialization methods`_ for Time columns when writing, as in the
following example:

.. doctest-skip::

    >>> from astropy.time import Time
    >>> from astropy.table import Table
    >>> from astropy.coordinates import EarthLocation
    >>> t = Table()
    >>> t['a'] = Time([100.0, 200.0], scale='tt', format='mjd')
    >>> t.write('my_table.fits', overwrite=True,
    ...         serialize_method={Time: 'formatted_value'})
    >>> tm = Table.read('my_table.fits')
    >>> tm['a']
    <Column name='a' dtype='float64' length=2>
    100.0
    200.0
    >>> all(tm['a'] == t['a'].value)
    True

By default, ``serialize_method`` for Time columns is equal to
``'jd1_jd2'``, that is, Time columns will be written in full precision.

.. note::

   The ``astropy`` `~astropy.time.Time` object does not precisely map to the
   FITS time standard.

   * FORMAT

     The FITS format considers only three formats: ISO 8601, JD, and MJD.
     ``astropy`` Time allows for many other formats like ``unix`` or ``cxcsec``
     for representing the values.

     Hence, the ``format`` attribute of Time is not stored. After reading from
     FITS the user must set the ``format`` as desired.

   * LOCATION

     In the FITS standard, the reference position for a time coordinate is a
     scalar expressed via keywords. However, vectorized reference position or
     location can be supported by the `Green Bank Keyword Convention
     <https://fits.gsfc.nasa.gov/registry/greenbank.html>`_ which is a
     Registered FITS Convention. In ``astropy`` Time, location can be an array
     which is broadcastable to the Time values.

     Hence, vectorized ``location`` attribute of Time is stored and read
     following this convention.

..
  EXAMPLE END

.. doctest-skip-all

.. _table_io_hdf5:

HDF5
----

.. _HDF5: https://www.hdfgroup.org/HDF5/
.. _h5py: http://www.h5py.org/

Reading/writing from/to HDF5_ files is supported with ``format='hdf5'`` (this
requires h5py_ to be installed). However, the ``.hdf5`` file extension is
automatically recognized when writing files, and HDF5 files are automatically
identified (even with a different extension) when reading in (using the first
few bytes of the file to identify the format), so in most cases you will not
need to explicitly specify ``format='hdf5'``.

Since HDF5 files can contain multiple tables, the full path to the table
should be specified via the ``path=`` argument when reading and writing.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading from and Writing to HDF5 Files

To read a table called ``data`` from an HDF5 file named ``observations.hdf5``,
you can do::

    >>> t = Table.read('observations.hdf5', path='data')

To read a table nested in a group in the HDF5 file, you can do::

    >>> t = Table.read('observations.hdf5', path='group/data')

To write a table to a new file, the path should also be specified::

    >>> t.write('new_file.hdf5', path='updated_data')

It is also possible to write a table to an existing file using ``append=True``::

    >>> t.write('observations.hdf5', path='updated_data', append=True)

As with other formats, the ``overwrite=True`` argument is supported for
overwriting existing files. To overwrite only a single table within an HDF5
file that has multiple datasets, use *both* the ``overwrite=True`` and
``append=True`` arguments.

Finally, when writing to HDF5 files, the ``compression=`` argument can be
used to ensure that the data is compressed on disk::

    >>> t.write('new_file.hdf5', path='updated_data', compression=True)

..
  EXAMPLE END

Metadata and Mixin Columns
^^^^^^^^^^^^^^^^^^^^^^^^^^

``astropy`` tables can contain metadata, both in the table ``meta`` attribute
(which is an ordered dictionary of arbitrary key/value pairs), and within the
columns, which each have attributes ``unit``, ``format``, ``description``,
and ``meta``.

By default, when writing a table to HDF5 the code will attempt to store each
key/value pair within the table ``meta`` as HDF5 attributes of the table
dataset. This will fail if the values within ``meta`` are not objects that can
be stored as HDF5 attributes. In addition, if the table columns being stored
have defined values for any of the above-listed column attributes, these
metadata will *not* be stored and a warning will be issued.

serialize_meta
~~~~~~~~~~~~~~

To enable storing all table and column metadata to the HDF5 file, call
the ``write()`` method with ``serialize_meta=True``. This will store metadata
in a separate HDF5 dataset, contained in the same file, which is named
``<path>.__table_column_meta__``. Here ``path`` is the argument provided in
the call to ``write()``::

    >>> t.write('observations.hdf5', path='data', serialize_meta=True)

The table metadata are stored as a dataset of strings by serializing the
metadata in YAML following the `ECSV header format
<https://github.com/astropy/astropy-APEs/blob/master/APE6.rst#header-details>`_
definition. Since there are YAML parsers for most common languages, one can
easily access and use the table metadata if reading the HDF5 in a non-astropy
application.

As of ``astropy`` 3.0, by specifying ``serialize_meta=True`` one can also store
to HDF5 tables that contain :ref:`mixin_columns` such as `~astropy.time.Time` or
`~astropy.coordinates.SkyCoord` columns.

.. _table_io_pandas:

Pandas
------

.. _pandas: https://pandas.pydata.org/pandas-docs/stable/index.html

``astropy`` `~astropy.table.Table` supports the ability to read or write tables
using some of the `I/O methods <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html>`_
available within pandas_. This interface thus provides convenient wrappers to
the following functions / methods:

.. csv-table::
    :header: "Format name", "Data Description", "Reader", "Writer"
    :widths: 25, 25, 25, 25
    :delim: ;

    ``pandas.csv``;`CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`__;`read_csv() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-read-csv-table>`_;`to_csv() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-store-in-csv>`_
    ``pandas.json``;`JSON <http://www.json.org/>`__;`read_json() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-json-reader>`_;`to_json() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-json-writer>`_
    ``pandas.html``;`HTML <https://en.wikipedia.org/wiki/HTML>`__;`read_html() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-read-html>`_;`to_html() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-html>`_
    ``pandas.fwf``;Fixed Width;`read_fwf() <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_fwf.html#pandas.read_fwf>`_;

**Notes**:

- There is no fixed-width writer in pandas_.
- Reading HTML requires `BeautifulSoup4 <https://pypi.org/project/beautifulsoup4/>`_ and
  `html5lib <https://pypi.org/project/html5lib/>`_ to be installed.

When reading or writing a table, any keyword arguments apart from the
``format`` and file name are passed through to pandas, for instance:

.. doctest-skip::

  >>> t.write('data.csv', format='pandas.csv', sep=' ', header=False)
  >>> t2 = Table.read('data.csv', format='pandas.csv', sep=' ', names=['a', 'b', 'c'])

.. _table_io_jsviewer:

JSViewer
--------

Provides an interactive HTML export of a Table, like the
:class:`~astropy.io.ascii.HTML` writer but using the DataTables_ library, which
allow to visualize interactively an HTML table (with columns sorting, search,
and pagination).

Example
^^^^^^^

..
  EXAMPLE START
  JSViewer to Provide an Interactive HTML Export of a Table

To write a table ``t`` to a new file::

    >>> t.write('new_table.html', format='jsviewer')

Several additional parameters can be used:

- *table_id*: the HTML ID of the ``<table>`` tag, defaults to ``'table{id}'``
  where ``id`` is the ID of the Table object.
- *max_lines*: maximum number of lines.
- *table_class*: HTML classes added to the ``<table>`` tag, can be useful to
  customize the style of the table.
- *jskwargs*: additional arguments passed to :class:`~astropy.table.JSViewer`.
- *css*: CSS style, default to ``astropy.table.jsviewer.DEFAULT_CSS``.
- *htmldict*: additional arguments passed to :class:`~astropy.io.ascii.HTML`.

.. _Datatables: https://www.datatables.net/

..
  EXAMPLE END

.. _table_io_votable:

VO Tables
---------

Reading/writing from/to `VO table <http://www.ivoa.net/documents/VOTable/>`_
files is supported with ``format='votable'``. In most cases, existing VO
tables should be automatically identified as such based on the header of the
file, but if not, or if writing to disk, then the format should be explicitly
specified.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading from and Writing to VO Tables

If a VO table file contains only a single table, then it can be read in with::

    >>> t = Table.read('aj285677t3_votable.xml')

If more than one table is present in the file, an error will be raised,
unless the table ID is specified via the ``table_id=`` argument::

    >>> t = Table.read('catalog.xml')
    Traceback (most recent call last):
    ...
    ValueError: Multiple tables found: table id should be set via the table_id= argument. The available tables are twomass, spitzer

    >>> t = Table.read('catalog.xml', table_id='twomass')

To write to a new file, the ID of the table should also be specified (unless
``t.meta['ID']`` is defined)::

    >>> t.write('new_catalog.xml', table_id='updated_table', format='votable')

When writing, the ``compression=True`` argument can be used to force
compression of the data on disk, and the ``overwrite=True`` argument can be
used to overwrite an existing file.

..
  EXAMPLE END

.. _table_serialization_methods:

Table Serialization Methods
===========================

``astropy`` supports fine-grained control of the way to write out (serialize)
the columns in a Table. For instance, if you are writing an ISO format
Time column to an ECSV ASCII table file, you may want to write this as a pair
of JD1/JD2 floating point values for full resolution (perfect "round-trip"),
or as a formatted ISO date string so that the values are easily readable by
your other applications.

The default method for serialization depends on the format (FITS, ECSV, HDF5).
For instance HDF5 is a binary format and so it would make sense to store a Time
object as JD1/JD2, while ECSV is a flat ASCII format and commonly you
would want to see the date in the same format as the Time object. The defaults
also reflect an attempt to minimize compatibility issues between ``astropy``
versions. For instance, it was possible to write Time columns to ECSV as
formatted strings in a version prior to the ability to write as JD1/JD2
pairs, so the current default for ECSV is to write as formatted strings.

The two classes which have configurable serialization methods are
`~astropy.time.Time` and `~astropy.table.MaskedColumn`. See the sections
on Time `Details`_ and `Masked columns`_, respectively, for additional
information. The defaults for each format are listed below:

====== ==================== ===============
Format    Time                MaskedColumn
====== ==================== ===============
FITS    ``jd1_jd2``          ``null_value``
ECSV    ``formatted_value``  ``null_value``
HDF5    ``jd1_jd2``          ``data_mask``
YAML    ``jd2_jd2``            ---
====== ==================== ===============

Examples
--------

..
  EXAMPLE START
  Table Serialization Methods in astropy.io

Start by making a table with a Time column and masked column:

  >>> import sys
  >>> from astropy.time import Time
  >>> from astropy.table import Table, MaskedColumn

  >>> t = Table(masked=True)
  >>> t['tm'] = Time(['2000-01-01', '2000-01-02'])
  >>> t['mc1'] = MaskedColumn([1.0, 2.0], mask=[True, False])
  >>> t['mc2'] = MaskedColumn([3.0, 4.0], mask=[False, True])
  >>> t
  <Table masked=True length=2>
             tm             mc1     mc2
           object         float64 float64
  ----------------------- ------- -------
  2000-01-01 00:00:00.000      --     3.0
  2000-01-02 00:00:00.000     2.0      --

Now specify that you want all `~astropy.time.Time` columns written as JD1/JD2
and the ``mc1`` column written as a data/mask pair and write to ECSV:

.. doctest-skip::

  >>> serialize_method = {Time: 'jd1_jd2', 'mc1': 'data_mask'}
  >>> t.write(sys.stdout, format='ascii.ecsv', serialize_method=serialize_method)
  # %ECSV 0.9
   ...
  # schema: astropy-2.0
   tm.jd1    tm.jd2  mc1  mc1.mask  mc2
  2451544.0    0.5   1.0   True     3.0
  2451546.0   -0.5   2.0   False     ""

(Spaces added for clarity)

Notice that the ``tm`` column has been replaced by the ``tm.jd1`` and ``tm.jd2``
columns, and likewise a new column ``mc1.mask`` has appeared and it explicitly
contains the mask values. When this table is read back with the ``ascii.ecsv``
reader then the original columns are reconstructed.

The ``serialize_method`` argument can be set in two different ways:

- As a single string like ``data_mask``. This value then applies to every
  column, and is a convenient strategy for a masked table with no Time columns.
- As a `dict`, where the key can be either a single column name or a class (as
  shown in the example above), and the value is the corresponding serialization
  method.

..
  EXAMPLE END
