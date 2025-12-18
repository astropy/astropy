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
as shown below. In this example we use a file that is installed with astropy::

    >>> from astropy.table import Table
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> chandra_events = get_pkg_data_filename('data/chandra_time.fits',
    ...                                        package='astropy.io.fits.tests')
    >>> t = Table.read(chandra_events)

A benefit of using the unified interface to read the table is that it will reconstruct
any :ref:`mixin_columns` that were written to that HDU.

If more than one table is present in the file, you can select the HDU
by index or by name::

    >>> t = Table.read(chandra_events, hdu="EVENTS")

In this case if the ``hdu`` argument is omitted, then the first table found
will be read in and a warning will be emitted.

You can also read a table from the HDUs of an in-memory FITS file. ::

    >>> from astropy.io import fits
    >>> with fits.open(chandra_events) as hdul:
    ...     t = Table.read(hdul["EVENTS"])

If a column contains unit information, it will have an associated
`astropy.units` object::

    >>> t["energy"].unit
    Unit("eV")

It is also possible to get directly a table with columns as
`~astropy.units.Quantity` objects by using the `~astropy.table.QTable` class::

    >>> from astropy.table import QTable
    >>> t2 = QTable.read(chandra_events, hdu=1)
    >>> t2['energy']
    <Quantity [7782.7305, 5926.725 ] eV>

Writing
^^^^^^^

To write a table ``t`` to a new file::

    >>> t.write('new_table.fits')

If the file already exists and you want to overwrite it, then set the
``overwrite`` keyword::

    >>> t.write('existing_table.fits', overwrite=True)

If you want to append a table to an existing file, set the ``append``
keyword::

    >>> t.write('existing_table.fits', append=True)


.. testcleanup::

    >>> import pathlib
    >>> pathlib.Path.unlink('new_table.fits')
    >>> pathlib.Path.unlink('existing_table.fits')

Alternatively, you can use the convenience function
:func:`~astropy.io.fits.table_to_hdu` to create a single
binary table HDU and insert or append that to an existing
:class:`~astropy.io.fits.HDUList`.

There is support for writing a table which contains :ref:`mixin_columns` such
as `~astropy.time.Time` or `~astropy.coordinates.SkyCoord`. This uses FITS
``COMMENT`` cards to capture additional information needed order to fully
reconstruct the mixin columns when reading back from FITS. The information is a
Python `dict` structure which is serialized using YAML.

Keywords
^^^^^^^^

The FITS keywords associated with an HDU table are represented in the ``meta``
ordered dictionary attribute of a :ref:`Table <astropy-table>`. After reading
a table you can view the available keywords in a readable format using:

  >>> for key, value in t.meta.items():
  ...     print(f'{key} = {value}')
  EXTNAME = EVENTS
  HDUNAME = EVENTS
  TLMIN2 = 0
  ...

This does not include the "internal" FITS keywords that are required to specify
the FITS table properties (e.g., ``NAXIS``, ``TTYPE1``). ``HISTORY`` and
``COMMENT`` keywords are treated specially and are returned as a list of

  >>> t.meta['MY_KEYWD'] = 'my value'
  >>> t.meta['COMMENT'] = ['First comment', 'Second comment', 'etc']
  >>> t.write('my_table.fits', overwrite=True)

.. testcleanup::
   >>> pathlib.Path.unlink('my_table.fits')

The keyword names (e.g., ``MY_KEYWD``) will be automatically capitalized prior
to writing.

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

The integer formats (decimal integer, binary, octal, hexadecimal) map to the
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

.. _unified_table_fits_masked_columns:

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
for problems in all three cases. It is possible to deactivate the masking with
``mask_invalid=False``.

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

::

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

::

  >>> t.write('data.fits', serialize_method='data_mask', overwrite=True)
  >>> Table.read('data.fits')
  <Table length=3>
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
incorporate the unit attribute of Quantity. For example::

    >>> from astropy.table import QTable
    >>> import astropy.units as u
    >>> t = QTable([[1, 2] * u.angstrom])
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

.. testcleanup::
   >>> pathlib.Path.unlink('my_table.fits')

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
avoiding any loss of precision and preserving information about the time
system if set in the fits header.

Example
~~~~~~~

..
  EXAMPLE START
  Writing and Reading Time Columns to/from FITS Tables

To read a FITS table into `~astropy.table.Table`:

    >>> from astropy.time import Time
    >>> from astropy.table import Table
    >>> from astropy.coordinates import EarthLocation
    >>> t = Table()
    >>> t['a'] = Time([100.0, 200.0], scale='tt', format='mjd',
    ...               location=EarthLocation(-2446354, 4237210, 4077985, unit='m'))
    >>> t.write('my_table.fits', overwrite=True)
    >>> tm = Table.read('my_table.fits', astropy_native=True)
    >>> tm['a']
    <Time object: scale='tt' format='jd' value=[2400100.5 2400200.5]>
    >>> tm['a'].location
    <EarthLocation (-2446354., 4237210., 4077985.) m>
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
supported by ``astropy``; a thorough examination of the FITS standard
time-related keywords is done and the time data is interpreted accordingly.

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

.. _unified_table_fits_details:

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
keywords.

Examples
~~~~~~~~

..
  EXAMPLE START
  Time Columns in FITS Files

Consider the following Time column::

    >>> t = Table()
    >>> t['a'] = Time([100.0, 200.0], scale='tt', format='mjd')

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
:ref:`table_serialization_methods` for Time columns when writing, as in the
following example::

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
