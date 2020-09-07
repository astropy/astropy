
.. currentmodule:: astropy.io.fits

Table Data
**********

In this chapter, we will discuss the data component in a table HDU. A table will
always be in an extension HDU, never in a primary HDU.

There are two kinds of tables in the FITS standard: binary tables and ASCII
tables. Binary tables are more economical in storage and faster in data access
and manipulation. ASCII tables store the data in a "human readable" form and
therefore take up more storage space as well as more processing time since the
ASCII text needs to be parsed into numerical values.

.. note::

    If you want to read or write a single table in FITS format then the
    most convenient method is often via the high-level :ref:`table_io`. In
    particular see the :ref:`Unified I/O FITS <table_io_fits>` section.

Table Data as a Record Array
============================


What is a Record Array?
-----------------------

A record array is an array which contains records (i.e., rows) of heterogeneous
data types. Record arrays are available through the records module in the NumPy
library.

Here is a sample record array::

    >>> import numpy as np
    >>> bright = np.rec.array([(1,'Sirius', -1.45, 'A1V'),
    ...                        (2,'Canopus', -0.73, 'F0Ib'),
    ...                        (3,'Rigil Kent', -0.1, 'G2V')],
    ...                       formats='int16,a20,float32,a10',
    ...                       names='order,name,mag,Sp')

In this example, there are three records (rows) and four fields (columns). The
first field is a short integer, the second a character string (of length 20),
the third a floating point number, and the fourth a character string (of length
10). Each record has the same (heterogeneous) data structure.

The underlying data structure used for FITS tables is a class called
:class:`FITS_rec` which is a specialized subclass of `numpy.recarray`. A
:class:`FITS_rec` can be instantiated directly using the same initialization
format presented for plain recarrays as in the example above. You may also
instantiate a new :class:`FITS_rec` from a list of `astropy.io.fits.Column`
objects using the :meth:`FITS_rec.from_columns` class method. This has the
exact same semantics as :meth:`BinTableHDU.from_columns` and
:meth:`TableHDU.from_columns`, except that it only returns an actual FITS_rec
array and not a whole HDU object.


Metadata of a Table
-------------------

The data in a FITS table HDU is basically a record array with added
attributes. The metadata (i.e., information about the table data) are stored in
the header. For example, the keyword TFORM1 contains the format of the first
field, TTYPE2 the name of the second field, etc. NAXIS2 gives the number of
records (rows) and TFIELDS gives the number of fields (columns). For FITS
tables, the maximum number of fields is 999. The data type specified in TFORM
is represented by letter codes for binary tables and a Fortran-like format
string for ASCII tables. Note that this is different from the format
specifications when constructing a record array.


Reading a FITS Table
--------------------

Like images, the ``.data`` attribute of a table HDU contains the data of the
table.

Example
^^^^^^^

..
  EXAMPLE START
  Reading a FITS Table with astropy.io.fits

To read a FITS Table::


    >>> from astropy.io import fits
    >>> fits_table_filename = fits.util.get_testdata_filepath('btable.fits')

    >>> hdul = fits.open(fits_table_filename)  # open a FITS file
    >>> data = hdul[1].data  # assume the first extension is a table
    >>> # show the first two rows
    >>> first_two_rows = data[:2]
    >>> first_two_rows  # doctest: +SKIP
    [(1, 'Sirius', -1.45000005, 'A1V') (2, 'Canopus', -0.73000002, 'F0Ib')]
    >>> # show the values in field "mag"
    >>> magnitudes = data['mag']
    >>> magnitudes  # doctest: +SKIP
    array([-1.45000005, -0.73000002, -0.1       ], dtype=float32)
    >>> # columns can be referenced by index too
    >>> names = data.field(1)
    >>> names.tolist() # doctest: +SKIP
    ['Sirius', 'Canopus', 'Rigil Kent']
    >>> hdul.close()

Note that in ``astropy``, when using the ``field()`` method, it is 0-indexed
while the suffixes in header keywords such as TFORM is 1-indexed. So,
``data.field(0)`` is the data in the column with the name specified in TTYPE1
and format in TFORM1.

.. warning::

    The FITS format allows table columns with a zero-width data format, such as
    ``'0D'``. This is probably intended as a space-saving measure on files in
    which that column contains no data. In such files, the zero-width columns
    are omitted when accessing the table data, so the indexes of fields might
    change when using the ``field()`` method. For this reason, if you expect
    to encounter files containing zero-width columns it is recommended to access
    fields by name rather than by index.

..
  EXAMPLE END


Table Operations
================


Selecting Records in a Table
----------------------------

Like image data, we can use the same "mask array" idea to pick out desired
records from a table and make a new table out of it.

Examples
^^^^^^^^

..
  EXAMPLE START
  Selecting Records in a Table Using a "Mask Array"

Assuming the table's second field as having the name 'magnitude', an output
table containing all the records of magnitude > -0.5 from the input table is
generated::

    >>> with fits.open(fits_table_filename) as hdul:
    ...     data = hdul[1].data
    ...     mask = data['mag'] > -0.5
    ...     newdata = data[mask]
    ...     hdu = fits.BinTableHDU(data=newdata)
    ...     hdu.writeto('newtable.fits')

It is also possible to update the data from the HDU object in-place::

    >>> with fits.open(fits_table_filename) as hdul:
    ...     hdu = hdul[1]
    ...     mask = hdu.data['mag'] > -0.5
    ...     hdu.data = hdu.data[mask]
    ...     hdu.writeto('newtable2.fits')

..
  EXAMPLE END

Merging Tables
--------------

Merging different tables is very convenient in ``astropy``.

Examples
^^^^^^^^

..
  EXAMPLE START
  Merging FITS Tables

To merge the column definitions of the input tables::

    >>> fits_other_table_filename = fits.util.get_testdata_filepath('table.fits')

    >>> with fits.open(fits_table_filename) as hdul1:
    ...     with fits.open(fits_other_table_filename) as hdul2:
    ...         new_columns = hdul1[1].columns + hdul2[1].columns
    ...         new_hdu = fits.BinTableHDU.from_columns(new_columns)
    >>> new_columns
    ColDefs(
            name = 'order'; format = 'I'
            name = 'name'; format = '20A'
            name = 'mag'; format = 'E'
            name = 'Sp'; format = '10A'
            name = 'target'; format = '20A'
            name = 'V_mag'; format = 'E'
        )

The number of fields in the output table will be the sum of numbers of fields
of the input tables. Users have to make sure the input tables do not share any
common field names. The number of records in the output table will be the
largest number of records of all input tables. The expanded slots for the
originally shorter table(s) will be zero (or blank) filled.

Another version of this example can be used to append a new column to a
table. Updating an existing table with a new column is generally more
difficult than it is worth, but you can "append" a column to a table by creating
a new table with columns from the existing table plus the new column(s)::

    >>> with fits.open(fits_table_filename) as hdul:
    ...     orig_table = hdul[1].data
    ...     orig_cols = orig_table.columns
    >>> new_cols = fits.ColDefs([
    ...     fits.Column(name='NEWCOL1', format='D',
    ...                 array=np.zeros(len(orig_table))),
    ...     fits.Column(name='NEWCOL2', format='D',
    ...                 array=np.zeros(len(orig_table)))])
    >>> hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

Now ``newtable.fits`` contains a new table with the original table, plus the
two new columns filled with zeros.

..
  EXAMPLE END

Appending Tables
----------------

Appending one table after another is slightly trickier, since the two tables
may have different field attributes.

Examples
^^^^^^^^

..
  EXAMPLE START
  Appending to FITS Tables

Here, the first example is to append by field indices, and the second one is to
append by field names. In both cases, the output table will inherit the column
attributes (name, format, etc.) of the first table::

    >>> with fits.open(fits_table_filename) as hdul1:
    ...     with fits.open(fits_table_filename) as hdul2:
    ...         nrows1 = hdul1[1].data.shape[0]
    ...         nrows2 = hdul2[1].data.shape[0]
    ...         nrows = nrows1 + nrows2
    ...         hdu = fits.BinTableHDU.from_columns(hdul1[1].columns, nrows=nrows)
    ...         for colname in hdul1[1].columns.names:
    ...             hdu.data[colname][nrows1:] = hdul2[1].data[colname]

..
  EXAMPLE END

Scaled Data in Tables
=====================

A table field's data, like an image, can also be scaled. Scaling in a table has
a more generalized meaning than in images. In images, the physical data is a
simple linear transformation from the storage data. The table fields do have
such a construct too, where BSCALE and BZERO are stored in the header as TSCALn
and TZEROn. In addition, boolean columns and ASCII tables' numeric fields are
also generalized "scaled" fields, but without TSCAL and TZERO.

All scaled fields, like the image case, will take extra memory space as well as
processing. So, if high performance is desired, try to minimize the use of
scaled fields.

All of the scalings are done for the user, so the user only sees the physical
data. Thus, there is no need to worry about scaling back and forth between the
physical and storage column values.


Creating a FITS Table
=====================

.. _column_creation:

Column Creation
---------------

To create a table from scratch, it is necessary to create individual columns
first. A :class:`Column` constructor needs the minimal information of column
name and format. Here is a summary of all allowed formats for a binary table:

.. parsed-literal::

    **FITS format code         Description                     8-bit bytes**

    L                        logical (Boolean)               1
    X                        bit                             \*
    B                        Unsigned byte                   1
    I                        16-bit integer                  2
    J                        32-bit integer                  4
    K                        64-bit integer                  8
    A                        character                       1
    E                        single precision floating point 4
    D                        double precision floating point 8
    C                        single precision complex        8
    M                        double precision complex        16
    P                        array descriptor                8
    Q                        array descriptor                16

We will concentrate on binary tables in this chapter. ASCII tables will be
discussed in a later chapter. The less frequently used X format (bit array) and
P format (used in variable length tables) will also be discussed in a later
chapter.

Besides the required name and format arguments in constructing a
:class:`Column`, there are many optional arguments which can be used in
creating a column. Here is a list of these arguments and their corresponding
header keywords and descriptions:

.. parsed-literal::

    **Argument        Corresponding         Description**
    **in Column()     header keyword**

    name            TTYPE                 column name
    format          TFORM                 column format
    unit            TUNIT                 unit
    null            TNULL                 null value (only for B, I, and J)
    bscale          TSCAL                 scaling factor for data
    bzero           TZERO                 zero point for data scaling
    disp            TDISP                 display format
    dim             TDIM                  multi-dimensional array spec
    start           TBCOL                 starting position for ASCII table
    coord_type      TCTYP                 coordinate/axis type
    coord_unit      TCUNI                 coordinate/axis unit
    coord_ref_point TCRPX                 pixel coordinate of the reference point
    coord_ref_value TCRVL                 coordinate value at reference point
    coord_inc       TCDLT                 coordinate increment at reference point
    time_ref_pos    TRPOS                 reference position for a time coordinate column
    ascii                                 specifies a column for an ASCII table
    array                                 the data of the column

Examples
^^^^^^^^

..
  EXAMPLE START
  Creating a FITS Table

Here are a few Columns using various combinations of the optional arguments::

    >>> counts = np.array([312, 334, 308, 317])
    >>> names = np.array(['NGC1', 'NGC2', 'NGC3', 'NGC4'])
    >>> values = np.arange(2*2*4).reshape(4, 2, 2)
    >>> col1 = fits.Column(name='target', format='10A', array=names)
    >>> col2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
    >>> col3 = fits.Column(name='notes', format='A10')
    >>> col4 = fits.Column(name='spectrum', format='10E')
    >>> col5 = fits.Column(name='flag', format='L', array=[True, False, True, True])
    >>> col6 = fits.Column(name='intarray', format='4I', dim='(2, 2)', array=values)

In this example, formats are specified with the FITS letter codes. When there
is a number (>1) preceding a (numeric type) letter code, it means each cell in
that field is a one-dimensional array. In the case of column "col4", each cell
is an array (a NumPy array) of 10 elements. And in the case of column "col6",
with the use of the "dim" argument, each cell is a multi-dimensional array of
2x2 elements.

For character string fields, the number should be to the *left* of the letter
'A' when creating binary tables, and should be to the *right* when creating
ASCII tables. However, as this is a common confusion, both formats are
understood when creating binary tables (note, however, that upon writing to a
file the correct format will be written in the header). So, for columns "col1"
and "col3", they both have 10 characters in each of their cells. For numeric
data type, the dimension number must be before the letter code, not after.

After the columns are constructed, the :meth:`BinTableHDU.from_columns` class
method can be used to construct a table HDU. We can either go through the
column definition object::

    >>> coldefs = fits.ColDefs([col1, col2, col3, col4, col5, col6])
    >>> hdu = fits.BinTableHDU.from_columns(coldefs)
    >>> coldefs
    ColDefs(
        name = 'target'; format = '10A'
        name = 'counts'; format = 'J'; unit = 'DN'
        name = 'notes'; format = '10A'
        name = 'spectrum'; format = '10E'
        name = 'flag'; format = 'L'
        name = 'intarray'; format = '4I'; dim = '(2, 2)'
    )

or directly use the :meth:`BinTableHDU.from_columns` method::

    >>> hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    >>> hdu.columns
    ColDefs(
        name = 'target'; format = '10A'
        name = 'counts'; format = 'J'; unit = 'DN'
        name = 'notes'; format = '10A'
        name = 'spectrum'; format = '10E'
        name = 'flag'; format = 'L'
        name = 'intarray'; format = '4I'; dim = '(2, 2)'
    )

.. note::

    Users familiar with older versions of ``astropy`` will wonder what
    happened to ``astropy.io.fits.new_table``. :meth:`BinTableHDU.from_columns`
    and its companion for ASCII tables :meth:`TableHDU.from_columns` are the
    same in the arguments they accept and their behavior, but make it
    more explicit as to what type of table HDU they create.

A look at the newly created HDU's header will show that relevant keywords are
properly populated::

    >>> hdu.header
    XTENSION= 'BINTABLE'           / binary table extension
    BITPIX  =                    8 / array data type
    NAXIS   =                    2 / number of array dimensions
    NAXIS1  =                   73 / length of dimension 1
    NAXIS2  =                    4 / length of dimension 2
    PCOUNT  =                    0 / number of group parameters
    GCOUNT  =                    1 / number of groups
    TFIELDS =                    6 / number of table fields
    TTYPE1  = 'target  '
    TFORM1  = '10A     '
    TTYPE2  = 'counts  '
    TFORM2  = 'J       '
    TUNIT2  = 'DN      '
    TTYPE3  = 'notes   '
    TFORM3  = '10A     '
    TTYPE4  = 'spectrum'
    TFORM4  = '10E     '
    TTYPE5  = 'flag    '
    TFORM5  = 'L       '
    TTYPE6  = 'intarray'
    TFORM6  = '4I      '
    TDIM6   = '(2, 2)  '

.. warning::

    It should be noted that when creating a new table with
    :meth:`BinTableHDU.from_columns`, an in-memory copy of all of the input
    column arrays is created. This is because it is not guaranteed that the
    columns are arranged contiguously in memory in row-major order (in fact,
    they are most likely not), so they have to be combined into a new array.

However, if the array data *is* already contiguous in memory, such as in an
existing record array, a kludge can be used to create a new table HDU without
any copying. First, create the Columns as before, but without using the
``array=`` argument::

    >>> col1 = fits.Column(name='target', format='10A')

Then call :meth:`BinTableHDU.from_columns`::

    >>> hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])

This will create a new table HDU as before, with the correct column
definitions, but an empty data section. Now you can assign your array directly
to the HDU's data attribute:

.. doctest-skip::

    >>> hdu.data = mydata

In a future version of ``astropy``, table creation will be simplified and this
process will not be necessary.

..
  EXAMPLE END

.. _fits_time_column:

FITS Tables with Time Columns
=============================

The `FITS Time standard paper
<https://ui.adsabs.harvard.edu/abs/2015A%26A...574A..36R/>`_ defines the formats
and keywords used to represent timing information in FITS files. The ``astropy``
FITS package provides support for reading and writing native
`~astropy.time.Time` columns and objects using this format. This is done
within the :ref:`table_io_fits` unified I/O interface and examples of usage can
be found in the :ref:`fits_astropy_native` section. The support is not
complete and only a subset of the full standard is implemented.

Example
-------

..
  EXAMPLE START
  FITS Tables with Time Columns

The following is an example of a Header extract of a binary table (event list)
with a time column:

.. parsed-literal::

    COMMENT      ---------- Globally valid key words ----------------
    TIMESYS = ’TT      ’          / Time system
    MJDREF  = 50814.000000000000  / MJD zero point for (native) TT (= 1998-01-01)
    MJD-OBS = 53516.257939301￼￼     / MJD for observation in (native) TT

    COMMENT      ---------- Time Column -----------------------
    TTYPE1  = ’Time    ’          / S/C TT corresponding to mid-exposure
    TFORM1  = ’2D      ’          / format of field
    TUNIT1  = ’s       ’
    TCTYP1  = ’TT      ’
    TCNAM1  = ’Terrestrial Time’  / This is TT
    TCUNI1  = ’s       ’

..
  EXAMPLE END

However, the FITS standard and the ``astropy`` Time object are not perfectly
mapped and some compromises must be made. To help the user understand how the
``astropy`` code deals with these situations, the following text describes the
approach that ``astropy`` takes in some detail.

To create FITS columns which adhere to the FITS Time standard, we have taken
into account the following important points stated in the `FITS Time paper
<https://ui.adsabs.harvard.edu/abs/2015A%26A...574A..36R/>`_.

The strategy used to store `~astropy.time.Time` columns in FITS tables is to
create a `~astropy.io.fits.Header` with the appropriate time coordinate
global reference keywords and the column-specific override keywords. The
module ``astropy.io.fits.fitstime`` deals with the reading and writing of
Time columns.

The following keywords set the Time Coordinate Frame:

* TIME SCALE

  The most important of all of the metadata is the time scale which is a
  specification for measuring time.

  .. parsed-literal::

      **TIMESYS** (string-valued)
      Time scale; default UTC

      **TCTYPn** (string-valued)
      Column-specific override keyword

  The global time scale may be overridden by a time scale recorded in the table
  equivalent keyword ``TCTYPn`` for time coordinates in FITS table columns.
  ``TCTYna`` is used for alternate coordinates.

* TIME REFERENCE

  The reference point in time to which all times in the HDU are relative.
  Since there are no context-specific reference times in case there are
  multiple time columns in the same table, we need to adjust the reference
  times for the columns using some other keywords.

  The reference point in time shall be specified through one of the three
  following keywords, which are listed in decreasing order of preference:

  .. parsed-literal::

      **MJDREF** (floating-valued)
      Reference time in MJD

      **JDREF** (floating-valued)
      Reference time in JD

      **DATEREF** (datetime-valued)
      Reference time in ISO-8601

  The time reference keywords (MJDREF, JDREF, DATEREF) are interpreted using the
  time scale specified in ``TIMESYS``.

  .. note::

     If none of the three keywords are present, there is no problem as long as
     all times in the HDU are expressed in ISO-8601 ``Datetime Strings`` format:
     ``CCYY-MM-DD[Thh:mm:ss[.s...]]`` (e.g., ``"2015-04-05T12:22:33.8"``);
     otherwise MJDREF = 0.0 must be assumed.

     The value of the reference time has global validity for all time values,
     but it does not have a particular time scale associated with it. Thus we
     need to use ``TCRVLn`` (time coordinate reference value) keyword to
     compensate for the time scale differences.

* TIME REFERENCE POSITION

  The reference position, specified by the keyword ``TREFPOS``, specifies the
  spatial location at which the time is valid, either where the observation was
  made or the point in space for which light-time corrections have been applied.
  This may be a standard location (such as ``GEOCENTER`` or ``TOPOCENTER``) or
  a point in space defined by specific coordinates.

  .. parsed-literal::

      **TREFPOS** (string-valued)
      Time reference position; default TOPOCENTER

      **TRPOSn** (string-valued)
      Column-specific override keyword

  .. note::

     For TOPOCENTER, we need to specify the observatory location
     (ITRS Cartesian coordinates or geodetic latitude/longitude/height) in the
     ``OBSGEO-*`` keywords.

* TIME REFERENCE DIRECTION

  If any pathlength corrections have been applied to the time stamps (i.e., if
  the reference position is not ``TOPOCENTER`` for observational data), the
  reference direction that is used in calculating the pathlength delay should
  be provided in order to maintain a proper analysis trail of the data.
  However, this is useful only if there is also information available on the
  location from where the observation was made (the observatory location).

  The reference direction is indicated through a reference to specific keywords.
  These keywords may explicitly hold the direction or indicate columns holding
  the coordinates.

  .. parsed-literal::

      **TREFDIR** (string-valued)
      Pointer to time reference direction

      **TRDIRn** (string-valued)
      Column-specific override keyword

* TIME UNIT

  The FITS standard recommends the time unit to be one of the allowed ones
  in the specification.

  .. parsed-literal::

      **TIMEUNIT** (string-valued)
      Time unit; default s

      **TCUNIn** (string-valued)
      Column-specific override

* TIME OFFSET

  It is sometimes convenient to be able to apply a uniform clock correction
  in bulk by putting that number in a single keyword. A second use
  for a time offset is to set a zero offset to a relative time series,
  allowing zero-relative times, or higher precision, in the time stamps.
  Its default value is zero.

  .. parsed-literal::

      **TIMEOFFS** (floating-valued)
      This has global validity

* The absolute, relative errors and time resolution, time binning can be used
  when needed.


The following keywords define the global time informational keywords:

* DATE and DATE-* keywords

  These define the date of HDU creation and observation in ISO-8601.
  ``DATE`` is in UTC if the file is constructed on the Earth’s surface
  and others are in the time scale given by ``TIMESYS``.

* MJD-* keywords

  These define the same as above, but in ``MJD`` (Modified Julian Date).

The implementation writes a subset of the above FITS keywords, which map
to the Time metadata. Time is intrinsically a coordinate and hence shares
keywords with the ``World Coordinate System`` specification for spatial
coordinates. Therefore, while reading FITS tables with time columns,
the verification that a coordinate column is indeed time is done using
the FITS WCS standard rules and suggestions.
