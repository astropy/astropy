.. doctest-skip-all

.. currentmodule:: astropy.io.fits

Table Data
----------

In this chapter, we'll discuss the data component in a table HDU. A table will
always be in an extension HDU, never in a primary HDU.

There are two kinds of table in the FITS standard: binary tables and ASCII
tables. Binary tables are more economical in storage and faster in data access
and manipulation. ASCII tables store the data in a "human readable" form and
therefore take up more storage space as well as more processing time since the
ASCII text needs to be parsed into numerical values.


Table Data as a Record Array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


What is a Record Array?
"""""""""""""""""""""""

A record array is an array which contains records (i.e. rows) of heterogeneous
data types. Record arrays are available through the records module in the numpy
library. Here is a simple example of record array::

    >>> from numpy import rec
    >>> bright = rec.array([(1,'Sirius', -1.45, 'A1V'),
    ...                     (2,'Canopus', -0.73, 'F0Ib'),
    ...                     (3,'Rigil Kent', -0.1, 'G2V')],
    ...                     formats='int16,a20,float32,a10',
    ...                     names='order,name,mag,Sp')

In this example, there are 3 records (rows) and 4 fields (columns). The first
field is a short integer, second a character string (of length 20), third a
floating point number, and fourth a character string (of length 10). Each
record has the same (heterogeneous) data structure.

The underlying data structure used for FITS tables is a class called
:class:`FITS_rec` which is a specialized subclass of `numpy.recarray`.  A
:class:`FITS_rec` can be instantiated directly using the same initialization
format presented for plain recarrays as in the example above.  One may also
instantiate a new :class:`FITS_rec` from a list of PyFITS `Column` objects
using the :meth:`FITS_rec.from_columns` class method.  This has the exact same
semantics as :meth:`BinTableHDU.from_columns` and
:meth:`TableHDU.from_columns`, except that it only returns an actual FITS_rec
array and not a whole HDU object.


Metadata of a Table
"""""""""""""""""""

The data in a FITS table HDU is basically a record array, with added
attributes. The metadata, i.e. information about the table data, are stored in
the header. For example, the keyword TFORM1 contains the format of the first
field, TTYPE2 the name of the second field, etc. NAXIS2 gives the number of
records(rows) and TFIELDS gives the number of fields (columns). For FITS
tables, the maximum number of fields is 999. The data type specified in TFORM
is represented by letter codes for binary tables and a FORTRAN-like format
string for ASCII tables. Note that this is different from the format
specifications when constructing a record array.


Reading a FITS Table
""""""""""""""""""""

Like images, the ``.data`` attribute of a table HDU contains the data of the
table.  To recap, the simple example in the Quick Tutorial::

    >>> f = fits.open('bright_stars.fits')  # open a FITS file
    >>> tbdata = f[1].data  # assume the first extension is a table
    >>> print tbdata[:2]  # show the first two rows
    [(1, 'Sirius', -1.4500000476837158, 'A1V'),
    (2, 'Canopus', -0.73000001907348633, 'F0Ib')]

    >>> print tbdata['mag']  # show the values in field "mag"
    [-1.45000005 -0.73000002 -0.1 ]
    >>> print tbdata.field(1)  # columns can be referenced by index too
    ['Sirius' 'Canopus' 'Rigil Kent']

Note that in Astropy, when using the ``field()`` method, it is 0-indexed while
the suffixes in header keywords, such as TFORM is 1-indexed. So,
``tbdata.field(0)`` is the data in the column with the name specified in TTYPE1
and format in TFORM1.

.. warning::

    The FITS format allows table columns with a zero-width data format, such as
    ``'0D'``.  This is probably intended as a space-saving measure on files in
    which that column contains no data.  In such files, the zero-width columns
    are ommitted when accessing the table data, so the indexes of fields might
    change when using the ``field()`` method.  For this reason, if you expect
    to encounter files containing zero-width columns it is recommended to access
    fields by name rather than by index.


Table Operations
^^^^^^^^^^^^^^^^


Selecting Records in a Table
""""""""""""""""""""""""""""

Like image data, we can use the same "mask array" idea to pick out desired
records from a table and make a new table out of it.

In the next example, assuming the table's second field having the name
'magnitude', an output table containing all the records of magnitude > 5 from
the input table is generated::

    >>> from astropy.io import fits
    >>> t = fits.open('table.fits')
    >>> tbdata = t[1].data
    >>> mask = tbdata.['magnitude'] > 5
    >>> newtbdata = tbdata[mask]
    >>> hdu = fits.BinTableHDU(data=newtbdata)
    >>> hdu.writeto('newtable.fits')


Merging Tables
""""""""""""""

Merging different tables is straightforward in Astropy. Simply merge the column
definitions of the input tables::

    >>> t1 = fits.open('table1.fits')
    >>> t2 = fits.open('table2.fits')
    >>> new_columns = t1[1].columns + t2[1].columns
    >>> hdu = fits.BinTableHDU.from_columns(new_columns)
    >>> hdu.writeto('newtable.fits')

The number of fields in the output table will be the sum of numbers of fields
of the input tables. Users have to make sure the input tables don't share any
common field names. The number of records in the output table will be the
largest number of records of all input tables. The expanded slots for the
originally shorter table(s) will be zero (or blank) filled.

A simpler version of this example can be used to append a new column to a
table.  Updating an existing table with a new column is generally more
difficult than it's worth, but one can "append" a column to a table by creating
a new table with columns from the existing table plus the new column(s)::

    >>> orig_table = fits.open('table.fits')[1].data
    >>> orig_cols = orig_table.columns
    >>> new_cols = fits.ColDefs([
    ...     fits.Column(name='NEWCOL1', format='D',
    ...                 array=np.zeros(len(orig_table))),
    ...     fits.Column(name='NEWCOL2', format='D',
    ...                 array=np.zeros(len(orig_table)))])
    >>> hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    >>> hdu.writeto('newtable.fits')

Now ``newtable.fits`` contains a new table with the original table, plus the
two new columns filled with zeros.


Appending Tables
""""""""""""""""

Appending one table after another is slightly trickier, since the two tables
may have different field attributes. Here are two examples. The first is to
append by field indices, the second one is to append by field names. In both
cases, the output table will inherit column attributes (name, format, etc.) of
the first table::

    >>> t1 = fits.open('table1.fits')
    >>> t2 = fits.open('table2.fits')
    >>> nrows1 = t1[1].data.shape[0]
    >>> nrows2 = t2[1].data.shape[0]
    >>> nrows = nrows1 + nrows2
    >>> hdu = fits.BinTableHDU.from_columns(t1[1].columns, nrows=nrows)
    >>> for colname in t1[1].columns.names:
    ...     hdu.data[colname][nrows1:] = t2[1].data[colname]
    >>> hdu.writeto('newtable.fits')


Scaled Data in Tables
^^^^^^^^^^^^^^^^^^^^^

A table field's data, like an image, can also be scaled. Scaling in a table has
a more generalized meaning than in images. In images, the physical data is a
simple linear transformation from the storage data. The table fields do have
such a construct too, where BSCALE and BZERO are stored in the header as TSCALn
and TZEROn. In addition, boolean columns and ASCII tables' numeric fields are
also generalized "scaled" fields, but without TSCAL and TZERO.

All scaled fields, like the image case, will take extra memory space as well as
processing. So, if high performance is desired, try to minimize the use of
scaled fields.

All the scalings are done for the user, so the user only sees the physical
data. Thus, this no need to worry about scaling back and forth between the
physical and storage column values.


Creating a FITS Table
^^^^^^^^^^^^^^^^^^^^^


Column Creation
"""""""""""""""

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
    K                        64-bit integer                  4
    A                        character                       1
    E                        single precision floating point 4
    D                        double precision floating point 8
    C                        single precision complex        8
    M                        double precision complex        16
    P                        array descriptor                8
    Q                        array descriptor                16

We'll concentrate on binary tables in this chapter. ASCII tables will be
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
    array                                 the data of the column


Here are a few Columns using various combination of these arguments:

    >>> import numpy as np
    >>> from fits import Column
    >>> counts = np.array([312, 334, 308, 317])
    >>> names = np.array(['NGC1', 'NGC2', 'NGC3', 'NGC4'])
    >>> c1 = Column(name='target', format='10A', array=names)
    >>> c2 = Column(name='counts', format='J', unit='DN', array=counts)
    >>> c3 = Column(name='notes', format='A10')
    >>> c4 = Column(name='spectrum', format='1000E')
    >>> c5 = Column(name='flag', format='L', array=[True, False, True, True])

In this example, formats are specified with the FITS letter codes. When there
is a number (>1) preceding a (numeric type) letter code, it means each cell in
that field is a one-dimensional array. In the case of column c4, each cell is
an array (a numpy array) of 1000 elements.

For character string fields, the number be to the *left* of the letter 'A' when
creating binary tables, and should be to the *right* when creating ASCII
tables.  However, as this is a common confusion both formats are understood
when creating binary tables (note, however, that upon writing to a file the
correct format will be written in the header).  So, for columns c1 and c3, they
both have 10 characters in each of their cells. For numeric data type, the
dimension number must be before the letter code, not after.

After the columns are constructed, the :meth:`BinTableHDU.from_columns` class
method can be used to construct a table HDU. We can either go through the
column definition object::

    >>> coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
    >>> tbhdu = fits.BinTableHDU.from_columns(coldefs)

or directly use the :meth:`BinTableHDU.from_columns` method::

    >>> tbhdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5])

.. note::

    Users familiar with older versions of PyFITS or Astropy will wonder what
    happened to :func:`~astropy.io.fits.new_table`.  It is still there, but is
    deprecated.  :meth:`BinTableHDU.from_columns` and its companion for ASCII
    tables :meth:`TableHDU.from_columns` are the same as
    :func:`~astropy.io.fits.new_table` in the arguments they accept and their
    behavior.  They just make it more explicit what type of table HDU they
    create.

A look of the newly created HDU's header will show that relevant keywords are
properly populated::

    >>> tbhdu.header
    XTENSION = 'BINTABLE'                      / binary table extension
    BITPIX   =                               8 / array data type
    NAXIS    =                               2 / number of array dimensions
    NAXIS1   =                            4025 / length of dimension 1
    NAXIS2   =                               4 / length of dimension 2
    PCOUNT   =                               0 / number of group parameters
    GCOUNT   =                               1 / number of groups
    TFIELDS  =                               5 / number of table fields
    TTYPE1   = 'target '
    TFORM1   = '10A '
    TTYPE2   = 'counts '
    TFORM2   = 'J '
    TUNIT2   = 'DN '
    TTYPE3   = 'notes '
    TFORM3   = '10A '
    TTYPE4   = 'spectrum'
    TFORM4   = '1000E '
    TTYPE5   = 'flag '
    TFORM5   = 'L '

.. warning::

    It should be noted that when creating a new table with
    :meth:`BinTableHDU.from_columns`, an in-memory copy of all of the input
    column arrays is created.  This is because it is not guaranteed that the
    columns are arranged contiguously in memory in row-major order (in fact,
    they are most likely not), so they have to be combined into a new array.

However, if the array data *is* already contiguous in memory, such as in an
existing record array, a kludge can be used to create a new table HDU without
any copying.  First, create the Columns as before, but without using the
``array=`` argument::

    >>> c1 = Column(name='target', format='10A')

Then call :meth:`BinTableHDU.from_columns`::

    >>> tbhdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5])

This will create a new table HDU as before, with the correct column
definitions, but an empty data section.  Now simply assign your array directly
to the HDU's data attribute::

    >>> tbhdu.data = mydata

In a future version of Astropy table creation will be simplified and this
process won't be necessary.
