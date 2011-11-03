.. currentmodule:: pyfits.core

**********
Table Data
**********

In this chapter, we'll discuss the data component in a table HDU. A table will
always be in an extension HDU, never in a primary HDU.

There are two kinds of table in the FITS standard: binary tables and ASCII
tables. Binary tables are more economical in storage and faster in data access
and manipulation. ASCII tables store the data in a "human readable" form and
therefore takes up more storage space as well as more processing time since the
ASCII text need to be parsed back into numerical values.


Table Data as a Record Array
============================


What is a Record Array?
-----------------------

A record array is an array which contains records (i.e. rows) of heterogeneous
data types. Record arrays are available through the records module in the numpy
library. Here is a simple example of record array:

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


Metadata of a Table
-------------------

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
--------------------

Like images, the .data attribute of a table HDU contains the data of the table.
To recap, the simple example in the Quick Tutorial:

    >>> f = pyfits.open('bright_stars.fits') # open a FITS file
    >>> tbdata = f[1].data # assume the first extension is a table
    >>> print tbdata[:2] # show the first two rows
    [(1, 'Sirius', -1.4500000476837158, 'A1V'),
    (2, 'Canopus', -0.73000001907348633, 'F0Ib')]
    
    >>> print tbdata.field('mag') # show the values in field "mag"
    [-1.45000005 -0.73000002 -0.1 ]
    >>> print tbdata.field(1) # field can be referred by index too
    ['Sirius' 'Canopus' 'Rigil Kent']
    >>> scidata[1,4] = 999 # update a pixel value
    >>> scidata[30:40, 10:20] = 0 # update values of a subsection
    >>> scidata[3] = scidata[2] # copy the 3rd row to the 4th row

Note that in PyFITS, when using the ``field()`` method, it is 0-indexed while
the suffixes in header keywords, such as TFORM is 1-indexed. So,
``tbdata.field(0)`` is the data in the column with the name specified in TTYPE1
and format in TFORM1.

**Warning:** The FITS format allows table columns with a zero-width data
format, such as '0D'.  This is probably intended as a space-saving measure on
files in which that column contains no data.  In such files, the zero-width
columns are ommitted when accessing the table data, so the indexes of fields
might change when using the ``field()`` method.  For this reason, if you expect
to encounter files containg zero-width columns it is recommended to access
fields by name rather than by index.


Table Operations
================


Selecting Records in a Table
----------------------------

Like image data, we can use the same "mask array" idea to pick out desired
records from a table and make a new table out of it.

In the next example, assuming the table's second field having the name
'magnitude', an output table containing all the records of magnitude > 5 from
the input table is generated:

    >>> import pyfits
    >>> t = pyfits.open('table.fits')
    >>> tbdata = t[1].data
    >>> mask = tbdata.field('magnitude') > 5
    >>> newtbdata = tbdata[mask]
    >>> hdu = pyfits.BinTableHDU(newtbdata)
    >>> hdu.writeto('newtable.fits')


Merging Tables
--------------

Merging different tables is straightforward in PyFITS. Simply merge the column
definitions of the input tables:

    >>> t1 = pyfits.open('table1.fits')
    >>> t2 = pyfits.open('table2.fits')
    # the column attribute is the column definitions
    >>> t = t1[1].columns + t2[1].columns
    >>> hdu = pyfits.new_table(t)
    >>> hdu.writeto('newtable.fits')

The number of fields in the output table will be the sum of numbers of fields
of the input tables. Users have to make sure the input tables don't share any
common field names. The number of records in the output table will be the
largest number of records of all input tables. The expanded slots for the
originally shorter table(s) will be zero (or blank) filled.


Appending Tables
----------------

Appending one table after another is slightly trickier, since the two tables
may have different field attributes. Here are two examples. The first is to
append by field indices, the second one is to append by field names. In both
cases, the output table will inherit column attributes (name, format, etc.) of
the first table.

    >>> t1 = pyfits.open('table1.fits')
    >>> t2 = pyfits.open('table2.fits')
    # one way to find the number of records
    >>> nrows1 = t1[1].data.shape[0]
    # another way to find the number of records
    >>> nrows2 = t2[1].header['naxis2']
    # total number of rows in the table to be generated
    >>> nrows = nrows1 + nrows2
    >>> hdu = pyfits.new_table(t1[1].columns, nrows=nrows)
    # first case, append by the order of fields
    >>> for i in range(len(t1[1].columns)):
    ... hdu.data.field(i)[nrows1:]=t2[1].data.field(i)
    # or, second case, append by the field names
    >>> for name in t1[1].columns.names:
    ... hdu.data.field(name)[nrows1:]=t2[1].data.field(name)
    # write the new table to a FITS file
    >>> hdu.writeto('newtable.fits')


Scaled Data in Tables
=====================

A table field's data, like an image, can also be scaled. Scaling in a table has
a more generalized meaning than in images. In images, the physical data is a
simple linear transformation from the storage data. The table fields do have
such construct too, where BSCALE and BZERO are stored in the header as TSCALn
and TZEROn. In addition, Boolean columns and ASCII tables' numeric fields are
also generalized "scaled" fields, but without TSCAL and TZERO.

All scaled fields, like the image case, will take extra memory space as well as
processing. So, if high performance is desired, try to minimize the use of
scaled fields.

All the scalings are done for the user, so the user only sees the physical
data. Thus, this no need to worry about scaling back and forth between the
physical and storage column values.


Creating a FITS Table
=====================


Column Creation
---------------

To create a table from scratch, it is necessary to create individual columns
first. A `Column` constructor needs the minimal information of column name and
format. Here is a summary of all allowed formats for a binary table:

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

We'll concentrate on binary tables in this chapter. ASCII tables will be
discussed in a later chapter. The less frequently used X format (bit array) and
P format (used in variable length tables) will also be discussed in a later
chapter.

Besides the required name and format arguments in constructing a `Column`,
there are many optional arguments which can be used in creating a column. Here
is a list of these arguments and their corresponding header keywords and
descriptions:

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
    >>> from pyfits import Column
    >>> counts = np.array([312, 334, 308, 317])
    >>> names = np.array(['NGC1', 'NGC2', 'NGC3', 'NGC4'])
    >>> c1 = Column(name='target', format='10A', array=names)
    >>> c2 = Column(name='counts', format='J', unit='DN', array=counts)
    >>> c3 = Column(name='notes', format='A10')
    >>> c4 = Column(name='spectrum', format='1000E')
    >>> c5 = Column(name='flag', format='L', array=[1, 0, 1, 1])

In this example, formats are specified with the FITS letter codes. When there
is a number (>1) preceding a (numeric type) letter code, it means each cell in
that field is a one-dimensional array. In the case of column c4, each cell is
an array (a numpy array) of 1000 elements.

For character string fields, the number can be either before or after the
letter 'A' and they will mean the same string size. So, for columns c1 and c3,
they both have 10 characters in each of their cells. For numeric data type, the
dimension number must be before the letter code, not after.

After the columns are constructed, the `new_table()` function can be used to
construct a table HDU. We can either go through the column definition object:

    >>> coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
    >>> tbhdu = pyfits.new_table(coldefs)

or directly use the `new_table()` function:

    >>> tbhdu = pyfits.new_table([c1, c2, c3, c4, c5])

A look of the newly created HDU's header will show that relevant keywords are
properly populated:

    >>> print tbhdu.header.ascardlist()
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

**Warning:** It should be noted that when creating a new table with
`new_table()`, an in-memory copy of all of the input column arrays is
created.  This is because it is not guaranteed that the columns are arranged
contiguously in memory in row-major order (in fact, they are most likely not),
so they have to be combined into a new array.

However, if the array data *is* already contiguous in memory, such as in an
existing record array, a kludge can be used to create a new table HDU without
any copying.  First, create the Columns as before, but without using the
``array=`` argument:

    >>> c1 = Column(name='target', format='10A')
    ...

Then call `new_table()`:

    >>> tbhdu = pyfits.new_table([c1, c2, c3, c4, c5])

This will create a new table HDU as before, with the correct column
definitions, but an empty data section.  Now simply assign your array directly
to the HDU's data attribute:

    >>> tbhdu.data = mydata

In a future version of PyFITS table creation will be simplified and this
process won't be necessary.
