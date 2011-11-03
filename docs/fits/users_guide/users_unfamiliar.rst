.. currentmodule:: pyfits.core

*********************
Less Familiar Objects
*********************

In this chapter, we'll discuss less frequently used FITS data structures. They
include ASCII tables, variable length tables, and random access group FITS
files.


ASCII Tables
============

FITS standard supports both binary and ASCII tables. In ASCII tables, all the
data are stored in a human readable text form, so it takes up more space and
extra processing to parse the text for numeric data.

In PyFITS, the interface for ASCII tables and binary tables is basically the
same, i.e. the data is in the ``.data`` attribute and the ``field()`` method
is used to refer to the columns and returns a numpy array. When reading the
table, PyFITS will automatically detect what kind of table it is.

    >>> hdus = pyfits.open('ascii_table.fits')
    >>> hdus[1].data[:1]
    FITS_rec(
    ... [(10.123000144958496, 37)],
    ... dtype=[('a', '>f4'),('b','>i4')])
    >>> hdus[1].data.field('a')
    array([ 10.12300014, 5.19999981, 15.60999966, 0. ,
    345. ], dtype=float32)
    >>> hdus[1].data.formats
    ['E10.4', 'I5']

Note that the formats in the record array refer to the raw data which are ASCII
strings (therefore 'a11' and 'a5'), but the .formats attribute of data retains
the original format specifications ('E10.4' and 'I5').


Creating an ASCII Table
-----------------------

Creating an ASCII table from scratch is similar to creating a binary table. The
difference is in the Column definitions. The columns/fields in an ASCII table
are more limited than in a binary table. It does not allow more than one
numerical value in a cell. Also, it only supports a subset of what allowed in a
binary table, namely character strings, integer, and (single and double
precision) floating point numbers. Boolean and complex numbers are not allowed.

The format syntax (the values of the TFORM keywords) is different from that of a
binary table, they are:

.. parsed-literal::

    Aw         Character string
    Iw         (Decimal) Integer
    Fw.d       Single precision real
    Ew.d       Single precision real, in exponential notation
    Dw.d       Double precision real, in exponential notation

where, w is the width, and d the number of digits after the decimal point. The
syntax difference between ASCII and binary tables can be confusing. For example,
a field of 3-character string is specified '3A' in a binary table and as 'A3' in
an ASCII table.

The other difference is the need to specify the table type when using either
`ColDef()` or `new_table()`.

The default value for tbtype is `BinTableHDU`.

     >>>
     # Define the columns
     >>> import numpy as np
     >>> import pyfits
     >>> a1 = np.array(['abcd', 'def'])
     >>> r1 = np.array([11., 12.])
     >>> c1 = pyfits.Column(name='abc', format='A3', array=a1)
     >>> c2 = pyfits.Column(name='def', format='E', array=r1, bscale=2.3,
     ...                    bzero=0.6)
     >>> c3 = pyfits.Column(name='t1', format='I', array=[91, 92, 93])
     # Create the table
     >>> x = pyfits.ColDefs([c1, c2, c3], tbtype='TableHDU')
     >>> hdu = pyfits.new_table(x, tbtype='TableHDU')
     # Or, simply,
     >>> hdu = pyfits.new_table([c1, c2, c3], tbtype='TableHDU')
     >>> hdu.writeto('ascii.fits')
     >>> hdu.data
     FITS_rec([('abcd', 11.0, 91), ('def', 12.0, 92), ('', 0.0, 93)],
              dtype=[('abc', '|S3'), ('def', '|S14'), ('t1', '|S10')])


Variable Length Array Tables
============================

The FITS standard also supports variable length array tables. The basic idea is
that sometimes it is desirable to have tables with cells in the same field
(column) that have the same data type but have different lengths/dimensions.
Compared with the standard table data structure, the variable length table can
save storage space if there is a large dynamic range of data lengths in
different cells.

A variable length array table can have one or more fields (columns) which are
variable length. The rest of the fields (columns) in the same table can still be
regular, fixed-length ones. PyFITS will automatically detect what kind of field
it is during reading; no special action is needed from the user. The data type
specification (i.e. the value of the TFORM keyword) uses an extra letter 'P' and
the format is

.. parsed-literal::

    rPt(max)

where r is 0, 1, or absent, t is one of the letter code for regular table data
type (L, B, X, I, J, etc. currently, the X format is not supported for variable
length array field in PyFITS), and max is the maximum number of elements. So,
for a variable length field of int32, The corresponding format spec is,
e.g. 'PJ(100)'.

    >>> f = pyfits.open('variable_length_table.fits')
    >>> print f[1].header['tform5']
    1PI(20)
    >>> print f[1].data.field(4)[:3]
    [array([1], dtype=int16) array([88, 2], dtype=int16)
    array([ 1, 88, 3], dtype=int16)]

The above example shows a variable length array field of data type int16 and its
first row has one element, second row has 2 elements etc. Accessing variable
length fields is almost identical to regular fields, except that operations on
the whole filed are usually not possible. A user has to process the field row by
row.


Creating a Variable Length Array Table
--------------------------------------

Creating a variable length table is almost identical to creating a regular
table. The only difference is in the creation of field definitions which are
variable length arrays. First, the data type specification will need the 'P'
letter, and secondly, the field data must be an objects array (as included in
the numpy module). Here is an example of creating a table with two fields,  one
is regular and the other variable length array.

    >>> import pyfits
    >>> import numpy as np
    >>> c1 = pyfits.Column(name='var', format='PJ()',
    ...                    array=np.array([[45., 56]
                                           [11, 12, 13]],
    ...                                   dtype=np.object))
    >>> c2 = pyfits.Column(name='xyz', format='2I', array=[[11, 3], [12, 4]])
    # the rest is the same as a regular table.
    # Create the table HDU
    >>> tbhdu = pyfits.new_table([c1, c2])
    >>> print tbhdu.data
    FITS_rec([(array([45, 56]), array([11,  3], dtype=int16)),
           (array([11, 12, 13]), array([12,  4], dtype=int16))], 
          dtype=[('var', '<i4', 2), ('xyz', '<i2', 2)])
    # write to a FITS file
    >>> tbhdu.writeto('var_table.fits')
    >>> hdu = pyfits.open('var_table.fits')
    # Note that heap info is taken care of (PCOUNT) when written to FITS file.
    >>> print hdu[1].header.ascardlist()
    XTENSION= 'BINTABLE'       / binary table extension
    BITPIX  =                8 / array data type
    NAXIS   =                2 / number of array dimensions
    NAXIS1  =               12 / length of dimension 1
    NAXIS2  =                2 / length of dimension 2
    PCOUNT  =               20 / number of group parameters
    GCOUNT  =                1 / number of groups
    TFIELDS =                2 / number of table fields
    TTYPE1  = 'var '
    TFORM1  = 'PJ(3) '
    TTYPE2  = 'xyz '
    TFORM2  = '2I '


Random Access Groups
====================

Another less familiar data structure supported by the FITS standard is the
random access group. This convention was established before the binary table
extension was introduced. In most cases its use can now be superseded by the
binary table. It is mostly used in radio interferometry.

Like Primary HDUs, a Random Access Group HDU is always the first HDU of a FITS
file. Its data has one or more groups. Each group may have any number
(including 0) of parameters, together with an image. The parameters and the
image have the same data type.

All groups in the same HDU have the same data structure, i.e. same data type
(specified by the keyword BITPIX, as in image HDU), same number of parameters
(specified by PCOUNT), and the same size and shape (specified by NAXISn
keywords) of the image data. The number of groups is specified by GCOUNT and
the keyword NAXIS1 is always 0. Thus the total data size for a Random Access
Group HDU is

.. parsed-literal::

    \|BITPIX\| \* GCOUNT \* (PCOUNT + NAXIS2 \* NAXIS3 \* ... \* NAXISn)


Header and Summary
------------------

Accessing the header of a Random Access Group HDU is no different from any
other HDU. Just use the .header attribute.

The content of the HDU can similarly be summarized by using the
`HDUList.info()` method:

    >>> f = pyfits.open('random_group.fits')
    >>> print f[0].header['groups']
    True
    >>> print f[0].header['gcount']
    7956
    >>> print f[0].header['pcount']
    6
    >>> f.info()
    Filename: random_group.fits
    No. Name Type Cards Dimensions Format
    0 AN GroupsHDU 158 (3, 4, 1, 1, 1) Float32 7956 Groups
    6 Parameters


Data: Group Parameters
----------------------

The data part of a random access group HDU is, like other HDUs, in the
``.data`` attribute. It includes both parameter(s) and image array(s).

1. show the data in 100th group, including parameters and data

    >>> print f[0].data[99]
    (-8.1987486677035799e-06, 1.2010923615889215e-05,
    -1.011189139244005e-05, 258.0, 2445728., 0.10, array([[[[[ 12.4308672 ,
    0.56860745, 3.99993873],
    [ 12.74043655, 0.31398511, 3.99993873],
    [ 0. , 0. , 3.99993873],
    [ 0. , 0. , 3.99993873]]]]], dtype=float32))

The data first lists all the parameters, then the image array, for the
specified group(s). As a reminder, the image data in this file has the shape of
(1,1,1,4,3) in Python or C convention, or (3,4,1,1,1) in IRAF or FORTRAN
convention.

To access the parameters, first find out what the parameter names are, with the
.parnames attribute:

    >>> f[0].data.parnames # get the parameter names
    ['uu--', 'vv--', 'ww--', 'baseline', 'date', 'date']

The group parameter can be accessed by the ``.par()`` method. Like the table
``field()`` method, the argument can be either index or name:

    >>> print f[0].data.par(0)[99] # Access group parameter by name or by index
    -8.1987486677035799e-06
    >>> print f[0].data.par('uu--')[99]
    -8.1987486677035799e-06

Note that the parameter name 'date' appears twice. This is a feature in the
random access group, and it means to add the values together. Thus:

    >>>
    # Duplicate group parameter name 'date' for 5th and 6th parameters
    >>> print f[0].data.par(4)[99]
    2445728.0
    >>> print f[0].data.par(5)[99]
    0.10
    # When accessed by name, it adds the values together if the name is shared
    # by more than one parameter
    >>> print f[0].data.par('date')[99]
    2445728.10

The ``.par()`` is a method for either the entire data object or one data item
(a group). So there are two possible ways to get a group parameter for a
certain group, this is similar to the situation in table data (with its
``field()`` method):

    >>>
    # Access group parameter by selecting the row (group) number last
    >>> print f[0].data.par(0)[99]
    -8.1987486677035799e-06
    # Access group parameter by selecting the row (group) number first
    >>> print f[0].data[99].par(0)
    -8.1987486677035799e-06

On the other hand, to modify a group parameter, we can either assign the new
value directly (if accessing the row/group number last) or use the ``setpar()``
method (if accessing the row/group number first). The method ``setpar()`` is
also needed for updating by name if the parameter is shared by more than one
parameters:

    >>>
    # Update group parameter when selecting the row (group) number last
    >>> f[0].data.par(0)[99] = 99.
    >>>
    # Update group parameter when selecting the row (group) number first
    >>> f[0].data[99].setpar(0, 99.) # or setpar('uu--', 99.)
    >>>
    # Update group parameter by name when the name is shared by more than
    # one parameters, the new value must be a tuple of constants or sequences
    >>> f[0].data[99].setpar('date', (2445729., 0.3))
    >>> f[0].data[:3].setpar('date', (2445729., [0.11, 0.22, 0.33]))
    >>> f[0].data[:3].par('date')
    array([ 2445729.11 , 2445729.22 , 2445729.33000001])


Data: Image Data
----------------

The image array of the data portion is accessible by the ``.data`` attribute of
the data object. A numpy array is returned:

    # image part of the data
    >>> print f[0].data.data[99]
    array([[[[[ 12.4308672 , 0.56860745, 3.99993873],
    [ 12.74043655, 0.31398511, 3.99993873],
    [ 0. , 0. , 3.99993873],
    [ 0. , 0. , 3.99993873]]]]], type=float32)


Creating a Random Access Group HDU
----------------------------------

To create a random access group HDU from scratch, use `GroupData()` to
encapsulate the data into the group data structure, and use `GroupsHDU()` to
create the HDU itself:

    >>>
    # Create the image arrays. The first dimension is the number of groups.
    >>> imdata = numpy.arange(100., shape=(10, 1, 1, 2, 5))
    # Next, create the group parameter data, we'll have two parameters.
    # Note that the size of each parameter's data is also the number of groups.
    # A parameter's data can also be a numeric constant.
    >>> pdata1 = numpy.arange(10) + 0.1
    >>> pdata2 = 42
    # Create the group data object, put parameter names and parameter data
    # in lists assigned to their corresponding arguments.
    # If the data type (bitpix) is not specified, the data type of the image
    # will be used.
    >>> x = pyfits.GroupData(imdata, parnames=['abc', 'xyz'],
    ...                      pardata=[pdata1, pdata2], bitpix=-32)
    # Now, create the GroupsHDU and write to a FITS file.
    >>> hdu = pyfits.GroupsHDU(x)
    >>> hdu.writeto('test_group.fits')
    >>> print hdu.header.ascardlist()[:]
    SIMPLE =            T / conforms to FITS standard
    BITPIX =          -32 / array data type
    NAXIS  =            5 / number of array dimensions
    NAXIS1 =            0
    NAXIS2 =            5
    NAXIS3 =            2
    NAXIS4 =            1
    NAXIS5 =            1
    EXTEND =            T
    GROUPS =            T / has groups
    PCOUNT =            2 / number of parameters
    GCOUNT =           10 / number of groups
    PTYPE1 = 'abc '
    PTYPE2 = 'xyz '
    >>> print hdu.data[:2]
    FITS_rec[
    (0.10000000149011612, 42.0, array([[[[ 0., 1., 2., 3., 4.],
    [ 5., 6., 7., 8., 9.]]]], dtype=float32)),
    (1.1000000238418579, 42.0, array([[[[ 10., 11., 12., 13., 14.],
    [ 15., 16., 17., 18., 19.]]]], dtype=float32))
    ]


Compressed Image Data
=====================

A general technique has been developed for storing compressed image data in
FITS binary tables.  The principle used in this convention is to first divide
the n-dimensional image into a rectangular grid of sub images or  'tiles'.
Each tile is then compressed as a continuous block of data, and the resulting
compressed byte stream is stored in a row of a variable  length column in a
FITS binary table.  Several commonly used algorithms for compressing image
tiles are supported.  These include, Gzip, Rice,  IRAF Pixel List (PLIO), and
Hcompress.

For more details, reference "A FITS Image Compression Proposal" from:

    http://www.adass.org/adass/proceedings/adass99/P2-42/

and "Registered FITS Convention, Tiled Image Compression Convention":

    http://fits.gsfc.nasa.gov/registry/tilecompression.html

Compressed image data is accessed, in PyFITS, using the optional
"pyfits.compression" module contained in a C shared library (compression.so).
If an attempt is made to access an HDU containing compressed image data when
the pyfitsComp module is not available, the user is notified of the  problem
and the HDU is treated like a standard binary table HDU.  This  notification
will only be made the first time compressed image data is encountered.  In this
way, the pyfitsComp module is not required in order for PyFITS to work.


Header and Summary
------------------

In PyFITS, the header of a compressed image HDU appears to the user like any
image header.  The actual header stored in the FITS file is that of a  binary
table HDU with a set of special keywords, defined by the convention, to
describe the structure of the compressed image.  The conversion between binary
table HDU header and image HDU header is all performed behind the scenes.
Since the HDU is actually a binary table, it may not appear as a primary HDU in
a FITS file.

The content of the HDU header may be accessed using the ``.header`` attribute:

    >>> f = pyfits.open('compressed_image.fits')
    >>> print f[1].header
    XTENSION= 'IMAGE   '           / extension type
    BITPIX  =                   16 / array data type
    NAXIS   =                    2 / number of array dimensions
    NAXIS1  =                  512 / length of data axis
    NAXIS2  =                  512 / length of data axis
    PCOUNT  =                    0 / number of parameters
    GCOUNT  =                    1 / one data group (required keyword)
    EXTNAME = 'COMPRESSED'         / name of this binary table extension

The contents of the corresponding binary table HDU may be accessed using the
hidden ``._header`` attribute.  However, all user interface with the HDU header
should be accomplished through the image header (the ``.header`` attribute).

    >>> f = pyfits.open('compressed_image.fits')
    >>> print f[1]._header
    XTENSION= 'BINTABLE'           / binary table extension
    BITPIX  =                    8 / 8-bit bytes
    NAXIS   =                    2 / 2-dimensional binary table
    NAXIS1  =                    8 / width of table in bytes
    NAXIS2  =                  512 / number of rows in table
    PCOUNT  =               157260 / size of special data area
    GCOUNT  =                    1 / one data group (required keyword)
    TFIELDS =                    1 / number of fields in each row
    TTYPE1  = 'COMPRESSED_DATA'    / label for field   1
    TFORM1  = '1PB(384)'           / data format of field: variable length array
    ZIMAGE  =                    T / extension contains compressed image
    ZBITPIX =                   16 / data type of original image
    ZNAXIS  =                    2 / dimension of original image
    ZNAXIS1 =                  512 / length of original image axis
    ZNAXIS2 =                  512 / length of original image axis
    ZTILE1  =                  512 / size of tiles to be compressed
    ZTILE2  =                    1 / size of tiles to be compressed
    ZCMPTYPE= 'RICE_1  '           / compression algorithm
    ZNAME1  = 'BLOCKSIZE'          / compression block size
    ZVAL1   =                   32 / pixels per block
    EXTNAME = 'COMPRESSED'         / name of this binary table extension

The contents of the HDU can be summarized by using either the `info()`
convenience function or method:

    >>> pyfits.info('compressed_image.fits')
    Filename: compressed_image.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU       6  ()            int16
    1    COMPRESSED  CompImageHDU    52  (512, 512)    int16
    >>>
    >>> f = pyfits.open('compressed_image.fits')
    >>> f.info()
    Filename: compressed_image.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU       6  ()            int16
    1    COMPRESSED  CompImageHDU    52  (512, 512)    int16
    >>>


Data
----

As with the header, the data of a compressed image HDU appears to the user as
standard uncompressed image data.  The actual data is stored in the fits file
as Binary Table data containing at least one column (COMPRESSED_DATA).  Each
row of this variable-length column contains the byte stream that was generated
as a result of compressing the corresponding image tile.  Several optional
columns may also appear.  These include, UNCOMPRESSED_DATA to hold the
uncompressed pixel values for tiles that cannot be compressed, ZSCALE and ZZERO
to hold the linear scale factor and zero point offset which may be needed to
transform the raw uncompressed values back to the original image pixel values,
and ZBLANK to hold the integer value used to represent undefined pixels (if
any) in the image.

The content of the HDU data may be accessed using the ``.data`` attribute:

    >>> f = pyfits.open('compressed_image.fits')
    >>> f[1].data
    array([[38, 43, 35, ..., 45, 43, 41],
           [36, 41, 37, ..., 42, 41, 39],
           [38, 45, 37, ..., 42, 35, 43],
           ...,
           [49, 52, 49, ..., 41, 35, 39],
           [57, 52, 49, ..., 40, 41, 43],
           [53, 57, 57, ..., 39, 35, 45]], dtype=int16)


Creating a Compressed Image HDU
-------------------------------

To create a compressed image HDU from scratch, simply construct a
`CompImageHDU` object from an uncompressed image data array and its associated
image header.  From there, the HDU can be treated just like any other image
HDU.

    >>> hdu = pyfits.CompImageHDU(imageData, imageHeader)
    >>> hdu.writeto('compressed_image.fits')
    >>>

The signature for the `CompImageHDU` initializer method describes the possible
options for constructing a `CompImageHDU` object:

.. parsed-literal::

    def __init__(self, data=None, header=None, name=None,
                 compressionType='RICE_1',
                 tileSize=None,
                 hcompScale=0.,
                 hcompSmooth=0
                 quantizeLevel=16.):
    """data:            data of the image
       header:          header to be associated with the image
       name:            the EXTNAME value; if this value is None, then
                         the name from the input image header will be
                         used; if there is no name in the input image
                         header then the default name 'COMPRESSED_IMAGE'
                         is used
       compressionType: compression algorithm 'RICE_1', 'PLIO_1',
                         'GZIP_1', 'HCOMPRESS_1'
       tileSize:        compression tile sizes default treats each row
                         of image as a tile
       hcompScale:      HCOMPRESS scale parameter
       hcompSmooth:     HCOMPRESS smooth parameter
       quantizeLevel:   floating point quantization level
    """

