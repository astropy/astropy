.. currentmodule:: astropy.io.fits

Less Familiar Objects
*********************

In this chapter, we will discuss less frequently used FITS data structures. They
include ASCII tables, variable length tables, and random access group FITS
files.


ASCII Tables
============

FITS standard supports both binary and ASCII tables. In ASCII tables, all of the
data are stored in a human-readable text form, so it takes up more space and
extra processing to parse the text for numeric data. Depending on how the
columns are formatted, floating point data may also lose precision.

In ``astropy``, the interface for ASCII tables and binary tables is basically
the same (i.e., the data is in the ``.data`` attribute and the ``field()``
method is used to refer to the columns and returns a ``numpy`` array). When
reading the table, ``astropy`` will automatically detect what kind of table it
is.

::

    >>> from astropy.io import fits
    >>> filename = fits.util.get_testdata_filepath('ascii.fits')
    >>> hdul = fits.open(filename)
    >>> hdul[1].data[:1]  # doctest: +FLOAT_CMP
    FITS_rec([(10.123, 37)],
             dtype=(numpy.record, {'names':['a','b'], 'formats':['S10','S5'], 'offsets':[0,11], 'itemsize':16}))
    >>> hdul[1].data['a']
    array([  10.123,    5.2  ,   15.61 ,    0.   ,  345.   ])
    >>> hdul[1].data.formats
    ['E10.4', 'I5']
    >>> hdul.close()

Note that the formats in the record array refer to the raw data which are ASCII
strings (therefore 'a11' and 'a5'), but the ``.formats`` attribute of data
retains the original format specifications ('E10.4' and 'I5').

.. _creating_ascii_table:

Creating an ASCII Table
-----------------------

Creating an ASCII table from scratch is similar to creating a binary table. The
difference is in the Column definitions. The columns/fields in an ASCII table
are more limited than in a binary table. It does not allow more than one
numerical value in a cell. Also, it only supports a subset of what is allowed
in a binary table, namely character strings, integer, and (single and double
precision) floating point numbers. Boolean and complex numbers are not allowed.

The format syntax (the values of the TFORM keywords) is different from that of a
binary table. They are:

.. parsed-literal::

    Aw         Character string
    Iw         (Decimal) Integer
    Fw.d       Double precision real
    Ew.d       Double precision real, in exponential notation
    Dw.d       Double precision real, in exponential notation

where w is the width, and d the number of digits after the decimal point. The
syntax difference between ASCII and binary tables can be confusing. For example,
a field of 3-character string is specified as '3A' in a binary table and as
'A3' in an ASCII table.

The other difference is the need to specify the table type when using the
:meth:`TableHDU.from_columns` method, and that `Column` should be provided the
``ascii=True`` argument in order to be unambiguous.

.. note::

    Although binary tables are more common in most FITS files, earlier versions
    of the FITS format only supported ASCII tables. That is why the class
    :class:`TableHDU` is used for representing ASCII tables specifically,
    whereas :class:`BinTableHDU` is more explicit that it represents a binary
    table. These names come from the value ``XTENSION`` keyword in the tables'
    headers, which is ``TABLE`` for ASCII tables and ``BINTABLE`` for binary
    tables.

:meth:`TableHDU.from_columns` can be used like so::

    >>> import numpy as np

    >>> a1 = np.array(['abcd', 'def'])
    >>> r1 = np.array([11., 12.])
    >>> col1 = fits.Column(name='abc', format='A3', array=a1, ascii=True)
    >>> col2 = fits.Column(name='def', format='E', array=r1, bscale=2.3,
    ...                    bzero=0.6, ascii=True)
    >>> col3 = fits.Column(name='t1', format='I', array=[91, 92, 93], ascii=True)
    >>> hdu = fits.TableHDU.from_columns([col1, col2, col3])
    >>> hdu.data
    FITS_rec([('abc', 11.0, 91), ('def', 12.0, 92), ('', 0.0, 93)],
             dtype=(numpy.record, [('abc', 'S3'), ('def', 'S15'), ('t1', 'S10')]))

It should be noted that when the formats of the columns are unambiguously
specific to ASCII tables it is not necessary to specify ``ascii=True`` in
the :class:`ColDefs` constructor. In this case there *is* ambiguity because
the format code ``'I'`` represents a 16-bit integer in binary tables, while in
ASCII tables it is not technically a valid format. ASCII table format codes
technically require a character width for each column, such as ``'I10'`` to
create a column that can hold integers up to 10 characters wide.

However, ``astropy`` allows the width specification to be omitted in some cases.
When it is omitted from ``'I'`` format columns the minimum width needed to
accurately represent all integers in the column is used. The only problem with
using this shortcut is its ambiguity with the binary table ``'I'`` format, so
specifying ``ascii=True`` is a good practice (though ``astropy`` will still
figure out what you meant in most cases).


Variable Length Array Tables
============================

The FITS standard also supports variable length array tables. The basic idea is
that sometimes it is desirable to have tables with cells in the same field
(column) that have the same data type but have different lengths/dimensions.
Compared with the standard table data structure, the variable length table can
save storage space if there is a large dynamic range of data lengths in
different cells.

A variable length array table can have one or more fields (columns) which are
variable length. The rest of the fields (columns) in the same table can still
be regular, fixed-length ones. ``astropy`` will automatically detect what kind
of field it is during reading; no special action is needed from the user. The
data type specification (i.e., the value of the TFORM keyword) uses an extra
letter 'P' and the format is:

.. parsed-literal::

    rPt(max)

where ``r`` may be 0 or 1 (typically omitted, as it is not applicable to
variable length arrays), ``t`` is one of the letter codes for basic data types
(L, B, I, J, etc.; currently, the X format is not supported for variable length
array field in ``astropy``), and ``max`` is the maximum number of elements of
any array in the column. So, for a variable length field of int16, the
corresponding format spec is, for example, 'PJ(100)'.

Example
-------

..
  EXAMPLE START
  Accessing Variable Length Array Columns in FITS Tables

This example shows a variable length array field of data type int16::

    >>> filename = fits.util.get_testdata_filepath('variable_length_table.fits')
    >>> hdul = fits.open(filename)
    >>> hdul[1].header['tform1']
    'PI(3)'
    >>> print(hdul[1].data.field(0))
    [array([45, 56], dtype=int16) array([11, 12, 13], dtype=int16)]
    >>> hdul.close()

In this field the first row has one element, the second row has two elements,
etc. Accessing variable length fields is almost identical to regular fields,
except that operations on the whole field simultaneously are usually not
possible. A user has to process the field row by row as though they are
independent arrays.

..
  EXAMPLE END


Creating a Variable Length Array Table
--------------------------------------

Creating a variable length table is almost identical to creating a regular
table. The only difference is in the creation of field definitions which are
variable length arrays. First, the data type specification will need the 'P'
letter, and secondly, the field data must be an objects array (as included in
the ``numpy`` module).

Example
^^^^^^^

..
  EXAMPLE START
  Creating a Variable Length Array Column in a FITS Table

Here is an example of creating a table with two fields; one is regular and the
other a variable length array::

    >>> col1 = fits.Column(
    ...    name='var', format='PI()',
    ...    array=np.array([[45, 56], [11, 12, 13]], dtype=np.object_))
    >>> col2 = fits.Column(name='xyz', format='2I', array=[[11, 3], [12, 4]])
    >>> hdu = fits.BinTableHDU.from_columns([col1, col2])
    >>> data = hdu.data
    >>> data  # doctest: +SKIP
    FITS_rec([([45, 56], [11,  3]), ([11, 12, 13], [12,  4])],
             dtype=(numpy.record, [('var', '<i4', (2,)), ('xyz', '<i2', (2,))]))
    >>> hdu.writeto('variable_length_table.fits')
    >>> with fits.open('variable_length_table.fits') as hdul:
    ...     print(repr(hdul[1].header))
    XTENSION= 'BINTABLE'           / binary table extension
    BITPIX  =                    8 / array data type
    NAXIS   =                    2 / number of array dimensions
    NAXIS1  =                   12 / length of dimension 1
    NAXIS2  =                    2 / length of dimension 2
    PCOUNT  =                   10 / number of group parameters
    GCOUNT  =                    1 / number of groups
    TFIELDS =                    2 / number of table fields
    TTYPE1  = 'var     '
    TFORM1  = 'PI(3)   '
    TTYPE2  = 'xyz     '
    TFORM2  = '2I      '

..
  EXAMPLE END

.. _random-groups:

Random Access Groups
====================

Another less familiar data structure supported by the FITS standard is the
random access group. This convention was established before the binary table
extension was introduced. In most cases its use can now be superseded by the
binary table. It is mostly used in radio interferometry.

Like primary HDUs, a Random Access Group HDU is always the first HDU of a FITS
file. Its data has one or more groups. Each group may have any number
(including 0) of parameters, together with an image. The parameters and the
image have the same data type.

All groups in the same HDU have the same data structure, that is, same data type
(specified by the keyword BITPIX, as in image HDU), same number of parameters
(specified by PCOUNT), and the same size and shape (specified by NAXISn
keywords) of the image data. The number of groups is specified by GCOUNT and
the keyword NAXIS1 is always 0. Thus the total data size for a Random Access
Group HDU is:

.. parsed-literal::

    \|BITPIX\| \* GCOUNT \* (PCOUNT + NAXIS2 \* NAXIS3 \* ... \* NAXISn)


Header and Summary
------------------

Accessing the header of a Random Access Group HDU is no different from any
other HDU; you can use the .header attribute.

The content of the HDU can similarly be summarized by using the
:meth:`HDUList.info` method::

    >>> filename = fits.util.get_testdata_filepath('group.fits')
    >>> hdul = fits.open(filename)
    >>> hdul[0].header['groups']
    True
    >>> hdul[0].header['gcount']
    10
    >>> hdul[0].header['pcount']
    3
    >>> hdul.info()
    Filename: ...group.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 GroupsHDU       15   (5, 3, 1, 1)   float32   10 Groups  3 Parameters


Data: Group Parameters
----------------------

The data part of a Random Access Group HDU is, like other HDUs, in the
``.data`` attribute. It includes both parameter(s) and image array(s).

Examples
^^^^^^^^

..
  EXAMPLE START
  Group Parameters in Random Access Group HDUs

To show the contents of the third group, including parameters and data::

    >>> hdul[0].data[2]  # doctest: +FLOAT_CMP
    (2.0999999, 42.0, 42.0, array([[[[30., 31., 32., 33., 34.],
             [35., 36., 37., 38., 39.],
             [40., 41., 42., 43., 44.]]]], dtype=float32))

The data first lists all of the parameters, then the image array, for the
specified group(s). As a reminder, the image data in this file has the shape of
(1,1,1,4,3) in Python or C convention, or (3,4,1,1,1) in IRAF or Fortran
convention.

To access the parameters, first find out what the parameter names are, with the
``.parnames`` attribute::

    >>> hdul[0].data.parnames # get the parameter names
    ['abc', 'xyz', 'xyz']

The group parameter can be accessed by the :meth:`~GroupData.par` method. Like
the table :meth:`~FITS_rec.field` method, the argument can be either index or
name::

    >>> hdul[0].data.par(0)[8]  # Access group parameter by name or by index  # doctest: +FLOAT_CMP
    8.1
    >>> hdul[0].data.par('abc')[8]  # doctest: +FLOAT_CMP
    8.1

Note that the parameter name 'xyz' appears twice. This is a feature in the
random access group, and it means to add the values together. Thus::

    >>> hdul[0].data.parnames  # get the parameter names
    ['abc', 'xyz', 'xyz']
    >>> hdul[0].data.par(1)[8]  # Duplicate parameter name 'xyz'
    42.0
    >>> hdul[0].data.par(2)[8]
    42.0
    >>> # When accessed by name, it adds the values together if the name is
    >>> # shared by more than one parameter
    >>> hdul[0].data.par('xyz')[8]
    84.0

The :meth:`~GroupData.par` is a method for either the entire data object or one
data item (a group). So there are two possible ways to get a group parameter
for a certain group, this is similar to the situation in table data (with its
:meth:`~FITS_rec.field` method)::

    >>> hdul[0].data.par(0)[8]  # doctest: +FLOAT_CMP
    8.1
    >>> hdul[0].data[8].par(0)  # doctest: +FLOAT_CMP
    8.1

On the other hand, to modify a group parameter, we can either assign the new
value directly (if accessing the row/group number last) or use the
:meth:`~Group.setpar` method (if accessing the row/group number first). The
method :meth:`~Group.setpar` is also needed for updating by name if the
parameter is shared by more than one parameters::

    >>> # Update group parameter when selecting the row (group) number last
    >>> hdul[0].data.par(0)[8] = 99.
    >>> # Update group parameter when selecting the row (group) number first
    >>> hdul[0].data[8].setpar(0, 99.)  # or:
    >>> hdul[0].data[8].setpar('abc', 99.)
    >>> # Update group parameter by name when the name is shared by more than
    >>> # one parameters, the new value must be a tuple of constants or
    >>> # sequences
    >>> hdul[0].data[8].setpar('xyz', (2445729., 0.3))
    >>> hdul[0].data[8:].par('xyz')  # doctest: +FLOAT_CMP
    array([2.44572930e+06, 8.40000000e+01])

..
  EXAMPLE END

Data: Image Data
----------------

The image array of the data portion is accessible by the
:attr:`~GroupData.data` attribute of the data object. A ``numpy`` array is
returned::

    >>> print(hdul[0].data.data[8])  # doctest: +FLOAT_CMP
    [[[[120. 121. 122. 123. 124.]
       [125. 126. 127. 128. 129.]
       [130. 131. 132. 133. 134.]]]]
    >>> hdul.close()


Creating a Random Access Group HDU
----------------------------------

To create a Random Access Group HDU from scratch, use :class:`GroupData` to
encapsulate the data into the group data structure, and use :class:`GroupsHDU`
to create the HDU itself.

Example
^^^^^^^

..
  EXAMPLE START
  Creating a Random Access Group HDU in a FITS File

To create a Random Access Group HDU::

    >>> # Create the image arrays. The first dimension is the number of groups.
    >>> imdata = np.arange(150.0).reshape(10, 1, 1, 3, 5)
    >>> # Next, create the group parameter data, we'll have two parameters.
    >>> # Note that the size of each parameter's data is also the number of
    >>> # groups.
    >>> # A parameter's data can also be a numeric constant.
    >>> pdata1 = np.arange(10) + 0.1
    >>> pdata2 = 42
    >>> # Create the group data object, put parameter names and parameter data
    >>> # in lists assigned to their corresponding arguments.
    >>> # If the data type (bitpix) is not specified, the data type of the
    >>> # image will be used.
    >>> x = fits.GroupData(imdata, bitpix=-32,
    ...                    parnames=['abc', 'xyz', 'xyz'],
    ...                    pardata=[pdata1, pdata2, pdata2])
    >>> # Now, create the GroupsHDU and write to a FITS file.
    >>> hdu = fits.GroupsHDU(x)
    >>> hdu.writeto('test_group.fits')
    >>> hdu.header
    SIMPLE  =                    T / conforms to FITS standard
    BITPIX  =                  -32 / array data type
    NAXIS   =                    5 / number of array dimensions
    NAXIS1  =                    0
    NAXIS2  =                    5
    NAXIS3  =                    3
    NAXIS4  =                    1
    NAXIS5  =                    1
    EXTEND  =                    T
    GROUPS  =                    T / has groups
    PCOUNT  =                    3 / number of parameters
    GCOUNT  =                   10 / number of groups
    PTYPE1  = 'abc     '
    PTYPE2  = 'xyz     '
    PTYPE3  = 'xyz     '
    >>> data = hdu.data
    >>> hdu.data  # doctest: +SKIP
    GroupData([ (0.1       , 42., 42., [[[[  0.,   1.,   2.,   3.,   4.], [  5.,   6.,   7.,   8.,   9.], [ 10.,  11.,  12.,  13.,  14.]]]]),
               (1.10000002, 42., 42., [[[[ 15.,  16.,  17.,  18.,  19.], [ 20.,  21.,  22.,  23.,  24.], [ 25.,  26.,  27.,  28.,  29.]]]]),
               (2.0999999 , 42., 42., [[[[ 30.,  31.,  32.,  33.,  34.], [ 35.,  36.,  37.,  38.,  39.], [ 40.,  41.,  42.,  43.,  44.]]]]),
               (3.0999999 , 42., 42., [[[[ 45.,  46.,  47.,  48.,  49.], [ 50.,  51.,  52.,  53.,  54.], [ 55.,  56.,  57.,  58.,  59.]]]]),
               (4.0999999 , 42., 42., [[[[ 60.,  61.,  62.,  63.,  64.], [ 65.,  66.,  67.,  68.,  69.], [ 70.,  71.,  72.,  73.,  74.]]]]),
               (5.0999999 , 42., 42., [[[[ 75.,  76.,  77.,  78.,  79.], [ 80.,  81.,  82.,  83.,  84.], [ 85.,  86.,  87.,  88.,  89.]]]]),
               (6.0999999 , 42., 42., [[[[ 90.,  91.,  92.,  93.,  94.], [ 95.,  96.,  97.,  98.,  99.], [100., 101., 102., 103., 104.]]]]),
               (7.0999999 , 42., 42., [[[[105., 106., 107., 108., 109.], [110., 111., 112., 113., 114.], [115., 116., 117., 118., 119.]]]]),
               (8.10000038, 42., 42., [[[[120., 121., 122., 123., 124.], [125., 126., 127., 128., 129.], [130., 131., 132., 133., 134.]]]]),
               (9.10000038, 42., 42., [[[[135., 136., 137., 138., 139.], [140., 141., 142., 143., 144.], [145., 146., 147., 148., 149.]]]])],
               dtype=(numpy.record, [('abc', '<f4'), ('xyz', '<f4'), ('_xyz', '<f4'), ('DATA', '<f4', (1, 1, 3, 5))]))

..
  EXAMPLE END

Compressed Image Data
=====================
.. _astropy-io-fits-compressedImageData:

A general technique has been developed for storing compressed image data in
FITS binary tables. The principle used in this convention is to first divide
the n-dimensional image into a rectangular grid of sub-images or 'tiles'.
Each tile is then compressed as a continuous block of data, and the resulting
compressed byte stream is stored in a row of a variable length column in a
FITS binary table. Several commonly used algorithms for compressing image
tiles are supported. These include Gzip, Rice, IRAF Pixel List (PLIO), and
Hcompress.

For more details, reference "A FITS Image Compression Proposal" from:

    https://www.adass.org/adass/proceedings/adass99/P2-42/

and "Registered FITS Convention, Tiled Image Compression Convention":

    https://fits.gsfc.nasa.gov/registry/tilecompression.html

Compressed image data is accessed, in ``astropy``, using the optional
``astropy.io.fits.compression`` module contained in a C shared library
(compression.so). If an attempt is made to access an HDU containing compressed
image data when the compression module is not available, the user is notified
of the problem and the HDU is treated like a standard binary table HDU. This
notification will only be made the first time compressed image data is
encountered. In this way, the compression module is not required in order for
``astropy`` to work.


Header and Summary
------------------

In ``astropy``, the header of a compressed image HDU appears to the user like
any image header. The actual header stored in the FITS file is that of a binary
table HDU with a set of special keywords, defined by the convention, to
describe the structure of the compressed image. The conversion between binary
table HDU header and image HDU header is all performed behind the scenes.
Since the HDU is actually a binary table, it may not appear as a primary HDU in
a FITS file.

Example
^^^^^^^

..
  EXAMPLE START
  Accessing Compressed FITS Image HDU Headers

The content of the HDU header may be accessed using the ``.header`` attribute::

    >>> filename = fits.util.get_testdata_filepath('compressed_image.fits')

    >>> hdul = fits.open(filename)
    >>> hdul[1].header
    XTENSION= 'IMAGE   '           / Image extension
    BITPIX  =                   16 / data type of original image
    NAXIS   =                    2 / dimension of original image
    NAXIS1  =                   10 / length of original image axis
    NAXIS2  =                   10 / length of original image axis
    PCOUNT  =                    0 / number of parameters
    GCOUNT  =                    1 / number of groups

The contents of the corresponding binary table HDU may be accessed using the
hidden ``._header`` attribute. However, all user interface with the HDU header
should be accomplished through the image header (the ``.header`` attribute)::

    >>> hdul[1]._header
    XTENSION= 'BINTABLE'           / binary table extension
    BITPIX  =                    8 / array data type
    NAXIS   =                    2 / number of array dimensions
    NAXIS1  =                    8 / width of table in bytes
    NAXIS2  =                   10 / number of rows in table
    PCOUNT  =                   60 / number of group parameters
    GCOUNT  =                    1 / number of groups
    TFIELDS =                    1 / number of fields in each row
    TTYPE1  = 'COMPRESSED_DATA'    / label for field 1
    TFORM1  = '1PB(6)  '           / data format of field: variable length array
    ZIMAGE  =                    T / extension contains compressed image
    ZTENSION= 'IMAGE   '           / Image extension
    ZBITPIX =                   16 / data type of original image
    ZNAXIS  =                    2 / dimension of original image
    ZNAXIS1 =                   10 / length of original image axis
    ZNAXIS2 =                   10 / length of original image axis
    ZPCOUNT =                    0 / number of parameters
    ZGCOUNT =                    1 / number of groups
    ZTILE1  =                   10 / size of tiles to be compressed
    ZTILE2  =                    1 / size of tiles to be compressed
    ZCMPTYPE= 'RICE_1  '           / compression algorithm
    ZNAME1  = 'BLOCKSIZE'          / compression block size
    ZVAL1   =                   32 / pixels per block
    ZNAME2  = 'BYTEPIX '           / bytes per pixel (1, 2, 4, or 8)
    ZVAL2   =                    2 / bytes per pixel (1, 2, 4, or 8)
    EXTNAME = 'COMPRESSED_IMAGE'   / name of this binary table extension

The contents of the HDU can be summarized by using either the :func:`info`
convenience function or method::

    >>> fits.info(filename)
    Filename: ...compressed_image.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU       4   ()
      1  COMPRESSED_IMAGE    1 CompImageHDU      7   (10, 10)   int16

    >>> hdul.info()
    Filename: ...compressed_image.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU       4   ()
      1  COMPRESSED_IMAGE    1 CompImageHDU      7   (10, 10)   int16

..
  EXAMPLE END

Data
----

As with the header, the data of a compressed image HDU appears to the user as
standard uncompressed image data. The actual data is stored in the FITS file
as Binary Table data containing at least one column (COMPRESSED_DATA). Each
row of this variable length column contains the byte stream that was generated
as a result of compressing the corresponding image tile. Several optional
columns may also appear. These include UNCOMPRESSED_DATA to hold the
uncompressed pixel values for tiles that cannot be compressed, ZSCALE and ZZERO
to hold the linear scale factor and zero point offset which may be needed to
transform the raw uncompressed values back to the original image pixel values,
and ZBLANK to hold the integer value used to represent undefined pixels (if
any) in the image.

Example
^^^^^^^

..
  EXAMPLE START
  Accessing Compressed FITS Image HDU Data

The contents of the uncompressed HDU data may be accessed using the ``.data``
attribute::

    >>> hdul[1].data
    array([[ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9],
           [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
           [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
           [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
           [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
           [50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
           [60, 61, 62, 63, 64, 65, 66, 67, 68, 69],
           [70, 71, 72, 73, 74, 75, 76, 77, 78, 79],
           [80, 81, 82, 83, 84, 85, 86, 87, 88, 89],
           [90, 91, 92, 93, 94, 95, 96, 97, 98, 99]], dtype=int16)
    >>> hdul.close()

The compressed data can be accessed via the ``.compressed_data`` attribute, but
this rarely needs be accessed directly. It may be useful for performing direct
copies of the compressed data without needing to decompress it first.

..
  EXAMPLE END


Creating a Compressed Image HDU
-------------------------------

To create a compressed image HDU from scratch, construct a
:class:`CompImageHDU` object from an uncompressed image data array and its
associated image header. From there, the HDU can be treated just like any
other image HDU.

Example
^^^^^^^

..
  EXAMPLE START
  Creating a Compressed FITS Image HDU

To create a compressed image HDU::

    >>> imageData = np.arange(100).astype('i2').reshape(10, 10)
    >>> imageHeader = fits.Header()
    >>> hdu = fits.CompImageHDU(imageData, imageHeader)
    >>> hdu.writeto('compressed_image.fits')

The API documentation for the :class:`CompImageHDU` initializer method
describes the possible options for constructing a :class:`CompImageHDU` object.

..
  EXAMPLE END
