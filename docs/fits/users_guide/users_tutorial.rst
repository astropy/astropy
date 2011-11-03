.. currentmodule:: pyfits.core

**************
Quick Tutorial
**************

This chapter provides a quick introduction of using PyFITS. The goal is to
demonstrate PyFITS's basic features without getting into too much detail. If
you are a first time user or an occasional PyFITS user, using only the most
basic functionality, this is where you should start. Otherwise, it is safe to
skip this chapter.

After installing numpy and PyFITS, start Python and load the PyFITS library.
Note that the module name is all lower case.

    >>> import pyfits


Reading and Updating Existing FITS Files
========================================


Opening a FITS file
-------------------

Once the PyFITS module is loaded, we can open an existing FITS file:

    >>> hdulist = pyfits.open('input.fits')

The open() function has several optional arguments which will be discussed in a
later chapter. The default mode, as in the above example, is "readonly".  The
open method returns a PyFITS object called an `HDUList` which is a Python-like
list, consisting of HDU objects. An HDU (Header Data Unit) is the highest level
component of the FITS file structure. So, after the above open call,
``hdulist[0]`` is the primary HDU, ``hdulist[1]``, if any, is the first
extension HDU, etc.  It should be noted that PyFITS is using zero-based
indexing when referring to HDUs and header cards, though the FITS standard
(which was designed with FORTRAN in mind) uses one-based indexing.

The `HDUList` has a useful method `HDUList.info()`, which summarizes the
content of the opened FITS file:

    >>> hdulist.info()
    Filename: test1.fits
    No. Name  Type       Cards Dimensions Format
    0 PRIMARY PrimaryHDU   220 ()         int16
    1 SCI     ImageHDU      61 (800, 800) float32
    2 SCI     ImageHDU      61 (800, 800) float32
    3 SCI     ImageHDU      61 (800, 800) float32
    4 SCI     ImageHDU      61 (800, 800) float32

After you are done with the opened file, close it with the `HDUList.close()`
method:

    >>> hdulist.close()

The headers will still be accessible after the HDUlist is closed. The data may
or may not be accessible depending on whether the data are touched and if they
are memory-mapped, see later chapters for detail.

Working with large files
^^^^^^^^^^^^^^^^^^^^^^^^

The `pyfits.open()` function supports a ``memmap=True`` argument that cause
the array data of each HDU to be accessed with mmap, rather than being read
into memory all at once.  This is particularly useful for working with very
large arrays that cannot fit entirely into physical memory.

This has minimal impact on smaller files as well, though some operations, such
as reading the array data sequentially, may incur some additional overhead.  On
32-bit systems arrays larger than 2-3 GB cannot be mmap'd (which is fine,
because by that point you're likely to run out of physical memory anyways), but
64-bit systems are much less limited in this respect.


Working With a FITS Header
--------------------------

As mentioned earlier, each element of an HDUList is an HDU object with
attributes of header and data, which can be used to access the header keywords
and the data.

The header attribute is a Header instance, another PyFITS object. To get the
value of a header keyword, simply do (a la Python dictionaries):

    >>> hdulist[0].header['targname']
    'NGC121'

to get the value of the keyword targname, which is a string 'NGC121'.

Although keyword names are always in upper case inside the FITS file,
specifying a keyword name with PyFITS is case-insensitive, for user's
convenience. If the specified keyword name does not exist, it will raise a
KeyError exception.

We can also get the keyword value by indexing (a la Python lists):

    >>> hdulist[0].header[27]
    96

This example returns the 28th (like Python lists, it is 0-indexed) keyword's
value, an integer, 96.

Similarly, it is easy to update a keyword's value in PyFITS, either through
keyword name or index:

    >>> prihdr = hdulist[0].header
    >>> prihdr['targname'] = 'NGC121-a'
    >>> prihdr[27] = 99

Use the above syntax if the keyword is already present in the header. If the
keyword might not exist and you want to add it if it doesn't, use the
`Header.update()` method:

    >>> prihdr.update('observer', 'Edwin Hubble')

Special methods must be used to add comment or history records:

    >>> prihdr.add_history('I updated this file 2/26/09')
    >>> prihdr.add_comment('Edwin Hubble really knew his stuff')

A header consists of `Card` objects (i.e. the 80-column card-images specified
in the FITS standard). Each Card normally has up to three parts: key, value,
and comment. To see the entire list of cardimages of an HDU, use the
`Header.ascardlist()` method :

    >>> print prihdr.ascardlist()[:3]
    SIMPLE  =                    T / file does conform to FITS standard
    BITPIX  =                   16 / number of bits per data pixel
    NAXIS   =                    0 / number of data axes

Only the first three cards are shown above.

To get a list of all keywords, use the `CardList.keys()` method of the card
list:

    >>> prihdr.ascardlist().keys()
    ['SIMPLE', 'BITPIX', 'NAXIS', ...]


Working With Image Data
-----------------------

If an HDU's data is an image, the data attribute of the HDU object will return
a numpy ndarray object. Refer to the numpy documentation for details on
manipulating these numerical arrays.

    >>> scidata = hdulist[1].data

Here, scidata points to the data object in the second HDU (the first HDU,
``hdulist[0]``, being the primary HDU) in ``hdulist``, which corresponds  to
the 'SCI' extension. Alternatively, you can access the extension by its
extension name (specified in the EXTNAME keyword):

    >>> scidata = hdulist['SCI'].data

If there is more than one extension with the same EXTNAME, EXTVER's value needs
to be specified as the second argument, e.g.:

    >>> scidata = hdulist['sci',2].data

The returned numpy object has many attributes and methods for a user to get
information about the array, e. g.:

    >>> scidata.shape
    (800, 800)
    >>> scidata.dtype.name
    'float32'

Since image data is a numpy object, we can slice it, view it, and perform
mathematical operations on it. To see the pixel value at x=5, y=2:

    >>> print scidata[1,4]

Note that, like C (and unlike FORTRAN), Python is 0-indexed and the indices
have the slowest axis first and fast axis last, i.e. for a 2-D image, the fast
axis (X-axis) which corresponds to the FITS NAXIS1 keyword, is the second
index. Similarly, the sub-section of x=11 to 20 (inclusive) and y=31 to 40
(inclusive) is:

    >>> scidata[30:40, 10:20]

To update the value of a pixel or a sub-section:

    >>> scidata[30:40,10:20] = scidata[1,4] = 999

This example changes the values of both the pixel \[1,4] and the sub-section
\[30:40,10:20] to the new value of 999.

The next example of array manipulation is to convert the image data from counts
to flux:

    >>> photflam = hdulist[1].header['photflam']
    >>> exptime = prihdr['exptime']
    >>> scidata *= photflam / exptime

This example performs the math on the array in-place, thereby keeping the
memory usage to a minimum.

If at this point you want to preserve all the changes you made and write it to
a new file, you can use the `HDUList.writeto()` method (see below).


Working With Table Data
-----------------------

If you are familiar with the record array in numpy, you will find the table
data is basically a record array with some extra properties. But familiarity
with record arrays is not a prerequisite for this Guide.

Like images, the data portion of a FITS table extension is in the ``.data``
attribute:

    >>> hdulist = pyfits.open('table.fits')
    >>> tbdata = hdulist[1].data # assuming the first extension is a table

To see the first row of the table:

    >>> print tbdata[0]
    (1, 'abc', 3.7000002861022949, 0)

Each row in the table is a `FITS_rec` object which looks like a (Python) tuple
containing elements of heterogeneous data types. In this example: an integer, a
string, a floating point number, and a Boolean value. So the table data are
just an array of such records. More commonly, a user is likely to access the
data in a column-wise way. This is accomplished by using the field() method. To
get the first column (or field) of the table, use:

    >>> tbdata.field(0)
    array([1, 2])

A numpy object with the data type of the specified field is returned.

Like header keywords, a field can be referred either by index, as above, or by
name:

    >>> tbdata.field('id')
    array([1, 2])

But how do we know what field names we've got? First, let's introduce another
attribute of the table HDU: the ``.columns`` attribute:

    >>> cols = hdulist[1].columns

This attribute is a `ColDefs` (column definitions) object. If we use the
`ColDefs.info()` method:

    >>> cols.info()
     name:
          ['c1', 'c2', 'c3', 'c4']
     format:
          ['1J', '3A', '1E', '1L']
     unit:
          ['', '', '', '']
     null:
          [-2147483647, '', '', '']
     bscale:
          ['', '', 3, '']
     bzero:
          ['', '', 0.40000000000000002, '']
     disp:
          ['I11', 'A3', 'G15.7', 'L6']
     start:
          ['', '', '', '']
     dim:
          ['', '', '', '']

it will show all its attributes, such as names, formats, bscales, bzeros, etc.
We can also get these properties individually, e.g.:

    >>> cols.names
    ['ID', 'name', 'mag', 'flag']

returns a (Python) list of field names.

Since each field is a numpy object, we'll have the entire arsenal of numpy
tools to use. We can reassign (update) the values:

    >>> tbdata.field('flag')[:] = 0


Save File Changes
-----------------

As mentioned earlier, after a user opened a file, made a few changes to either
header or data, the user can use `HDUList.writeto()` to save the changes. This
takes the version of headers and data in memory and writes them to a new FITS
file on disk. Subsequent operations can be performed to the data in memory and
written out to yet another different file, all without recopying the original
data to (more) memory.

    >>> hdulist.writeto('newimage.fits')

will write the current content of ``hdulist`` to a new disk file newfile.fits.
If a file was opened with the update mode, the `HDUList.flush()` method can
also be used to write all the changes made since ``open()``, back to the
original file. The ``close()`` method will do the same for a FITS file opened
with update mode.

    >>> f = pyfits.open('original.fits', mode='update')
    ... # making changes in data and/or header
    >>> f.flush() # changes are written back to original.fits


Creating a New FITS File
========================


Creating a New Image File
-------------------------

So far we have demonstrated how to read and update an existing FITS file. But
how about creating a new FITS file from scratch? Such task is very easy in
PyFITS for an image HDU. We'll first demonstrate how to create a FITS file
consisting only the primary HDU with image data.

First, we create a numpy object for the data part:

    >>> import numpy as np
    >>> n = np.arange(100) # a simple sequence from 0 to 99

Next, we create a `PrimaryHDU` object to encapsulate the data:

    >>> hdu = pyfits.PrimaryHDU(n)

We then create a HDUList to contain the newly created primary HDU, and write to
a new file:

    >>> hdulist = pyfits.HDUList([hdu])
    >>> hdulist.writeto('new.fits')

That's it! In fact, PyFITS even provides a short cut for the last two lines to
accomplish the same behavior:

    >>> hdu.writeto('new.fits')


Creating a New Table File
-------------------------

To create a table HDU is a little more involved than image HDU, because a
table's structure needs more information. First of all, tables can only be an
extension HDU, not a primary. There are two kinds of FITS table extensions:
ASCII and binary. We'll use binary table examples here.

To create a table from scratch, we need to define columns first, by
constructing the `Column` objects and their data. Suppose we have two columns,
the first containing strings, and the second containing floating point numbers:

    >>> import pyfits
    >>> import numpy as np
    >>> a1 = np.array(['NGC1001', 'NGC1002', 'NGC1003'])
    >>> a2 = np.array([11.1, 12.3, 15.2])
    >>> col1 = pyfits.Column(name='target', format='20A', array=a1)
    >>> col2 = pyfits.Column(name='V_mag', format='E', array=a2)

Next, create a `ColDefs` (column-definitions) object for all columns:

    >>> cols = pyfits.ColDefs([col1, col2])

Now, create a new binary table HDU object by using the PyFITS function
`new_table()`:

    >>> tbhdu = pyfits.new_table(cols)

This function returns (in this case) a `BinTableHDU`.

Of course, you can do this more concisely:

    >>> tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='target',
    ...                                                        format='20A',
    ...                                                        array=a1),
    ...                                          pyfits.Column(name='V_mag',
    ...                                                        format='E',
    ...                                                        array=a2)]
    ...                                        ))

As before, we create a `PrimaryHDU` object to encapsulate the data:

    >>> hdu = pyfits.PrimaryHDU(n)

We then create a HDUList containing both the primary HDU and the newly created
table extension, and write to a new file:

    >>> thdulist = pyfits.HDUList([hdu, tbhdu])
    >>> thdulist.writeto('table.fits')

If this will be the only extension of the new FITS file and you only have a
minimal primary HDU with no data, PyFITS again provides a short cut:

    >>> tbhdu.writeto('table.fits')

Alternatively, you can append it to the hdulist we have already created from
the image file section:

    >>> hdulist.append(tbhdu)

So far, we have covered the most basic features of PyFITS. In the following
chapters we'll show more advanced examples and explain options in each class
and method.


Convenience Functions
=====================

PyFITS also provides several high level ("convenience") functions. Such a
convenience function is a "canned" operation to achieve one simple task. By
using these "convenience" functions, a user does not have to worry about
opening or closing a file, all the housekeeping is done implicitly.

The first of these functions is `getheader()`, to get the header of an HDU.
Here are several examples of getting the header. Only the file name is required
for this function. The rest of the arguments are optional and flexible to
specify which HDU the user wants to get:

    >>> from pyfits import getheader
    >>> getheader('in.fits') # get default HDU (=0), i.e. primary HDU's header
    >>> getheader('in.fits', 0) # get primary HDU's header
    >>> getheader('in.fits', 2) # the second extension
    # the HDU with EXTNAME='sci' (if there is only 1)
    >>> getheader('in.fits', 'sci')
    # the HDU with EXTNAME='sci' and EXTVER=2
    >>> getheader('in.fits', 'sci', 2)
    >>> getheader('in.fits', ('sci', 2)) # use a tuple to do the same
    >>> getheader('in.fits', ext=2) # the second extension
    # the 'sci' extension, if there is only 1
    >>> getheader('in.fits', extname='sci')
    # the HDU with EXTNAME='sci' and EXTVER=2
    >>> getheader('in.fits', extname='sci', extver=2)
    # ambiguous specifications will raise an exception, DON'T DO IT!!
    >>> getheader('in.fits', ext=('sci',1), extname='err', extver=2)

After you get the header, you can access the information in it, such as getting
and modifying a keyword value:

    >>> from pyfits import getheader
    >>> hdr = getheader('in.fits', 1) # get first extension's header
    >>> filter = hdr['filter'] # get the value of the keyword "filter'
    >>> val = hdr[10] # get the 11th keyword's value
    >>> hdr['filter'] = 'FW555' # change the keyword value

For the header keywords, the header is like a dictionary, as well as a list.
The user can access the keywords either by name or by numeric index, as
explained earlier in this chapter.

If a user only needs to read one keyword, the `getval()` function can further
simplify to just one call, instead of two as shown in the above examples:

    >>> from pyfits import getval
    >>> flt = getval('in.fits', 'filter', 1) # get 1st extension's keyword
                                             # FILTER's value
    >>> val = getval('in.fits', 10, 'sci', 2) # get the 2nd sci extension's
                                              # 11th keyword's value

The function `getdata()` gets the data of an HDU. Similar to `getheader()`, it
only requires the input FITS file name while the extension is specified through
the optional arguments. It does have one extra optional argument header. If
header is set to True, this function will return both data and header,
otherwise only data is returned.

    >>> from pyfits import getdata
    >>> dat = getdata('in.fits', 'sci', 3) # get 3rd sci extension's data
    # get 1st extension's data and header
    >>> data, hdr = getdata('in.fits', 1, header=True)

The functions introduced above are for reading. The next few functions
demonstrate convenience functions for writing:

    >>> pyfits.writeto('out.fits', data, header)

The `writeto()` function uses the provided data and an optional header to write
to an output FITS file.

    >>> pyfits.append('out.fits', data, header)

The `append()` function will use the provided data and the optional header to
append to an existing FITS file. If the specified output file does not exist,
it will create one.

    >>> from pyfits import update
    >>> update(file, dat, hdr, 'sci') # update the 'sci' extension
    >>> update(file, dat, 3) # update the 3rd extension
    >>> update(file, dat, hdr, 3) # update the 3rd extension
    >>> update(file, dat, 'sci', 2) # update the 2nd SCI extension
    >>> update(file, dat, 3, header=hdr) # update the 3rd extension
    >>> update(file, dat, header=hdr, ext=5) # update the 5th extension

The `update()` function will update the specified extension with the input
data/header. The 3rd argument can be the header associated with the data. If
the 3rd argument is not a header, it (and other positional arguments) are
assumed to be the extension specification(s). Header and extension specs can
also be keyword arguments.

Finally, the `info()` function will print out information of the specified FITS
file:

    >>> pyfits.info('test0.fits')
    Filename: test0.fits
    No. Name Type Cards Dimensions Format
    0 PRIMARY PrimaryHDU 138 () Int16
    1 SCI ImageHDU 61 (400, 400) Int16
    2 SCI ImageHDU 61 (400, 400) Int16
    3 SCI ImageHDU 61 (400, 400) Int16
    4 SCI ImageHDU 61 (400, 400) Int16
