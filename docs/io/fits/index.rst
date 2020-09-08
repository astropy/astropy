.. currentmodule:: astropy.io.fits

.. _astropy-io-fits:

**************************************
FITS File Handling (`astropy.io.fits`)
**************************************

Introduction
============

The :mod:`astropy.io.fits` package provides access to FITS files. FITS
(Flexible Image Transport System) is a portable file standard widely used in
the astronomy community to store images and tables.

.. _tutorial:

Getting Started
===============

This section provides a quick introduction of using :mod:`astropy.io.fits`. The
goal is to demonstrate the package's basic features without getting into too
much detail. If you are a first time user or have never used ``astropy`` or
PyFITS, this is where you should start. See also the :ref:`FAQ <io-fits-faq>`
for answers to common questions and issues.

.. note::

    If you want to read or write a single table in FITS format, the
    recommended method is via the high-level :ref:`table_io`. In particular
    see the :ref:`Unified I/O FITS <table_io_fits>` section.

Reading and Updating Existing FITS Files
----------------------------------------

Opening a FITS File
^^^^^^^^^^^^^^^^^^^

.. note::

    The ``astropy.io.fits.util.get_testdata_filepath()`` function,
    used in the examples here, is for accessing data shipped with ``astropy``.
    To work with your own data instead, please use :func:`astropy.io.fits.open`,
    which takes either the relative or absolute path.

Once the `astropy.io.fits` package is loaded using the standard convention
[#f1]_, we can open an existing FITS file::

    >>> from astropy.io import fits
    >>> fits_image_filename = fits.util.get_testdata_filepath('test0.fits')

    >>> hdul = fits.open(fits_image_filename)

The :func:`open` function has several optional arguments which will be
discussed in a later chapter. The default mode, as in the above example, is
"readonly". The open function returns an object called an :class:`HDUList`
which is a `list`-like collection of HDU objects. An HDU (Header Data Unit) is
the highest level component of the FITS file structure, consisting of a header
and (typically) a data array or table.

After the above open call, ``hdul[0]`` is the primary HDU, ``hdul[1]`` is
the first extension HDU, etc. (if there are any extensions), and so on. It
should be noted that ``astropy`` uses zero-based indexing when referring to
HDUs and header cards, though the FITS standard (which was designed with
Fortran in mind) uses one-based indexing.

The :class:`HDUList` has a useful method :meth:`HDUList.info`, which
summarizes the content of the opened FITS file:

    >>> hdul.info()
    Filename: ...test0.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     138   ()
      1  SCI           1 ImageHDU        61   (40, 40)   int16
      2  SCI           2 ImageHDU        61   (40, 40)   int16
      3  SCI           3 ImageHDU        61   (40, 40)   int16
      4  SCI           4 ImageHDU        61   (40, 40)   int16

After you are done with the opened file, close it with the
:meth:`HDUList.close` method:

    >>> hdul.close()

You can avoid closing the file manually by using :func:`open` as context
manager::

    >>> with fits.open(fits_image_filename) as hdul:
    ...     hdul.info()
    Filename: ...test0.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     138   ()
      1  SCI           1 ImageHDU        61   (40, 40)   int16
      2  SCI           2 ImageHDU        61   (40, 40)   int16
      3  SCI           3 ImageHDU        61   (40, 40)   int16
      4  SCI           4 ImageHDU        61   (40, 40)   int16

After exiting the ``with`` scope the file will be closed automatically. That is
(generally) the preferred way to open a file in Python, because it will close
the file even if an exception happens.

If the file is opened with ``lazy_load_hdus=False``, all of the headers will
still be accessible after the HDUList is closed. The headers and data may or
may not be accessible depending on whether the data are touched and if they
are memory-mapped; see later chapters for detail.

.. _fits-large-files:

Working with large files
""""""""""""""""""""""""

The :func:`open` function supports a ``memmap=True`` argument that allows the
array data of each HDU to be accessed with mmap, rather than being read into
memory all at once. This is particularly useful for working with very large
arrays that cannot fit entirely into physical memory. Here ``memmap=True`` by
default, and this value is obtained from the configuration item ``astropy.io.fits.Conf.use_memmap``.

This has minimal impact on smaller files as well, though some operations, such
as reading the array data sequentially, may incur some additional overhead. On
32-bit systems, arrays larger than 2 to 3 GB cannot be mmap'd (which is fine,
because by that point you are likely to run out of physical memory anyways), but
64-bit systems are much less limited in this respect.

.. warning::
    When opening a file with ``memmap=True``, because of how mmap works this
    means that when the HDU data is accessed (i.e., ``hdul[0].data``) another
    handle to the FITS file is opened by mmap. This means that even after
    calling ``hdul.close()`` the mmap still holds an open handle to the data so
    that it can still be accessed by unwary programs that were built with the
    assumption that the .data attribute has all of the data in-memory.

    In order to force the mmap to close, either wait for the containing
    ``HDUList`` object to go out of scope, or manually call
    ``del hdul[0].data``. (This works so long as there are no other references
    held to the data array.)

Unsigned integers
"""""""""""""""""

Due to the FITS format's Fortran origins, FITS does not natively support
unsigned integer data in images or tables. However, there is a common
convention to store unsigned integers as signed integers, along with a
*shift* instruction (a ``BZERO`` keyword with value ``2 ** (BITPIX - 1)``) to
shift up all signed integers to unsigned integers. For example, when writing
the value ``0`` as an unsigned 32-bit integer, it is stored in the FITS
file as ``-32768``, along with the header keyword ``BZERO = 32768``.

``astropy`` recognizes and applies this convention by default, so that all data
that looks like it should be interpreted as unsigned integers is automatically
converted (this applies to both images and tables). In ``astropy`` versions
prior to v1.1.0 this was *not* applied automatically, and it is necessary to
pass the argument ``uint=True`` to :func:`open`. In v1.1.0 or later this is the
default.

Even with ``uint=False``, the ``BZERO`` shift is still applied, but the
returned array is of "float64" type. To disable scaling/shifting entirely, use
``do_not_scale_image_data=True`` (see :ref:`fits-scaled-data-faq` in the FAQ
for more details).

Working with compressed files
"""""""""""""""""""""""""""""

.. note::

    Files that use compressed HDUs within the FITS file are discussed
    in :ref:`Compressed Image Data <astropy-io-fits-compressedImageData>`.


The :func:`open` function will seamlessly open FITS files that have been
compressed with gzip, bzip2 or pkzip. Note that in this context we are talking
about a FITS file that has been compressed with one of these utilities (e.g., a
.fits.gz file).

There are some limitations when working with compressed files. For example,
with Zip files that contain multiple compressed files, only the first file will
be accessible. Also bzip2 does not support the append or update access modes.

When writing a file (e.g., with the :func:`writeto` function), compression will
be determined based on the filename extension given, or the compression used in
a pre-existing file that is being written to.


Working with non-standard files
"""""""""""""""""""""""""""""""
When `astropy.io.fits` reads a FITS file which does not conform to the FITS standard it will try
to make an educated interpretation of non-compliant fields. This may not always
succeed and may trigger warnings when accessing headers or exceptions when writing
to file. Verification of fields written to an output file can be controlled with
the ``output_verify`` parameter of :func:`open`. Files opened for reading can be
verified and fixed with method ``HDUList.verify``. This method is typically invoked
after opening the file but before accessing any headers or data::

    >>> with fits.open(fits_image_filename) as hdul:
    ...    hdul.verify('fix')
    ...    data = hdul[1].data

In the above example, the call to ``hdul.verify("fix")`` requests that `astropy.io.fits`
fix non-compliant fields and print informative messages. Other options in addition to ``"fix"``
are described under FITS :ref:`fits_io_verification`

.. seealso:: FITS :ref:`fits_io_verification`.

Working with FITS Headers
^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned earlier, each element of an :class:`HDUList` is an HDU object with
``.header`` and ``.data`` attributes, which can be used to access the header
and data portions of the HDU.

For those unfamiliar with FITS headers, they consist of a list of 80 byte
"cards", where a card contains a keyword, a value, and a comment. The keyword
and comment must both be strings, whereas the value can be a string or an
integer, floating point number, complex number, or ``True``/``False``. Keywords
are usually unique within a header, except in a few special cases.

The header attribute is a Header instance, another ``astropy`` object. To get
the value associated with a header keyword, do (à la Python dicts)::

    >>> hdul = fits.open(fits_image_filename)
    >>> hdul[0].header['DATE']
    '01/04/99'

to get the value of the keyword "DATE", which is a string '01/04/99'.

Although keyword names are always in upper case inside the FITS file,
specifying a keyword name with ``astropy`` is case-insensitive for the user's
convenience. If the specified keyword name does not exist, it will raise a
`KeyError` exception.

We can also get the keyword value by indexing (à la Python lists)::

    >>> hdul[0].header[7]
    32768.0

This example returns the eighth (like Python lists, it is 0-indexed) keyword's
value — a float — 32768.0.

Similarly, it is possible to update a keyword's value in ``astropy``, either
through keyword name or index::

    >>> hdr = hdul[0].header
    >>> hdr['targname'] = 'NGC121-a'
    >>> hdr[27] = 99

Please note however that almost all application code should update header
values via their keyword name and not via their positional index. This is
because most FITS keywords may appear at any position in the header.

It is also possible to update both the value and comment associated with a
keyword by assigning them as a tuple::

    >>> hdr = hdul[0].header
    >>> hdr['targname'] = ('NGC121-a', 'the observation target')
    >>> hdr['targname']
    'NGC121-a'
    >>> hdr.comments['targname']
    'the observation target'

Like a dict, you may also use the above syntax to add a new keyword/value pair
(and optionally a comment as well). In this case the new card is appended to
the end of the header (unless it is a commentary keyword such as COMMENT or
HISTORY, in which case it is appended after the last card with that keyword).

Another way to either update an existing card or append a new one is to use the
:meth:`Header.set` method::

    >>> hdr.set('observer', 'Edwin Hubble')

Comment or history records are added like normal cards, though in their case a
new card is always created, rather than updating an existing HISTORY or COMMENT
card::

    >>> hdr['history'] = 'I updated this file 2/26/09'
    >>> hdr['comment'] = 'Edwin Hubble really knew his stuff'
    >>> hdr['comment'] = 'I like using HST observations'
    >>> hdr['history']
    I updated this file 2/26/09
    >>> hdr['comment']
    Edwin Hubble really knew his stuff
    I like using HST observations

Note: Be careful not to confuse COMMENT cards with the comment value for normal
cards.

To update existing COMMENT or HISTORY cards, reference them by index::

    >>> hdr['history'][0] = 'I updated this file on 2/27/09'
    >>> hdr['history']
    I updated this file on 2/27/09
    >>> hdr['comment'][1] = 'I like using JWST observations'
    >>> hdr['comment']
    Edwin Hubble really knew his stuff
    I like using JWST observations


To see the entire header as it appears in the FITS file (with the END card and
padding stripped), enter the header object by itself, or
``print(repr(hdr))``::

    >>> hdr  # doctest: +ELLIPSIS
    SIMPLE  =                    T / file does conform to FITS standard
    BITPIX  =                   16 / number of bits per data pixel
    NAXIS   =                    0 / number of data axes
    ...
    >>> print(repr(hdr))  # doctest: +ELLIPSIS
    SIMPLE  =                    T / file does conform to FITS standard
    BITPIX  =                   16 / number of bits per data pixel
    NAXIS   =                    0 / number of data axes
    ...

Entering only ``print(hdr)`` will also work, but may not be very legible
on most displays, as this displays the header as it is written in the FITS file
itself, which means there are no line breaks between cards. This is a common
source of confusion for new users.

It is also possible to view a slice of the header::

   >>> hdr[:2]
   SIMPLE  =                    T / file does conform to FITS standard
   BITPIX  =                   16 / number of bits per data pixel

Only the first two cards are shown above.

To get a list of all keywords, use the :meth:`Header.keys` method just as you
would with a dict::

    >>> list(hdr.keys())  # doctest: +ELLIPSIS
    ['SIMPLE', 'BITPIX', 'NAXIS', ...]

.. topic:: Examples:

    See also :ref:`sphx_glr_generated_examples_io_modify-fits-header.py`.

Working with Image Data
^^^^^^^^^^^^^^^^^^^^^^^

If an HDU's data is an image, the data attribute of the HDU object will return
a ``numpy`` `~numpy.ndarray` object. Refer to the ``numpy`` documentation for
details on manipulating these numerical arrays::

    >>> data = hdul[1].data

Here, ``data`` points to the data object in the second HDU (the first HDU,
``hdul[0]``, being the primary HDU) which corresponds to the 'SCI'
extension. Alternatively, you can access the extension by its extension name
(specified in the EXTNAME keyword)::

    >>> data = hdul['SCI'].data

If there is more than one extension with the same EXTNAME, the EXTVER value
needs to be specified along with the EXTNAME as a tuple; for example::

    >>> data = hdul['sci',2].data

Note that the EXTNAME is also case-insensitive.

The returned ``numpy`` object has many attributes and methods for a user to get
information about the array, for example::

    >>> data.shape
    (40, 40)
    >>> data.dtype.name
    'int16'

Since image data is a ``numpy`` object, we can slice it, view it, and perform
mathematical operations on it. To see the pixel value at x=5, y=2::

    >>> print(data[1, 4])
    348

Note that, like C (and unlike Fortran), Python is 0-indexed and the indices
have the slowest axis first and fastest changing axis last; that is, for a 2D
image, the fast axis (X-axis) which corresponds to the FITS NAXIS1 keyword, is
the second index. Similarly, the 1-indexed subsection of x=11 to 20
(inclusive) and y=31 to 40 (inclusive) would be given in Python as::

    >>> data[30:40, 10:20]
    array([[350, 349, 349, 348, 349, 348, 349, 347, 350, 348],
           [348, 348, 348, 349, 348, 349, 347, 348, 348, 349],
           [348, 348, 347, 349, 348, 348, 349, 349, 349, 349],
           [349, 348, 349, 349, 350, 349, 349, 347, 348, 348],
           [348, 348, 348, 348, 349, 348, 350, 349, 348, 349],
           [348, 347, 349, 349, 350, 348, 349, 348, 349, 347],
           [347, 348, 347, 348, 349, 349, 350, 349, 348, 348],
           [349, 349, 350, 348, 350, 347, 349, 349, 349, 348],
           [349, 348, 348, 348, 348, 348, 349, 347, 349, 348],
           [349, 349, 349, 348, 350, 349, 349, 350, 348, 350]], dtype=int16)

To update the value of a pixel or a subsection::

    >>> data[30:40, 10:20] = data[1, 4] = 999

This example changes the values of both the pixel \[1, 4] and the subsection
\[30:40, 10:20] to the new value of 999. See the `Numpy documentation`_ for
more details on Python-style array indexing and slicing.

The next example of array manipulation is to convert the image data from counts
to flux::

    >>> photflam = hdul[1].header['photflam']
    >>> exptime = hdr['exptime']
    >>> data = data * photflam / exptime
    >>> hdul.close()

Note that performing an operation like this on an entire image requires holding
the entire image in memory. This example performs the multiplication in-place
so that no copies are made, but the original image must first be able to fit in
main memory. For most observations this should not be an issue on modern
personal computers.

If at this point you want to preserve all of the changes you made and write it
to a new file, you can use the :meth:`HDUList.writeto` method (see below).

.. _Numpy documentation: https://numpy.org/doc/stable/reference/arrays.indexing.html

.. topic:: Examples:

    See also :ref:`sphx_glr_generated_examples_io_plot_fits-image.py`.

Working with Table Data
^^^^^^^^^^^^^^^^^^^^^^^

This section describes reading and writing table data in the FITS format using
the `~astropy.io.fits` package directly. For some cases, however, the
high-level :ref:`table_io` will often suffice and is somewhat more convenient
to use. See the :ref:`Unified I/O FITS <table_io_fits>` section for details.

Like images, the data portion of a FITS table extension is in the ``.data``
attribute::

    >>> fits_table_filename = fits.util.get_testdata_filepath('tb.fits')
    >>> hdul = fits.open(fits_table_filename)
    >>> data = hdul[1].data # assuming the first extension is a table

If you are familiar with ``numpy`` `~numpy.recarray` (record array) objects, you
will find the table data is basically a record array with some extra
properties. But familiarity with record arrays is not a prerequisite for this
guide.

To see the first row of the table::

    >>> print(data[0])
    (1, 'abc', 3.7000000715255736, False)

Each row in the table is a :class:`FITS_record` object which looks like a
(Python) tuple containing elements of heterogeneous data types. In this
example: an integer, a string, a floating point number, and a Boolean value. So
the table data are just an array of such records. More commonly, a user is
likely to access the data in a column-wise way. This is accomplished by using
the :meth:`~FITS_rec.field` method. To get the first column (or "field" in
NumPy parlance — it is used here interchangeably with "column") of the table,
use::

    >>> data.field(0)
    array([1, 2]...)

A ``numpy`` object with the data type of the specified field is returned.

Like header keywords, a column can be referred either by index, as above, or by
name::

    >>> data.field('c1')
    array([1, 2]...)

When accessing a column by name, dict-like access is also possible (and even
preferable)::

    >>> data['c1']
    array([1, 2]...)

In most cases it is preferable to access columns by their name, as the column
name is entirely independent of its physical order in the table. As with
header keywords, column names are case-insensitive.

But how do we know what columns we have in a table? First, we will introduce
another attribute of the table HDU: the :attr:`~BinTableHDU.columns`
attribute::

    >>> cols = hdul[1].columns

This attribute is a :class:`ColDefs` (column definitions) object. If we use the
:meth:`ColDefs.info` method from the interactive prompt::

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
        ['', '', 0.4, '']
    disp:
        ['I11', 'A3', 'G15.7', 'L6']
    start:
        ['', '', '', '']
    dim:
        ['', '', '', '']
    coord_type:
        ['', '', '', '']
    coord_unit:
        ['', '', '', '']
    coord_ref_point:
        ['', '', '', '']
    coord_ref_value:
        ['', '', '', '']
    coord_inc:
        ['', '', '', '']
    time_ref_pos:
        ['', '', '', '']

it will show the attributes of all columns in the table, such as their names,
formats, bscales, bzeros, etc. A similar output that will display the column
names and their formats can be printed from within a script with::

    >>> hdul[1].columns
    ColDefs(
        name = 'c1'; format = '1J'; null = -2147483647; disp = 'I11'
        name = 'c2'; format = '3A'; disp = 'A3'
        name = 'c3'; format = '1E'; bscale = 3; bzero = 0.4; disp = 'G15.7'
        name = 'c4'; format = '1L'; disp = 'L6'
    )

We can also get these properties individually; for example::

    >>> cols.names
    ['c1', 'c2', 'c3', 'c4']

returns a (Python) list of field names.

Since each field is a ``numpy`` object, we will have the entire arsenal of
``numpy`` tools to use. We can reassign (update) the values::

    >>> data['c4'][:] = 0

take the mean of a column::

    >>> data['c3'].mean()  # doctest: +FLOAT_CMP
    5.19999989271164

and so on.

.. topic:: Examples:

    See also :ref:`sphx_glr_generated_examples_io_fits-tables.py`.

Save File Changes
^^^^^^^^^^^^^^^^^

As mentioned earlier, after a user opened a file, made a few changes to either
header or data, the user can use :meth:`HDUList.writeto` to save the changes.
This takes the version of headers and data in memory and writes them to a new
FITS file on disk. Subsequent operations can be performed to the data in memory
and written out to yet another different file, all without recopying the
original data to (more) memory:

.. code:: python

    hdul.writeto('newtable.fits')

will write the current content of ``hdulist`` to a new disk file newfile.fits.
If a file was opened with the update mode, the :meth:`HDUList.flush` method can
also be used to write all of the changes made since :func:`open`, back to the
original file. The :meth:`~HDUList.close` method will do the same for a FITS
file opened with update mode:

.. code:: python

    with fits.open('original.fits', mode='update') as hdul:
        # Change something in hdul.
        hdul.flush()  # changes are written back to original.fits

    # closing the file will also flush any changes and prevent further writing


Creating a New FITS File
------------------------

Creating a New Image File
^^^^^^^^^^^^^^^^^^^^^^^^^

So far we have demonstrated how to read and update an existing FITS file. But
how about creating a new FITS file from scratch? Such tasks are very convenient
in ``astropy`` for an image HDU. We will first demonstrate how to create a FITS
file consisting of only the primary HDU with image data.

First, we create a ``numpy`` object for the data part::

    >>> import numpy as np
    >>> n = np.arange(100.0) # a simple sequence of floats from 0.0 to 99.9

Next, we create a :class:`PrimaryHDU` object to encapsulate the data::

    >>> hdu = fits.PrimaryHDU(n)

We then create an :class:`HDUList` to contain the newly created primary HDU, and write to
a new file::

    >>> hdul = fits.HDUList([hdu])
    >>> hdul.writeto('new1.fits')

That is it! In fact, ``astropy`` even provides a shortcut for the last two
lines to accomplish the same behavior::

    >>> hdu.writeto('new2.fits')

This will write a single HDU to a FITS file without having to manually
encapsulate it in an :class:`HDUList` object first.


Creating a New Table File
^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    If you want to create a **binary** FITS table with no other HDUs,
    you can use :class:`~astropy.table.Table` instead and then write to FITS.
    This is less complicated than "lower-level" FITS interface::

    >>> from astropy.table import Table
    >>> t = Table([[1, 2], [4, 5], [7, 8]], names=('a', 'b', 'c'))
    >>> t.write('table1.fits', format='fits')

    The equivalent code using ``astropy.io.fits`` would look like this:

    >>> from astropy.io import fits
    >>> import numpy as np
    >>> c1 = fits.Column(name='a', array=np.array([1, 2]), format='K')
    >>> c2 = fits.Column(name='b', array=np.array([4, 5]), format='K')
    >>> c3 = fits.Column(name='c', array=np.array([7, 8]), format='K')
    >>> t = fits.BinTableHDU.from_columns([c1, c2, c3])
    >>> t.writeto('table2.fits')

To create a table HDU is a little more involved than an image HDU, because a
table's structure needs more information. First of all, tables can only be an
extension HDU, not a primary. There are two kinds of FITS table extensions:
ASCII and binary. We will use binary table examples here.

To create a table from scratch, we need to define columns first, by
constructing the :class:`Column` objects and their data. Suppose we have two
columns, the first containing strings, and the second containing floating point
numbers::

    >>> import numpy as np
    >>> a1 = np.array(['NGC1001', 'NGC1002', 'NGC1003'])
    >>> a2 = np.array([11.1, 12.3, 15.2])
    >>> col1 = fits.Column(name='target', format='20A', array=a1)
    >>> col2 = fits.Column(name='V_mag', format='E', array=a2)

.. note::

    It is not necessary to create a :class:`Column` object explicitly
    if the data is stored in a
    `structured array <https://numpy.org/doc/stable/user/basics.rec.html>`_.

Next, create a :class:`ColDefs` (column-definitions) object for all columns::

    >>> cols = fits.ColDefs([col1, col2])

Now, create a new binary table HDU object by using the
:func:`BinTableHDU.from_columns` function::

    >>> hdu = fits.BinTableHDU.from_columns(cols)

This function returns (in this case) a :class:`BinTableHDU`.

The data structure used to represent FITS tables is called a :class:`FITS_rec`
and is derived from the :class:`numpy.recarray` interface. When creating
a new table HDU the individual column arrays will be assembled into a single
:class:`FITS_rec` array.

You can create a :class:`BinTableHDU` more concisely without creating intermediate
variables for the individual columns and without manually creating a
:class:`ColDefs` object::

    >>> hdu = fits.BinTableHDU.from_columns(
    ...     [fits.Column(name='target', format='20A', array=a1),
    ...      fits.Column(name='V_mag', format='E', array=a2)])

Now you may write this new table HDU directly to a FITS file like so::

    >>> hdu.writeto('table3.fits')

This shortcut will automatically create a minimal primary HDU with no data and
prepend it to the table HDU to create a valid FITS file. If you require
additional data or header keywords in the primary HDU you may still create a
:class:`PrimaryHDU` object and build up the FITS file manually using an
:class:`HDUList`, as described in the next section.

Creating a File with Multiple Extensions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the previous examples we created files with a single meaningful extension (a
:class:`PrimaryHDU` or :class:`BinTableHDU`). To create a file with multiple
extensions we need to create extension HDUs and append them to an :class:`HDUList`.

First, we create some data for Image extensions::

    >>> import numpy as np
    >>> n = np.ones((3, 3))
    >>> n2 = np.ones((100, 100))
    >>> n3 = np.ones((10, 10, 10))

Note that the data shapes of the different extensions do not need to be the same.
Next, place the data into separate :class:`PrimaryHDU` and :class:`ImageHDU`
objects::

    >>> primary_hdu = fits.PrimaryHDU(n)
    >>> image_hdu = fits.ImageHDU(n2)
    >>> image_hdu2 = fits.ImageHDU(n3)

A multi-extension FITS file is not constrained to be only imaging or table data, we
can mix them. To show this we'll use the example from the previous section to make a
:class:`BinTableHDU`::

    >>> c1 = fits.Column(name='a', array=np.array([1, 2]), format='K')
    >>> c2 = fits.Column(name='b', array=np.array([4, 5]), format='K')
    >>> c3 = fits.Column(name='c', array=np.array([7, 8]), format='K')
    >>> table_hdu = fits.BinTableHDU.from_columns([c1, c2, c3])

Now when we create the :class:`HDUList` we list all extensions we want to
include::

    >>> hdul = fits.HDUList([primary_hdu, image_hdu, table_hdu])

Because :class:`HDUList` acts like a :class:`list` we can also append, for example,
an :class:`ImageHDU` to an already existing :class:`HDUList`::

    >>> hdul.append(image_hdu2)

Multi-extension :class:`HDUList` are treated just like those with only a
:class:`PrimaryHDU`, so to save the file use :func:`HDUList.writeto` as shown above.

.. note::

    The FITS standard enforces all files to have exactly one :class:`PrimaryHDU` that
    is the first HDU present in the file. This standard is enforced during the call to
    :func:`HDUList.writeto` and an error will be raised if it is not met. See the
    ``output_verify`` option in :func:`HDUList.writeto` for ways to fix or ignore these
    warnings.

In the previous example the :class:`PrimaryHDU` contained actual data. In some cases it
is desirable to have a minimal :class:`PrimaryHDU` with only basic header information.
To do this, first create a new :class:`Header` object to encapsulate any keywords you want
to include in the primary HDU, then as before create a :class:`PrimaryHDU`::

    >>> hdr = fits.Header()
    >>> hdr['OBSERVER'] = 'Edwin Hubble'
    >>> hdr['COMMENT'] = "Here's some commentary about this FITS file."
    >>> empty_primary = fits.PrimaryHDU(header=hdr)

When we create a new primary HDU with a custom header as in the above example,
this will automatically include any additional header keywords that are
*required* by the FITS format (keywords such as ``SIMPLE`` and ``NAXIS`` for
example). In general, users should not have to manually manage such keywords,
and should only create and modify observation-specific informational keywords.

We then create an HDUList containing both the primary HDU and any other HDUs want::

    >>> hdul = fits.HDUList([empty_primary, image_hdu2, table_hdu])

.. topic:: Examples:

    See also :ref:`sphx_glr_generated_examples_io_create-mef.py`.

Convenience Functions
---------------------

`astropy.io.fits` also provides several high-level ("convenience") functions.
Such a convenience function is a "canned" operation to achieve one task.
By using these "convenience" functions, a user does not have to worry about
opening or closing a file; all of the housekeeping is done implicitly.

.. warning::

    These functions are useful for interactive Python sessions and less complex
    analysis scripts, but should not be used for application code, as they
    are highly inefficient. For example, each call to :func:`getval`
    requires re-parsing the entire FITS file. Code that makes repeated use
    of these functions should instead open the file with :func:`open`
    and access the data structures directly.

The first of these functions is :func:`getheader`, to get the header of an HDU.
Here are several examples of getting the header. Only the file name is required
for this function. The rest of the arguments are optional and flexible to
specify which HDU the user wants to access::

    >>> from astropy.io.fits import getheader
    >>> hdr = getheader(fits_image_filename)  # get default HDU (=0), i.e. primary HDU's header
    >>> hdr = getheader(fits_image_filename, 0)  # get primary HDU's header
    >>> hdr = getheader(fits_image_filename, 2)  # the second extension
    >>> hdr = getheader(fits_image_filename, 'sci')  # the first HDU with EXTNAME='SCI'
    >>> hdr = getheader(fits_image_filename, 'sci', 2)  # HDU with EXTNAME='SCI' and EXTVER=2
    >>> hdr = getheader(fits_image_filename, ('sci', 2))  # use a tuple to do the same
    >>> hdr = getheader(fits_image_filename, ext=2)  # the second extension
    >>> hdr = getheader(fits_image_filename, extname='sci')  # first HDU with EXTNAME='SCI'
    >>> hdr = getheader(fits_image_filename, extname='sci', extver=2)

Ambiguous specifications will raise an exception::

    >>> getheader(fits_image_filename, ext=('sci', 1), extname='err', extver=2)
    Traceback (most recent call last):
        ...
    TypeError: Redundant/conflicting extension arguments(s): ...

After you get the header, you can access the information in it, such as getting
and modifying a keyword value::

    >>> fits_image_2_filename = fits.util.get_testdata_filepath('o4sp040b0_raw.fits')
    >>> hdr = getheader(fits_image_2_filename, 0)    # get primary hdu's header
    >>> filter = hdr['filter']                       # get the value of the keyword "filter'
    >>> val = hdr[10]                                # get the 11th keyword's value
    >>> hdr['filter'] = 'FW555'                      # change the keyword value

For the header keywords, the header is like a dictionary, as well as a list.
The user can access the keywords either by name or by numeric index, as
explained earlier in this chapter.

If a user only needs to read one keyword, the  :func:`getval` function can
further simplify to just one call, instead of two as shown in the above
examples::

    >>> from astropy.io.fits import getval
    >>> # get 0th extension's keyword FILTER's value
    >>> flt = getval(fits_image_2_filename, 'filter', 0)
    >>> flt
    'Clear'

    >>> # get the 2nd sci extension's 11th keyword's value
    >>> val = getval(fits_image_2_filename, 10, 'sci', 2)
    >>> val
    False

The function :func:`getdata` gets the data of an HDU. Similar to
:func:`getheader`, it only requires the input FITS file name while the
extension is specified through the optional arguments. It does have one extra
optional argument header. If header is set to True, this function will return
both data and header, otherwise only data is returned::

    >>> from astropy.io.fits import getdata
    >>> # get 3rd sci extension's data:
    >>> data = getdata(fits_image_filename, 'sci', 3)
    >>> # get 1st extension's data AND header:
    >>> data, hdr = getdata(fits_image_filename, 1, header=True)

The functions introduced above are for reading. The next few functions
demonstrate convenience functions for writing::

    >>> fits.writeto('out.fits', data, hdr)

The :func:`writeto` function uses the provided data and an optional header to
write to an output FITS file.

::

    >>> fits.append('out.fits', data, hdr)

The :func:`append` function will use the provided data and the optional header
to append to an existing FITS file. If the specified output file does not
exist, it will create one.

.. code:: python

    from astropy.io.fits import update
    update(filename, dat, hdr, 'sci')         # update the 'sci' extension
    update(filename, dat, 3)                  # update the 3rd extension
    update(filename, dat, hdr, 3)             # update the 3rd extension
    update(filename, dat, 'sci', 2)           # update the 2nd SCI extension
    update(filename, dat, 3, header=hdr)      # update the 3rd extension
    update(filename, dat, header=hdr, ext=5)  # update the 5th extension

The :func:`update` function will update the specified extension with the input
data/header. The third argument can be the header associated with the data. If
the third argument is not a header, it (and other positional arguments) are
assumed to be the extension specification(s). Header and extension specs can
also be keyword arguments.

The :func:`printdiff` function will print a difference report of two FITS files,
including headers and data. The first two arguments must be two FITS
filenames or FITS file objects with matching data types (i.e., if using strings
to specify filenames, both inputs must be strings). The third
argument is an optional extension specification, with the same call format
of :func:`getheader` and :func:`getdata`. In addition you can add any keywords
accepted by the :class:`FITSDiff` class.

.. code:: python

    from astropy.io.fits import printdiff
    # get a difference report of ext 2 of inA and inB
    printdiff('inA.fits', 'inB.fits', ext=2)
    # ignore HISTORY and COMMMENT keywords
    printdiff('inA.fits', 'inB.fits', ignore_keywords=('HISTORY','COMMENT')

Finally, the :func:`info` function will print out information of the specified
FITS file::

    >>> fits.info(fits_image_filename)
    Filename: ...test0.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     138   ()
      1  SCI           1 ImageHDU        61   (40, 40)   int16
      2  SCI           2 ImageHDU        61   (40, 40)   int16
      3  SCI           3 ImageHDU        61   (40, 40)   int16
      4  SCI           4 ImageHDU        61   (40, 40)   int16

This is one of the most useful convenience functions for getting an overview of
what a given file contains without looking at any of the details.


Using `astropy.io.fits`
=======================
.. toctree::
   :maxdepth: 2

   usage/headers
   usage/image
   usage/table
   usage/verification
   usage/unfamiliar
   usage/scripts
   usage/misc

Command-Line Utilities
======================

For convenience, several of ``astropy``'s sub-packages install utility programs
on your system which allow common tasks to be performed without having
to open a Python interpreter. These utilities include:

- `~astropy.io.fits.scripts.fitsheader`: prints the headers of a FITS file.

- `~astropy.io.fits.scripts.fitscheck`: verifies and optionally rewrites
  the CHECKSUM and DATASUM keywords of a FITS file.

- :ref:`fitsdiff`: compares two FITS files and reports the differences.

- :ref:`fits2bitmap`: converts FITS images to bitmaps, including scaling and
  stretching.

- :ref:`wcslint <wcslint>`: checks the :ref:`WCS <astropy-wcs>` keywords in a
  FITS file for compliance against the standards.

Other Information
=================

.. toctree::
    :maxdepth: 1

    appendix/faq
    appendix/header_transition
    appendix/history

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. automodule:: astropy.io.fits

.. toctree::
    :maxdepth: 3

    api/files.rst
    api/hdulists.rst
    api/hdus.rst
    api/headers.rst
    api/cards.rst
    api/tables.rst
    api/images.rst
    api/diff.rst
    api/verification.rst

.. rubric:: Footnotes

.. [#f1]  For legacy code only that already depends on PyFITS, it's acceptable to continue using "from astropy.io import fits as pyfits".
