.. doctest-skip-all

.. _io-fits-faq:

astropy.io.fits FAQ
-------------------

.. contents::

General Questions
^^^^^^^^^^^^^^^^^

What is PyFITS and how does it relate to Astropy?
"""""""""""""""""""""""""""""""""""""""""""""""""

PyFITS_ is a library written in, and for use with the Python_ programming
language for reading, writing, and manipulating FITS_ formatted files.  It
includes a high-level interface to FITS headers with the ability for high and
low-level manipulation of headers, and it supports reading image and table
data as Numpy_ arrays.  It also supports more obscure and non-standard formats
found in some FITS files.

The `astropy.io.fits` module is identical to PyFITS but with the names changed.
When development began on Astropy it was clear that one of the core
requirements would be a FITS reader.  Rather than starting from scratch,
PyFITS--being the most flexible FITS reader available for Python--was ported
into Astropy.  There are plans to gradually phase out PyFITS as a stand-alone
module and deprecate it in favor of `astropy.io.fits`.  See more about that in
the next question.

Although PyFITS is written mostly in Python, it includes an optional module
written in C that's required to read/write compressed image data.  However,
the rest of PyFITS functions without this extension module.

.. _PyFITS: http://www.stsci.edu/institute/software_hardware/pyfits
.. _Python: http://www.python.org
.. _FITS: http://fits.gsfc.nasa.gov/
.. _Numpy: http://numpy.scipy.org/


What is the development status of PyFITS?
"""""""""""""""""""""""""""""""""""""""""

PyFITS is written and maintained by the Science Software Branch at the `Space
Telescope Science Institute`_, and is licensed by AURA_ under a `3-clause BSD
license`_ (see `LICENSE.txt`_ in the PyFITS source code).

It is now primarily developed as primarily as a component of Astropy
(`astropy.io.fits`) rather than as a stand-alone module.  There are a few
reasons for this: The first is simply to reduce development effort; the
overhead of maintaining both PyFITS *and* `astropy.io.fits` in separate code
bases is non-trivial.  The second is that there are many features of Astropy
(units, tables, etc.) from which the `astropy.io.fits` module can benefit
greatly.  Since PyFITS is already integrated into Astropy, it makes more sense
to continue development there rather than make Astropy a dependency of PyFITS.

PyFITS' current primary developer and active maintainer is `Erik Bray`_, though
patch submissions are welcome from anyone.  PyFITS is now primarily developed
in a Git repository for ease of merging to and from Astropy.  Patches and issue
reports can be posted to the `GitHub project`_ for PyFITS, or for Astropy.
There is also a legacy `Trac site`_ with some older issue reports still open,
but new issues should be submitted via GitHub if possible.  An `SVN mirror`_ of
the repository is still maintained as well.

The current stable release series is 3.3.x.  Each 3.3.x release tries to
contain only bug fixes, and to not introduce any significant behavioral or API
changes (though this isn't guaranteed to be perfect).  Patch releases for older
release series may be released upon request.  Older versions of PyFITS (2.4 and
earlier) are no longer actively supported.

.. _Space Telescope Science Institute: http://www.stsci.edu/
.. _AURA: http://www.aura-astronomy.org/
.. _3-clause BSD license: http://en.wikipedia.org/wiki/BSD_licenses#3-clause_license_.28.22New_BSD_License.22_or_.22Modified_BSD_License.22.29
.. _LICENSE.txt: https://aeon.stsci.edu/ssb/trac/pyfits/browser/trunk/LICENSE.txt
.. _Erik Bray: mailto:embray@stsci.edu
.. _Trac site: https://aeon.stsci.edu/ssb/trac/pyfits/
.. _SVN mirror: https://aeon.stsci.edu/ssb/svn/pyfits/
.. _GitHub project: https://github.com/spacetelescope/PyFITS


Usage Questions
^^^^^^^^^^^^^^^

Something didn't work as I expected.  Did I do something wrong?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Possibly.  But if you followed the documentation and things still did not work
as expected, it is entirely possible that there is a mistake in the
documentation, a bug in the code, or both.  So feel free to report it as a bug.
There are also many, many corner cases in FITS files, with new ones discovered
almost every week.  `astropy.io.fits` is always improving, but does not support
all cases perfectly.  There are some features of the FITS format (scaled data,
for example) that are difficult to support correctly and can sometimes cause
unexpected behavior.

For the most common cases, however, such as reading and updating FITS headers,
images, and tables, `astropy.io.fits`. is very stable and well-tested.  Before
every Astropy/PyFITS release it is ensured that all its tests pass on a variety
of platforms, and those tests cover the majority of use-cases (until new corner
cases are discovered).


Astropy crashed and output a long string of code.  What do I do?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This listing of code is what is knows as a `stack trace`_ (or in Python
parlance a "traceback").  When an unhandled exception occurs in the code,
causing the program to end, this is a way of displaying where the exception
occurred and the path through the code that led to it.

As Astropy is meant to be used as a piece in other software projects, some
exceptions raised by Astropy are by design.  For example, one of the most
common exceptions is a `~.exceptions.KeyError` when an attempt is made to read
the value of a non-existent keyword in a header::

    >>> from astropy.io import fits
    >>> h = fits.Header()
    >>> h['NAXIS']
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/path/to/astropy/io/fits/header.py", line 125, in __getitem__
        return self._cards[self._cardindex(key)].value
      File "/path/to/astropy/io/fits/header.py", line 1535, in _cardindex
        raise KeyError("Keyword %r not found." % keyword)
    KeyError: "Keyword 'NAXIS' not found."

This indicates that something was looking for a keyword called "NAXIS" that
does not exist.  If an error like this occurs in some other software that uses
Astropy, it may indicate a bug in that software, in that it expected to find a
keyword that didn't exist in a file.

Most "expected" exceptions will output a message at the end of the traceback
giving some idea of why the exception occurred and what to do about it.  The
more vague and mysterious the error message in an exception appears, the more
likely that it was caused by a bug in Astropy.  So if you're getting an
exception and you really don't know why or what to do about it, feel free to
report it as a bug.

.. _stack trace: http://en.wikipedia.org/wiki/Stack_trace


Why does opening a file work in CFITSIO, ds9, etc. but not in Astropy?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

As mentioned elsewhere in this FAQ, there are many unusual corner cases when
dealing with FITS files.  It's possible that a file should work, but isn't
support due to a bug.  Sometimes it's even possible for a file to work in an
older version of Astropy or PyFITS, but not a newer version due to a regression
that isn't tested for yet.

Another problem with the FITS format is that, as old as it is, there are many
conventions that appear in files from certain sources that do not meet the FITS
standard.  And yet they are so common-place that it is necessary to support
them in any FITS readers.  CONTINUE cards are one such example.  There are
non-standard conventions supported by Astropy/PyFITS that are not supported by
CFITSIO and possibly vice-versa.  You may have hit one of those cases.

If Astropy is having trouble opening a file, a good way to rule out whether not
the problem is with Astropy is to run the file through the `fitsverify`_
program.  For smaller files you can even use the `online FITS verifier`_.
These use CFITSIO under the hood, and should give a good indication of whether
or not there is something erroneous about the file.  If the file is
malformatted, fitsverify will output errors and warnings.

If fitsverify confirms no problems with a file, and Astropy is still having
trouble opening it (especially if it produces a traceback) then it's possible
there is a bug in Astropy.

.. _fitsverify: http://heasarc.gsfc.nasa.gov/docs/software/ftools/fitsverify/
.. _online FITS verifier: http://fits.gsfc.nasa.gov/fits_verify.html


How do I turn off the warning messages Astropy keeps outputting to my console?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Astropy uses Python's built-in `warnings`_ subsystem for informing about
exceptional conditions in the code that are recoverable, but that the user may
want to be informed of.  One of the most common warnings in `astropy.io.fits`
occurs when updating a header value in such a way that the comment must be
truncated to preserve space::

    Card is too long, comment is truncated.

Any console output generated by Astropy can be assumed to be from the warnings
subsystem.  See Astropy's documentation on the :ref:`python-warnings` for more
information on how to control and quiet warnings.

.. _warnings: http://docs.python.org/library/warnings.html


What convention does Astropy use for indexing, such as of image coordinates?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

All arrays and sequences in Astropy use a zero-based indexing scheme.  For
example, the first keyword in a header is ``header[0]``, not ``header[1]``.
This is in accordance with Python itself, as well as C, on which Python is
based.

This may come as a surprise to veteran FITS users coming from IRAF, where
1-based indexing is typically used, due to its origins in FORTRAN.

Likewise, the top-left pixel in an N x N array is ``data[0,0]``.  The indices
for 2-dimensional arrays are row-major order, in that the first index is the
row number, and the second index is the column number.  Or put in terms of
axes, the first axis is the y-axis, and the second axis is the x-axis.  This is
the opposite of column-major order, which is used by FORTRAN and hence FITS.
For example, the second index refers to the axis specified by NAXIS1 in the
FITS header.

In general, for N-dimensional arrays, row-major orders means that the
right-most axis is the one that varies the fastest while moving over the
array data linearly.  For example, the 3-dimensional array::

    [[[1, 2],
      [3, 4]],
     [[5, 6],
      [7, 8]]]

is represented linearly in row-major order as::

    [1, 2, 3, 4, 5, 6, 7, 8]

Since 2 immediately follows 1, you can see that the right-most (or inner-most)
axis is the one that varies the fastest.

The discrepancy in axis-ordering may take some getting used to, but it is a
necessary evil.  Since most other Python and C software assumes row-major
ordering, trying to enforce column-major ordering in arrays returned by Astropy
is likely to cause more difficulties than it's worth.


How do I open a very large image that won't fit in memory?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In PyFITS, prior to version 3.1, when the data portion of an HDU is accessed,
the data is read into memory in its entirety.  For example::

    >>> hdul = pyfits.open('myimage.fits')
    >>> hdul[0].data
    ...

reads the entire image array from disk into memory.  For very large images or
tables this is clearly undesirable, if not impossible given the available
resources.

However, `astropy.io.fits.open` has an option to access the data portion of an
HDU by memory mapping using `mmap`_.  In both Astropy and newer versions of
PyFITS this is used by *default*.

What this means is that accessing the data as in the example above only reads
portions of the data into memory on demand.  For example, if I request just a
slice of the image, such as ``hdul[0].data[100:200]``, then just rows 100-200
will be read into memory.  This happens transparently, as though the entire
image were already in memory.  This works the same way for tables.  For most
cases this is your best bet for working with large files.

To ensure use of memory mapping, just add the ``memmap=True`` argument to
`fits.open <astropy.io.fits.open>`.  Likewise, using ``memmap=False`` will
force data to be read entirely into memory.


The default can also be controlled through a configuration option called
``USE_MEMMAP``.  Setting this to ``0`` will disable mmap by default.

Unfortunately, memory mapping does not currently work as well with scaled
image data, where BSCALE and BZERO factors need to be applied to the data to
yield physical values.  Currently this requires enough memory to hold the
entire array, though this is an area that will see improvement in the future.

An alternative, which currently only works for image data (that is, non-tables)
is the sections interface.  It is largely replaced by the better support for
mmap, but may still be useful on systems with more limited virtual-memory
space, such as on 32-bit systems.  Support for scaled image data is flakey with
sections too, though that will be fixed.  See the documentation on :ref:`image
sections <data-sections>` for more details on using this interface.

.. _mmap: http://en.wikipedia.org/wiki/Mmap


How can I create a very large FITS file from scratch?
"""""""""""""""""""""""""""""""""""""""""""""""""""""

This is a very common issue, but unfortunately Astropy does not come with any
built-in facilities for creating large files (larger than will fit in memory)
from scratch (though it may in the future).

Normally to create a single image FITS file one would do something like::

    >>> import numpy
    >>> from astropy.io import fits
    >> data = numpy.zeros((40000, 40000), dtype=numpy.float64)
    >> hdu = fits.PrimaryHDU(data=data)
    >> hdu.writeto('large.fits')

However, a 40000 x 40000 array of doubles is nearly twelve gigabytes!  Most
systems won't be able to create that in memory just to write out to disk.  In
order to create such a large file efficiently requires a little extra work,
and a few assumptions.

First, it is helpful to anticipate about how large (as in, how many keywords)
the header will have in it.  FITS headers must be written in 2880 byte
blocks--large enough for 36 keywords per block (including the END keyword in
the final block).  Typical headers have somewhere between 1 and 4 blocks,
though sometimes more.

Since the first thing we write to a FITS file is the header, we want to write
enough header blocks so that there is plenty of padding in which to add new
keywords without having to resize the whole file.  Say you want the header to
use 4 blocks by default.  Then, excluding the END card which Astropy will add
automatically, create the header and pad it out to 36 * 4 cards like so::

    >>> data = numpy.zeros((100, 100), dtype=numpy.float64)
    # This is a stub array that we'll be using the initialize the HDU; its
    # exact size is irrelevant, as long as it has the desired number of
    # dimensions
    >>> hdu = fits.PrimaryHDU(data=data)
    >>> header = hdu.header
    >>> while len(header) < (36 * 4 - 1):
    ...     header.append()  # Adds a blank card to the end

Now adjust the NAXISn keywords to the desired size of the array, and write
*only* the header out to a file.  Using the ``hdu.writeto()`` method will
cause Astropy to "helpfully" reset the NAXISn keywords to match the size of the
dummy array.  That is because it works hard to ensure that only valid FITS
files are written.  Instead, we can write *just* the header to a file using
the `Header.tofile <astropy.io.fits.Header.tofile>` method::

    >>> header['NAXIS1'] = 40000
    >>> header['NAXIS2'] = 40000
    >>> header.tofile('large.fits')

Finally, we need to grow out the end of the file to match the length of the
data (plus the length of the header).  This can be done very efficiently on
most systems by seeking past the end of the file and writing a single byte,
like so::

    >>> with open('large.fits', 'rb+') as fobj:
    ...     # Seek past the length of the header, plus the length of the
    ...     # Data we want to write.
    ...     # The -1 is to account for the final byte taht we are about to
    ...     # write:
    ...     fobj.seek(len(header.tostring()) + (40000 * 40000 * 8) - 1)
    ...     fobj.write('\0')

On modern operating systems this will cause the file (past the header) to be
filled with zeros out to the ~12GB needed to hold a 40000 x 40000 image.  On
filesystems that support sparse file creation (most Linux filesystems, but not
the HFS+ filesystem used by most Macs) this is a very fast, efficient
operation.  On other systems your mileage may vary.

This isn't the only way to build up a large file, but probably one of the
safest.  This method can also be used to create large multi-extension FITS
files, with a little care.

For creating very large tables, this method may also be used.  Though it can be
difficult to determine ahead of time how many rows a table will need.  In
general, use of the `astropy.io.fits` module is currently discouraged for the
creation and manipulation of large tables.  The FITS format itself is not
designed for efficient on-disk or in-memory manipulation of table structures.
For large, heavy-duty table data it might be better too look into using `HDF5`_
through the `PyTables`_ library.  The :ref:`Astropy Table <astropy-table>`
interface can provide an abstraction layer between different on-disk table
formats as well (for example for converting a table between FITS and HDF5).

PyTables makes use of Numpy under the hood, and can be used to write binary
table data to disk in the same format required by FITS.  It is then possible
to serialize your table to the FITS format for distribution.  At some point
this FAQ might provide an example of how to do this.

.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _PyTables: http://www.pytables.org/moin


How do I create a multi-extension FITS file from scratch?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When you open a FITS file with `astropy.io.fits.open`, an
`~astropy.io.fits.HDUList` object is returned, which holds all the HDUs in the
file.  This ``HDUList`` class is a subclass of Python's builtin `list`, and can
be created from scratch and used as such::

    >>> from astropy.io import fits
    >>> new_hdul = fits.HDUList()
    >>> new_hdul.append(fits.ImageHDU())
    >>> new_hdul.append(fits.ImageHDU())
    >>> new_hdul.writeto('test.fits')

Or the HDU instances can be created first (or read from an existing FITS file)
and the HDUList instantiated like so::

    >>> hdu1 = fits.PrimaryHDU()
    >>> hdu2 = fits.ImageHDU()
    >>> new_hdul = fits.HDUList([hdu1, hdu2])
    >>> new_hdul.writeto('test.fits')

That will create a new multi-extension FITS file with two empty IMAGE
extensions (a default PRIMARY HDU is prepended automatically if one was not
provided manually).


.. _fits-scaled-data-faq:

Why is an image containing integer data being converted unexpectedly to floats?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If the header for your image contains non-trivial values for the optional
BSCALE and/or BZERO keywords (that is, BSCALE != 1 and/or BZERO != 0), then
the raw data in the file must be rescaled to its physical values according to
the formula::

    physical_value = BZERO + BSCALE * array_value

As BZERO and BSCALE are floating point values, the resulting value must be a
float as well.  If the original values were 16-bit integers, the resulting
values are single-precision (32-bit) floats.  If the original values were
32-bit integers the resulting values are double-precision (64-bit floats).

This automatic scaling can easily catch you of guard if you're not expecting
it, because it doesn't happen until the data portion of the HDU is accessed
(to allow things like updating the header without rescaling the data).  For
example::

    >>> hdul = fits.open('scaled.fits')
    >>> image = hdul['SCI', 1]
    >>> image.header['BITPIX']
    32
    >>> image.header['BSCALE']
    2.0
    >>> data = image.data  # Read the data into memory
    >>> data.dtype
    dtype('float64')  # Got float64 despite BITPIX = 32 (32-bit int)
    >>> image.header['BITPIX']  # The BITPIX will automatically update too
    -64
    >>> 'BSCALE' in image.header  # And the BSCALE keyword removed
    False

The reason for this is that once a user accesses the data they may also
manipulate it and perform calculations on it.  If the data were forced to
remain as integers, a great deal of precision is lost.  So it is best to err
on the side of not losing data, at the cost of causing some confusion at
first.

If the data must be returned to integers before saving, use the `ImageHDU.scale
<astropy.io.fits.hdu.image.ImageHDU.scale>` method::

    >>> image.scale('int32')
    >>> image.header['BITPIX']
    32

Alternatively, if a file is opened with ``mode='update'`` along with the
``scale_back=True`` argument, the original BSCALE and BZERO scaling will
be automatically re-applied to the data before saving.  Usually this is
not desirable, especially when converting from floating point back to
unsigned integer values.  But this may be useful in cases where the raw
data needs to be modified corresponding to changes in the physical values.

To prevent rescaling from occurring at all (good for updating headers--even if
you don't intend for the code to access the data, it's good to err on the side
of caution here), use the ``do_not_scale_image_data`` argument when opening
the file::

    >>> hdul = fits.open('scaled.fits', do_not_scale_image_data=True)
    >>> image = hdul['SCI', 1]
    >>> image.data.dtype
    dtype('int32')


Why am I losing precision when I assign floating point values in the header?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The FITS standard allows two formats for storing floating-point numbers in a
header value.  The "fixed" format requires the ASCII representation of the
number to be in bytes 11 through 30 of the header card, and to be
right-justified.  This leaves a standard number of characters for any comment
string.

The fixed format is not wide enough to represent the full range of values that
can be stored in a 64-bit float with full precision.  So FITS also supports a
"free" format in which the ASCII representation can be stored anywhere, using
the full 70 bytes of the card (after the keyword).

Currently Astropy/PyFITS only supports writing fixed format (it can read both
formats), so all floating point values assigned to a header are stored in the
fixed format.  There are plans to add support for more flexible formatting.

In the meantime it is possible to add or update cards by manually formatting
the card image from a string, as it should appear in the FITS file::

    >>> c = fits.Card.fromstring('FOO     = 1234567890.123456789')
    >>> h = fits.Header()
    >>> h.append(c)
    >>> h
    FOO     = 1234567890.123456789

As long as you don't assign new values to 'FOO' via ``h['FOO'] = 123``, will
maintain the header value exactly as you formatted it (as long as it is valid
according to the FITS standard).


Why is reading rows out of a FITS table so slow?
""""""""""""""""""""""""""""""""""""""""""""""""

Underlying every table data array returned by `astropy.io.fits` is a Numpy
`~numpy.recarray` which is a Numpy array type specifically for representing
structured array data (i.e. a table).  As with normal image arrays, Astropy
accesses the underlying binary data from the FITS file via mmap (see the
question "`What performance differences are there between astropy.io.fits and
fitsio?`_" for a deeper explanation of this).  The underlying mmap is then
exposed as a `~numpy.recarray` and in general this is a very efficient way to
read the data.

However, for many (if not most) FITS tables it isn't all that simple.  For
many columns there are conversions that have to take place between the actual
data that's "on disk" (in the FITS file) and the data values that are returned
to the user.  For example FITS binary tables represent boolean values
differently from how Numpy expects them to be represented, "Logical" columns
need to be converted on the fly to a format Numpy (and hence the user) can
understand.  This issue also applies to data that is linearly scaled via the
``TSCALn`` and ``TZEROn`` header keywords.

Supporting all of these "FITS-isms" introduces a lot of overhead that might
not be necessary for all tables, but are still common nonetheless.  That's
not to say it can't be faster even while supporting the peculiarities of
FITS--CFITSIO for example supports all the same features but is orders of
magnitude faster.  Astropy could do much better here too, and there are many
known issues causing slowdown.  There are plenty of opportunities for speedups,
and patches are welcome.  In the meantime for high-performance applications
with FITS tables some users might find the ``fitsio`` library more to their
liking.


I'm opening many FITS files in a loop and getting OSError: Too many open files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Say you have some code like:

.. doctest-skip::

    >>> from astropy.io import fits
    >>> for filename in filenames:
    ...     hdul = fits.open(filename)
    ...     for hdu in hdul:
    ...         hdu_data = hdul.data
    ...         # Do some stuff with the data
    ...     hdul.close()
    ...

The details may differ, but the qualitative point is that the data to many
HDUs and/or FITS files are being accessed in a loop.  This may result in
an exception like::

    Traceback (most recent call last):
      File "<stdin>", line 2, in <module>
    OSError: [Errno 24] Too many open files: 'my_data.fits'

As explained in the :ref:`note on working with large files <fits-large-files>`,
because Astropy uses mmap by default to read the data in a FITS file, even if
you correctly close a file with `HDUList.close <astropy.io.fits.HDUList.close>`
a handle is kept open to that file so that the memory-mapped data array can
still be continued to be read transparently.

The way Numpy supports mmap is such that the file mapping is not closed until
the overlying `~numpy.ndarray` object has no references to it and is freed
memory.  However, when looping over a large number of files (or even just HDUs)
rapidly, this may not happen immediately.  Or in some cases if the HDU object
persists, the data array attached to it may persist too.  The easiest
workaround is to *manually* delete the ``.data`` attribute on the HDU object so
that the `~numpy.ndarray` reference is freed and the mmap can be closed:

.. doctest-skip::

    >>> from astropy.io import fits
    >>> for filename in filenames:
    ...     hdul = fits.open(filename)
    ...     for hdu in hdul:
    ...         hdu_data = hdul.data
    ...         # Do some stuff with the data
    ...         # ...
    ...         # Don't need the data anymore; delete all references to it
    ...         # so that it can be garbage collected
    ...         del hdu_data
    ...         del hdu.data
    ...     hdul.close()
    ...

In some extreme cases files are opened and closed fast enough that Python's
garbage collector does not free them (and hence free the file handles) often
enough.  To mitigate this your code can manually force a garbage collection
by calling :func:`gc.collect` at the end of the loop.

In a future release it will be easier to automatically perform this sort of
cleanup when closing FITS files, where needed.


Comparison with Other FITS Readers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

What is the difference between astropy.io.fits and fitsio?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The `astropy.io.fits` module (originally PyFITS) is a "pure Python" FITS
reader in that all the code for parsing the FITS file format is in Python,
though Numpy is used to provide access to the FITS data via the
`~numpy.ndarray` interface.  `astropy.io.fits` currently also accesses the
`CFITSIO <http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html>`_ to support the
FITS Tile Compression convention, but this feature is optional.  It does not
use CFITSIO outside of reading compressed images.

`fitsio <https://github.com/esheldon/fitsio>`_, on the other hand, is a Python
wrapper for the CFITSIO library.  All the heavy lifting of reading the FITS
format is handled by CFITSIO, while ``fitsio`` provides an easier to use
object-oriented API including providing a Numpy interface to FITS files read
from CFITSIO.  Much of it is written in C (to provide the interface between
Python and CFITSIO), and the rest is in Python.  The Python end mostly
provides the documentation and user-level API.

Because ``fitsio`` wraps CFITSIO it inherits most of its strengths and
weaknesses, though it has an added strength of providing an easier to use
API than if one were to use CFITSIO directly.


Why did Astropy adopt PyFITS as its FITS reader instead of fitsio?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When the Astropy project was first started it was clear from the start that
one of its core components should be a submodule for reading and writing FITS
files, as many other components would be likely to depend on this
functionality.  At the time, the ``fitsio`` package was in its infancy (it
goes back to roughly 2011) while PyFITS had already been established going
back to before the year 2000).  It was already a mature package with support
for the vast majority of FITS files found in the wild, including outdated
formats such as "Random Groups" FITS files still used extensively in the
radio astronomy community.

Although many aspects of PyFITS' interface have evolved over the years, much
of it has also remained the same, and is already familiar to astronomers
working with FITS files in Python.  Most of not all existing training
materials were also based around PyFITS.  PyFITS was developed at STScI, which
also put forward significant resources to develop Astropy, with an eye toward
integrating Astropy into STScI's own software stacks.  As most of the Python
software at STScI uses PyFITS it was the only practical choice for making that
transition.

Finally, although CFITSIO (and by extension ``fitsio``) can read any FITS files
that conform to the FITS standard, it does not support all of the non-standard
conventions that have been added to FITS files in the wild.  It does have some
support for some of these conventions (such as CONTINUE cards and, to a limited
extent, HIERARCH cards), it is not easy to add support for other conventions
to a large and complex C codebase.

PyFITS' object-oriented design makes supporting non-standard conventions
somewhat easier in most cases, and as such PyFITS can be more flexible in the
types of FITS files it can read and return *useful* data from.  This includes
better support for files that fail to meet the FITS standard, but still contain
useful data that should still be readable at least well-enough to correct any
violations of the FITS standard.  For example, a common error in non-English-
speaking regions is to insert non-ASCII characters into FITS headers.  This
is not a valid FITS file, but still should be readable in some sense.
Supporting structural errors such as this is more difficult in CFITSIO which
assumes a more rigid structure.


What performance differences are there between astropy.io.fits and fitsio?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

There are two main performance areas to look at: reading/parsing FITS headers
and reading FITS data (image-like arrays as well as tables).

In the area of headers ``fitsio`` is significantly faster in most cases.  This
is due in large part to the (almost) pure C implementation (due to the use of
CFITSIO), but also due to fact that it is more rigid and does not support as
many local conventions and other special cases as `astropy.io.fits` tries to
support in its pure Python implementation.

That said the difference is small, and only likely to be a bottleneck either
when opening files containing thousands of HDUs, or reading the headers out
of thousands of FITS files in succession (in either case the difference is
not even an order of magnitude).

Where data is concerned the situation is a little more complicated, and
requires some understanding of how PyFITS is implemented versus CFITSIO and
``fitsio``.  First it's important to understand how they differ in terms of
memory management.

`astropy.io.fits`/PyFITS uses mmap, by default, to provide access to the raw
binary data in FITS files.  Mmap is a system call (or in most cases these days
a wrapper in your libc for a lower-level system call) which allows user-space
applications to essentially do the same thing your OS is doing when it uses a
pagefile (swap space) for virtual memory:  It allows data in a file on disk to
be paged into physical memory one page (or in practice usually several pages)
at a time on an as-needed basis.  These cached pages of the file are also
accessible from all processes on the system, so multiple processes can read
from the same file with little additional overhead.  In the case of reading
over all the data in the file the performance difference between using mmap
versus reading the entire data into physical memory at once can vary widely
between systems, hardware, and depending on what else is happening on the
system at the moment, but mmap almost always going to be better.

In principle it requires more overhead since accessing each page will result in
a page fault, and the system requires more requests to the disk.  But in
practice the OS will optimize this pretty aggressively, especially for the most
common case of sequential access--also in reality reading the entire thing into
memory is still going to result in a whole lot of page faults too.  For random
access having all the data in physical memory is always going to be best,
though with mmap it's usually going to be pretty good too (one doesn't normally
access all the data in a file in totally random order--usually a few sections
of it will be accessed most frequently, the OS will keep those pages in
physical memory as best it can).  So for the most general case of reading FITS
files (or most large data on disk) this is the best choice, especially for
casual users, and is hence enabled by default.

CFITSIO/``fitsio``, on the other hand, doesn't assume the existence of
technologies like mmap and page caching.  Thus it implements its own LRU cache
of I/O buffers that store sections of FITS files read from disk in memory in
FITS' famous 2880 byte chunk size.  The I/O buffers are used heavily in
particular for keeping the headers in memory.  Though for large data reads (for
example reading an entire image from a file) it *does* bypass the cache and
instead does a read directly from disk into a user-provided memory buffer.

However, even when CFITSIO reads direct from the file, this is still largely
less efficient than using mmap:  Normally when your OS reads a file from disk,
it caches as much of that read as it can in physical memory (in its page cache)
so that subsequent access to those same pages does not require a subsequent
expensive disk read.  This happens when using mmap too, since the data has to
be copied from disk into RAM at some point.  The difference is that when using
mmap to access the data, the program is able to read that data *directly* out
of the OS's page cache (so long as it's only being read).  On the other hand
when reading data from a file into a local buffer such as with fread(), the
data is first read into the page cache (if not already present) and then copied
from the page cache into the local buffer.  So every read performs at least one
additional memory copy per page read (requiring twice as much physical memory,
and possibly lots of paging if the file is large and pages need to dropped from
the cache).

The user API for CFITSIO usually works by having the user allocate a memory
buffer large enough to hold the image/table they want to read (or at least the
section they're interested in).  There are some helper functions for
determining the appropriate amount of space to allocate.  Then you just pass it
a pointer to your buffer and CFITSIO handles all the reading (usually using the
process described above), and copies the results into your user buffer.  For
large reads it reads directly from the file into your buffer.  Though if the
data needs to be scaled it makes a stop in CFITSIO's own buffer first, then
writes the rescaled values out to the user buffer (if rescaling has been
requested).  Regardless, this means that if your program wishes to hold an
entire image in memory at once it will use as much RAM as the size of the
data.  For most applications it's better (and sufficient) to write it work on
smaller sections of the data, but this requires extra complexity.  Using mmap
on the other hand makes managing this complexity simpler and more efficient.

A very simple and informal test demonstrates this difference.  This test was
performed on four simple FITS images (one of which is a cube) of dimensions
256x256, 1024x1024, 4096x4096, and 256x1024x1024.  Each image was generated
before the test and filled with randomized 64-bit floating point values.  A
similar test was performed using both `astropy.io.fits` and ``fitsio``:  A
handle to the FITS file is opened using each library's basic semantics, and
then the entire data array of the files is copied into a temporary array in
memory (for example if we were blitting the image to a video buffer).  For
Astropy the test is written:

.. code:: python

    def read_test_pyfits(filename):
        with fits.open(filename, memmap=True) as hdul:
            data = hdul[0].data
            c = data.copy()

The test was timed in IPython on a Linux system with kernel version 2.6.32, a
6-core Intel Xeon X5650 CPU clocked at 2.67 GHz per core, and 11.6 GB of RAM
using:

.. code:: python

    for filename in filenames:
        print(filename)
        %timeit read_test_pyfits(filename)

where ``filenames`` is just a list of the aforementioned generated sample
files.  The results were::

    256x256.fits
    1000 loops, best of 3: 1.28 ms per loop
    1024x1024.fits
    100 loops, best of 3: 4.24 ms per loop
    4096x4096.fits
    10 loops, best of 3: 60.6 ms per loop
    256x1024x1024.fits
    1 loops, best of 3: 1.15 s per loop

For ``fitsio`` the test was:

.. code:: python

    def read_test_fitsio(filename):
        with fitsio.FITS(filename) as f:
            data = f[0].read()
            c = data.copy()

This was also run in a loop over all the sample files, producing the results::

    256x256.fits
    1000 loops, best of 3: 476 µs per loop
    1024x1024.fits
    100 loops, best of 3: 12.2 ms per loop
    4096x4096.fits
    10 loops, best of 3: 136 ms per loop
    256x1024x1024.fits
    1 loops, best of 3: 3.65 s per loop

It should be made clear that the sample files were rewritten with new random
data between the Astropy test and the fitsio test, so they were not reading
the same data from the OS's page cache.  Fitsio was much faster on the small
(256x256) image because in that case the time is dominated by parsing the
headers.  As already explained this is much faster in CFITSIO.  However, as
the data size goes up and the header parsing no longer dominates the time,
`astropy.io.fits` using mmap is roughly twice as fast.  This discrepancy would
be almost entirely due to it requiring roughly half as many in-memory copies
to read the data, as explained earlier.  That said, more extensive benchmarking
could be very interesting.

This is also not to say that `astropy.io.fits` does better in all cases.  There
are some cases where it is currently blown away by fitsio.  See the subsequent
question.


Why is fitsio so much faster than Astropy at reading tables?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In many cases it isn't--there is either no difference, or it may be a little
faster in Astropy depending on what you're trying to do with the table and
what types of columns or how many columns the table has.  There are some
cases, however, where ``fitsio`` can be radically faster, mostly for reasons
explained above in "`Why is reading rows out of a FITS table so slow?`_"

In principle a table is no different from, say, an array of pixels.  But
instead of pixels each element of the array is some kind of record structure
(for example two floats, a boolean, and a 20 character string field).  Just as
a 64-bit float is an 8 byte record in an array, a row in such a table can be
thought of as a 37 byte (in the case of the previous example) record in a 1-D
array of rows.  So in principle everything that was explained in the answer to
the question "`What performance differences are there between astropy.io.fits
and fitsio?`_" applies just as well to tables as it does to any other array.

However, FITS tables have many additional complexities that sometimes preclude
streaming the data directly from disk, and instead require transformation from
the on-disk FITS format to a format more immediately useful to the user.  A
common example is how FITS represents boolean values in binary tables.
Another, significantly more complicated example, is variable length arrays.

As explained in "`Why is reading rows out of a FITS table so slow?`_",
`astropy.io.fits`/PyFITS does not currently handle some of these cases as
efficiently as it could, in particular in cases where a user only wishes to
read a few rows out of a table.  Fitsio, on the other hand, has a better
interface for copying one row at a time out of a table and performing the
necessary transformations on that row *only*, rather than on the entire column
or columns that the row is taken from.  As such, for many cases ``fitsio`` gets
much better performance and should be preferred for many performance-critical
table operations.

Fitsio also exposes a microlanguage (implemented in CFITSIO) for making
efficient SQL-like queries of tables (single tables only though--no joins or
anything like that).  This format, described in the `CFITSIO documentation
<http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node97.html>`_ can
in some cases perform more efficient selections of rows than might be possible
with Numpy alone, which requires creating an intermediate mask array in order
to perform row selection.
