.. currentmodule:: astropy.io.fits

Image Data
**********

In this chapter, we will discuss the data component in an image HDU.


Image Data as an Array
======================

A FITS primary HDU or an image extension HDU may contain image data. The
following discussions apply to both of these HDU classes. For most cases in
``astropy``, it is a ``numpy`` array, having the shape specified by the NAXIS
keywords and the data type specified by the BITPIX keyword — unless the data is
scaled, in which case see the next section. Here is a quick cross reference
between allowed BITPIX values in FITS images and the ``numpy`` data types:

.. parsed-literal::

    **BITPIX**    **Numpy Data Type**
    8         numpy.uint8 (note it is UNsigned integer)
    16        numpy.int16
    32        numpy.int32
    64        numpy.int64
    -32       numpy.float32
    -64       numpy.float64

To recap, in ``numpy`` the arrays are 0-indexed and the axes are
ordered from slow to fast. So, if a FITS image has NAXIS1=300 and NAXIS2=400,
the ``numpy`` array of its data will have the shape of (400, 300).

Examples
--------

..
  EXAMPLE START
  Image Data as an Array in astropy.io.fits

Here is a summary of reading and updating image data values::

    >>> from astropy.io import fits
    >>> fits_image_filename = fits.util.get_testdata_filepath('test0.fits')

    >>> with fits.open(fits_image_filename) as hdul:  # open a FITS file
    ...     data = hdul[1].data  # assume the first extension is an image
    >>> print(data[1, 4])   # get the pixel value at x=5, y=2
    313
    >>> # get values of the subsection from x=11 to 20, y=31 to 40 (inclusive)
    >>> data[30:40, 10:20]
    array([[314, 314, 313, 312, 313, 313, 313, 313, 313, 312],
           [314, 314, 312, 313, 313, 311, 313, 312, 312, 314],
           [314, 315, 313, 313, 313, 313, 315, 312, 314, 312],
           [314, 313, 313, 314, 311, 313, 313, 313, 313, 313],
           [313, 314, 312, 314, 312, 314, 314, 315, 313, 313],
           [312, 311, 311, 312, 312, 312, 312, 313, 311, 312],
           [314, 314, 314, 314, 312, 313, 314, 314, 314, 311],
           [314, 313, 312, 313, 313, 314, 312, 312, 311, 314],
           [313, 313, 313, 314, 313, 313, 315, 313, 312, 313],
           [314, 313, 313, 314, 313, 312, 312, 314, 310, 314]], dtype=int16)
    >>> data[1,4] = 999  # update a pixel value
    >>> data[30:40, 10:20] = 0  # update values of a subsection
    >>> data[3] = data[2]    # copy the 3rd row to the 4th row

Here are some more complicated examples by using the concept of the "mask
array." The first example is to change all negative pixel values in ``data`` to
zero. The second one is to take logarithm of the pixel values which are
positive::

    >>> data[data < 0] = 0
    >>> import numpy as np
    >>> data[data > 0] = np.log(data[data > 0])

These examples show the concise nature of ``numpy`` array operations.

..
  EXAMPLE END


Scaled Data
===========

Sometimes an image is scaled; that is, the data stored in the file is not the
image's physical (true) values, but linearly transformed according to the
equation:

.. parsed-literal::

    physical value = BSCALE \* (storage value) + BZERO

BSCALE and BZERO are stored as keywords of the same names in the header of the
same HDU. The most common use of a scaled image is to store unsigned 16-bit
integer data because the FITS standard does not allow it. In this case, the
stored data is signed 16-bit integer (BITPIX=16) with BZERO=32768
(:math:`2^{15}`), BSCALE=1.


Reading Scaled Image Data
-------------------------

Images are scaled only when either of the BSCALE/BZERO keywords are present in
the header and either of their values is not the default value (BSCALE=1,
BZERO=0).

For unscaled data, the data attribute of an HDU in ``astropy`` is a ``numpy``
array of the same data type specified by the BITPIX keyword. For a scaled
image, the ``.data`` attribute will be the physical data (i.e., already
transformed from the storage data and may not be the same data type as
prescribed in BITPIX). This means an extra step of copying is needed and thus
the corresponding memory requirement. This also means that the advantage of
memory mapping is reduced for scaled data.

For floating point storage data, the scaled data will have the same data type.
For integer data type, the scaled data will always be single precision floating
point (``numpy.float32``).

Example
^^^^^^^

..
  EXAMPLE START
  Reading Scaled Image Data with astropy.io.fits

Here is an example of what happens to scaled data, before and after the data is
touched::

    >>> fits_scaledimage_filename = fits.util.get_testdata_filepath('scale.fits')

    >>> hdul = fits.open(fits_scaledimage_filename)
    >>> hdu = hdul[0]
    >>> hdu.header['bitpix']
    16
    >>> hdu.header['bzero']
    1500.0
    >>> hdu.data[0, 0]  # once data is touched, it is scaled  #  doctest: +FLOAT_CMP
    557.7563
    >>> hdu.data.dtype.name
    'float32'
    >>> hdu.header['bitpix']  # BITPIX is also updated
    -32
    >>> # BZERO and BSCALE are removed after the scaling
    >>> hdu.header['bzero']
    Traceback (most recent call last):
        ...
    KeyError: "Keyword 'BZERO' not found."

.. warning::

    An important caveat to be aware of when dealing with scaled data in
    ``astropy``, is that when accessing the data via the ``.data`` attribute,
    the data is automatically scaled with the BZERO and BSCALE parameters. If
    the file was opened in "update" mode, it will be saved with the rescaled
    data. This surprising behavior is a compromise to err on the side of not
    losing data: if some floating point calculations were made on the data,
    rescaling it when saving could result in a loss of information.

    To prevent this automatic scaling, open the file with the
    ``do_not_scale_image_data=True`` argument to ``fits.open()``. This is
    especially useful for updating some header values, while ensuring that the
    data is not modified.

    You may also manually reapply scale parameters by using ``hdu.scale()``
    (see below). Alternately, you may open files with the ``scale_back=True``
    argument. This assures that the original scaling is preserved when saving
    even when the physical values are updated. In other words, it reapplies
    the scaling to the new physical values upon saving.

..
  EXAMPLE END


Writing Scaled Image Data
-------------------------

With the extra processing and memory requirement, we discourage the use of
scaled data as much as possible. However, ``astropy`` does provide ways to
write scaled data with the `~ImageHDU.scale` method.

Examples
^^^^^^^^

..
  EXAMPLE START
  Writing Scaled Image Data in astropy.io.fits

To write scaled data with the `~ImageHDU.scale` method::

    >>> # scale the data to Int16 with user specified bscale/bzero
    >>> hdu.scale('int16', bzero=32768)
    >>> # scale the data to Int32 with the min/max of the data range, emits
    >>> # RuntimeWarning: overflow encountered in short_scalars
    >>> hdu.scale('int32', 'minmax')  # doctest: +SKIP
    >>> # scale the data, using the original BSCALE/BZERO, emits
    >>> # RuntimeWarning: invalid value encountered in add
    >>> hdu.scale('int32', 'old')  # doctest: +SKIP
    >>> hdul.close()

The first example above shows how to store an unsigned short integer array.

Caution must be exercised when using the :meth:`~ImageHDU.scale` method.
The :attr:`~ImageHDU.data` attribute of an image HDU, after the
:meth:`~ImageHDU.scale` call, will become the storage values, not the physical
values. So, only call :meth:`~ImageHDU.scale` just before writing out to FITS
files (i.e., calls of :meth:`~HDUList.writeto`, :meth:`~HDUList.flush`, or
:meth:`~HDUList.close`). No further use of the data should be exercised. Here is
an example of what happens to the :attr:`~ImageHDU.data` attribute after the
:meth:`~ImageHDU.scale` call::

    >>> hdu = fits.PrimaryHDU(np.array([0., 1, 2, 3]))
    >>> print(hdu.data)  # doctest: +FLOAT_CMP
    [0. 1. 2. 3.]
    >>> hdu.scale('int16', bzero=32768)
    >>> print(hdu.data)  # now the data has storage values
    [-32768 -32767 -32766 -32765]
    >>> hdu.writeto('new.fits')

..
  EXAMPLE END

.. _data-sections:

Data Sections
=============

When a FITS image HDU's :attr:`~ImageHDU.data` is accessed, either the whole
data is copied into memory (in cases of NOT using memory mapping or if the data
is scaled) or a virtual memory space equivalent to the data size is allocated
(in the case of memory mapping of non-scaled data). If there are several very
large image HDUs being accessed at the same time, the system may run out of
memory.

If a user does not need the entire image(s) at the same time (e.g., processing
the images(s) ten rows at a time), the :attr:`~ImageHDU.section` attribute of an
HDU can be used to alleviate such memory problems.

With ``astropy``'s improved support for memory-mapping, the sections feature is
not as necessary as it used to be for handling very large images. However, if
the image's data is scaled with non-trivial BSCALE/BZERO values, accessing the
data in sections may still be necessary under the current implementation.
Memmap is also insufficient for loading images larger than 2 to 4 GB on a 32-bit
system — in such cases it may be necessary to use sections.

Example
-------

..
  EXAMPLE START
  Data Sections in astropy.io.fits

Here is an example of getting the median image from three input images of the
size 5000x5000.

.. code:: python

    hdul1 = fits.open('file1.fits')
    hdul2 = fits.open('file2.fits')
    hdul3 = fits.open('file3.fits')
    output = np.zeros((5000, 5000))
    for i in range(50):
        j = i * 100
        k = j + 100
        x1 = hdul1[0].section[j:k,:]
        x2 = hdul2[0].section[j:k,:]
        x3 = hdul3[0].section[j:k,:]
        output[j:k, :] = np.median([x1, x2, x3], axis=0)

Data in each :attr:`~ImageHDU.section` does not need to be contiguous for
memory savings to be possible. ``astropy`` will do its best to join together
discontiguous sections of the array while reading as little as possible into
main memory.

Sections cannot currently be assigned. Any modifications made to a data
section are not saved back to the original file.

..
  EXAMPLE END
