.. _io_unified_image:

Image Data
==========

Reading and writing image data in the unified I/O interface is supported
though the `~astropy.nddata.CCDData` class using FITS file format:

.. doctest-skip::

    >>> # Read CCD image
    >>> ccd = CCDData.read('image.fits')

.. doctest-skip::

    >>> # Write back CCD image
    >>> ccd.write('new_image.fits')

Note that the unit is stored in the ``BUNIT`` keyword in the header on saving,
and is read from the header if it is present.

Detailed help on the available keyword arguments for reading and writing
can be obtained via the ``help()`` method as follows:

.. doctest-skip::

    >>> CCDData.read.help('fits')  # Get help on the CCDData FITS reader
    >>> CCDData.writer.help('fits')  # Get help on the CCDData FITS writer
