I/O mixin
=========

The I/O mixin, `~astropy.nddata.NDIOMixin`, adds read and write methods that
us the astropy I/O registry.

At the moment, the only available format is FITS.

To use the I/O mixin, start by creating a class that includes it::

    >>> from astropy.nddata import NDData, NDIOMixin
    >>> class NDDataIO(NDIOMixin, NDData): pass

Once you create an instance, you can write it::

    >>> ndd = NDDataIO([1, 2, 3])
    >>> ndd.write('my_data.fits')

You can also create an instance by reading::

    >>> ndd2 = NDDataIO.read('my_data.fits')

A couple important notes about I/O to FITS:

+ FITS is not very tolerant of long keywords or values, though by default
  `astropy.io.fits` will attempt to use conventions that work around those
  limitations.
+ If you create an object by reading from a FITS file, as in the case of
  ``ndd2`` above, its ``meta`` attribute will be a FITS header.
