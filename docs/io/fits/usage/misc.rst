.. currentmodule:: astropy.io.fits

Miscellaneous Features
----------------------

This section describes some of the miscellaneous features of :mod:`astropy.io.fits`.

Differs
^^^^^^^

The :mod:`astropy.io.fits.diff` module contains several facilities for
generating and reporting the differences between two FITS files, or two
components of a FITS file.

The :class:`FITSDiff` class can be used to generate and represent the
differences between either two FITS files on disk, or two existing
:class:`HDUList` objects (or some combination thereof).

Likewise, the :class:`HeaderDiff` class can be used to find the differences
just between two :class:`Header` objects.  Other available differs include
:class:`HDUDiff`, :class:`ImageDataDiff`, :class:`TableDataDiff`, and
:class:`RawDataDiff`.

Each of these classes are instantiated with two instances of the objects that
they diff.  The returned diff instance has a number of attributes starting with
``.diff_`` that describe differences between the two objects.

For example the :class:`HeaderDiff` class cam be used to find the differences
between two :class:`Header` objects like so::

    >>> from astropy.io import fits
    >>> header1 = fits.Header([('KEY_A', 1), ('KEY_B', 2)])
    >>> header2 = fits.Header([('KEY_A', 3), ('KEY_C', 4)])
    >>> diff = fits.diff.HeaderDiff(header1, header2)
    >>> diff.identical
    False
    >>> diff.diff_keywords
    (['KEY_B'], ['KEY_C'])
    >>> diff.diff_keyword_values
    defaultdict(<function <lambda> at ...>, {'KEY_A': [(1, 3)]})

See the API documentation for details on the different differ classes.
