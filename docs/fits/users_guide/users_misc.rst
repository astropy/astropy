.. currentmodule:: astropy.io.fits

**********************
Miscellaneous Features
**********************

In this chapter, we'll describe some of the miscellaneous features of Astropy.


Warning Messages
================

Astropy uses the Python warnings module to issue warning messages.  The user can
suppress the warnings using the python command line argument ``-W"ignore"``
when starting an interactive python session.  For example:

.. parsed-literal::

     python -W"ignore"

The user may also use the command line argument when running a python script as
follows:

.. parsed-literal::

     python -W"ignore" myscript.py

It is also possible to suppress warnings from within a python script.  For
instance, the warnings issued from a single call to the writeto convenience
function may be suppressed from within a python script as follows:

.. parsed-literal::

     import warnings
     from astropy.io import fits

     # ...

     warnings.resetwarnings()
     warnings.filterwarnings('ignore', category=UserWarning, append=True)
     fits.writeto(file, im, clobber=True)
     warnings.resetwarnings()
     warnings.filterwarnings('always', category=UserWarning, append=True)

     # ...

Astropy also issues warnings when deprecated API features are used.  In Python
2.7 and up deprecation warnings are ignored by default.  To run Python with
deprecation warnings enabled, either start Python with the ``-Wall`` argument,
or you can enable deprecation warnings specifically with ``-Wd``.

In Python versions below 2.7, if you wish to *squelch* deprecation warnings,
you can start Python with ``-Wi::Deprecation``.  This sets all deprecation
warnings to ignored.  See
http://docs.python.org/using/cmdline.html#cmdoption-unittest-discover-W
for more information on the -W argument.

Differs
=======

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
``.diff_`` that describe differences between the two objects.  See the API
documentation for details on the different differ classes.
