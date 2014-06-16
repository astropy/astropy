****************************
Getting Started with Astropy
****************************

Importing Astropy
=================

In order to encourage consistency amongst users in importing and using Astropy
functionality, we have put together the following guidelines.

Since most of the functionality in Astropy resides in sub-packages, importing
astropy as::

    >>> import astropy

is not very useful. Instead, it is best to import the desired sub-package
with the syntax::

    >>> from astropy import subpackage  # doctest: +SKIP

For example, to access the FITS-related functionality, you can import
`astropy.io.fits` with::

    >>> from astropy.io import fits
    >>> hdulist = fits.open('data.fits')  # doctest: +SKIP

In specific cases, we have recommended shortcuts in the documentation for
specific sub-packages, for example::

    >>> from astropy import units as u
    >>> from astropy import coordinates as coord
    >>> coord.SkyCoord(ra=10.68458*u.deg, dec=41.26917*u.deg, frame='icrs')
    <SkyCoord (ICRS): ra=10.68458 deg, dec=41.26917 deg>

Finally, in some cases, most of the required functionality is contained in a
single class (or a few classes). In those cases, the class can be directly
imported::

    >>> from astropy.cosmology import WMAP7
    >>> from astropy.table import Table
    >>> from astropy.wcs import WCS

Note that for clarity, and to avoid any issues, we recommend to **never**
import any Astropy functionality using ``*``, for example::

    >>> from astropy.io.fits import *  # NOT recommended

Some components of Astropy started off as standalone packages (e.g. PyFITS, PyWCS),
so in cases where Astropy needs to be used as a drop-in replacement, the following
syntax is also acceptable::

    >>> from astropy.io import fits as pyfits

Getting started with subpackages
================================

Because different subpackages have very different functionality, further
suggestions for getting started are in the documentation for the subpackages,
which you can reach by browsing the sections listed in the :ref:`user-docs`.

Or, if you want to dive right in, you can either look at docstrings for
particular a package or object, or access their documentation using the
`~astropy.utils.misc.find_api_page` function. For example, doing this::

    >>> from astropy import find_api_page
    >>> from astropy.units import Quantity
    >>> find_api_page(Quantity)  # doctest: +SKIP

Will bring up the documentation for the `~astropy.units.Quantity` class
in your browser.

Command-line utilities
======================

For convenience, several of Astropy's subpackages install utility programs
on your system which allow common tasks to be performed without having
to open a Python interpreter. These utilities include:

- `~astropy.io.fits.scripts.fitsheader`: prints the headers of a FITS file.

- `~astropy.io.fits.scripts.fitscheck`: verifies and optionally re-writes
  the CHECKSUM and DATASUM keywords of a FITS file.

- :ref:`fitsdiff`: compares two FITS files and reports the differences.

- :ref:`samp_hub <vo-samp-example_hub>`: starts a :ref:`SAMP <vo-samp>` hub.

- ``volint``: checks a :ref:`VOTable <astropy-io-votable>`
  file for compliance against the standards.

- :ref:`wcslint <wcslint>`: checks the :ref:`WCS <astropy-wcs>` keywords in a
  FITS file for compliance against the standards.
