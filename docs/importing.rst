*****************
Importing Astropy
*****************

In order to encourage consistency amongst users in importing and using Astropy
functionality, we have put together the following guidelines.

Since most of the functionality in Astropy resides in sub-packages, importing
astropy as::

    >>> import astropy

is not very useful. Instead, it is best to import the desired sub-pacakge
with the syntax::

    >>> from astropy import subpackage

For example, to access the FITS-related functionality, you can import
`astropy.io.fits` with::

    >>> from astropy.io import fits
    >>> hdulist = fits.open('data.fits')

In specific cases, we have recommended shortcuts in the documentation for
specific sub-packages, for example::

    >>> from astropy import units as u
    >>> from astropy import coordinates as coord
    >>> coord.ICRSCoordinates(ra=10.68458, dec=41.26917, unit=(u.degree, u.degree)))
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg>

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
syntax is acceptable::

    >>> from astropy.io import fits as pyfits
