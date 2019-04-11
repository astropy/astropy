**************************************
Importing ``astropy`` and Sub-packages
**************************************

In order to encourage consistency among users in importing and using Astropy
functionality, we have put together the following guidelines.

Since most of the functionality in Astropy resides in sub-packages, importing
``astropy`` as::

    >>> import astropy

is not very useful. Instead, it's best to import the desired sub-package
with the syntax::

    >>> from astropy import subpackage  # doctest: +SKIP

For example, to access the FITS-related functionality, you can import
`astropy.io.fits` with::

    >>> from astropy.io import fits
    >>> hdulist = fits.open('data.fits')  # doctest: +SKIP

In specific cases, we have recommended shortcuts in the documentation for
specific sub-packages. For example::

    >>> from astropy import units as u
    >>> from astropy import coordinates as coord
    >>> coord.SkyCoord(ra=10.68458*u.deg, dec=41.26917*u.deg, frame='icrs')  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        ( 10.68458,  41.26917)>

Finally, in some cases, most of the required functionality is contained in a
single class (or a few classes). In those cases, the class can be directly
imported::

    >>> from astropy.cosmology import WMAP7
    >>> from astropy.table import Table
    >>> from astropy.wcs import WCS

Note that for clarity, and to avoid any issues, we recommend **never**
importing any Astropy functionality using ``*``, for example::

    >>> from astropy.io.fits import *  # NOT recommended

Some components of Astropy started off as standalone packages (e.g. PyFITS,
PyWCS), so in cases where Astropy needs to be used as a drop-in replacement,
the following syntax is also acceptable::

    >>> from astropy.io import fits as pyfits

*********************************
Getting Started with Sub-packages
*********************************

Because different sub-packages have very different functionalities, each
sub-package has its own getting started guide. These can be found by browsing
the sections listed in the :ref:`user-docs`.

You can also look at docstrings for a particular package or object, or access
their documentation using the `~astropy.utils.misc.find_api_page` function. For
example, ::

    >>> from astropy import find_api_page
    >>> from astropy.units import Quantity
    >>> find_api_page(Quantity)  # doctest: +SKIP

will bring up the documentation for the `~astropy.units.Quantity` class
in your browser.
