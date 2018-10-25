****************************************************
Shared Python interface for World Coordinate Systems
****************************************************

Background
==========

The :class:`~astropy.wcs.WCS` class implements what is considered the
most common 'standard' for representing world coordinate systems in
FITS files, but it cannot represent arbitrarily complex transformations
and there is no agreement on how to use the standard beyond FITS files.
Therefore, other world coordinate system transformation approaches exist,
such as the `gwcs <http://gwcs.readthedocs.io/>`_ package being developed
for the James Webb Space Telescope (which is also applicable to other data).

Since one of the goals of the Astropy Project is to improve interoperability
between packages, we have collaboratively defined a standardized application
programming interface (API) for world coordinate system objects to be used
in Python. This API is described in the Astropy Proposal for Enhancements (APE) 14:
`A shared Python interface for World Coordinate Systems
<https://doi.org/10.5281/zenodo.1188874>`_.

The core astropy package provides base classes that define the low- and high-
level APIs described in APE 14 in the :mod:`astropy.wcs.wcsapi` module, and
these are listed in the `Reference/API`_ section below.

Overview
========
While the full  details and motivation for the API are detailed in APE 14, 
this documentation summarizes the elements that are implemented directly in the `astropy` core package.  The high-level interface is likely of most interest to the average user.  In particular, the most important methods are the `world_to_pixel` and `pixel_to_world` methods. These provide the essential elements of WCS: mapping to and from world coordinates. The remainder generally provide information about the *kind* of world coordinates or similar information about the structure of the WCS.

In a bit more detail, the key classes implemented here are a high-level that provides the main user interface (`~astropy.wcs.wcsapi.BaseHighLevelWCS` and subclasses), and a lower-level interface (`~astropy.wcs.wcsapi.BaseLowLevelWCS` and subclasses).  These can be distinct objects *or* the same one.  For FITS-WCS, the `~astropy.wcs.WCS` object meant for FITS-WCS follows both interfaces, allowing immediate use of this API with files that already contain FITS-WCS.  More concrete examples are outlined below.
Pixel conventions and definitions
=================

This API assumes that integer pixel values fall at the center of pixels (as
assumed in the FITS-WCS standard, see Section 2.1.4 of `Greisen et al., 2002,
A&A 446, 747 <https://doi.org/10.1051/0004-6361:20053818>`_), while at the same
time matching the Python 0-index philosophy.  That is, the first pixel is
considered pixel ``0``, but pixel coordinates ``(0, 0)`` are the *center* of
that pixel.  Hence the first pixel spans pixel values ``-0.5`` to ``0.5``.

The order of the pixel coordinates (``(x, y)`` vs ``(row, column)``) depends on
the method or property used, and this can normally be determined from the
property or method name. Properties and methods containing ``pixel`` assume
``(x, y)`` ordering, while properties and methods containing ``array`` assume
``(row, column)`` ordering.

Basic usage
===========

Let's take a look at the shared Python interface for WCS by using a spectral
cube (two celestial axes and one spectral axis)::

    >>> from astropy.wcs import WCS
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> filename = get_pkg_data_filename('l1448/l1448_13co.fits')  # doctest: +REMOTE_DATA
    >>> wcs = WCS(filename)  # doctest: +REMOTE_DATA
    >>> wcs  # doctest: +REMOTE_DATA
    WCS Keywords
    Number of WCS axes: 3
    CTYPE : 'RA---SFL'  'DEC--SFL'  'VOPT'
    CRVAL : 57.6599999999  0.0  -9959.44378305
    CRPIX : -799.0  -4741.913  -187.0
    PC1_1 PC1_2 PC1_3  : 1.0  0.0  0.0
    PC2_1 PC2_2 PC2_3  : 0.0  1.0  0.0
    PC3_1 PC3_2 PC3_3  : 0.0  0.0  1.0
    CDELT : -0.006388889  0.006388889  66.42361
    NAXIS : 105  105  53

We can check how many pixel and world axes are in the transformation as well
as the shape of the data the WCS applies to::

    >>> wcs.pixel_n_dim  # doctest: +REMOTE_DATA
    3
    >>> wcs.world_n_dim  # doctest: +REMOTE_DATA
    3
    >>> wcs.array_shape  # doctest: +REMOTE_DATA
    [53, 105, 105]

Let's now check what the physical type of each axis is::

    >>> wcs.world_axis_physical_types  # doctest: +REMOTE_DATA
    ['pos.eq.ra', 'pos.eq.dec', 'spect.dopplerVeloc.opt']

This is indeed a spectral cube, with RA/Dec and a velocity axis.

The main part of the new interface defines standard methods for transforming
coordinates. The most convenient way is to use the high-level methods
:meth:`~astropy.wcs.wcsapi.BaseHighLevelWCS.pixel_to_world` and
:meth:`~astropy.wcs.wcsapi.BaseHighLevelWCS.world_to_pixel`, which can
transform directly to astropy objects::

    >>> celestial, spectral = wcs.pixel_to_world([1, 2], [4, 3], [2, 3])  # doctest: +REMOTE_DATA
    >>> celestial  # doctest: +REMOTE_DATA
    <SkyCoord (ICRS): (ra, dec) in deg
        [(51.73115731, 30.32750025), (51.72414268, 30.32111136)]>
    >>> spectral  # doctest: +REMOTE_DATA
    <Quantity [2661.04211695, 2727.46572695] m / s>

Similarly, we can transform astropy objects back::

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> coord = SkyCoord('03h26m36.4901s +30d45m22.2012s')
    >>> pixels = wcs.world_to_pixel(coord, 3000 * u.m / u.s)  # doctest: +REMOTE_DATA
    >>> pixels  # doctest: +REMOTE_DATA
    [array(8.11341207), array(71.0956641), array(7.10297292)]

If you are looking to index the original data using these pixel coordinates,
be sure to instead use
:meth:`~astropy.wcs.wcsapi.BaseHighLevelWCS.world_to_array_index` which returns
the coordinates in the correct order to index Numpy arrays, and also rounds to
the nearest integer values::

    >>> index = wcs.world_to_array_index(coord, 3000 * u.m / u.s)  # doctest: +REMOTE_DATA
    >>> index  # doctest: +REMOTE_DATA
    (7, 71, 8)
    >>> from astropy.io import fits
    >>> data = fits.getdata(filename)  # doctest: +REMOTE_DATA
    >>> data[index]  # doctest: +REMOTE_DATA +FLOAT_CMP
    0.22262384

If you are interested in converting to/from world values as simple Python scalars
or Numpy arrays without using high-level astropy objects, there are methods
such as :meth:`~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_to_world_values` to
do this - see `Reference/API`_ for more details.

Extending the physical types in FITS-WCS
========================================

As shown above, the :attr:`~astropy.wcs.WCS.world_axis_physical_types` property
returns the list of physical types for each axis. For FITS-WCS, this is
determined from the CTYPE values in the header. In cases where the physical
type is not known, `None` is returned. However, it is possible to override the
physical types returned by using the
:class:`~astropy.wcs.wcsapi.fitswcs.custom_ctype_to_ucd_mapping` context
manager. Consider a WCS with the following CTYPE::

    >>> from astropy.wcs import WCS
    >>> wcs = WCS(naxis=1)
    >>> wcs.wcs.ctype = ['SPAM']
    >>> wcs.world_axis_physical_types
    [None]

We can specify that for this CTYPE, the physical type should be
``'food.spam'``::

    >>> from astropy.wcs.wcsapi.fitswcs import custom_ctype_to_ucd_mapping
    >>> with custom_ctype_to_ucd_mapping({'SPAM': 'food.spam'}):
    ...     wcs.world_axis_physical_types
    ['food.spam']

Reference/API
=============

.. automodapi:: astropy.wcs.wcsapi

.. automodapi:: astropy.wcs.wcsapi.fitswcs
