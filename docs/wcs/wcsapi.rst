.. _wcsapi:

Shared Python Interface for World Coordinate Systems
****************************************************

Background
^^^^^^^^^^

The :class:`~astropy.wcs.WCS` class implements what is considered the
most common 'standard' for representing world coordinate systems in
FITS files, but it cannot represent arbitrarily complex transformations
and there is no agreement on how to use the standard beyond FITS files.
Therefore, other world coordinate system transformation approaches exist,
such as the `gwcs <https://gwcs.readthedocs.io/>`_ package being developed
for the James Webb Space Telescope (which is also applicable to other data).

Since one of the goals of the Astropy Project is to improve interoperability
between packages, we have collaboratively defined a standardized application
programming interface (API) for world coordinate system objects to be used
in Python. This API is described in the Astropy Proposal for Enhancements (APE) 14:
`A shared Python interface for World Coordinate Systems
<https://doi.org/10.5281/zenodo.1188874>`_.

The core astropy package provides base classes that define the low- and
high-level APIs described in APE 14 in the :mod:`astropy.wcs.wcsapi` module, and
these are listed in the :ref:`wcs-reference-api` section below.

Overview
^^^^^^^^

While the full  details and motivation for the API are detailed in APE 14,  this
documentation summarizes the elements that are implemented directly in the
astropy core package.  The high-level interface is likely of most interest to
the average user.  In particular, the most important methods are the
:meth:`~astropy.wcs.wcsapi.BaseHighLevelWCS.pixel_to_world` and
:meth:`~astropy.wcs.wcsapi.BaseHighLevelWCS.world_to_pixel` methods. These
provide the essential elements of WCS: mapping to and from world coordinates.
The remainder generally provide information about the *kind* of world
coordinates or similar information about the structure of the WCS.

In a bit more detail, the key classes implemented here are a high-level that
provides the main user interface (:class:`~astropy.wcs.wcsapi.BaseHighLevelWCS` and
subclasses), and a lower-level interface (:class:`~astropy.wcs.wcsapi.BaseLowLevelWCS`
and subclasses).  These can be distinct objects *or* the same one.  For
FITS-WCS, the `~astropy.wcs.WCS` object meant for FITS-WCS follows both
interfaces, allowing immediate use of this API with files that already contain
FITS-WCS. More concrete examples are outlined below.

Basic usage
^^^^^^^^^^^

Let's start off by looking at the shared Python interface for WCS by using a
simple image with two celestial axes (Right Ascension and Declination)::

    >>> from astropy.wcs import WCS
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> from astropy.io import fits
    >>> filename = get_pkg_data_filename('galactic_center/gc_2mass_k.fits')  # doctest: +REMOTE_DATA
    >>> hdulist = fits.open(filename)  # doctest: +REMOTE_DATA
    >>> hdu = hdulist[0]  # doctest: +REMOTE_DATA
    >>> wcs = WCS(hdu.header)  # doctest: +REMOTE_DATA
    >>> wcs  # doctest: +REMOTE_DATA
    WCS Keywords
    <BLANKLINE>
    Number of WCS axes: 2
    CTYPE : 'RA---TAN'  'DEC--TAN'
    CRVAL : 266.4  -28.93333
    CRPIX : 361.0  360.5
    NAXIS : 721  720

We can check how many pixel and world axes are in the transformation as well
as the shape of the data the WCS applies to::

    >>> wcs.pixel_n_dim  # doctest: +REMOTE_DATA
    2
    >>> wcs.world_n_dim  # doctest: +REMOTE_DATA
    2
    >>> wcs.array_shape  # doctest: +REMOTE_DATA
    (720, 721)

Note that the array shape should match that of the data::

    >>> hdu.data.shape  # doctest: +REMOTE_DATA
    (720, 721)

As mentioned in :ref:`pixel_conventions`, what would normally be
considered the 'y-axis' of the image (when looking at it visually) is the first
dimension, while the 'x-axis' of the image is the second dimension. Thus
:attr:`~astropy.wcs.WCS.array_shape` returns the shape in the *opposite* order
to the NAXIS keywords in the FITS header (in the case of FITS-WCS). If you are
interested in the data shape in the reverse order (which would match the NAXIS
order in the case of FITS-WCS), then you can use
:attr:`~astropy.wcs.WCS.pixel_shape`::

    >>> wcs.pixel_shape  # doctest: +REMOTE_DATA
    (721, 720)

Let's now check what the physical type of each axis is::

    >>> wcs.world_axis_physical_types  # doctest: +REMOTE_DATA
    ['pos.eq.ra', 'pos.eq.dec']

This is indeed an image with two celestial axes.

The main part of the new interface defines standard methods for transforming
coordinates. The most convenient way is to use the high-level methods
:meth:`~astropy.wcs.wcsapi.BaseHighLevelWCS.pixel_to_world` and
:meth:`~astropy.wcs.wcsapi.BaseHighLevelWCS.world_to_pixel`, which can
transform directly to astropy objects::

    >>> coord = wcs.pixel_to_world([1, 2], [4, 3])  # doctest: +REMOTE_DATA
    >>> coord  # doctest: +REMOTE_DATA
    <SkyCoord (FK5: equinox=2000.0): (ra, dec) in deg
        [(266.97242993, -29.42584415), (266.97084321, -29.42723968)]>

Similarly, we can transform astropy objects back - we can test this by creating
Galactic coordinates and these will automatically be converted::

    >>> from astropy.coordinates import SkyCoord
    >>> coord = SkyCoord('00h00m00s +00d00m00s', frame='galactic')
    >>> pixels = wcs.world_to_pixel(coord)  # doctest: +REMOTE_DATA
    >>> pixels  # doctest: +REMOTE_DATA
    (array(356.85179997), array(357.45340331))

If you are looking to index the original data using these pixel coordinates,
be sure to instead use
:meth:`~astropy.wcs.wcsapi.BaseHighLevelWCS.world_to_array_index` which returns
the coordinates in the correct order to index Numpy arrays, and also rounds to
the nearest integer values::

    >>> index = wcs.world_to_array_index(coord)  # doctest: +REMOTE_DATA
    >>> index  # doctest: +REMOTE_DATA
    (357, 357)
    >>> hdu.data[index]  # doctest: +REMOTE_DATA +FLOAT_CMP
    np.float32(563.7532)
    >>> hdulist.close()  # doctest: +REMOTE_DATA

Advanced usage
^^^^^^^^^^^^^^

Let's now take a look at a WCS for a spectral cube (two celestial axes and one
spectral axis)::

    >>> filename = get_pkg_data_filename('l1448/l1448_13co.fits')  # doctest: +REMOTE_DATA
    >>> hdulist = fits.open(filename)  # doctest: +REMOTE_DATA
    >>> hdu = hdulist[0]  # doctest: +REMOTE_DATA
    >>> wcs = WCS(hdu.header)  # doctest: +REMOTE_DATA
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

As before we can check how many pixel and world axes are in the transformation
as well as the shape of the data the WCS applies to, as well as the physical
types of each axis::

    >>> wcs.pixel_n_dim  # doctest: +REMOTE_DATA
    3
    >>> wcs.world_n_dim  # doctest: +REMOTE_DATA
    3
    >>> wcs.array_shape  # doctest: +REMOTE_DATA
    (53, 105, 105)
    >>> wcs.world_axis_physical_types  # doctest: +REMOTE_DATA
    ['pos.eq.ra', 'pos.eq.dec', 'spect.dopplerVeloc.opt']

This is indeed a spectral cube, with RA/Dec and a velocity axis.

As before, we can convert between pixels and high-level Astropy objects::

    >>> celestial, spectral = wcs.pixel_to_world([1, 2], [4, 3], [2, 3])  # doctest: +REMOTE_DATA
    >>> celestial  # doctest: +REMOTE_DATA
    <SkyCoord (ICRS): (ra, dec) in deg
        [(51.73115731, 30.32750025), (51.72414268, 30.32111136)]>
    >>> spectral  # doctest: +REMOTE_DATA
    <SpectralCoord
       (target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
                    (57.66, 0., 1000.)
                 (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
                    (0., 0., 0.)>)
      [2661.04211695, 2727.46572695] m / s>

and back::

    >>> from astropy import units as u
    >>> coord = SkyCoord('03h26m36.4901s +30d45m22.2012s')
    >>> pixels = wcs.world_to_pixel(coord, 3000 * u.m / u.s)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    >>> pixels  # doctest: +REMOTE_DATA
    (array(8.11341207), array(71.0956641), array(7.10297292))

And as before we can index array values using::

    >>> index = wcs.world_to_array_index(coord, 3000 * u.m / u.s)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    >>> index  # doctest: +REMOTE_DATA
    (7, 71, 8)
    >>> hdu.data[index]  # doctest: +REMOTE_DATA +FLOAT_CMP
    np.float32(0.22262384)
    >>> hdulist.close()  # doctest: +REMOTE_DATA

If you are interested in converting to/from world values as simple Python scalars
or Numpy arrays without using high-level astropy objects, there are methods
such as :meth:`~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_to_world_values` to
do this - see :ref:`wcs-reference-api` section for more details.

Extending the physical types in FITS-WCS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Slicing of WCS objects
^^^^^^^^^^^^^^^^^^^^^^

A common operation when dealing with data with WCS information attached is to
slice the WCS - this can be either to extract the WCS for a sub-region of the
data, preserving the overall number of dimensions (e.g. a cutout from an image)
or it can be reducing the dimensionality of the data and associated WCS (e.g.
extracting a slice from a spectral cube).

The :class:`~astropy.wcs.wcsapi.SlicedLowLevelWCS` class can be used to slice
any WCS object that conforms to the :class:`~astropy.wcs.wcsapi.BaseLowLevelWCS`
API. To demonstrate this, let's start off by reading in a spectral cube file::

    >>> filename = get_pkg_data_filename('l1448/l1448_13co.fits')  # doctest: +REMOTE_DATA
    >>> wcs = WCS(fits.getheader(filename, ext=0))  # doctest: +REMOTE_DATA

The ``wcs`` object is an instance of :class:`~astropy.wcs.WCS` which conforms to the
:class:`~astropy.wcs.wcsapi.BaseLowLevelWCS` API. We can then use the
:class:`~astropy.wcs.wcsapi.SlicedLowLevelWCS` class to slice the cube::

    >>> from astropy.wcs.wcsapi import SlicedLowLevelWCS
    >>> slices = [10, slice(30, 100), slice(30, 100)]  # doctest: +REMOTE_DATA
    >>> subwcs = SlicedLowLevelWCS(wcs, slices=slices)  # doctest: +REMOTE_DATA

The ``slices`` argument takes any combination of slices, integer values, and
ellipsis which would normally slice a Numpy array. In the above case, we are
extracting a spectral slice, and in that slice we are extracting a sub-region
on the sky.

If you are implementing your own WCS class, you could choose to implement
``__getitem__`` and have it internally use
:class:`~astropy.wcs.wcsapi.SlicedLowLevelWCS`. In fact, the
:class:`~astropy.wcs.WCS` class does this - the example above can be written
more succinctly as::

    >>> wcs[10, 30:100, 30:100]  # doctest: +REMOTE_DATA +ELLIPSIS
    <...>
    SlicedFITSWCS Transformation
    <BLANKLINE>
    This transformation has 2 pixel and 2 world dimensions
    <BLANKLINE>
    Array shape (Numpy order): (70, 70)
    <BLANKLINE>
    Pixel Dim  Axis Name  Data size  Bounds
            0  None              70  None
            1  None              70  None
    <BLANKLINE>
    World Dim  Axis Name  Physical Type  Units
            0  None       pos.eq.ra      deg
            1  None       pos.eq.dec     deg
    <BLANKLINE>
    Correlation between pixel and world axes:
    <BLANKLINE>
               Pixel Dim
    World Dim    0    1
            0  yes  yes
            1  yes  yes

This slicing infrastructure is able to deal with slicing of WCS objects which
have correlated axes - in this case, you may end up with a WCS that has a
different number of pixel and world coordinates. For example, if we slice
a spectral cube to extract a 1D dataset corresponding to a row in the
image plane of a spectral slice, the final WCS will have one pixel dimension
and two world dimensions (since both RA/Dec vary over the extracted 1D slice)::

    >>> wcs[10, 40, :]  # doctest: +REMOTE_DATA +ELLIPSIS
    <...>
    SlicedFITSWCS Transformation
    <BLANKLINE>
    This transformation has 1 pixel and 2 world dimensions
    <BLANKLINE>
    Array shape (Numpy order): (105,)
    <BLANKLINE>
    Pixel Dim  Axis Name  Data size  Bounds
            0  None             105  None
    <BLANKLINE>
    World Dim  Axis Name  Physical Type  Units
            0  None       pos.eq.ra      deg
            1  None       pos.eq.dec     deg
    <BLANKLINE>
    Correlation between pixel and world axes:
    <BLANKLINE>
               Pixel Dim
    World Dim    0
            0  yes
            1  yes
