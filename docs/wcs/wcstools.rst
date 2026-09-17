.. _wcstools:

Subsetting and Pixel Scales
^^^^^^^^^^^^^^^^^^^^^^^^^^^

WCS objects can be broken apart into their constituent axes using the
`~astropy.wcs.WCS.sub` function.  There is also a `~astropy.wcs.WCS.celestial`
convenience function that will return a WCS object with only the celestial axes
included.

The pixel scales of a celestial image or the pixel dimensions of a non-celestial
image can be extracted with the utility functions
`~astropy.wcs.utils.proj_plane_pixel_scales` and
`~astropy.wcs.utils.non_celestial_pixel_scales`. Likewise, celestial pixel
area can be extracted with the utility function
`~astropy.wcs.utils.proj_plane_pixel_area`.

Matplotlib plots with correct WCS projection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`WCSAxes <wcsaxes>` framework, previously a standalone package, allows
the :class:`~astropy.wcs.WCS` to be used to define projections in Matplotlib.
More information on using WCSAxes can be found :ref:`here <wcsaxes>`.

.. plot::
    :context: reset
    :include-source:
    :align: center

    import warnings
    from matplotlib import pyplot as plt
    from astropy.io import fits
    from astropy.wcs import WCS, FITSFixedWarning
    from astropy.utils.data import get_pkg_data_filename

    filename = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')

    hdu = fits.open(filename)[0]
    with warnings.catch_warnings():
        # Ignore a warning on using DATE-OBS in place of MJD-OBS
        warnings.filterwarnings('ignore', message="'datfix' made the change",
                                category=FITSFixedWarning)
        wcs = WCS(hdu.header)

    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.imshow(hdu.data, origin='lower', cmap='viridis')
    ax.set(xlabel='RA', ylabel='Dec')


Fitting a WCS from matched pixel and sky coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a fiducial point and a celestial projection type,
`~astropy.wcs.utils.fit_wcs_from_points` constructs a FITS WCS from two
matched lists: detector pixel positions and their corresponding celestial
coordinates. This is accomplished by fitting the WCS parameters (``CRPIX``,
``CD`` matrix, and optionally SIP distortion coefficients) to the provided
matched points.

The fiducial point of the spherical projection can be specified as
a `~astropy.coordinates.SkyCoord`; by default, it is set to the mean of the
input sky coordinates (``proj_point='center'``). The projection type can be
specified either as a three-letter projection code (for example, ``'TAN'``
for the gnomonic projection) or as a WCS object with a defined projection
type. If not provided, the projection defaults to ``'TAN'``.

Pixel coordinates must follow the FITS convention: the center of the
bottom-left pixel is ``(1, 1)``.  Units of the celestial coordinates of the
returned WCS are always degrees.

.. doctest-requires:: scipy

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy.wcs.utils import fit_wcs_from_points
    >>> x, y = np.meshgrid([5.0, 10.0, 15.0], [2.0, 4.0, 6.0])
    >>> x, y = x.ravel(), y.ravel()
    >>> world = SkyCoord(
    ...     (10.0 + x * 0.01) * u.deg,
    ...     (20.0 + y * 0.01) * u.deg,
    ...     frame="icrs",
    ... )
    >>> xy = (x, y)
    >>> wcs = fit_wcs_from_points(xy, world, projection="TAN")
    >>> list(wcs.wcs.ctype)
    WCS Keywords

    Number of WCS axes: 2
    CTYPE : 'RA---TAN' 'DEC--TAN'
    CUNIT : 'deg' 'deg'
    CRVAL : 10.100006366283209 20.0400070233876
    CRPIX : 10.000636628442498 4.0002341129981005
    CD1_1 CD1_2  : 0.00939453813176144 3.807555317347415e-10
    CD2_1 CD2_2  : -3.577033090976085e-10 0.010000002348464697
    NAXIS : 14  5

See :func:`~astropy.wcs.utils.fit_wcs_from_points` for the full argument
list, including ``sip_degree`` and ``projection``.
