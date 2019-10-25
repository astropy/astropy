# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from astropy import units as u

from .wcs import WCS, WCSSUB_LONGITUDE, WCSSUB_LATITUDE, WCSSUB_CELESTIAL

__doctest_skip__ = ['wcs_to_celestial_frame', 'celestial_frame_to_wcs']

__all__ = ['add_stokes_axis_to_wcs', 'celestial_frame_to_wcs',
           'wcs_to_celestial_frame', 'proj_plane_pixel_scales',
           'proj_plane_pixel_area', 'is_proj_plane_distorted',
           'non_celestial_pixel_scales', 'skycoord_to_pixel',
           'pixel_to_skycoord', 'custom_wcs_to_frame_mappings',
           'custom_frame_to_wcs_mappings', 'fit_wcs_from_points']


def add_stokes_axis_to_wcs(wcs, add_before_ind):
    """
    Add a new Stokes axis that is uncorrelated with any other axes.
    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        The WCS to add to
    add_before_ind : int
        Index of the WCS to insert the new Stokes axis in front of.
        To add at the end, do add_before_ind = wcs.wcs.naxis
        The beginning is at position 0.
    Returns
    -------
    A new `~astropy.wcs.WCS` instance with an additional axis
    """

    inds = [i + 1 for i in range(wcs.wcs.naxis)]
    inds.insert(add_before_ind, 0)
    newwcs = wcs.sub(inds)
    newwcs.wcs.ctype[add_before_ind] = 'STOKES'
    newwcs.wcs.cname[add_before_ind] = 'STOKES'
    return newwcs


def _wcs_to_celestial_frame_builtin(wcs):

    # Import astropy.coordinates here to avoid circular imports
    from astropy.coordinates import FK4, FK4NoETerms, FK5, ICRS, ITRS, Galactic

    # Import astropy.time here otherwise setup.py fails before extensions are compiled
    from astropy.time import Time

    if wcs.wcs.lng == -1 or wcs.wcs.lat == -1:
        return None

    radesys = wcs.wcs.radesys

    if np.isnan(wcs.wcs.equinox):
        equinox = None
    else:
        equinox = wcs.wcs.equinox

    xcoord = wcs.wcs.ctype[wcs.wcs.lng][:4]
    ycoord = wcs.wcs.ctype[wcs.wcs.lat][:4]

    # Apply logic from FITS standard to determine the default radesys
    if radesys == '' and xcoord == 'RA--' and ycoord == 'DEC-':
        if equinox is None:
            radesys = "ICRS"
        elif equinox < 1984.:
            radesys = "FK4"
        else:
            radesys = "FK5"

    if radesys == 'FK4':
        if equinox is not None:
            equinox = Time(equinox, format='byear')
        frame = FK4(equinox=equinox)
    elif radesys == 'FK4-NO-E':
        if equinox is not None:
            equinox = Time(equinox, format='byear')
        frame = FK4NoETerms(equinox=equinox)
    elif radesys == 'FK5':
        if equinox is not None:
            equinox = Time(equinox, format='jyear')
        frame = FK5(equinox=equinox)
    elif radesys == 'ICRS':
        frame = ICRS()
    else:
        if xcoord == 'GLON' and ycoord == 'GLAT':
            frame = Galactic()
        elif xcoord == 'TLON' and ycoord == 'TLAT':
            frame = ITRS(obstime=wcs.wcs.dateobs or None)
        else:
            frame = None

    return frame


def _celestial_frame_to_wcs_builtin(frame, projection='TAN'):

    # Import astropy.coordinates here to avoid circular imports
    from astropy.coordinates import BaseRADecFrame, FK4, FK4NoETerms, FK5, ICRS, ITRS, Galactic

    # Create a 2-dimensional WCS
    wcs = WCS(naxis=2)

    if isinstance(frame, BaseRADecFrame):

        xcoord = 'RA--'
        ycoord = 'DEC-'
        if isinstance(frame, ICRS):
            wcs.wcs.radesys = 'ICRS'
        elif isinstance(frame, FK4NoETerms):
            wcs.wcs.radesys = 'FK4-NO-E'
            wcs.wcs.equinox = frame.equinox.byear
        elif isinstance(frame, FK4):
            wcs.wcs.radesys = 'FK4'
            wcs.wcs.equinox = frame.equinox.byear
        elif isinstance(frame, FK5):
            wcs.wcs.radesys = 'FK5'
            wcs.wcs.equinox = frame.equinox.jyear
        else:
            return None
    elif isinstance(frame, Galactic):
        xcoord = 'GLON'
        ycoord = 'GLAT'
    elif isinstance(frame, ITRS):
        xcoord = 'TLON'
        ycoord = 'TLAT'
        wcs.wcs.radesys = 'ITRS'
        wcs.wcs.dateobs = frame.obstime.utc.isot
    else:
        return None

    wcs.wcs.ctype = [xcoord + '-' + projection, ycoord + '-' + projection]

    return wcs


WCS_FRAME_MAPPINGS = [[_wcs_to_celestial_frame_builtin]]
FRAME_WCS_MAPPINGS = [[_celestial_frame_to_wcs_builtin]]


class custom_wcs_to_frame_mappings:
    def __init__(self, mappings=[]):
        if hasattr(mappings, '__call__'):
            mappings = [mappings]
        WCS_FRAME_MAPPINGS.append(mappings)

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        WCS_FRAME_MAPPINGS.pop()


# Backward-compatibility
custom_frame_mappings = custom_wcs_to_frame_mappings


class custom_frame_to_wcs_mappings:
    def __init__(self, mappings=[]):
        if hasattr(mappings, '__call__'):
            mappings = [mappings]
        FRAME_WCS_MAPPINGS.append(mappings)

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        FRAME_WCS_MAPPINGS.pop()


def wcs_to_celestial_frame(wcs):
    """
    For a given WCS, return the coordinate frame that matches the celestial
    component of the WCS.
    Parameters
    ----------
    wcs : :class:`~astropy.wcs.WCS` instance
        The WCS to find the frame for
    Returns
    -------
    frame : :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame` subclass instance
        An instance of a :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame`
        subclass instance that best matches the specified WCS.
    Notes
    -----
    To extend this function to frames not defined in astropy.coordinates, you
    can write your own function which should take a :class:`~astropy.wcs.WCS`
    instance and should return either an instance of a frame, or `None` if no
    matching frame was found. You can register this function temporarily with::
        >>> from astropy.wcs.utils import wcs_to_celestial_frame, custom_wcs_to_frame_mappings
        >>> with custom_wcs_to_frame_mappings(my_function):
        ...     wcs_to_celestial_frame(...)
    """
    for mapping_set in WCS_FRAME_MAPPINGS:
        for func in mapping_set:
            frame = func(wcs)
            if frame is not None:
                return frame
    raise ValueError("Could not determine celestial frame corresponding to "
                     "the specified WCS object")


def celestial_frame_to_wcs(frame, projection='TAN'):
    """
    For a given coordinate frame, return the corresponding WCS object.
    Note that the returned WCS object has only the elements corresponding to
    coordinate frames set (e.g. ctype, equinox, radesys).
    Parameters
    ----------
    frame : :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame` subclass instance
        An instance of a :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame`
        subclass instance for which to find the WCS
    projection : str
        Projection code to use in ctype, if applicable
    Returns
    -------
    wcs : :class:`~astropy.wcs.WCS` instance
        The corresponding WCS object
    Examples
    --------
    ::
        >>> from astropy.wcs.utils import celestial_frame_to_wcs
        >>> from astropy.coordinates import FK5
        >>> frame = FK5(equinox='J2010')
        >>> wcs = celestial_frame_to_wcs(frame)
        >>> wcs.to_header()
        WCSAXES =                    2 / Number of coordinate axes
        CRPIX1  =                  0.0 / Pixel coordinate of reference point
        CRPIX2  =                  0.0 / Pixel coordinate of reference point
        CDELT1  =                  1.0 / [deg] Coordinate increment at reference point
        CDELT2  =                  1.0 / [deg] Coordinate increment at reference point
        CUNIT1  = 'deg'                / Units of coordinate increment and value
        CUNIT2  = 'deg'                / Units of coordinate increment and value
        CTYPE1  = 'RA---TAN'           / Right ascension, gnomonic projection
        CTYPE2  = 'DEC--TAN'           / Declination, gnomonic projection
        CRVAL1  =                  0.0 / [deg] Coordinate value at reference point
        CRVAL2  =                  0.0 / [deg] Coordinate value at reference point
        LONPOLE =                180.0 / [deg] Native longitude of celestial pole
        LATPOLE =                  0.0 / [deg] Native latitude of celestial pole
        RADESYS = 'FK5'                / Equatorial coordinate system
        EQUINOX =               2010.0 / [yr] Equinox of equatorial coordinates
    Notes
    -----
    To extend this function to frames not defined in astropy.coordinates, you
    can write your own function which should take a
    :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame` subclass
    instance and a projection (given as a string) and should return either a WCS
    instance, or `None` if the WCS could not be determined. You can register
    this function temporarily with::
        >>> from astropy.wcs.utils import celestial_frame_to_wcs, custom_frame_to_wcs_mappings
        >>> with custom_frame_to_wcs_mappings(my_function):
        ...     celestial_frame_to_wcs(...)
    """
    for mapping_set in FRAME_WCS_MAPPINGS:
        for func in mapping_set:
            wcs = func(frame, projection=projection)
            if wcs is not None:
                return wcs
    raise ValueError("Could not determine WCS corresponding to the specified "
                     "coordinate frame.")


def proj_plane_pixel_scales(wcs):
    """
    For a WCS returns pixel scales along each axis of the image pixel at
    the ``CRPIX`` location once it is projected onto the
    "plane of intermediate world coordinates" as defined in
    `Greisen & Calabretta 2002, A&A, 395, 1061 <http://adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_.
    .. note::
        This function is concerned **only** about the transformation
        "image plane"->"projection plane" and **not** about the
        transformation "celestial sphere"->"projection plane"->"image plane".
        Therefore, this function ignores distortions arising due to
        non-linear nature of most projections.
    .. note::
        In order to compute the scales corresponding to celestial axes only,
        make sure that the input `~astropy.wcs.WCS` object contains
        celestial axes only, e.g., by passing in the
        `~astropy.wcs.WCS.celestial` WCS object.
    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        A world coordinate system object.
    Returns
    -------
    scale : `~numpy.ndarray`
        A vector (`~numpy.ndarray`) of projection plane increments
        corresponding to each pixel side (axis). The units of the returned
        results are the same as the units of `~astropy.wcs.Wcsprm.cdelt`,
        `~astropy.wcs.Wcsprm.crval`, and `~astropy.wcs.Wcsprm.cd` for
        the celestial WCS and can be obtained by inquiring the value
        of `~astropy.wcs.Wcsprm.cunit` property of the input
        `~astropy.wcs.WCS` WCS object.
    See Also
    --------
    astropy.wcs.utils.proj_plane_pixel_area
    """
    return np.sqrt((wcs.pixel_scale_matrix**2).sum(axis=0, dtype=float))


def proj_plane_pixel_area(wcs):
    """
    For a **celestial** WCS (see `astropy.wcs.WCS.celestial`) returns pixel
    area of the image pixel at the ``CRPIX`` location once it is projected
    onto the "plane of intermediate world coordinates" as defined in
    `Greisen & Calabretta 2002, A&A, 395, 1061 <http://adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_.
    .. note::
        This function is concerned **only** about the transformation
        "image plane"->"projection plane" and **not** about the
        transformation "celestial sphere"->"projection plane"->"image plane".
        Therefore, this function ignores distortions arising due to
        non-linear nature of most projections.
    .. note::
        In order to compute the area of pixels corresponding to celestial
        axes only, this function uses the `~astropy.wcs.WCS.celestial` WCS
        object of the input ``wcs``.  This is different from the
        `~astropy.wcs.utils.proj_plane_pixel_scales` function
        that computes the scales for the axes of the input WCS itself.
    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        A world coordinate system object.
    Returns
    -------
    area : float
        Area (in the projection plane) of the pixel at ``CRPIX`` location.
        The units of the returned result are the same as the units of
        the `~astropy.wcs.Wcsprm.cdelt`, `~astropy.wcs.Wcsprm.crval`,
        and `~astropy.wcs.Wcsprm.cd` for the celestial WCS and can be
        obtained by inquiring the value of `~astropy.wcs.Wcsprm.cunit`
        property of the `~astropy.wcs.WCS.celestial` WCS object.
    Raises
    ------
    ValueError
        Pixel area is defined only for 2D pixels. Most likely the
        `~astropy.wcs.Wcsprm.cd` matrix of the `~astropy.wcs.WCS.celestial`
        WCS is not a square matrix of second order.
    Notes
    -----
    Depending on the application, square root of the pixel area can be used to
    represent a single pixel scale of an equivalent square pixel
    whose area is equal to the area of a generally non-square pixel.
    See Also
    --------
    astropy.wcs.utils.proj_plane_pixel_scales
    """
    psm = wcs.celestial.pixel_scale_matrix
    if psm.shape != (2, 2):
        raise ValueError("Pixel area is defined only for 2D pixels.")
    return np.abs(np.linalg.det(psm))


def is_proj_plane_distorted(wcs, maxerr=1.0e-5):
    r"""
    For a WCS returns `False` if square image (detector) pixels stay square
    when projected onto the "plane of intermediate world coordinates"
    as defined in
    `Greisen & Calabretta 2002, A&A, 395, 1061 <http://adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_.
    It will return `True` if transformation from image (detector) coordinates
    to the focal plane coordinates is non-orthogonal or if WCS contains
    non-linear (e.g., SIP) distortions.
    .. note::
        Since this function is concerned **only** about the transformation
        "image plane"->"focal plane" and **not** about the transformation
        "celestial sphere"->"focal plane"->"image plane",
        this function ignores distortions arising due to non-linear nature
        of most projections.
    Let's denote by *C* either the original or the reconstructed
    (from ``PC`` and ``CDELT``) CD matrix. `is_proj_plane_distorted`
    verifies that the transformation from image (detector) coordinates
    to the focal plane coordinates is orthogonal using the following
    check:
    .. math::
        \left \| \frac{C \cdot C^{\mathrm{T}}}
        {| det(C)|} - I \right \|_{\mathrm{max}} < \epsilon .
    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        World coordinate system object
    maxerr : float, optional
        Accuracy to which the CD matrix, **normalized** such
        that :math:`|det(CD)|=1`, should be close to being an
        orthogonal matrix as described in the above equation
        (see :math:`\epsilon`).
    Returns
    -------
    distorted : bool
        Returns `True` if focal (projection) plane is distorted and `False`
        otherwise.
    """
    cwcs = wcs.celestial
    return (not _is_cd_orthogonal(cwcs.pixel_scale_matrix, maxerr) or
            _has_distortion(cwcs))


def _is_cd_orthogonal(cd, maxerr):
    shape = cd.shape
    if not (len(shape) == 2 and shape[0] == shape[1]):
        raise ValueError("CD (or PC) matrix must be a 2D square matrix.")

    pixarea = np.abs(np.linalg.det(cd))
    if (pixarea == 0.0):
        raise ValueError("CD (or PC) matrix is singular.")

    # NOTE: Technically, below we should use np.dot(cd, np.conjugate(cd.T))
    # However, I am not aware of complex CD/PC matrices...
    I = np.dot(cd, cd.T) / pixarea
    cd_unitary_err = np.amax(np.abs(I - np.eye(shape[0])))

    return (cd_unitary_err < maxerr)


def non_celestial_pixel_scales(inwcs):
    """
    Calculate the pixel scale along each axis of a non-celestial WCS,
    for example one with mixed spectral and spatial axes.
    Parameters
    ----------
    inwcs : `~astropy.wcs.WCS`
        The world coordinate system object.
    Returns
    -------
    scale : `numpy.ndarray`
        The pixel scale along each axis.
    """

    if inwcs.is_celestial:
        raise ValueError("WCS is celestial, use celestial_pixel_scales instead")

    pccd = inwcs.pixel_scale_matrix

    if np.allclose(np.extract(1-np.eye(*pccd.shape), pccd), 0):
        return np.abs(np.diagonal(pccd))*u.deg
    else:
        raise ValueError("WCS is rotated, cannot determine consistent pixel scales")


def _has_distortion(wcs):
    """
    `True` if contains any SIP or image distortion components.
    """
    return any(getattr(wcs, dist_attr) is not None
               for dist_attr in ['cpdis1', 'cpdis2', 'det2im1', 'det2im2', 'sip'])


# TODO: in future, we should think about how the following two functions can be
# integrated better into the WCS class.

def skycoord_to_pixel(coords, wcs, origin=0, mode='all'):
    """
    Convert a set of SkyCoord coordinates into pixels.
    Parameters
    ----------
    coords : `~astropy.coordinates.SkyCoord`
        The coordinates to convert.
    wcs : `~astropy.wcs.WCS`
        The WCS transformation to use.
    origin : int
        Whether to return 0 or 1-based pixel coordinates.
    mode : 'all' or 'wcs'
        Whether to do the transformation including distortions (``'all'``) or
        only including only the core WCS transformation (``'wcs'``).
    Returns
    -------
    xp, yp : `numpy.ndarray`
        The pixel coordinates
    See Also
    --------
    astropy.coordinates.SkyCoord.from_pixel
    """

    if _has_distortion(wcs) and wcs.naxis != 2:
        raise ValueError("Can only handle WCS with distortions for 2-dimensional WCS")

    # Keep only the celestial part of the axes, also re-orders lon/lat
    wcs = wcs.sub([WCSSUB_LONGITUDE, WCSSUB_LATITUDE])

    if wcs.naxis != 2:
        raise ValueError("WCS should contain celestial component")

    # Check which frame the WCS uses
    frame = wcs_to_celestial_frame(wcs)

    # Check what unit the WCS needs
    xw_unit = u.Unit(wcs.wcs.cunit[0])
    yw_unit = u.Unit(wcs.wcs.cunit[1])

    # Convert positions to frame
    coords = coords.transform_to(frame)

    # Extract longitude and latitude. We first try and use lon/lat directly,
    # but if the representation is not spherical or unit spherical this will
    # fail. We should then force the use of the unit spherical
    # representation. We don't do that directly to make sure that we preserve
    # custom lon/lat representations if available.
    try:
        lon = coords.data.lon.to(xw_unit)
        lat = coords.data.lat.to(yw_unit)
    except AttributeError:
        lon = coords.spherical.lon.to(xw_unit)
        lat = coords.spherical.lat.to(yw_unit)

    # Convert to pixel coordinates
    if mode == 'all':
        xp, yp = wcs.all_world2pix(lon.value, lat.value, origin)
    elif mode == 'wcs':
        xp, yp = wcs.wcs_world2pix(lon.value, lat.value, origin)
    else:
        raise ValueError("mode should be either 'all' or 'wcs'")

    return xp, yp


def pixel_to_skycoord(xp, yp, wcs, origin=0, mode='all', cls=None):
    """
    Convert a set of pixel coordinates into a `~astropy.coordinates.SkyCoord`
    coordinate.
    Parameters
    ----------
    xp, yp : float or `numpy.ndarray`
        The coordinates to convert.
    wcs : `~astropy.wcs.WCS`
        The WCS transformation to use.
    origin : int
        Whether to return 0 or 1-based pixel coordinates.
    mode : 'all' or 'wcs'
        Whether to do the transformation including distortions (``'all'``) or
        only including only the core WCS transformation (``'wcs'``).
    cls : class or None
        The class of object to create.  Should be a
        `~astropy.coordinates.SkyCoord` subclass.  If None, defaults to
        `~astropy.coordinates.SkyCoord`.
    Returns
    -------
    coords : Whatever ``cls`` is (a subclass of `~astropy.coordinates.SkyCoord`)
        The celestial coordinates
    See Also
    --------
    astropy.coordinates.SkyCoord.from_pixel
    """

    # Import astropy.coordinates here to avoid circular imports
    from astropy.coordinates import SkyCoord, UnitSphericalRepresentation

    # we have to do this instead of actually setting the default to SkyCoord
    # because importing SkyCoord at the module-level leads to circular
    # dependencies.
    if cls is None:
        cls = SkyCoord

    if _has_distortion(wcs) and wcs.naxis != 2:
        raise ValueError("Can only handle WCS with distortions for 2-dimensional WCS")

    # Keep only the celestial part of the axes, also re-orders lon/lat
    wcs = wcs.sub([WCSSUB_LONGITUDE, WCSSUB_LATITUDE])

    if wcs.naxis != 2:
        raise ValueError("WCS should contain celestial component")

    # Check which frame the WCS uses
    frame = wcs_to_celestial_frame(wcs)

    # Check what unit the WCS gives
    lon_unit = u.Unit(wcs.wcs.cunit[0])
    lat_unit = u.Unit(wcs.wcs.cunit[1])

    # Convert pixel coordinates to celestial coordinates
    if mode == 'all':
        lon, lat = wcs.all_pix2world(xp, yp, origin)
    elif mode == 'wcs':
        lon, lat = wcs.wcs_pix2world(xp, yp, origin)
    else:
        raise ValueError("mode should be either 'all' or 'wcs'")

    # Add units to longitude/latitude
    lon = lon * lon_unit
    lat = lat * lat_unit

    # Create a SkyCoord-like object
    data = UnitSphericalRepresentation(lon=lon, lat=lat)
    coords = cls(frame.realize_frame(data))

    return coords


def _linear_wcs_fit(params, lon, lat, x, y, w_obj):
    """Objective function for fitting linear terms.
    Parameters
    ----------
    params : array
        6 element array. First 4 elements are PC matrix, last 2 are CRPIX.
    lon, lat: array
        Sky coordinates.
    x, y: array
        Pixel coordinates
    w_obj: `~astropy.wcs.WCS`
        WCS object
        """
    cd = params[0:4]
    crpix = params[4:6]

    w_obj.wcs.cd = ((cd[0], cd[1]), (cd[2], cd[3]))
    w_obj.wcs.crpix = crpix
    lon2, lat2 = w_obj.wcs_pix2world(x, y, 1)
    resids = np.concatenate((lon-lon2, lat-lat2))

    return resids
    

def _sip_fit(params, lon, lat, u, v, w_obj, order, coeff_names):

    """ Objective function for fitting SIP.
     Parameters
    -----------
    params : array
        6 element array. First 4 elements are PC matrix, last 2 are CRPIX.
    lon, lat: array
        Sky coordinates.
    u, v: array
        Pixel coordinates
    w_obj: `~astropy.wcs.WCS`
        WCS object
    """

    from ..modeling.models import SIP, InverseSIP   # here to avoid circular import

    # unpack params
    crpix = params[0:2]
    cdx = params[2:6].reshape((2, 2))
    a_params = params[6:6+len(coeff_names)]
    b_params = params[6+len(coeff_names):]

    # assign to wcs, used for transfomations in this function
    w_obj.wcs.cd = cdx
    w_obj.wcs.crpix = crpix

    a_coeff, b_coeff = {}, {}
    for i in range(len(coeff_names)):
        a_coeff['A_' + coeff_names[i]] = a_params[i]
        b_coeff['B_' + coeff_names[i]] = b_params[i]

    sip = SIP(crpix=crpix, a_order=order, b_order=order,
              a_coeff=a_coeff, b_coeff=b_coeff)
    fuv, guv = sip(u, v)

    xo, yo = np.dot(cdx, np.array([u+fuv-crpix[0], v+guv-crpix[1]]))

    # use all pix2world in case `projection` contains distortion table
    x, y = w_obj.all_world2pix(lon, lat, 1) 
    x, y = np.dot(w_obj.wcs.cd, (x-w_obj.wcs.crpix[0], y-w_obj.wcs.crpix[1]))

    resids = np.concatenate((x-xo, y-yo))

    return resids


def fit_wcs_from_points(xy, world_coords, proj_point='center', 
                        projection='TAN', sip_distortion=True, degree=4):

    """ 
    Given two matching sets of coordinates on detector and sky,
    compute the WCS. Optionally, a SIP can be fit to account for geometric
    distortion. Returns an `~astropy.wcs.WCS` object with the best fit
    parameters for mapping between input pixel and sky coordinates.

    The projection type (default 'TAN') can passed in as a string, one of 
    the valid three-letter projection codes - or as a WCS object with
    projection keywords already set. Note that if an input WCS has any 
    non-polynomial distortion, this will be applied and reflected in the 
    fit terms and coefficients. Passing in a WCS object in this way essentially 
    allows it to be refit based on the matched input coordinates and projection 
    point, but take care when using this option as non-projection related 
    keywords in the input might cause unexpected behavior. 

    Notes
    ------
    - The fiducial point for the spherical projection can be set to 'center'
      to use the mean position of input sky coordinates, or as an
      `~astropy.coordinates.SkyCoord` object.
    - All output will be in degrees.
    - If the coordinate frame differs between `~astropy.coordinates.SkyCoord`
      objects passed in for `world_coords` and `proj_point`, the frame for
      `world_coords`  will override.
    - If `sip_distortion` is False, `degree` will be ignored. 

    Parameters
    ----------
    xy : tuple of two `numpy.ndarray`
        x & y pixel coordinates.
    world_coords : `~astropy.coordinates.SkyCoord`
        Skycoord object with world coordinates.
    proj_point : 'center' or ~astropy.coordinates.SkyCoord`
        Defaults to 'center', in which the geometric center of input world
        coordinates will be used as the projection point. To specify an exact
        point for the projection, a Skycoord object with a coordinate pair can
        be passed in. For consistency, the units and frame of these coordinates
        will be transformed to match 'world_coords' if they don't.
    projection : str or `~astropy.wcs.WCS`
        Three letter projection code, of any of standard projections defined 
        in the FITS WCS standard. Optionally, a WCS object with projection 
        keywords set may be passed in. 
    sip_distortion : bool
        If True, will fit SIP of degree `degree` to points to account for 
        geometric distortion. If False, only linear terms are fit. Defaults to 
        True.  
    degree : int
        Degree of polynomial to fit. Only used if `sip_distortion` is True. 
        Defaults to 4.
    """

    from scipy.optimize import least_squares
    from .wcs import Sip
    import copy

    xp, yp = xy
    lon, lat = world_coords.data.lon.deg, world_coords.data.lat.deg

    # verify input
    if (proj_point != 'center') and (type(proj_point) != type(world_coords)):
        raise ValueError("proj_point must be set to 'center', or an" +
                         "`~astropy.coordinates.SkyCoord` object with " +
                         "a pair of points.")
    if proj_point != 'center':
        assert proj_point.size == 1

    proj_codes = [
        'AZP', 'SZP', 'TAN', 'STG', 'SIN', 'ARC', 'ZEA', 'AIR', 'CYP',
        'CEA', 'CAR', 'MER', 'SFL', 'PAR', 'MOL', 'AIT', 'COP', 'COE',
        'COD', 'COO', 'BON', 'PCO', 'TSC', 'CSC', 'QSC', 'HPX', 'XPH'
    ]
    if type(projection) == str:
        if projection not in proj_codes:
            raise ValueError("Must specify valid projection code from list of " 
                             + "supported types: ", ', '.join(proj_codes))
        # empty wcs to fill in with fit values
        wcs = celestial_frame_to_wcs(frame=world_coords.frame,
                                     projection=projection)
    else: #if projection is not string, should be wcs object. use as template.
        wcs = copy.deepcopy(projection)
        wcs.cdelt = (1., 1.) # make sure cdelt is 1
        wcs.sip = None

    # Change PC to CD, since cdelt will be set to 1
    if wcs.wcs.has_pc():
        wcs.wcs.cd = wcs.wcs.pc
        wcs.wcs.__delattr__('pc')

    if type(sip_distortion) != bool:
       raise ValueError("sip_distortion must be set to True or False.")
    if (sip_distortion is True) & (type(degree) != int):
        raise ValueError("If sip_distorion is True, an integer value for " +
                         "the polynomial order must be provided.")

    # set pixel_shape to span of input points
    wcs.pixel_shape = (max(xp)-min(xp), max(yp)-min(yp))

    # determine CRVAL from input
    close = lambda l, p: p[np.where(np.abs(l) == min(np.abs(l)))[0][0]]
    if str(proj_point) == 'center':  # use center of input points
        wcs.wcs.crval = ((max(lon)+min(lon))/2., (max(lat)+min(lat))/2.)
        wcs.wcs.crpix = ((max(xp)+min(xp))/2., (max(yp)+min(yp))/2.)
    elif proj_point is not None:  # convert units, initial guess for crpix
        proj_point.transform_to(world_coords)
        wcs.wcs.crval = (proj_point.data.lon.deg, proj_point.data.lat.deg)
        wcs.wcs.crpix = (close(lon-wcs.wcs.crval[0], xp),
                         close(lon-wcs.wcs.crval[1], yp))

    # fit linear terms, assign to wcs
    # use (1, 0, 0, 1) as initial guess, in case input wcs was passed in
    # and cd terms are way off. 
    p0 = np.concatenate([[1., 0., 0., 1], wcs.wcs.crpix.flatten()])
    fit = least_squares(_linear_wcs_fit, p0,
                        args=(lon, lat, xp, yp, wcs))
    wcs.wcs.crpix = np.array(fit.x[4:6])
    wcs.wcs.cd = np.array(fit.x[0:4].reshape((2, 2)))

    # fit SIP, if specified. Only fit forward coefficients
    if sip_distortion:
        if '-SIP' not in wcs.wcs.ctype[0]:
            wcs.wcs.ctype = [x + '-SIP' for x in wcs.wcs.ctype]

        coef_names = ['{0}_{1}'.format(i, j) for i in range(degree+1)
                      for j in range(degree+1) if (i+j) < (degree+1) and
                      (i+j) > 1]
        p0 = np.concatenate((np.array(wcs.wcs.crpix), wcs.wcs.cd.flatten(),
                             np.zeros(2*len(coef_names))))

        fit = least_squares(_sip_fit, p0,
                            args=(lon, lat, xp, yp, wcs, degree, coef_names))
        coef_fit = (list(fit.x[6:6+len(coef_names)]),
                    list(fit.x[6+len(coef_names):]))

        # put fit values in wcs
        wcs.wcs.cd = fit.x[2:6].reshape((2, 2))
        wcs.wcs.crpix = fit.x[0:2]

        a_vals = np.zeros((degree+1, degree+1))
        b_vals = np.zeros((degree+1, degree+1))

        for coef_name in coef_names:
            a_vals[int(coef_name[0])][int(coef_name[2])] = coef_fit[0].pop(0)
            b_vals[int(coef_name[0])][int(coef_name[2])] = coef_fit[1].pop(0)

        wcs.sip = Sip(a_vals, b_vals, np.zeros((degree+1, degree+1)), 
                      np.zeros((degree+1, degree+1)), wcs.wcs.crpix) 

    return wcs