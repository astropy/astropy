import numpy as np

from .. import units as u

from .wcs import WCS, WCSSUB_LONGITUDE, WCSSUB_LATITUDE

__doctest_skip__ = ['wcs_to_celestial_frame', 'celestial_frame_to_wcs']

__all__ = ['add_stokes_axis_to_wcs', 'celestial_frame_to_wcs',
           'wcs_to_celestial_frame', 'proj_plane_pixel_scales',
           'proj_plane_pixel_area', 'is_proj_plane_distorted',
           'non_celestial_pixel_scales', 'skycoord_to_pixel',
           'pixel_to_skycoord', 'custom_wcs_to_frame_mappings',
           'custom_frame_to_wcs_mappings','pixel_sky_to_wcs']

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
    from ..coordinates import FK4, FK4NoETerms, FK5, ICRS, Galactic

    # Import astropy.time here otherwise setup.py fails before extensions are compiled
    from ..time import Time

    # Keep only the celestial part of the axes
    wcs = wcs.sub([WCSSUB_LONGITUDE, WCSSUB_LATITUDE])

    if wcs.wcs.lng == -1 or wcs.wcs.lat == -1:
        return None

    radesys = wcs.wcs.radesys

    if np.isnan(wcs.wcs.equinox):
        equinox = None
    else:
        equinox = wcs.wcs.equinox

    xcoord = wcs.wcs.ctype[0][:4]
    ycoord = wcs.wcs.ctype[1][:4]

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
        else:
            frame = None

    return frame


def _celestial_frame_to_wcs_builtin(frame, projection='TAN'):

    # Import astropy.coordinates here to avoid circular imports
    from ..coordinates import BaseRADecFrame, FK4, FK4NoETerms, FK5, ICRS, Galactic

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
    from ..coordinates import SkyCoord, UnitSphericalRepresentation

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

def _linear_transformation_fit(params, ra, dec, x, y, w_obj, proj):
	"""
	Objective function for fitting linear terms
	"""	   
	pc = params[0:4]
	crpix = params[4:6]
	
	w_obj.wcs.pc = ((pc[0],pc[1]),(pc[2],pc[3]))
	w_obj.wcs.crpix = crpix
	lon, lat = w_obj.wcs_pix2world(x,y,1)
	resids = np.concatenate((ra-lon, dec-lat)) 
	return resids
	

def _sip_fit(params, u, v, x, y, coeff_names):
	"""
	Objective function for fitting SIP coefficients
	"""

	crpix = params[0:2]
	cdx = params[2:6].reshape((2,2))
	a_params = params[6:6+len(coeff_names)]
	b_params = params[6+len(coeff_names):]
		
	a_coeff = {}
	b_coeff = {}
	for i in range(len(coeff_names)):
		a_coeff['A_' + coeff_names[i]] = a_params[i]
		b_coeff['B_' + coeff_names[i]] = b_params[i]
		
	sip = SIP(crpix=crpix, a_order=4, b_order=4, a_coeff=a_coeff, b_coeff=b_coeff)
	
	fuv, guv = sip(u,v)
	xo, yo = np.dot(cdx, np.array([u+fuv-crpix[0], v+guv-crpix[1]]))
	resids = np.concatenate((x-xo, y-yo)) 
	return resids	

def _new_sip_fit(params, lon, lat, u, v, w_obj, coeff_names):	
		from ..modeling.models import SIP #here instead of top to avoid circular import
		crpix = params[0:2]
		cdx = params[2:6].reshape((2,2))
		a_params, b_params = params[6:6+len(coeff_names)],params[6+len(coeff_names):]
		w_obj.wcs.pc = cdx
		w_obj.wcs.crpix = crpix
		x,y = w_obj.wcs_world2pix(lon,lat,1) #'intermediate world coordinates', x & y 
		x,y = np.dot(w_obj.wcs.pc,(x - w_obj.wcs.crpix[0],y - w_obj.wcs.crpix[1]))
		
		a_params, b_params = params[6:6+len(coeff_names)], params[6+len(coeff_names):]
		a_coeff, b_coeff = {}, {}

		for i in range(len(coeff_names)):
			a_coeff['A_' + coeff_names[i]] = a_params[i]
			b_coeff['B_' + coeff_names[i]] = b_params[i]
			
		sip = SIP(crpix=crpix, a_order=4, b_order=4, a_coeff=a_coeff, b_coeff=b_coeff)
	
		fuv, guv = sip(u,v)
		xo, yo = np.dot(cdx, np.array([u+fuv-crpix[0], v+guv-crpix[1]]))
		resids = np.concatenate((x-xo, y-yo))
		return resids	
	
def pixel_sky_to_wcs(xp, yp, coords, projection='TAN', proj_point=None, 
					 mode='sip', order=None, inwcs=None):	 
	"""
	Given a set of matched x,y pixel positions and celestial coordinates, solves for
	WCS parameters and returns a WCS with these best fit values and other keywords based 
	on projection type, frame, and units. Optionally, a SIP may be fit. 
	
	Along with the matched coordinate pairs, users must provide the projection type 
	(e.g. 'TAN'), celestial coordinate pair for the projection point (or 'center' to use 
	the center of input coordinates), mode ('wcs' to fit only linear terms, or 'sip' ). 
	If a coordinate pair is passed to 'proj_point', it is assumed to be in the 
	same units respectivley as the input (lon, lat) coordinates. Additionaly, if mode is 
	set to 'sip', the polynomial order must be provided.
	
	Optionally, an existing ~astropy.wcs.WCS object with some may be passed in. For 
	example, this is useful if the user wishes to refit an existing WCS with better
	astrometry, or are using a projection type with non-standard keywords. If any overlap 
	between keyword arguments passed to the function and values in the input WCS exist, 
	the keyword argument will override, including the CD/PC convention. 
		
	NAXIS1 and NAXIS2 will returned as (0, 0). These can be set in the WCS output from 
	this function (wcs._naxis1, wcs._naxis2). 
	
	All output will be in degrees. 
	
	Parameters
	----------
	xp, yp : `numpy.ndarray`
		X and Y pixel coordinates, respectivley. 
	coords : `~astropy.coordinates.SkyCoord`
		Skycoord object.
	inwcs : `~astropy.wcs.WCS`
		Optional input WCS object. Populated keyword values will be used. For any 
		overlap between keyword arguments passed to the function and values in 'inwcs', 
		the keyword argument will be used and replace the value in 'inwcs'.
	projection : str
		Three-letter projection code. Defaults to TAN.
	proj_point: None, str, or tuple of int or float
		Celestial coordinates of projection point (lat, lon). If None, the geometric
		center of input pixel and sky coordinates will be used for the projection.
	mode : 'all' or 'wcs'
		Whether to do the transformation including distortions (``'all'``) or
		only including only the core WCS transformation (``'wcs'``).
	order : int
		Order for SIP polynomial to be fit. 

	Returns
	-------
	wcs : `~astropy.wcs.WCS`
		WCS object with best fit coefficients given the matched set of X, Y pixel and 
		sky positions of sources.
		
	Notes
    -----
    If passed an input wcs, please note this object will be modified. 
    
    Please pay careful attention to the logic that applies when both an input WCS with
    keywords populated is passed in. For example, if no 'CUNIT' is set in this input
    WCS, the units of the CRVAL etc. are assumed to be the same as the input Skycoord.
    Additionally, if 'RADESYS' is not set in the input WCS, this will be taken from the
    Skycoord as well. 

	"""
	from scipy.optimize import least_squares
	from .wcs import Sip
	
	lon, lat = coords.data.lon.deg, coords.data.lat.deg

	if mode not in ("wcs", "sip"): raise ValueError("mode must be 'wcs' or 'sip'")
	if mode == 'sip':
		if (order is None) or (type(order) != int):
			raise ValueError("Must provide integer order for SIP if mode = 'sip'")

	if projection is None:
		if (inwcs is None) or ('' in inwcs.wcs.ctype):
			raise ValueError("Must provide projection type or input WCS with CTYPE.")
		projection = inwcs.wcs.ctype[0].replace('-SIP','').split('-')[-1] 
		
	wcs = celestial_frame_to_wcs(frame = coords.frame, projection = projection)

	close = lambda l, p: p[np.where(np.abs(l) == min(np.abs(l)))[0][0]] 
	if str(proj_point) == "center":
		wcs.wcs.crval = ((max(lon) + min(lon))/2.,(max(lat) + min(lat))/2.)
		wcs.wcs.crpix = ((max(xp) + min(xp))/ 2., (max(yp) + min(yp))/2.) #initial guess
	elif (proj_point is None) and (inwcs is None):
		raise ValueError("Must give proj_point as argument or as CRVAL in input wcs.")
	elif proj_point is not None: #convert units + rough initial guess for crpix for fitter
		lon_u, lat_u = u.Unit(coords.data.lon.unit), u.Unit(coords.data.lat.unit)
		wcs.wcs.crval = (proj_point[0]*lon_u.to(u.deg), proj_point[1]*lat_u.to(u.deg))
		wcs.wcs.crpix = (close(lon-wcs.wcs.crval[0],xp),close(lon-wcs.wcs.crval[0],yp))

	if inwcs is not None: 
		if inwcs.wcs.radesys == '': #no frame specified, use Skycoord's frame (warn???)
			inwcs.wcs.radesys = coords.frame.name
		wcs.wcs.radesys = inwcs.wcs.radesys
		coords.transform_to(inwcs.wcs.radesys.lower()) #work in inwcs units
		if proj_point is None: #crval, crpix from wcs. crpix used for initial guess. 
			wcs.wcs.crval, wcs.wcs.crpix= inwcs.wcs.crval, inwcs.wcs.crpix 
			if wcs.wcs.crpix[0] == wcs.wcs.crpix[1] == 0: #assume 0 wasn't intentional 
				wcs.wcs.crpix = (close(lon-wcs.wcs.crval[0],xp), \
								 close(lon-wcs.wcs.crval[0],yp))
		if inwcs.wcs.has_pc():
			wcs.wcs.pc = inwcs.wcs.pc
		if inwcs.wcs.has_cd():
			wcs.wcs.pc = inwcs.wcs.cd
		wcs_dict = dict(wcs.to_header(relax=False))
		in_wcs_dict = dict(inwcs.to_header(relax=False)) 
		wcs = WCS({**in_wcs_dict,**wcs_dict})

	p0 = np.array((wcs.wcs.pc[0][0], wcs.wcs.pc[0][1], wcs.wcs.pc[1][0], 
				   wcs.wcs.pc[1][1], wcs.wcs.crpix[0], wcs.wcs.crpix[1]))
	fit = least_squares(_linear_transformation_fit, p0, 
						args = (lon, lat, xp, yp, wcs, projection))
	wcs.wcs.crpix = np.array(fit.x[4:6])		
	wcs.wcs.pc = np.array(fit.x[0:4].reshape((2,2)))
	
	if mode == "sip": 
		wcs.wcs.ctype = [x + '-SIP' for x in wcs.wcs.ctype]

		#x,y = wcs.wcs_world2pix(lon,lat,1) #'intermediate world coordinates', x & y 
		#x,y = np.dot(wcs.wcs.pc,(x - wcs.wcs.crpix[0],y - wcs.wcs.crpix[1]))
		
		coef_names = ['{0}_{1}'.format(i,j) for i in range(order+1) \
						for j in range(order+1) if (i+j) < 5 and (i+j) > 1]
	
		p0 = np.concatenate((np.array(wcs.wcs.crpix), wcs.wcs.pc.flatten(),\
							 np.zeros(len(coef_names)*2)))
		#fit = least_squares(_sip_fit, p0, args = (xp, yp, x, y, coef_names))
		fit = least_squares(_new_sip_fit, p0, args = (lon, lat, xp, yp, wcs, coef_names))
		wcs.wcs.pc = fit.x[2:6].reshape((2,2))
		wcs.wcs.crpix = fit.x[0:2]

		coef_fit = (list(fit.x[6:6+len(coef_names)]),list(fit.x[6+len(coef_names):]))
		a_vals, b_vals = np.zeros((order+1,order+1)), np.zeros((order+1,order+1))
		for coef_name in coef_names:
			a_vals[int(coef_name[0])][int(coef_name[2])] = coef_fit[0].pop(0)
			b_vals[int(coef_name[0])][int(coef_name[2])] = coef_fit[1].pop(0)

		wcs.sip = Sip(a_vals, b_vals, a_vals * -1., b_vals * -1., wcs.wcs.crpix)
		
		if (inwcs is not None):
			if inwcs.wcs.has_cd():
				wcs.wcs.cd = wcs.wcs.pc
				wcs.wcs.__delattr__('pc')
	return wcs