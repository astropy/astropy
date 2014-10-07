# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import warnings
from .. import units as u
from ..utils.exceptions import AstropyUserWarning

__doctest_skip__ = ['wcs_to_celestial_frame']

__all__ = ['add_stokes_axis_to_wcs',
           'wcs_to_celestial_frame', 'custom_frame_mappings',
           'wcs_to_celestial_frame', 'celestial_pixel_scale',
           'non_celestial_pixel_scales', 'skycoord_to_pixel',
           'pixel_to_skycoord']


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

    from ..coordinates import FK4, FK4NoETerms, FK5, ICRS, Galactic
    from ..time import Time
    from . import WCSSUB_CELESTIAL

    # Keep only the celestial part of the axes
    wcs = wcs.sub([WCSSUB_CELESTIAL])

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


WCS_FRAME_MAPPINGS = [[_wcs_to_celestial_frame_builtin]]


class custom_frame_mappings(object):
    def __init__(self, mappings=[]):
        if hasattr(mappings, '__call__'):
            mappings = [mappings]
        WCS_FRAME_MAPPINGS.append(mappings)

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        WCS_FRAME_MAPPINGS.pop()


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

        >>> from astropy.wcs.utils import wcs_to_celestial_frame, custom_frame_mappings
        >>> with custom_frame_mappings(my_function):
        ...     wcs_to_celestial_frame(...)

    """
    for mapping_set in WCS_FRAME_MAPPINGS:
        for func in mapping_set:
            frame = func(wcs)
            if frame is not None:
                return frame
    raise ValueError("Could not determine celestial frame corresponding to "
                     "the specified WCS object")


def celestial_pixel_scale(inwcs, allow_nonsquare=False):
    """
    For a WCS, if the pixels are square, return the pixel scale in the spatial
    dimensions

    Parameters
    ----------
    inwcs: `astropy.wcs.WCS`
        The world coordinate system object
    allow_nonsquare : bool
        Return the average of the X and Y scales if True

    Returns
    -------
    scale : float
        The square pixel scale

    Raises
    ------
    ValueError if the pixels are nonsquare and ``allow_nonsquare==False``
    """
    cwcs = inwcs.celestial
    if cwcs.wcs.ctype[0][-3:] != 'CAR':
        warnings.warn("Pixel sizes may very over the image for "
                      "projection class {0}".format(cwcs.wcs.ctype[0][-3:]),
                      AstropyUserWarning)
    scale = (cwcs.pixel_scale_matrix**2).sum(axis=0)**0.5
    if not np.allclose(scale[0],scale[1]):
        if allow_nonsquare:
            warnings.warn("Pixels are not square, using an average pixel scale")
            return np.mean(scale)*u.deg
        else:
            raise ValueError("Pixels are not square: 'pixel scale' is ambiguous")
    # Return a quantity: WCS always stores in degrees
    return scale[0]*u.deg


def non_celestial_pixel_scales(inwcs):
    """
    For a non-celestial WCS, e.g. one with mixed spectral and spatial axes, it
    is still sometimes possible to define a pixel scale.

    Parameters
    ----------
    inwcs: `astropy.wcs.WCS`
        The world coordinate system object

    Returns
    -------
    scale : `numpy.ndarray`
        The pixel scale along each axis
    """

    if inwcs.is_celestial:
        raise ValueError("WCS is celestial, use celestial_pixel_scale instead")

    pccd = inwcs.pixel_scale_matrix

    if np.allclose(np.extract(1-np.eye(*pccd.shape), pccd), 0):
        return np.abs(np.diagonal(pccd))*u.deg
    else:
        raise ValueError("WCS is rotated, cannot determine consistent pixel scales")

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
    """

    from .. import units as u
    from . import WCSSUB_CELESTIAL

    if wcs.has_distortion and wcs.naxis != 2:
        raise ValueError("Can only handle WCS with distortions for 2-dimensional WCS")

    # Keep only the celestial part of the axes, also re-orders lon/lat
    wcs = wcs.sub([WCSSUB_CELESTIAL])

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


def pixel_to_skycoord(xp, yp, wcs, origin=0, mode='all'):
    """
    Convert a set of pixel coordinates into a SkyCoord coordinate.

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

    Returns
    -------
    coords : `~astropy.coordinates.SkyCoord`
        The celestial coordinates
    """

    from .. import units as u
    from . import WCSSUB_CELESTIAL
    from ..coordinates import SkyCoord, UnitSphericalRepresentation

    if wcs.has_distortion and wcs.naxis != 2:
        raise ValueError("Can only handle WCS with distortions for 2-dimensional WCS")

    # Keep only the celestial part of the axes, also re-orders lon/lat
    wcs = wcs.sub([WCSSUB_CELESTIAL])

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

    # Create SkyCoord object
    data = UnitSphericalRepresentation(lon=lon, lat=lat)
    coords = SkyCoord(frame.realize_frame(data))

    return coords
