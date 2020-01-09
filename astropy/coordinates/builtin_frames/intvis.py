from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import numpy as np

import astropy.units as u

from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                         RepresentationMapping, frame_transform_graph)

from astropy.coordinates.attributes import (CoordinateAttribute,
                         TimeAttribute, EarthLocationAttribute)

from astropy.coordinates.transformations import FunctionTransform
from astropy.coordinates.representation import (SphericalRepresentation,
                              UnitSphericalRepresentation,CartesianRepresentation)
from . import (ITRS,ICRS,AltAz)

class InterferometricVisibility(BaseCoordinateFrame):
    """
    A coordinate or frame in the interferometric visibilities (a.k.a. "uvw"),
    as is common in radio astronomy.
    This system are the local East, North, and Up unit vectors projected
    into the plane of the sky, whose origin is fixed to a direction ``phase``.
    These coordinates are used to measure the interferometric baselines in 
    in radio interferometry.

    This frame has the following frame attributes, which are necessary for
    transforming from InterferometricVisibility to some other system:
    * ``obstime``
        The time at which the observation is taken.  Used for determining the
        position and orientation of the Earth.
    * ``location``
        The location of the center of the radio interferometer on the Earth,
        typically the centroid of all antenna locations, though it can be
        fixed to any particular antenna.  This can be specified either as an
        `~astropy.coordinates.EarthLocation` object or as anything that can be
        transformed to an `~astropy.coordinates.ITRS` frame.
    * ``phase``
        The phase tracking center of the frame.  This can be specified either
        as an (ra,dec) `~astropy.units.Qunatity` or as anything that can be
        transformed to an `~astropy.coordinates.ICRS` frame.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other
        keywords)
    u : :class:`~astropy.units.Quantity`, optional, must be keyword
        The u coordinate for this object (``v`` and ``w`` must also be given
        and ``representation`` must be None).
    v : :class:`~astropy.units.Quantity`, optional, must be keyword
        The v coordinate for this object (``u`` and ``w`` must also be given
        and ``representation`` must be None).
    w : :class:`~astropy.units.Quantity`, optional, must be keyword
        The w coordinate for this object (``u`` and ``v`` must also be given
        and ``representation`` must be None).

    Notes
    -----
    This is useful for radio astronomers.
    It enables simple transformation between various other frames and the
    observers frame during the radio synthesis.
    It enables easily calculating elevation of a field with e.g.
    #>>> coord = astropy.coordinates.SkyCoord(ra=..., dec=...,
                frame=InterferometricVisibility(location=..., obstime=..., phase=...))
    #>>> elevation = coord.transform_to(astropy.coordinates.AltAz).alt
    The local East-North-Up frame coordinates can easily be calculated by 
    transforming to a InterferometricVisibility with ``phase`` pointing at zenith.
    """

    frame_specific_representation_info = {
        'cartesian': [RepresentationMapping('x', 'u'),
                      RepresentationMapping('y', 'v'),
                     RepresentationMapping('z','w')],
    }
    default_representation = CartesianRepresentation

    obstime = TimeAttribute(default=None)
    location = EarthLocationAttribute(default=None)
    phase = CoordinateAttribute(ICRS,default=None)

    def __init__(self, *args, **kwargs):
        super(InterferometricVisibility, self).__init__(*args, **kwargs)

@frame_transform_graph.transform(FunctionTransform, ITRS, InterferometricVisibility)
def itrs_to_uvw(itrs_coo, uvw_frame):
    '''Defines the transformation between ITRS and the InterferometricVisibility frame.'''
    is_unitspherical = (isinstance(itrs_coo.data, UnitSphericalRepresentation) or
                        itrs_coo.cartesian.x.unit == u.one)
    lon, lat, height = uvw_frame.location.to_geodetic('WGS84')
    lst = uvw_frame.obstime.sidereal_time('mean',longitude=lon).to(u.radian).value
    ha = (lst - uvw_frame.phase.ra).to(u.radian).value
    dec = uvw_frame.phase.dec.to(u.radian).value
    lonrad = lon.to(u.radian).value - ha
    latrad = dec
    sinlat = np.sin(latrad)
    coslat = np.cos(latrad)
    sinlon = np.sin(lonrad)
    coslon = np.cos(lonrad)
    north = [-sinlat*coslon,
                      -sinlat*sinlon,
                      coslat]
    east = [-sinlon,coslon,0]
    up = [coslat*coslon,coslat*sinlon,sinlat]
    R = np.array([east,north,up])
    if is_unitspherical:
        p = itrs_coo.cartesian.xyz.value
        diff = p
        penu = R.dot(diff)
        rep = CartesianRepresentation(x = u.Quantity(penu[0],u.one,copy=False),
                                     y = u.Quantity(penu[1],u.one,copy=False),
                                     z = u.Quantity(penu[2],u.one,copy=False),
                                     copy=False)
    else:
        p = itrs_coo.cartesian.xyz
        p0 = ITRS(*uvw_frame.location.geocentric,obstime=uvw_frame.obstime).cartesian.xyz
        diff = (p.T-p0).T
        penu = R.dot(diff)
        rep = CartesianRepresentation(x = penu[0],
                y = penu[1],
                z = penu[2],
                copy=False)

    return uvw_frame.realize_frame(rep)

@frame_transform_graph.transform(FunctionTransform, InterferometricVisibility, ITRS)
def uvw_to_itrs(uvw_coo, itrs_frame):
    lon, lat, height = uvw_coo.location.to_geodetic('WGS84')
    lst = uvw_coo.obstime.sidereal_time('mean',longitude=lon).to(u.radian).value
    ha = (lst - uvw_coo.phase.ra).to(u.radian).value
    dec = uvw_coo.phase.dec.to(u.radian).value
    lonrad = lon.to(u.radian).value - ha
    latrad = dec
    sinlat = np.sin(latrad)
    coslat = np.cos(latrad)
    sinlon = np.sin(lonrad)
    coslon = np.cos(lonrad)
    north = [-sinlat*coslon,
                      -sinlat*sinlon,
                      coslat]
    east = [-sinlon,coslon,0]
    up = [coslat*coslon,coslat*sinlon,sinlat]
    R = np.array([east,north,up])
    if isinstance(uvw_coo.data, UnitSphericalRepresentation) or uvw_coo.cartesian.x.unit == u.one:
        diff = R.T.dot(uvw_coo.cartesian.xyz)
        p = diff
        rep = CartesianRepresentation(x = u.Quantity(p[0],u.one,copy=False),
                                     y = u.Quantity(p[1],u.one,copy=False),
                                     z = u.Quantity(p[2],u.one,copy=False),
                                     copy=False)
    else:
        diff = R.T.dot(uvw_coo.cartesian.xyz)
        p0 = ITRS(*uvw_coo.location.geocentric,obstime=uvw_coo.obstime).cartesian.xyz
        p = (diff.T + p0).T
        rep = CartesianRepresentation(x = p[0],
                y = p[1],
                z = p[2],
                copy=False)
    return itrs_frame.realize_frame(rep)

@frame_transform_graph.transform(FunctionTransform, InterferometricVisibility, InterferometricVisibility)
def uvw_to_uvw(from_coo, to_frame):
    return from_coo.transform_to(ITRS(obstime=from_coo.obstime)).transform_to(to_frame)
