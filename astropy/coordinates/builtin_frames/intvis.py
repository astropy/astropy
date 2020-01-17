from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import numpy as np

import astropy.units as u

from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                         RepresentationMapping, frame_transform_graph)

from astropy.coordinates.attributes import (CoordinateAttribute, EarthLocationAttribute,
                                            TimeAttribute, QuantityAttribute)

from astropy.coordinates.earth import EarthLocation

from astropy.coordinates.transformations import FunctionTransform
from astropy.coordinates.representation import (SphericalRepresentation,
                              UnitSphericalRepresentation,CartesianRepresentation)

from . import (ITRS, ICRS)
               
class InterferometricVisibility(BaseCoordinateFrame):
    """
    A coordinate or frame in the interferometric visibilities (a.k.a. "uvw"),
    as is common in radio astronomy.
    This system are the local East, North, and Up unit vectors projected
    into the plane of the sky, whose origin is fixed to a direction ``phase``.
    These coordinates are used to measure the projected antenna position
    that can be differenced against others for a baseline measurement.

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
        as an (ra,dec) `~astropy.units.Quantity` or as anything that can be
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

    def __init__(self, *args, **kwargs):
        super(InterferometricVisibility, self).__init__(*args, **kwargs)

@frame_transform_graph.transform(FunctionTransform, ICRS, InterferometricVisibility)
def icrs_to_uvw(icrs_loc, uvw_frame):
    '''Defines the transformation between ICRS and the InterferometricVisibility frame.'''

    lon, lat, height = uvw_frame.location.to_geodetic('WGS84')
    lst = uvw_frame.obstime.sidereal_time('mean', longitude=lon).to(u.radian).value
    ha = lst - icrs_loc.ra.to(u.radian).value
    dec = icrs_loc.dec.to(u.rad)

    sinha = np.sin(ha)
    cosha = np.cos(ha)
    sindec = np.sin(dec)
    cosdec = np.cos(dec)
    east = [sinha, cosha, 0]
    north = [-sindec*cosha,
             sindec*sinha,
             cosdec]
    up = [cosdec*cosha, -cosdec*sinha, sindec]
    R = np.array([east,north,up])
    xyz = [uvw_frame.location.x.value, uvw_frame.location.y.value, uvw_frame.location.z.value]

    uvw = R.dot(xyz)
    rep = CartesianRepresentation(x = uvw[0],
                                  y = uvw[1],
                                  z = uvw[2],
                                  copy=False)

    return uvw_frame.realize_frame(rep)
