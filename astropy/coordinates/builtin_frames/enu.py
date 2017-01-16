# An East-North-Up right handed coordinate system
# Written by Joshua G. Albert - albert@strw.leidenuniv.nl


from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import numpy as np

import astropy.units as u

from astropy.coordinates.baseframe import (BaseCoordinateFrame, FrameAttribute,
                         TimeFrameAttribute, QuantityFrameAttribute,
                         RepresentationMapping, EarthLocationAttribute,
                              frame_transform_graph)

from astropy.coordinates.transformations import FunctionTransform
from astropy.coordinates.representation import (SphericalRepresentation,
                              UnitSphericalRepresentation,CartesianRepresentation)
from astropy.coordinates import ITRS

class ENU(BaseCoordinateFrame):
    """
    Written by Joshua G. Albert - albert@strw.leidenuniv.nl
    A coordinate or frame in the East-North-Up (ENU) system.  

    This frame has the following frame attributes, which are necessary for
    transforming from ENU to some other system:

    * ``obstime``
        The time at which the observation is taken.  Used for determining the
        position and orientation of the Earth.
    * ``location``
        The location on the Earth.  This can be specified either as an
        `~astropy.coordinates.EarthLocation` object or as anything that can be
        transformed to an `~astropy.coordinates.ITRS` frame.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    east : :class:`~astropy.units.Quantity`, optional, must be keyword
        The east coordinate for this object (``north`` and ``up`` must also be given and
        ``representation`` must be None).
    north : :class:`~astropy.units.Quantity`, optional, must be keyword
        The east coordinate for this object (``east`` and ``up`` must also be given and
        ``representation`` must be None).
    up : :class:`~astropy.units.Quantity`, optional, must be keyword
        The east coordinate for this object (``north`` and ``east`` must also be given and
        ``representation`` must be None).

    Notes
    -----
    This is useful as an intermediate frame between ITRS and UVW for radio astronomy

    """

    frame_specific_representation_info = {
        'cartesian': [RepresentationMapping('x', 'east'),
                      RepresentationMapping('y', 'north'),
                     RepresentationMapping('z','up')],
    }
    
    default_representation = CartesianRepresentation

    obstime = TimeFrameAttribute(default=None)
    location = EarthLocationAttribute(default=None)
    #pressure = QuantityFrameAttribute(default=0, unit=u.hPa)
    #temperature = QuantityFrameAttribute(default=0, unit=u.deg_C)
    #relative_humidity = FrameAttribute(default=0)
    #obswl = QuantityFrameAttribute(default=1*u.micron, unit=u.micron)

    def __init__(self, *args, **kwargs):
        super(ENU, self).__init__(*args, **kwargs)

    @property
    def elevation(self):
        """
        Elevation above the horizon of the direction, in radians
        """
        return np.arctan2(self.up,np.sqrt(self.north**2 + self.east**2))

@frame_transform_graph.transform(FunctionTransform, ITRS, ENU)
def itrs_to_enu(itrs_coo, enu_frame):
    '''Defines the transformation between ITRS and the ENU frame.
    ITRS usually has units attached but ENU does not require units 
    if it specifies a direction.'''
    
    if np.any(itrs_coo.obstime != enu_frame.obstime):
        itrs_coo = itrs_coo.transform_to(ITRS(obstime=enu_frame.obstime))
        
    # if the data are UnitSphericalRepresentation, we can skip the distance calculations
    is_unitspherical = (isinstance(itrs_coo.data, UnitSphericalRepresentation) or
                        itrs_coo.cartesian.x.unit == u.one)
    
    lon, lat, height = enu_frame.location.to_geodetic('WGS84')
    lonrad = lon.to(u.radian).value
    latrad = lat.to(u.radian).value
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
        #don't need to do distance calculation
        p = itrs_coo.cartesian.xyz.value
        diff = p
        penu = R.dot(diff)
    
        rep = CartesianRepresentation(x = u.Quantity(penu[0],u.one,copy=False),
                                     y = u.Quantity(penu[1],u.one,copy=False),
                                     z = u.Quantity(penu[2],u.one,copy=False),
                                     copy=False)
    else:
        p = itrs_coo.cartesian.xyz
        p0 = ITRS(*enu_frame.location.geocentric,obstime=enu_frame.obstime).cartesian.xyz
        diff = (p.T-p0).T
        penu = R.dot(diff)
      
        rep = CartesianRepresentation(x = penu[0],#u.Quantity(penu[0],u.m,copy=False),
                                     y = penu[1],#u.Quantity(penu[1],u.m,copy=False),
                                     z = penu[2],#u.Quantity(penu[2],u.m,copy=False),
                                     copy=False)

    return enu_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransform, ENU, ITRS)
def enu_to_itrs(enu_coo, itrs_frame):
    #p = itrs_frame.cartesian.xyz.to(u.m).value
    #p0 = np.array(enu_coo.location.to(u.m).value)
    #p = np.array(itrs_frame.location.to(u.m).value)
    
    
    lon, lat, height = enu_coo.location.to_geodetic('WGS84')
    sinlat = np.sin(lat.to(u.radian).value)
    coslat = np.cos(lat.to(u.radian).value)
    sinlon = np.sin(lon.to(u.radian).value)
    coslon = np.cos(lon.to(u.radian).value)
    north = [-sinlat*coslon,
                      -sinlat*sinlon,
                      coslat]
    east = [-sinlon,coslon,0]
    up = [coslat*coslon,coslat*sinlon,sinlat]
    R = np.array([east,north,up])
    
    if isinstance(enu_coo.data, UnitSphericalRepresentation) or enu_coo.cartesian.x.unit == u.one:
        diff = R.T.dot(enu_coo.cartesian.xyz)
        p = diff
        rep = CartesianRepresentation(x = u.Quantity(p[0],u.one,copy=False),
                                     y = u.Quantity(p[1],u.one,copy=False),
                                     z = u.Quantity(p[2],u.one,copy=False),
                                     copy=False)
    else:
        diff = R.T.dot(enu_coo.cartesian.xyz)
        p0 = ITRS(*enu_coo.location.geocentric,obstime=enu_coo.obstime).cartesian.xyz
        #print (R,diff)
        p = (diff.T + p0).T
        #print (p)
        rep = CartesianRepresentation(x = p[0],#u.Quantity(p[0],u.m,copy=False),
                                     y = p[1],#u.Quantity(p[1],u.m,copy=False),
                                     z = p[2],#u.Quantity(p[2],u.m,copy=False),
                                     copy=False)

    return itrs_frame.realize_frame(rep)
    
    #return ITRS(*p*u.m,obstime=enu_coo.obstime).transform_to(itrs_frame)
    
@frame_transform_graph.transform(FunctionTransform, ENU, ENU)
def enu_to_enu(from_coo, to_frame):
    # for now we just implement this through ITRS to make sure we get everything
    # covered
    return from_coo.transform_to(ITRS(obstime=from_coo.obstime,location=from_coo.location)).transform_to(to_frame)

if __name__ == '__main__':
    import astropy.coordinates as ac
    import astropy.time as at
    
    #test to see if vector input works both ways
    
    loc = ac.EarthLocation(x=6731*u.km,y=1*u.km,z=1*u.km)
    time = at.Time(1,format='gps')
    enu = ENU(obstime=time,location=loc)
    print("With dim test:")
    enucoords = ac.SkyCoord(east = np.array([0,1])*u.m,
                            north=np.array([0,1])*u.m,
                            up=np.array([0,1])*u.m,frame=enu)
    print (enucoords)
    print(enucoords.transform_to('itrs'))
    print(enucoords.transform_to('itrs').transform_to(enu))
    print("Without dim test:")
    enucoords = ac.SkyCoord(east = np.array([0,1]),
                            north=np.array([0,1]),
                            up=np.array([0,1]),frame=enu)
    print (enucoords)
    print(enucoords.transform_to('itrs'))
    print(enucoords.transform_to('itrs').transform_to(enu))
