# Create a planetary WCS structure

from astropy import units as u
from astropy.coordinates import BaseBodycentricRepresentation, BaseCoordinateFrame
from astropy.wcs.utils import celestial_frame_to_wcs


class MARSCustomBodycentricRepresentation(BaseBodycentricRepresentation):
    _equatorial_radius = 3399190.0 * u.m
    _flattening = 0.5886 * u.percent


class MARSCustomBodyFrame(BaseCoordinateFrame):
    name = "Mars"


frame = MARSCustomBodyFrame()
frame.representation_type = MARSCustomBodycentricRepresentation
mywcs = celestial_frame_to_wcs(frame, projection="CAR")
print(mywcs.wcs.ctype)
print(mywcs.wcs.name)
print(mywcs.wcs.aux.a_radius)
print(mywcs.wcs.aux.c_radius)
