# Create a planetary WCS structure

from astropy.coordinates import BaseCoordinateFrame
from astropy.wcs.utils import celestial_frame_to_wcs

class IAUMARS2000BodycentricRepresentation(BaseBodycentricRepresentation):
    _equatorial_radius = 3396190.0 * u.m
    _flattening = 0.5886007555512007 * u.percent

class IAUMARS2000BodyFrame(BaseCoordinateFrame):
    name = "Mars"

frame = IAUMARS2000BodyFrame()
frame.representation_type = IAUMARS2000BodycentricRepresentation
mywcs = celestial_frame_to_wcs(frame, projection="CAR")
print(mywcs.wcs.ctype)
print(mywcs.wcs.name)
print(mywcs.wcs.aux.a_radius)
print(mywcs.wcs.aux.c_radius)
