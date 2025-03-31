from astropy.coordinates import SkyCoord, set_inf_dist
import astropy.units as u

# Enable infinite distance assumption
with set_inf_dist(True):
    coord = SkyCoord(0*u.deg, 0*u.deg, frame='icrs')  # No distance provided
    print("ICRS -> HCRS:", coord.hcrs)  # Now works
