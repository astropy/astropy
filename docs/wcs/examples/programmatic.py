# Set the WCS information manually by setting properties of the WCS
# object.

import numpy as np
from astropy import wcs
from astropy.io import fits

# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
w.wcs.crpix = [-234.75, 8.3393]
w.wcs.cdelt = np.array([-0.066667, 0.066667])
w.wcs.crval = [0, -90]
w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
w.wcs.set_pv([(2, 1, 45.0)])

# Three pixel coordinates of interest.
# The pixel coordinates are pairs of [X, Y].
# The "origin" argument indicates whether the input coordinates
# are 0-based (as in Numpy arrays) or
# 1-based (as in the FITS convention, for example coordinates
# coming from DS9).
pixcrd = np.array([[0, 0], [24, 38], [45, 98]], dtype=np.float64)

# Convert pixel coordinates to world coordinates.
# The second argument is "origin" -- in this case we're declaring we
# have 0-based (Numpy-like) coordinates.
world = w.wcs_pix2world(pixcrd, 0)
print(world)

# Convert the same coordinates back to pixel coordinates.
pixcrd2 = w.wcs_world2pix(world, 0)
print(pixcrd2)

# These should be the same as the original pixel coordinates, modulo
# some floating-point error.
assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6

# The example below illustrates the use of "origin" to convert between
# 0- and 1- based coordinates when executing the forward and backward
# WCS transform.
x = 0
y = 0
origin = 0
assert (w.wcs_pix2world(x, y, origin) ==
        w.wcs_pix2world(x + 1, y + 1, origin + 1))

# Now, write out the WCS object as a FITS header
header = w.to_header()

# header is an astropy.io.fits.Header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
hdu = fits.PrimaryHDU(header=header)
# Save to FITS file
# hdu.writeto('test.fits')
