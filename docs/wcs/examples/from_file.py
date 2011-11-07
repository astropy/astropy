# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.

from __future__ import division # confidence high

import numpy
from astropy import wcs
import pyfits
import sys

# Load the FITS hdulist using pyfits
hdulist = pyfits.open(sys.argv[-1])

# Parse the WCS keywords in the primary HDU
w = wcs.WCS(hdulist[0].header)

# Print out the "name" of the WCS, as defined in the FITS header
print w.wcs.name

# Print out all of the settings that were parsed from the header
w.wcs.print_contents()

# Some pixel coordinates of interest.
pixcrd = numpy.array([[0,0],[24,38],[45,98]], numpy.float_)

# Convert pixel coordinates to world coordinates
# The second argument is "origin" -- in this case we're declaring we
# have 1-based (Fortran-like) coordinates.
sky = w.wcs_pix2sky(pixcrd, 1)
print sky

# Convert the same coordinates back to pixel coordinates.
pixcrd2 = w.wcs_sky2pix(sky, 1)
print pixcrd2

# These should be the same as the original pixel coordinates, modulo
# some floating-point error.
assert numpy.max(numpy.abs(pixcrd - pixcrd2)) < 1e-6
