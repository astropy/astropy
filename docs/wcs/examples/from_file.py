# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.

import numpy as np
from astropy import wcs
from astropy.io import fits
import sys


def load_wcs_from_file(filename):
    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(filename)

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)

    # Print out the "name" of the WCS, as defined in the FITS header
    print(w.wcs.name)

    # Print out all of the settings that were parsed from the header
    w.wcs.print_contents()

    # Three pixel coordinates of interest.
    # Note we've silently assumed a NAXIS=2 image here.
    # Note also that the pixel coordinates are pairs of [X, Y],
    # and since WCS built from FITS header automatically has
    # origin set to 1, [0, 0] is actually not inside the image.
    pixcrd = np.array([[0, 0], [24, 38], [45, 98]], dtype=np.float64)

    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)
    print(world)

    # Convert the same coordinates back to pixel coordinates.
    pixcrd2 = w.wcs_world2pix(world, 1)
    print(pixcrd2)

    # These should be the same as the original pixel coordinates, modulo
    # some floating-point error.
    assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6

    # As illustrated below, origin defines whether the origin pixel position
    # is 0- or 1-indexed.
    x_origin = 0
    y_origin = 0
    origin = 0
    assert (w.wcs_pix2world(x_origin, y_origin, origin) ==
            w.wcs_pix2world(x_origin + 1, y_origin + 1, origin + 1))

    # Following the origin logic from above, X=0 and Y=0 are actually
    # not within the image when origin=1.
    # This example passes in X=0 and Y=0 for origin=1, and converts the result
    # back to pixels for origin=0. The final result of (-1, -1) shows that
    # pixels are out-of-bounds.
    print(w.wcs_world2pix(
        *w.wcs_pix2world(x_origin, y_origin, origin + 1), origin))


if __name__ == '__main__':
    load_wcs_from_file(sys.argv[-1])
