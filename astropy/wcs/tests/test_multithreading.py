# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This module is specifically for tests that check that the WCS functionality
# works correctly when used in a multi-threading context.

import numpy as np
from numpy.testing import assert_allclose
from astropy.wcs import WCS
from multiprocessing.pool import ThreadPool


def test_no_inplace_modify_pixel():

    # Regression test for a bug when converting coordinates in a multi-threaded
    # context - the astropy WCSLIB wrapper would modify the input pixel array
    # inplace to account for the 0-based vs 1-based pixel coordinates, but this
    # caused results to be wrong when using multi-threading because each thread
    # would re-apply the offset.

    wcs = WCS(naxis=2)

    N = 100_000

    # The input array needs to have shape (N, 2) and be of the correct type
    # (float) to avoid any copies being made.
    pixel = np.random.randint(-1000, 1000, N * 2).reshape((N, 2)).astype(float)

    n_threads = 4
    pool = ThreadPool(n_threads)

    def round_trip_transform(pixel):
        world = wcs.wcs_pix2world(pixel, 0)
        pixel = wcs.wcs_world2pix(world, 0)
        return pixel

    # Note that this is a deliberately inefficient example where each thread
    # modifies all coordinates - if the input pixels were split up sensibly the
    # issue wouldn't happen, but we can't assume users will never accidentally
    # convert the same coordinates in different threads.
    results = pool.map(round_trip_transform, (pixel,) * n_threads)

    for pixel2 in results:
        assert_allclose(pixel, pixel2)

    pool.close()
