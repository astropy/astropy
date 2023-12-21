# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np
from matplotlib.lines import Path

from astropy.coordinates import angular_separation

# Tolerance for WCS round-tripping, relative to the scale size
ROUND_TRIP_RTOL = 1.0

# Tolerance for discontinuities relative to the median
DISCONT_FACTOR = 10.0


def get_lon_lat_path(lon_lat, pixel, lon_lat_check):
    """
    Draw a curve, taking into account discontinuities.

    Parameters
    ----------
    lon_lat : ndarray
        The longitude and latitude values along the curve, given as a (n,2)
        array.
    pixel : ndarray
        The pixel coordinates corresponding to ``lon_lat``
    lon_lat_check : ndarray
        The world coordinates derived from converting from ``pixel``, which is
        used to ensure round-tripping.
    """
    # In some spherical projections, some parts of the curve are 'behind' or
    # 'in front of' the plane of the image, so we find those by reversing the
    # transformation and finding points where the result is not consistent.

    sep = angular_separation(
        np.radians(lon_lat[:, 0]),
        np.radians(lon_lat[:, 1]),
        np.radians(lon_lat_check[:, 0]),
        np.radians(lon_lat_check[:, 1]),
    )

    # Define the relevant scale size using the separation between the first two points
    scale_size = angular_separation(
        *np.radians(lon_lat[0, :]), *np.radians(lon_lat[1, :])
    )

    with np.errstate(invalid="ignore"):
        sep[sep > np.pi] -= 2.0 * np.pi

        mask = np.abs(sep > ROUND_TRIP_RTOL * scale_size)

    # Mask values with invalid pixel positions
    mask = mask | np.isnan(pixel[:, 0]) | np.isnan(pixel[:, 1])

    # We can now start to set up the codes for the Path.
    codes = np.zeros(lon_lat.shape[0], dtype=np.uint8)
    codes[:] = Path.LINETO
    codes[0] = Path.MOVETO
    codes[mask] = Path.MOVETO

    # Also need to move to point *after* a hidden value
    codes[1:][mask[:-1]] = Path.MOVETO

    # We now go through and search for discontinuities in the curve that would
    # be due to the curve going outside the field of view, invalid WCS values,
    # or due to discontinuities in the projection.

    # We start off by pre-computing the step in pixel coordinates from one
    # point to the next. The idea is to look for large jumps that might indicate
    # discontinuities.
    step = np.sqrt(
        (pixel[1:, 0] - pixel[:-1, 0]) ** 2 + (pixel[1:, 1] - pixel[:-1, 1]) ** 2
    )

    # We search for discontinuities by looking for places where the step
    # is larger by more than a given factor compared to the median
    # discontinuous = step > DISCONT_FACTOR * np.median(step)
    discontinuous = step[1:] > DISCONT_FACTOR * step[:-1]

    # Skip over discontinuities
    codes[2:][discontinuous] = Path.MOVETO

    # The above missed the first step, so check that too
    if step[0] > DISCONT_FACTOR * step[1]:
        codes[1] = Path.MOVETO

    # Create the path
    return Path(pixel, codes=codes)


def get_gridline_path(world, pixel):
    """
    Draw a grid line.

    Parameters
    ----------
    world : ndarray
        The longitude and latitude values along the curve, given as a (n,2)
        array.
    pixel : ndarray
        The pixel coordinates corresponding to ``lon_lat``
    """
    # Mask values with invalid pixel positions
    mask = np.isnan(pixel[:, 0]) | np.isnan(pixel[:, 1])

    # We can now start to set up the codes for the Path.
    codes = np.zeros(world.shape[0], dtype=np.uint8)
    codes[:] = Path.LINETO
    codes[0] = Path.MOVETO
    codes[mask] = Path.MOVETO

    # Also need to move to point *after* a hidden value
    codes[1:][mask[:-1]] = Path.MOVETO

    # We now go through and search for discontinuities in the curve that would
    # be due to the curve going outside the field of view, invalid WCS values,
    # or due to discontinuities in the projection.

    # Create the path
    return Path(pixel, codes=codes)
